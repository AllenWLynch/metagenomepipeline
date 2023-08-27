
from scipy.special import xlogy, gammaln
from scipy.optimize import minimize, Bounds
import numpy as np


def model_log_likelihood(*,y, _lambda, theta):
        return xlogy(y, _lambda) - (y + theta)*np.log(_lambda + theta)\
                + gammaln(y + theta) - gammaln(theta) + theta*np.log(theta) - gammaln(y + 1)


def predict(*, exposure, features, beta):
    return np.exp( np.log(exposure)[:,np.newaxis] + features @ beta[np.newaxis,:].T )\
            .ravel()


def score(*, y, exposure, features, beta, theta):
    
    return model_log_likelihood(
        y = y, 
        _lambda = predict(
            exposure = exposure, 
            features = features, 
            beta = beta), 
        theta = theta).sum()



def _fisher_update_NB_weights(*, y, exposure, features, theta,
                     init_beta = None
):
    
    y = y[:,np.newaxis]
    exposure = exposure[:,np.newaxis]
    n_features = features.shape[1]

    def _objective_jac(params):

        beta = params[np.newaxis,:]
        
        log_lambda = np.log(exposure) + features @ beta.T
        _lambda = np.exp(log_lambda)
        
        # negative binomial likelihood without regularizers
        obj_val = y.T @ log_lambda - (y + theta).T @ np.log(_lambda + theta)
        
        # jacobian
        error = theta/(_lambda + theta) * (y - _lambda)
        
        dL_dbeta = error.T @ features #- 2*beta
        
        jac = np.squeeze(dL_dbeta)
        
        return -obj_val, -jac
    
        
    def _hess(params):
        
        beta = params[np.newaxis,:]
        
        log_lambda = np.log(exposure) + features @ beta.T
        _lambda = np.exp(log_lambda)
        
        w = -theta * _lambda * (y + theta)/np.square(_lambda + theta)
        
        hess = (w * features).T @ features #- 2
        
        return -hess
    
    
    if init_beta is None:
        init_beta = [0.]*n_features
    
    res = minimize(
            _objective_jac,
            init_beta,
            jac = True,
            method = 'newton-cg',
            hess = _hess,
        )
    
    return res.x



def _method_of_moments_update_dispersion(*,
    y, exposure, features, beta, init_theta = None,
    max_iter = 100, tol = 1e-8,
):

    n_samples, n_features = features.shape
    _lambda = predict(exposure = exposure, features = features, beta = beta)

    squared_residuals = np.square(y - _lambda)

    # The following functions work with the inverse of the dispersion parameter, alpha
    # Minimize the pearson residuals given the current estimate of alpha
    def _update_alpha(alpha_hat):

        precision = 1/( _lambda/alpha_hat + np.square(_lambda) )
        
        alpha_hat = 1/(n_samples - n_features) * np.dot(squared_residuals.T, precision).item()

        return alpha_hat

    if init_theta is None:
        alpha = 1/10.
    else:
        alpha = 1/init_theta

    for i in range(max_iter):
        alpha_new = _update_alpha(alpha)
        if np.abs(alpha_new - alpha) <= tol:
            break
        alpha = alpha_new

    return 1/alpha



def fit_NB_regression(
    y, exposure, features, init_beta = None, init_theta = None,
    max_iter_dispersion = 100, tol_dispersion = 1e-8,
    em_iters = 100, loglike_tol = 1e-8,
):

    kwargs = dict(
        y = y,
        exposure = exposure,
        features = features,
    )

    if init_beta is None:
        beta = np.array([0.]*features.shape[1])
    else:
        beta = np.array(init_beta)

    if init_theta is None:
        theta = 10.
    else:
        theta = init_theta

    scores = [score(**kwargs, beta = beta, theta = theta)]
    
    for i in range(em_iters):
        
        beta = _fisher_update_NB_weights(**kwargs, init_beta = beta, theta = theta )
        
        theta = _method_of_moments_update_dispersion(
            **kwargs, beta = beta, init_theta = init_theta,
            max_iter = max_iter_dispersion, tol = tol_dispersion,
        )
        
        curr_score = score(**kwargs, beta = beta, theta = theta)
        
        if curr_score - scores[-1] < loglike_tol:
            break
        
        scores.append(curr_score)

    return (beta, theta), scores
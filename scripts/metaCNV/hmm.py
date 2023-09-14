
import hmmlearn.base
import numpy as np
from nb_regression import fit_NB_regression, predict, model_log_likelihood
import re
from numpy import log2,e
from patsy import dmatrix


def split_features(X, design_info):
    colkeys = dict(zip(design_info.column_names, range(len(design_info.column_names))))
    y = X[:,colkeys['coverage']]
    #exposure = X[:,colkeys['exposure']]

    features = X[:, [colkeys[col] for col in design_info.column_names if not col in ['coverage', 'exposure']]]

    return y, features



class CNVHMM(hmmlearn.base.BaseHMM):

    def __init__(self, n_components=1, 
                 beta = None, 
                 theta = None,
                 startprob_prior=1.0, 
                 transmat_prior=1.0,
                 algorithm="viterbi", random_state=None,
                 n_iter=10, tol=1e-2, verbose=False,
                 params='st',
                 init_params='st',
                 implementation="log",
                 design_info = None):

        super().__init__(n_components,
                         startprob_prior=startprob_prior,
                         transmat_prior=transmat_prior,
                         algorithm=algorithm,
                         random_state=random_state,
                         n_iter=n_iter, tol=tol, verbose=verbose,
                         params=params, init_params=init_params,
                         implementation=implementation)

        self.beta = beta
        self.theta = theta
        self.design_info = design_info
        

    def fit(self, X, lengths = None):
        self.design_info = X.design_info
        return super().fit(X, lengths = lengths)


    def get_PTR(self):
        
        ptr_cols = [ 
                (re.search(r'\[(.*?)\]', col).group(1), i)  
                for i, col in enumerate([col for col in self.design_info.column_names if not col in ['coverage', 'exposure']]) 
                if col.endswith('ori_distance') 
            ]

        return {
            config : -log2(e) * self.beta[i]
            for config, i in ptr_cols
        }


    def split_features(self, X):
        return split_features(X, self.design_info)
    

    @staticmethod
    def get_trasmat_prior(n_components, pseudo_observations = 1000):
        return np.ones((n_components, n_components)) + \
            np.eye(n_components)*(pseudo_observations - n_components)


    def _get_ploidies(self):  
        return np.array(
            [0.001] + [2**p for p in range(-2, self.n_components-3)]
        )
    

    def predict_coverage(self, X, lengths = None):

        _, states = self.decode(X, lengths=lengths)
        y, features = self.split_features(X)

        ploidies = self._get_ploidies()[states]
        
        return ploidies, predict(
                exposure = ploidies,
                features = features,
                beta = self.beta
            )


    def _compute_log_likelihood(self, X):
        
        y, features = self.split_features(X)

        ploidies = self._get_ploidies()[None,:] # 1 x n_states

        mu = predict(
                exposure = np.ones(len(ploidies)), 
                features = features, 
                beta = self.beta
            )[:,None]               # n_samples x 1
        
        _lambda = mu * ploidies  # n_samples x n_states

        return model_log_likelihood(
                y = y[:,None],
                _lambda = _lambda,
                theta = self.theta,
            )
    

    def _generate_sample_from_state(self, state, random_state=None):
        pass

    def _get_n_fit_scalars_per_param(self):
        nc = self.n_components
        nf = self.n_features - 2
        return {
            "s": nc - 1,
            "t": nc * (nc - 1),
        }
    


def get_hmm_model(*,
    beta, theta, 
    training = True,
    transmat = None,
    startprob = None,
    n_components = 7,
    **kw,
):

        init_params = ''
        if not training:
            assert transmat is not None
            assert startprob is not None
        else:
            if transmat is None:
                init_params += 't'
            if startprob is None:
                init_params += 's'


        init_kw = dict(
            beta = beta,
            theta = theta,
            transmat_prior = CNVHMM.get_trasmat_prior(n_components),
            algorithm='map',
            params = 'st' if training else '',
            init_params = init_params,
            n_components=n_components, 
            **kw,
        )

        model = CNVHMM(**init_kw)

        if not transmat is None:
            model.transmat_ = transmat

        if not startprob is None:
            model.startprob_ = startprob

        return model


def regions_to_matrix(regions):

    assert all([col in regions.columns for col in ['chr', 'start', 'end', 'coverage', 'gc', 'ori_distance']])
    
    X = dmatrix(
        'coverage + standardize(gc) + chr:ori_distance + chr', # + ' -1' if n_contigs == 1 else '',
        regions
    )
    contig_lengths = list(regions.groupby('chr').size().values)

    return X, contig_lengths



def run_CMV_EM(
    regions, # dataframe with columns 
    random_state = None,
    n_components = 7,
    n_iters = 20,
    tol = 1e-3,
):
    
    hmm_options = dict(n_components=n_components, random_state = random_state)

    def m_step_markov_chain(X, lengths,**params):

        model = get_hmm_model(**params, **hmm_options, training = True)

        model.fit(X, lengths = lengths)

        return model
    

    def e_step_markov_chain(X, lengths, model):
        return model.predict_proba(X, lengths)
    

    def m_step_observation_model(
            X, lengths,*,
            ploidies,
            posterior_probs,
            **params,                         
    ):
        
        y, X = split_features(X, X.design_info)

        expected_ploidy = (posterior_probs @ ploidies[:,None]).ravel()

        return fit_NB_regression(
            y = y,
            features = X,
            exposure = expected_ploidy,
            init_beta = params['beta'],
            init_theta = params['theta'],
        )[0]
    

    def em_iteration(
        X, lengths,**params,  
    ):
        
        model = m_step_markov_chain(
            X, lengths,**params,
        )
        
        posterior_probs = e_step_markov_chain(X, lengths, model)
        
        beta, theta = m_step_observation_model(
            X, lengths,
            ploidies = model._get_ploidies(),
            posterior_probs = posterior_probs,
            **params,
        )

        params = {
            'beta' : beta,
            'theta' : theta,
            'transmat' : model.transmat_,
            'startprob' : model.startprob_,
            }

        score = get_hmm_model(
                    **params, **hmm_options, training = False, design_info = X.design_info,
                ).score(X, lengths = lengths)

        return params, score


    X, contig_lengths = regions_to_matrix(regions)
    design_info = X.design_info

    y, features = split_features(X, X.design_info)
    nonzero_mask = y > 0
    (beta, theta),_ = fit_NB_regression(
            y = y[nonzero_mask],
            features = features[nonzero_mask],
            exposure = np.ones(nonzero_mask.sum()),
        )

    params = {
        'beta' : beta,
        'theta' : theta,
        'transmat' : None,
        'startprob' : None,
    }
    scores = []
    
    for i in range(n_iters):

        params, score = em_iteration(X, contig_lengths, **params)

        print(f'Iteration {i} complete, score: {score}')
        print(f"Params:\n\tBeta: {', '.join(map(str, params['beta']))}\n\tTheta: {str(params['theta'])}")
        
        if i > 1 and score - scores[-1] <= tol:
            break

        scores.append(score)


    model = get_hmm_model(
                **params, **hmm_options, 
                design_info = design_info, 
                training = False,
            )
    
    return model, scores

        
        
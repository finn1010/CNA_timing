#run tetraploidy simulations

import numpy as np
import matplotlib.pyplot as plt
from state_evolve.diploid_evolution import run_simulation_diploid, diploid_beta_vals, diploid_prob_matrix, ss_initialisation, state_initialisation
from state_evolve.tetraploidy_evolution import run_simulation_tetraploidy, tetraploidy_beta_vals, tetraploidy_prob_matrix
from state_evolve.tetraploidy_event import tetraploidy_event, tetraploidy_event_prob
from plot import hist_plot
from helpers.colours import pallet_dip, pallet_tet
from scipy import stats, linalg
    

def beta_convert_params(mu, kappa):
    """
    Convert mean/dispersion parameterization of a beta distribution to the ones
    scipy supports
    """

    if np.any(kappa <= 0):
        raise Exception("kappa must be greater than 0")
    elif np.any(mu <= 0) or np.any(mu >= 1):
        raise Exception("mu must be between 0 and 1")
    
    alpha = kappa * mu 
    beta = kappa * (1- mu)

    return alpha, beta

def beta_rvs(mean, kappa, **kwargs):
    """
    Generate random samples from a beta distribution with mean/dispersion
    specified
    """
    alpha, beta = beta_convert_params(mean, kappa)

    return stats.beta.rvs(alpha, beta, **kwargs)

def rescale_beta(beta, delta, eta):
    """
    Linear transform of beta values from between 0 and 1 to between delta and 
    eta
    """
    return (eta - delta) * beta + delta

def add_noise(beta, delta, eta, kappa):
    """
    Rescale distribution to lie between delta and eta and add beta distributed 
    noise
    """
    beta_rescale = rescale_beta(beta, delta, eta)
 
    return beta_rvs(beta_rescale, kappa)
          


def diploid_to_tetraploidy(mu, gamma, init_fn, num_sites, event_time, patient_age):
    final_diploid_states = run_simulation_diploid(mu, gamma,init_fn, num_sites, start_evoln=0, end_evoln=event_time, initial_state=None)
    beta_vals_before = []
    beta_vals_after = []
    beta_vals_before.append(diploid_beta_vals(final_diploid_states))
    beta_vals_before1 = np.array(beta_vals_before)

    for beta_val in beta_vals_before1:
        noisy_beta_before = add_noise(beta_val, 0.05,0.95,30)
        
    tetraploidy_state_list = []
    for state in final_diploid_states:
        tetraploidy_initial_state = tetraploidy_event(state)
        tetraploidy_states = run_simulation_tetraploidy(mu, gamma, tetraploidy_initial_state, start_evoln=0, end_evoln=patient_age-event_time)
        tetraploidy_state_list.extend(tetraploidy_states)
    beta_vals_after.append(tetraploidy_beta_vals(tetraploidy_state_list))
    beta_vals_after1 = np.array(beta_vals_after)
    for beta_val in beta_vals_after1:
        noisy_beta_after = add_noise(beta_val, 0.05,0.95,30)
    # hist_plot(noisy_beta_before, noisy_beta_after,'Trisomy', event_time, patient_age-event_time, 'e.png')

    return noisy_beta_before, noisy_beta_after

# mu = 0.02                          
# gamma = 0.02 
# diploid_to_tetraploidy(mu, gamma, ss_initialisation, 1000, 10, 60)
#run cnloh evolution

import numpy as np
from scipy import stats, linalg
import matplotlib.pyplot as plt
from state_evolve.diploid_evolution import run_simulation_diploid, diploid_beta_vals, diploid_prob_matrix, state_initialisation, ss_initialisation, ss_init_prob
from state_evolve.cnLOH_event import cnLOH_event, cnLOH_event_prob
from plot import hist_plot
from helpers.colours import pallet_dip

def beta_convert_params(mu, kappa):
    if np.any(kappa <= 0):
        raise Exception("kappa must be greater than 0")
    elif np.any(mu <= 0) or np.any(mu >= 1):
        raise Exception("mu must be between 0 and 1")
    
    alpha = kappa * mu 
    beta = kappa * (1- mu)

    return alpha, beta

def beta_rvs(mean, kappa, **kwargs):
    alpha, beta = beta_convert_params(mean, kappa)
    return stats.beta.rvs(alpha, beta, **kwargs)

def rescale_beta(beta, delta, eta):
    return (eta - delta) * beta + delta

def add_noise(beta, delta, eta, kappa):
    beta_rescale = rescale_beta(beta, delta, eta)
    return beta_rvs(beta_rescale, kappa)
          


def diploid_to_cnLOH_hist(mu, gamma, init_fn, num_sites, event_time, patient_age,fig_name):
    final_diploid_states = run_simulation_diploid(mu, gamma,init_fn, num_sites, start_evoln=0, end_evoln=event_time, initial_state=None)
    beta_vals_before = []
    beta_vals_after = []
    beta_vals_before.append(diploid_beta_vals(final_diploid_states))
    beta_vals_before1 = np.array(beta_vals_before)
    for beta_val in beta_vals_before1:
        noisy_beta_before = add_noise(beta_val, 0.05,0.95,30)
    
    cnLOH_state_list = []
    for state in final_diploid_states:
        cnLOH_initial_state = cnLOH_event(state)
        cnLOH_states = run_simulation_diploid(mu, gamma, num_sites, start_evoln=event_time, end_evoln=patient_age, initial_state=cnLOH_initial_state)
        cnLOH_state_list.extend(cnLOH_states)
    beta_vals_after.append(diploid_beta_vals(cnLOH_state_list))
    beta_vals_after1 = np.array(beta_vals_after)
    
    for beta_val in beta_vals_after1:
        noisy_beta_after = add_noise(beta_val,0.05,0.95,30)
    # hist_plot(beta_vals_before, beta_vals_after,'cnLOH', event_time, patient_age-event_time,fig_name)

    hist_plot(noisy_beta_before, noisy_beta_after,'cnLOH', event_time, patient_age-event_time,fig_name)

mu = 0.01                       
gamma = 0.01
event_time = 5
# diploid_to_cnLOH_hist(mu, gamma, state_initialisation, 5000, event_time, 60,f'/Users/finnkane/Desktop/ICR/plots/cnLOH/Hist/hist_t={event_time}')
# event_time = 35
# diploid_to_cnLOH_hist(mu, gamma, state_initialisation, 5000, event_time, 60,f'/Users/finnkane/Desktop/ICR/plots/cnLOH/Hist/hist_t={event_time}')
# event_time = 55
# diploid_to_cnLOH_hist(mu, gamma, state_initialisation, 5000, event_time, 60,f'/Users/finnkane/Desktop/ICR/plots/cnLOH/Hist/hist_t={event_time}')
# diploid_to_cnLOH_hist(mu, gamma, state_initialisation, 1000, 1, 60,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Hist/tau=1')
# diploid_to_cnLOH_hist(mu, gamma, state_initialisation, 1000, 30, 60,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Hist/tau=30')
# diploid_to_cnLOH_hist(mu, gamma, state_initialisation, 1000, 59, 60,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Hist/tau=59')


def diploid_to_cnLOH_prob_dist(initial_state, mu, gamma, event_time, evoln_time, fig_name):

    diploid_evoln_time = np.linspace(0,event_time)
    diploid_probs = diploid_prob_matrix(initial_state, mu, gamma, diploid_evoln_time)
    initial_cnLOH_probs = cnLOH_event_prob(diploid_probs)

    cnLOH_evoln_time = np.linspace(0,evoln_time-event_time)
    cnLOH_probs = diploid_prob_matrix(initial_cnLOH_probs, mu, gamma, cnLOH_evoln_time)
    initial_cnLOH_probs = np.array([cnLOH_event_prob(diploid_probs)])

    methylated_dip = [2, 1, 0]  

    for i in range(diploid_probs.shape[1]):
        plt.plot(diploid_evoln_time, diploid_probs[:, i], label=f'Diploid: {methylated_dip[i]} Methylated Alleles', color=pallet_dip[i])
 
    for i in range(cnLOH_probs.shape[1]):
        plt.plot(cnLOH_evoln_time + event_time, cnLOH_probs[:, i], label=f'Diploid: {methylated_dip[i]} Methylated Alleles', color=pallet_dip[i])

    plt.axvline(x=event_time, color='gray', linestyle='--', label=f'τ={event_time}')

    plt.xlabel('Time (years)',fontsize=13)
    plt.ylabel('Probability',fontsize=13)
    plt.legend(loc='lower right', fontsize=11)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid()
    plt.savefig(f'{fig_name}.pdf', format='pdf', dpi=300)
    plt.show()

# initial_state = ss_init_prob(0.02,0.02)
# print(initial_state)
# event_time = 10
# evoln_time = 60
# diploid_to_cnLOH_prob_dist(initial_state, mu, gamma, event_time, evoln_time,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Prob/cnLOH_prob_tau=10')
# event_time = 50
# evoln_time = 60
# diploid_to_cnLOH_prob_dist(initial_state, mu, gamma, event_time, evoln_time,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Prob/cnLOH_prob_tau=50')

# initial_state = state_initialisation()
# mu=0.01
# gamma=0.01
# event_time = 25
# evoln_time = 50
# diploid_to_cnLOH_prob_dist(initial_state, mu, gamma, event_time, evoln_time,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Prob/cnLOH_prob_tau=10')
# event_time = 50
# evoln_time = 60
# diploid_to_cnLOH_prob_dist(initial_state, mu, gamma, event_time, evoln_time,'/Users/finnkane/Desktop/ICR/plots/cnLOH/Prob/cnLOH_prob_tau=50')


def sim1(initial_state, mu, gamma, t):

    RateMatrix = np.array([
        [-2 * mu, gamma, 0], 
        [2 * mu, -(mu + gamma), 2 * gamma], 
        [0, mu, -2 * gamma]
    ])

    Probabilities = linalg.expm(RateMatrix * t) @ initial_state

    return Probabilities


def diploid_to_cnLOH_prob_dist1(initial_state, kappa, mu, gamma, event_time, evoln_time, n):
    beta_vals_after = []
    diploid_probs = sim1(initial_state, mu, gamma, event_time)
    initial_cnLOH_probs = cnLOH_event_prob(diploid_probs)
    cnLOH_probs = sim1(initial_cnLOH_probs, mu, gamma, (evoln_time - event_time))
    rng = np.random.default_rng()
    sampled_states = rng.multinomial(1, cnLOH_probs, size=n)
    beta_vals_after.append(diploid_beta_vals(sampled_states))
    beta_vals_after = np.array(beta_vals_after)
    for beta_val in beta_vals_after:
        noisy_beta_after = add_noise(beta_val,0.05,0.95,kappa)
    return noisy_beta_after




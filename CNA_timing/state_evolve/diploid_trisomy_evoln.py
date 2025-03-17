#run trisomy simulations
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, linalg

from state_evolve.diploid_evolution import run_simulation_diploid, diploid_beta_vals, diploid_prob_matrix, state_initialisation, ss_initialisation, ss_init_prob
from state_evolve.trisomy_evolution import run_simulation_trisomy, trisomy_beta_vals, trisomy_prob_matrix
from state_evolve.trisomy_event import trisomy_event, trisomy_event_prob
from plot import hist_plot, plot_prob_dist
from helpers.colours import pallet_dip, pallet_tri
              
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

def diploid_to_trisomy_hist(mu, gamma, init_fn, num_sites, event_time, patient_age,fig_name):
    final_diploid_states = run_simulation_diploid(mu, gamma,init_fn, num_sites, start_evoln=0, end_evoln=event_time, initial_state=None)
    beta_vals_before = []
    beta_vals_after = []
    beta_vals_before.append(diploid_beta_vals(final_diploid_states))
    beta_vals_before1 = np.array(beta_vals_before)

    for beta_val in beta_vals_before1:
        noisy_beta_before = add_noise(beta_val, 0.05,0.95,30)

    trisomy_state_list = []
    for state in final_diploid_states:
        trisomy_initial_state = trisomy_event(state)
        trisomy_states = run_simulation_trisomy(mu, gamma, trisomy_initial_state, start_evoln=0, end_evoln=patient_age-event_time)
        trisomy_state_list.extend(trisomy_states)
    beta_vals_after.append(trisomy_beta_vals(trisomy_state_list))
    beta_vals_after1 = np.array(beta_vals_after)
    for beta_val in beta_vals_after1:
        noisy_beta_after = add_noise(beta_val, 0.05,0.95,30)
 
    hist_plot(noisy_beta_before, noisy_beta_after,'Trisomy', event_time, patient_age-event_time,fig_name)

mu = 0.01                          
gamma = 0.01  
# diploid_to_trisomy_hist(mu, gamma, state_initialisation, 1000, 1, 60,'/Users/finnkane/Desktop/ICR/plots/Trisomy/Hist/tau=1')

# diploid_to_trisomy_hist(mu, gamma, state_initialisation, 1000, 30, 60,'/Users/finnkane/Desktop/ICR/plots/Trisomy/Hist/tau=30')
# diploid_to_trisomy_hist(mu, gamma, state_initialisation, 1000, 59, 60,'/Users/finnkane/Desktop/ICR/plots/Trisomy/Hist/tau=59')

def diploid_to_trisomy_prob_dist(initial_state, mu, gamma, event_time, evoln_time, fig_name):

    diploid_evoln_time = np.linspace(0,event_time)
    diploid_probs = diploid_prob_matrix(initial_state, mu, gamma, diploid_evoln_time)
    initial_trisomy_probs = trisomy_event_prob(diploid_probs)

    trisomy_evoln_time = np.linspace(0,evoln_time-event_time)
    trisomy_probs = trisomy_prob_matrix(initial_trisomy_probs, mu, gamma, trisomy_evoln_time)
    initial_trisomy_probs = np.array([trisomy_event_prob(diploid_probs)])
    methylated_dip = [2, 1, 0]  
    methylated_tri = [3, 2, 1, 0]  
    for i in range(diploid_probs.shape[1]):
        plt.plot(diploid_evoln_time, diploid_probs[:, i], label=f'Diploid: {methylated_dip[i]} Methylated Alleles', color=pallet_dip[i])

    for i in range(trisomy_probs.shape[1]):
        plt.plot(trisomy_evoln_time + event_time, trisomy_probs[:, i], label=f'Trisomy: {methylated_tri[i]} Methylated Alleles', color=pallet_tri[i])

    plt.axvline(x=event_time, color='gray', linestyle='--', label=f'τ={event_time}')

    plt.xlabel('Time (years)')
    plt.ylabel('Probability')
    plt.title('Diploid to Trisomy Evolution')   
    plt.legend(loc='upper right', fontsize=9)
    plt.grid()
    plt.savefig(f'{fig_name}.pdf', format='pdf', dpi=300)
    plt.show()
    
# initial_state = ss_init_prob(0.02,0.02)
# event_time = 10
# evoln_time = 60
# diploid_to_trisomy_prob_dist(initial_state, mu, gamma, event_time, evoln_time,'/Users/finnkane/Desktop/ICR/plots/Trisomy/Prob/tri_prob_tau=10')
# event_time = 50
# evoln_time = 60
# diploid_to_trisomy_prob_dist(initial_state, mu, gamma, event_time, evoln_time,'/Users/finnkane/Desktop/ICR/plots/Trisomy/Prob/tri_prob_tau=50')

def sim1(initial_state, mu, gamma, t):

    RateMatrix = np.array([
        [-2 * mu, gamma, 0], 
        [2 * mu, -(mu + gamma), 2 * gamma], 
        [0, mu, -2 * gamma]
    ])

    Probabilities = linalg.expm(RateMatrix * t) @ initial_state

    return Probabilities

def sim2(initial_state, mu, gamma, t):

    RateMatrix = np.array([[-3*mu, gamma, 0, 0], 
                                [3*mu, -(gamma+2*mu), 2*gamma, 0], 
                                [0, 2*mu, -(2*gamma+mu), 3*gamma],
                                [0, 0, mu, -3*gamma]])

    Probabilities = linalg.expm(RateMatrix * t) @ initial_state

    return Probabilities


def diploid_to_trisomy_prob_dist1(initial_state, mu, gamma, event_time, evoln_time, n):
    beta_vals_after = []
    diploid_probs = sim1(initial_state, mu, gamma, event_time)
    rng = np.random.default_rng()
    initial_trisomy_probs = trisomy_event_prob(diploid_probs)
    trisomy_probs = sim2(initial_trisomy_probs, mu, gamma, (evoln_time - event_time))
    x = rng.multinomial(1, trisomy_probs, size=n)
    beta_vals_after.append(trisomy_beta_vals(x))
    beta_vals_after = np.array(beta_vals_after)
    for beta_val in beta_vals_after:
        noisy_beta_after = add_noise(beta_val,0.05,0.95,30)
    return noisy_beta_after

# mu = 0.02
# gamma = 0.02

# s = diploid_to_trisomy_prob_dist1([0.5,0,0.5], 0.02, 0.02, 10, 40, 1000).flatten()
# print(np.array(s).mean())
# event_times = [10,25,10,15,40,20,10,26,45,18,35,30,23,36,16,31,20,45,25,56]
# patient_ages=[30,40,55,74,70,30,60,72,60,35,75,60,44,65,52,40,28,60,72,68]
# vals = []
# num_sites = [1000]*len(patient_ages)
# for i in range(len(patient_ages)):
#     noisy_beta_after = diploid_to_trisomy_prob_dist1([0.5,0,0.5], mu, gamma, event_times[i], patient_ages[i], num_sites[i])
#     vals.append(noisy_beta_after)

# print(np.array(vals).mean())
# plt.hist(np.array(vals).flatten(),bins=20)
# plt.show()
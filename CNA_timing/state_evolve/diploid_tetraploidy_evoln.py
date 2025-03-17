#run tetraploidy simualtions
import numpy as np
import matplotlib.pyplot as plt
from state_evolve.diploid_evolution import run_simulation_diploid, diploid_beta_vals, diploid_prob_matrix, ss_initialisation, state_initialisation
from state_evolve.tetraploidy_evolution import run_simulation_tetraploidy, tetraploidy_beta_vals, tetraploidy_prob_matrix
from state_evolve.tetraploidy_event import tetraploidy_event, tetraploidy_event_prob
from plot import hist_plot
from helpers.colours import pallet_dip, pallet_tet
from scipy import stats, linalg
   
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
          

def diploid_to_tetraploidy_hist(mu, gamma, init_fn, num_sites, event_time, patient_age,fig_name):
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
        noisy_beta_after = add_noise(beta_val, 0.05,0.95,50)
 
    hist_plot(noisy_beta_before, noisy_beta_after,'tetraploidy', event_time, patient_age-event_time,fig_name)

mu = 0.01                       
gamma = 0.01  
# diploid_to_tetraploidy_hist(mu, gamma, ss_initialisation, 1000, 10, 60,'2')
# diploid_to_tetraploidy_hist(mu, gamma, state_initialisation, 1000, 1, 60,'/Users/finnkane/Desktop/ICR/plots/Tetraploidy/Hist/tau=1')

# diploid_to_tetraploidy_hist(mu, gamma, state_initialisation, 1000, 30, 60,'/Users/finnkane/Desktop/ICR/plots/Tetraploidy/Hist/tau=30')
# diploid_to_tetraploidy_hist(mu, gamma, state_initialisation, 1000, 59, 60,'/Users/finnkane/Desktop/ICR/plots/Tetraploidy/Hist/tau=59')


# def diploid_to_tetraploidy_prob_dist(initial_state, mu, gamma, event_time, evoln_time, fig_name):

#     diploid_evoln_time = np.linspace(0,event_time)
#     diploid_probs = diploid_prob_matrix(initial_state, mu, gamma, diploid_evoln_time)
#     initial_tetraploidy_probs = tetraploidy_event_prob(diploid_probs)

#     tetraploidy_evoln_time = np.linspace(0,evoln_time-event_time)
#     tetraploidy_probs = tetraploidy_prob_matrix(initial_tetraploidy_probs, mu, gamma, tetraploidy_evoln_time)
#     initial_tetraploidy_probs = np.array([tetraploidy_event_prob(diploid_probs)])
    
#     methylated_dip = [2,1,0]
#     for i in range(diploid_probs.shape[1]):
#         plt.plot(diploid_evoln_time, diploid_probs[:, i],label=f'Diploid: {methylated_dip[i]} Methylated Alleles', color=pallet_dip[i])

#     methylated_tet = [4,3,2,1,0]
#     for i in range(tetraploidy_probs.shape[1]):
#         plt.plot(tetraploidy_evoln_time + event_time, tetraploidy_probs[:, i], label=f'Tetraploidy: {methylated_tet[i]} Methylated Alleles', color=pallet_tet[i])

#     plt.axvline(x=event_time, color='gray', linestyle='--', label=f'τ={event_time}')

#     plt.xlabel('Time (years)', fontsize = 14)
#     plt.ylabel('Probability', fontsize = 14)
#     plt.legend(loc='upper right', fontsize=12)
#     plt.grid()
#     plt.savefig(f'{fig_name}.pdf', format='pdf', dpi=300)
#     plt.show()
    
# initial_state = ss_initialisation(mu,gamma)
# mu = 0.1
# gamma = mu
# event_time = 10
# evoln_time = 60
# diploid_to_tetraploidy_prob_dist([1,0,0], mu, gamma, event_time, evoln_time, '/Users/finnkane/Desktop/ICR/plots/Tetraploidy/Prob/tau=10')
# event_time = 30
# evoln_time = 60
# diploid_to_tetraploidy_prob_dist([1,0,0], mu, gamma, event_time, evoln_time, '/Users/finnkane/Desktop/ICR/plots/Tetraploidy/Prob/tau=10')
# event_time = 50
# evoln_time = 60
# diploid_to_tetraploidy_prob_dist([1,0,0], mu, gamma, event_time, evoln_time, '/Users/finnkane/Desktop/ICR/plots/Tetraploidy/Prob/tau=50')
def sim1(initial_state, mu, gamma, t):

    RateMatrix = np.array([
        [-2 * mu, gamma, 0], 
        [2 * mu, -(mu + gamma), 2 * gamma], 
        [0, mu, -2 * gamma]
    ])

    Probabilities = linalg.expm(RateMatrix * t) @ initial_state

    return Probabilities

def sim2(initial_state, mu, gamma, t):

    RateMatrix = np.array([
        [-4*mu, gamma, 0, 0, 0],                
        [4*mu, -(3 * mu + gamma), 2 * gamma, 0, 0], 
        [0, 3 * mu, -(2 * mu + 2*gamma), 3*gamma, 0], 
        [0, 0, 2*mu, -(3*gamma+mu), 4*gamma],
        [0, 0, 0, mu, -4*gamma]           
    ])

    Probabilities = linalg.expm(RateMatrix * t) @ initial_state

    return Probabilities


def diploid_to_tetraploidy_prob_dist1(initial_state, kappa, mu, gamma, event_time, evoln_time, n):
    beta_vals_after = []
    diploid_probs = sim1(initial_state, mu, gamma, event_time)
    rng = np.random.default_rng()
    initial_tetraploidy_probs = tetraploidy_event_prob(diploid_probs)
    tetraploidy_probs = sim2(initial_tetraploidy_probs, mu, gamma, (evoln_time - event_time))
    x = rng.multinomial(1, tetraploidy_probs, size=n)
    beta_vals_after.append(tetraploidy_beta_vals(x))
    beta_vals_after = np.array(beta_vals_after)
    for beta_val in beta_vals_after:
        noisy_beta_after = add_noise(beta_val,0.05,0.95,kappa)
    return noisy_beta_after

# s = diploid_to_tetraploidy_prob_dist1([0.5,0,0.5], 0.02, 0.02, 10, 40, 1000)
# print(s)
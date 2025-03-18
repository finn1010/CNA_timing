#diploid evolution functions

from scipy import linalg
import numpy as np

def state_initialisation():
    rng = np.random.default_rng()
    m = [1,0,0]
    w = [0,0,1]
    if rng.random() > 0.5:
        state = m
    else:
        state = w
    return state

def ss_initialisation(mu, gamma):
    prob = [
    gamma ** 2 / (mu + gamma) ** 2, 
    2 * mu * gamma / (mu + gamma) ** 2, 
    mu ** 2 / (mu + gamma) ** 2
    ]

    states = [
    [1, 0, 0], 
    [0, 1, 0],  
    [0, 0, 1]   
    ]
    
    rng = np.random.default_rng()
    state = rng.choice(states, p=prob)
    
    return state

def ss_init_prob(mu,gamma):
    return [
    gamma ** 2 / (mu + gamma) ** 2, 
    2 * mu * gamma / (mu + gamma) ** 2, 
    mu ** 2 / (mu + gamma) ** 2
    ]



def diploid_prob_matrix(initial_state, mu, gamma, time_points):
    RateMatrix = np.array([[-2*mu, gamma, 0], 
                                [2*mu, -(mu + gamma), 2*gamma], 
                                [0, mu, -2*gamma]])

    Probabilities = np.array([linalg.expm(RateMatrix * t) @ initial_state / np.sum(initial_state) 
                            for t in time_points])
    return Probabilities


def calc_dt_max_diploid(mu, gamma):
    dt_max = 0.1 / np.max((
    2*mu, 
    2*gamma)
    )
    return dt_max 

def diploid_dt(start_evoln, end_evoln, dt_max):
    n = max([int((end_evoln-start_evoln) / dt_max), 5])
    t = np.linspace(start_evoln, end_evoln, n) 
    dt = t[1] - t[0]
    return dt



def state_simulation(initial_state, mu, gamma, dt):
    rng = np.random.default_rng()
    m,k,w = initial_state
    
    p_m_to_k = 2 * mu * dt 
    p_k_to_m = gamma * dt        
    p_k_to_w = mu * dt      
    p_w_to_k = 2 * gamma * dt   

    if m == 1:
        if rng.random() < p_m_to_k:
            m, k, w = 0, 1, 0  

    elif k == 1:
        rand_val = rng.random()
        if rand_val < p_k_to_m:
            m, k, w = 1, 0, 0  
            
        elif rand_val < p_k_to_w + p_k_to_m:
            m, k, w = 0, 0, 1  

        else:
            m, k, w = 0, 1, 0

    elif w == 1:
        if rng.random() < p_w_to_k:
            m, k, w = 0, 1, 0  

    return [m, k, w]

def run_simulation_diploid(mu, gamma, init_fn, num_sites=100, start_evoln=0, end_evoln=10, initial_state=None):

    c = []
    states = []
    final_states = []
    dt_max = calc_dt_max_diploid(mu, gamma)
    dt = diploid_dt(start_evoln, end_evoln, dt_max)
    output_file="simulation_results.txt"
    with open(output_file, "w") as f:
        if initial_state is None:
            for _ in range(num_sites):
                current_state = init_fn(mu, gamma) if init_fn.__code__.co_argcount > 0 else init_fn()
                counter = 0
                for _ in range(int((end_evoln - start_evoln) / dt) + 1):
                    current_state = state_simulation(current_state, mu, gamma, dt)
                    counter += dt
                    
                    states.append(current_state)
                final_states.append(states[-1])
                # f.write(f"Site {_}: Final State: {states[-1]}\n")

        else:
            current_state = initial_state
            counter = 0
            for _ in range(int((end_evoln - start_evoln) / dt) + 1):
                current_state = state_simulation(current_state, mu, gamma, dt)
                states.append(current_state)
                counter += dt
            final_states.append(states[-1])
            # f.write(f"Initial Site: Final State: {states[-1]}\n")

        f.write(f"{counter}")
    return final_states


def cnLOH_event(mkw):
    rng = np.random.default_rng()  
    m_cancer, k_cancer, w_cancer = mkw
    if k_cancer == 1:
        if rng.random() < 0.5:
            m_cancer += k_cancer
            k_cancer = 0
        else:
            w_cancer += k_cancer
        k_cancer = 0
    
    return np.array([m_cancer, k_cancer, w_cancer])


def diploid_beta_vals(states):
    beta_vals = [(state[1] + 2 * state[2]) / 2 for state in states]
    return beta_vals


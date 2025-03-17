# Trisomy evolution functions


from scipy import linalg
import numpy as np

def state_initialisation():
    rng = np.random.default_rng()
    m = [1,0,0,0]
    w = [0,0,0,1]
    if rng.random() > 0.5:
        state = m
    else:
        state = w
    return state

def trisomy_prob_matrix(initial_state, mu, gamma, time_points):
    RateMatrix = np.array([[-3*mu, gamma, 0, 0], 
                                [3*mu, -(gamma+2*mu), 2*gamma, 0], 
                                [0, 2*mu, -(2*gamma+mu), 3*gamma],
                                [0, 0, mu, -3*gamma]])

    Probabilities = np.zeros((len(time_points), len(initial_state)))
    for i, t in enumerate(time_points):
        ProbStates = linalg.expm(RateMatrix * t) @ initial_state
        Probabilities[i] = ProbStates / np.sum(ProbStates) 
    return Probabilities


def calc_dt_trisomy(mu, gamma):
    dt_max = 0.1 / np.max((
    3*mu, 
    3*gamma)
    )
    return dt_max 
  

def tri_dt(start_evoln, end_evoln, dt_max):
    n = max([int((end_evoln-start_evoln) / dt_max), 5])
    t = np.linspace(start_evoln, end_evoln, n) 
    dt = t[1] - t[0]
    return dt

def state_simulation(initial_state, mu, gamma, dt):
    rng = np.random.default_rng()
    m,k,d,w = initial_state
    
    p_m_to_k = 3 * mu * dt if m == 1 else 0  
    p_k_to_m = gamma * dt if k == 1 else 0         
    p_k_to_d = 2*mu * dt if k == 1 else 0      
    p_d_to_k = 2*gamma * dt if d == 1 else 0    
    p_d_to_w = mu * dt if d == 1 else 0      
    p_w_to_d = 3*gamma * dt if w == 1 else 0  


    if m == 1:
        if rng.random() < p_m_to_k:
            m, k, d, w = 0, 1, 0, 0 

    elif k == 1:
        rand_val = rng.random()
        if rand_val < p_k_to_m:
            m, k, d, w = 1, 0, 0, 0
            
        elif rand_val < p_k_to_d + p_k_to_m:
            m, k, d, w = 0, 0, 1, 0

        else:
            m, k, d, w = 0, 1, 0, 0

    elif d == 1:
        rand_val = rng.random()
        if rand_val < p_d_to_k:
            m, k, d, w = 0, 1, 0, 0
            
        elif rand_val < p_d_to_w + p_d_to_k:
            m, k, d, w = 0, 0, 0, 1

        else:
            m, k, d, w = 0, 0, 1, 0
        
    elif w == 1:
        if rng.random() < p_w_to_d:
            m, k, d, w = 0, 0, 1, 0

    return [m, k, d, w]

def run_simulation_trisomy(mu, gamma, initial_state, start_evoln, end_evoln):
    states = [] 
    final_states = []
    dt_max = calc_dt_trisomy(mu, gamma)
    dt = tri_dt(start_evoln, end_evoln, dt_max)
    current_state = initial_state
    for _ in range(int(end_evoln/dt)+1):
        current_state = state_simulation(current_state, mu, gamma,dt)
        states.append(current_state)
    final_states.append(states[-1])
    
    return final_states

def trisomy_beta_vals(states):
    beta_vals = [(state[1] + 2 * state[2] + 3 * state[3]) / 3 for state in states]
    return beta_vals
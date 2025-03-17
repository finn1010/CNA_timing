#Trisomy event functions
import numpy as np

def trisomy_event(mkw):
    rng = np.random.default_rng()  
    m_cancer, k_cancer, w_cancer = mkw
    u_cancer, d_cancer = 0, 0
    if k_cancer == 1:
        if rng.random() < 0.5:
            u_cancer = 1
        else:
            d_cancer = 1
    mudw = np.array([m_cancer, u_cancer, d_cancer, w_cancer])

    return mudw
# def trisomy_event(mkw):
#     rng = np.random.default_rng()  
#     m_cancer, k_cancer, w_cancer = mkw
#     u_cancer = rng.binomial(n = k_cancer, p = 0.5)
#     d_cancer = k_cancer - u_cancer

#     mudw = np.array([m_cancer, u_cancer, d_cancer, w_cancer])

#     return mudw

def trisomy_event_prob(probabilities):
    final_diploid_probs = probabilities
    m_prob, k_prob, w_prob = final_diploid_probs
    u_prob = d_prob = k_prob / 2
    initial_trisomy_probs = [m_prob, u_prob, d_prob, w_prob]
    return initial_trisomy_probs 


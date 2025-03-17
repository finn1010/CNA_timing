#cnloh event functions

import numpy as np

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


def cnLOH_event_prob(probabilities):
    final_diploid_probs = probabilities
    m_prob, k_prob, w_prob = final_diploid_probs
    m_prob = m_prob + 0.5 * k_prob
    w_prob = w_prob + 0.5 * k_prob
    k_prob = 0
    initial_cnLOH_probs = [m_prob, k_prob, w_prob]
    return initial_cnLOH_probs 
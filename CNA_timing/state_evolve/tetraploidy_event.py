#tetraploidy event functions

import numpy as np

def tetraploidy_event(mkw):
    m_cancer, k_cancer, w_cancer = mkw
    u_cancer, d_cancer = 0, 0
    mukdw = np.array([m_cancer, u_cancer, k_cancer, d_cancer, w_cancer])

    return mukdw

def tetraploidy_event_prob(probabilities):
    final_diploid_probs = probabilities
    m_prob,d_prob, w_prob = final_diploid_probs
    k_prob = v_prob = 0
    initial_tetraploidy_probs = [m_prob, k_prob, d_prob, v_prob, w_prob]
    return initial_tetraploidy_probs 
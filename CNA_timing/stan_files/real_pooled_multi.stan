\\pooled model for real data with multiple CNA types

functions {
    vector ss_init_prob(real mu, real gamma) {
        vector[3] probs;
        real x = (mu + gamma)^2;

        probs[1] = gamma^2 / x;
        probs[2] = 2 * mu * gamma / x;
        probs[3] = mu^2 / x;

        return probs;
    }

    matrix dip_rate_matrix(real gamma, real mu) {
        matrix[3, 3] R;
        R[1, 1] = -2 * mu;
        R[1, 2] = gamma;
        R[1, 3] = 0;

        R[2, 1] = 2 * mu;
        R[2, 2] = -(gamma + mu);
        R[2, 3] = 2 * gamma;
        
        R[3, 1] = 0;
        R[3, 2] = mu;
        R[3, 3] = -2 * gamma;

        return R;  
    }
    
    matrix tri_rate_matrix(real mu, real gamma) {
        matrix[4, 4] R;
        R[1, 1] = -3 * mu;
        R[1, 2] = gamma;
        R[1, 3] = 0;
        R[1, 4] = 0;

        R[2, 1] = 3 * mu;
        R[2, 2] = -(gamma + 2 * mu);
        R[2, 3] = 2 * gamma;
        R[2, 4] = 0;

        R[3, 1] = 0;
        R[3, 2] = 2 * mu;
        R[3, 3] = -(2 * gamma + mu);
        R[3, 4] = 3 * gamma;

        R[4, 1] = 0;
        R[4, 2] = 0;
        R[4, 3] = mu;
        R[4, 4] = -3 * gamma;
        
        return R;
    }

    matrix tet_rate_matrix(real mu, real gamma) {
        matrix[5, 5] R;

        R[1, 1] = -4 * mu;
        R[1, 2] = gamma;
        R[1, 3] = 0;
        R[1, 4] = 0;
        R[1, 5] = 0;

        R[2, 1] = 4 * mu;
        R[2, 2] = -(3 * mu + gamma);
        R[2, 3] = 2 * gamma;
        R[2, 4] = 0;
        R[2, 5] = 0;

        R[3, 1] = 0;
        R[3, 2] = 3 * mu;
        R[3, 3] = -(2 * mu + 2 * gamma);
        R[3, 4] = 3 * gamma;
        R[3, 5] = 0;

        R[4, 1] = 0;
        R[4, 2] = 0;
        R[4, 3] = 2 * mu;
        R[4, 4] = -(3 * gamma + mu);
        R[4, 5] = 4 * gamma;

        R[5, 1] = 0;
        R[5, 2] = 0;
        R[5, 3] = 0;
        R[5, 4] = mu;
        R[5, 5] = -4 * gamma;

        return R;
     }

    vector diploid_prob(real t, vector initial_state, real mu, real gamma) {
        matrix[3, 3] R = dip_rate_matrix(mu, gamma);
        matrix[3, 1] initial_state_matrix;
        initial_state_matrix[, 1] = initial_state;
        matrix[3, 1] F = scale_matrix_exp_multiply(t, R, initial_state_matrix);
        return to_vector(F);
    }


    vector tri_prob(real t, vector initial_state, real mu, real gamma) {
        matrix[4, 4] R = tri_rate_matrix(mu, gamma);  
        matrix[4, 1] initial_state_matrix;
        initial_state_matrix[, 1] = initial_state;
        matrix[4, 1] F = scale_matrix_exp_multiply(t, R, initial_state_matrix);
        return to_vector(F);
    }

    vector tet_prob(real t, vector initial_state, real mu, real gamma) {
        matrix[5, 5] R = tet_rate_matrix(mu, gamma); 
        matrix[5, 1] initial_state_matrix;
        initial_state_matrix[, 1] = initial_state;
        matrix[5, 1] F = scale_matrix_exp_multiply(t, R, initial_state_matrix);
        return to_vector(F);
    }

    vector cnLOH_event_prob(vector probabilities) {
        real m_prob = probabilities[1];
        real k_prob = probabilities[2];
        real w_prob = probabilities[3];

        m_prob += 0.5 * k_prob;
        w_prob += 0.5 * k_prob;
        k_prob = 0;

        vector[3] initial_cnLOH_probs;
        initial_cnLOH_probs[1] = m_prob;
        initial_cnLOH_probs[2] = k_prob;
        initial_cnLOH_probs[3] = w_prob;

        return initial_cnLOH_probs;
        }

    vector tri_event_prob(vector probabilities) {
        real m_prob = probabilities[1];
        real k_prob = probabilities[2];
        real w_prob = probabilities[3];

        real u_prob = 0.5 * k_prob;  // Up event
        real d_prob = 0.5 * k_prob;  // Down event

        vector[4] initial_trisomy_probs;
        initial_trisomy_probs[1] = m_prob;
        initial_trisomy_probs[2] = u_prob;
        initial_trisomy_probs[3] = d_prob;
        initial_trisomy_probs[4] = w_prob;

        return initial_trisomy_probs;
    }
    vector tet_event_prob(vector probabilities) {
        real m_prob = probabilities[1];  
        real d_prob = probabilities[2]; 
        real w_prob = probabilities[3];  

        real k_prob = 0;  
        real v_prob = 0; 

        vector[5] initial_tetraploidy_probs;
        initial_tetraploidy_probs[1] = m_prob;
        initial_tetraploidy_probs[2] = k_prob;
        initial_tetraploidy_probs[3] = d_prob;
        initial_tetraploidy_probs[4] = v_prob;
        initial_tetraploidy_probs[5] = w_prob;

        return initial_tetraploidy_probs;
    }

    vector ss_cnloh_evln(real t, real age, real mu, real gamma) {

        vector[3] initial_state = ss_init_prob(mu, gamma);
        vector[3] post_cnLOH_state_probs = cnLOH_event_prob(initial_state);
        vector[3] altered_state_probs = diploid_prob(age - t, post_cnLOH_state_probs, mu, gamma);

        return altered_state_probs;
    }
    vector ss_tri_evln(real t, real age, real mu, real gamma) {

        vector[3] initial_state = ss_init_prob(mu, gamma);
        vector[4] post_tri_state_probs = tri_event_prob(initial_state);
        vector[4] altered_state_probs = tri_prob(age - t, post_tri_state_probs, mu, gamma);

        return altered_state_probs;
    }

    vector ss_tet_evln(real t, real age, real mu, real gamma) {

        vector[3] initial_state = ss_init_prob(mu, gamma);
        vector[5] post_tet_state_probs = tet_event_prob(initial_state);
        vector[5] altered_state_probs = tet_prob(age - t, post_tet_state_probs, mu, gamma);

        return altered_state_probs;
    }

    vector dip_cnloh_evln(real t, real age, real mu, real gamma, vector initial_state) {
       
        vector[3] state_probs = diploid_prob(t, initial_state, mu, gamma);
        vector[3] post_cnloh_state_probs = cnLOH_event_prob(state_probs);
        vector[3] altered_state_probs = diploid_prob(age - t, post_cnloh_state_probs, mu, gamma);

        return altered_state_probs;
    }

    vector dip_tri_evln(real t, real age, real mu, real gamma, vector initial_state) {

        vector[3] state_probs = diploid_prob((t), initial_state, mu, gamma);
        vector[4] post_tri_state_probs = tri_event_prob(state_probs);
        vector[4] altered_state_probs = tri_prob(age - t, post_tri_state_probs, mu, gamma);

        return altered_state_probs;
    }

    vector dip_tet_evln(real t, real age, real mu, real gamma, vector initial_state) {

        vector[3] state_probs = diploid_prob((t), initial_state, mu, gamma);
        vector[5] post_tet_state_probs = tet_event_prob(state_probs);
        vector[5] altered_state_probs = tet_prob(age - t, post_tet_state_probs, mu, gamma);

        return altered_state_probs;
    }
}

data{
    int<lower=0> P;                             // number of patients
    array[P] int<lower=0> J;                    // number of CNAs per patient
    array[sum(J)] int<lower=0> n;               // number of sites per CNA
    array[sum(n)] real<lower=0, upper=1> y;     // observed beta values
    array[sum(J)] int<lower=0> patient_idx;
    array[sum(J)] int<lower=1> K;               // number of mixture components
    array[sum(J)] int<lower=1> type;            // type of CNA
    array[sum(J)] real<lower=0> restrict_ages;  //max cna time
    array[sum(J)] real<lower=0> age;            // patient age

}

transformed data {
    array[sum(J) + 1] int start_idx; 

    start_idx[1] = 1;
    for (j in 1:sum(J)) {
        start_idx[j + 1] = start_idx[j] + n[j];
    }
}


parameters{
    array[sum(J)] real<lower=0> kappa;                //standard deviation of peak
    real<lower=0> mu;
    array[sum(J)] real<lower=0,upper=1> eta;          //Beta dist parameter
    array[sum(J)] real<lower=0,upper=1> delta;        //Beta dist parameter
    real<lower=0> gamma_raw;
    array[sum(J)] real<lower=0, upper=1> t_raw;  

    

}

transformed parameters {
    real<lower=0> gamma;
    gamma = gamma_raw * mu;

    array[sum(J)] real t;

    for (j in 1:sum(J)) {
        t[j] = t_raw[j] * restrict_ages[j]; 
    }
}
    


model {
    mu ~ normal(0,0.01);
    gamma_raw ~ normal(1,0.3);

    for (j in 1:sum(J)) {
        kappa[j] ~ lognormal(3.4, 0.5);
        eta[j]   ~ beta(95, 5);
        delta[j] ~ beta(5, 95);
    }
    
    for (j in 1:sum(J)) {
        int K_val = K[j];

        vector[3] cache_theta_loh;
        vector[4] cache_theta_tri;
        vector[5] cache_theta_tet;

        if (type[j] == 4) {
            cache_theta_loh = dip_cnloh_evln(t[j], age[j], mu, gamma, [0.5,0,0.5]');
            cache_theta_loh /= sum(cache_theta_loh);  

        } else if (type[j] == 5) {
            cache_theta_tri = dip_tri_evln(t[j], age[j], mu, gamma, [0.5,0,0.5]');
            cache_theta_tri /= sum(cache_theta_tri); 

        } else if (type[j] == 6) {
            cache_theta_tet = dip_tet_evln(t[j], age[j], mu, gamma, [0.5,0,0.5]');
            cache_theta_tet /= sum(cache_theta_tet);  
        }

        for (i in start_idx[j] : (start_idx[j + 1] - 1)) {
            if (type[j] == 4) {
                vector[3] log_lik_loh;
                for (k in 1:3) {
                    real pos_obs = (eta[j] - delta[j]) * (k - 1) / 2.0 + delta[j];
                    real a = kappa[j] * pos_obs;
                    real b = kappa[j] * (1 - pos_obs);
                    log_lik_loh[k] = log(cache_theta_loh[k]) + beta_lpdf(y[i] | a, b);
                }
                target += log_sum_exp(log_lik_loh);

            } else if (type[j] == 5) {
                vector[4] log_lik_tri;
                for (k in 1:4) {
                    real pos_obs = (eta[j] - delta[j]) * (k - 1) / 3.0 + delta[j];
                    real a = kappa[j] * pos_obs;
                    real b = kappa[j] * (1 - pos_obs);
                    log_lik_tri[k] = log(cache_theta_tri[k]) + beta_lpdf(y[i] | a, b);
                }
                target += log_sum_exp(log_lik_tri);

            } else if (type[j] == 6) {
                vector[5] log_lik_tet;
                for (k in 1:5) {
                    real pos_obs = (eta[j] - delta[j]) * (k - 1) / 4.0 + delta[j];
                    real a = kappa[j] * pos_obs;
                    real b = kappa[j] * (1 - pos_obs);
                    log_lik_tet[k] = log(cache_theta_tet[k]) + beta_lpdf(y[i] | a, b);
                }
                target += log_sum_exp(log_lik_tet);
            }
        }
    }
}

generated quantities {
    // Posterior predictive samples
    array[sum(J)] real t_pred;
    array[sum(n)] real<lower=0, upper=1> y_rep;
    
    // Separate log-likelihoods for different CNA types
    vector[sum(n)] log_lik_loh;
    vector[sum(n)] log_lik_tri;
    vector[sum(n)] log_lik_tet;

    for (j in 1:sum(J)) {
        t_pred[j] = t_raw[j] * restrict_ages[j];
    }

    for (j in 1:sum(J)) {
        vector[3] cache_theta_loh;
        vector[4] cache_theta_tri;
        vector[5] cache_theta_tet;

        if (type[j] == 4) {
            cache_theta_loh = dip_cnloh_evln(t_pred[j], age[j], mu, gamma, [0.5, 0, 0.5]');
            cache_theta_loh /= sum(cache_theta_loh);
        } else if (type[j] == 5) {
            cache_theta_tri = dip_tri_evln(t_pred[j], age[j], mu, gamma, [0.5, 0, 0.5]');
            cache_theta_tri /= sum(cache_theta_tri);
        } else if (type[j] == 6) {
            cache_theta_tet = dip_tet_evln(t_pred[j], age[j], mu, gamma, [0.5, 0, 0.5]');
            cache_theta_tet /= sum(cache_theta_tet);
        }
        
        for (i in start_idx[j]:(start_idx[j+1] - 1)) {
            int k_sample;
            real pos_obs;
            real beta_a;
            real beta_b;

            if (type[j] == 4) {
                k_sample = categorical_rng(softmax(cache_theta_loh));
                pos_obs = (eta[j] - delta[j]) * (k_sample - 1) / 2.0 + delta[j];
            } else if (type[j] == 5) {
                k_sample = categorical_rng(softmax(cache_theta_tri));
                pos_obs = (eta[j] - delta[j]) * (k_sample - 1) / 3.0 + delta[j];
            } else if (type[j] == 6) {
                k_sample = categorical_rng(softmax(cache_theta_tet));
                pos_obs = (eta[j] - delta[j]) * (k_sample - 1) / 4.0 + delta[j];
            }

            beta_a = kappa[j] * pos_obs;
            beta_b = kappa[j] * (1 - pos_obs);
            y_rep[i] = beta_rng(beta_a, beta_b);

            // Compute separate log-likelihoods
            if (type[j] == 4) {
                vector[3] log_likelihood_loh = log(cache_theta_loh);
                for (k in 1:3) {
                    real pos_k = (eta[j] - delta[j]) * (k - 1) / 2.0 + delta[j];
                    log_likelihood_loh[k] += beta_lpdf(y_rep[i] | kappa[j] * pos_k, kappa[j] * (1 - pos_k));
                }
                log_lik_loh[i] = log_sum_exp(log_likelihood_loh);
                log_lik_tri[i] = not_a_number();
                log_lik_tet[i] = not_a_number();
            } else if (type[j] == 5) {
                vector[4] log_likelihood_tri = log(cache_theta_tri);
                for (k in 1:4) {
                    real pos_k = (eta[j] - delta[j]) * (k - 1) / 3.0 + delta[j];
                    log_likelihood_tri[k] += beta_lpdf(y_rep[i] | kappa[j] * pos_k, kappa[j] * (1 - pos_k));
                }
                log_lik_loh[i] = not_a_number();
                log_lik_tri[i] = log_sum_exp(log_likelihood_tri);
                log_lik_tet[i] = not_a_number();
            } else if (type[j] == 6) {
                vector[5] log_likelihood_tet = log(cache_theta_tet);
                for (k in 1:5) {
                    real pos_k = (eta[j] - delta[j]) * (k - 1) / 4.0 + delta[j];
                    log_likelihood_tet[k] += beta_lpdf(y_rep[i] | kappa[j] * pos_k, kappa[j] * (1 - pos_k));
                }
                log_lik_loh[i] = not_a_number();
                log_lik_tri[i] = not_a_number();
                log_lik_tet[i] = log_sum_exp(log_likelihood_tet);
            }
        }
    }
}


import pickle
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plot import hist_plot
from scipy.special import softmax

with open('dict_loh_cll.pkl', 'rb') as f:
    dict_tri = pickle.load(f)

# with open('dict_tri_mcl.pkl', 'rb') as f:
#     dict_tri = pickle.load(f)
type = 5
selected_dict = dict_tri


num_sites = []
vals = []   
patient_ages = []
patient_idx = []
cna_per_patient = []
restrict_ages = []
cna_types = []
k_vals = []

def extract_data_from_dict_tri(dict_tri):
    cna_per_patient = []
    patient_idx = []
    num_sites = []
    vals = []
    patient_ages = []
    restrict_ages = []
    epi_rates = []
    cna_types = []
    k_vals = []
    filtered_patient_ids = []
    patient_counter = {}

    for patient_id, data in dict_tri.items():
        chr_data = data.get('chr_data', {})
        for chr_values in chr_data.values():
            cna_type = chr_values.get('CNA.type', None)
            # if cna_type == 5:
            beta_values = chr_values.get('beta', [])
            num_sites.append(len(beta_values))
            vals.extend(beta_values)
            patient_ages.append(chr_values.get('age', None))
            restrict_ages.append(chr_values.get('restrictAge', None))
            epi_rates.append(chr_values.get('epiRate', None))
            cna_types.append(cna_type)
            filtered_patient_ids.append(patient_id)
            if patient_id not in patient_counter:
                patient_counter[patient_id] = len(patient_counter) + 1
            patient_idx.append(patient_counter[patient_id])

    nested_key_counts = {
        sample: dict_tri[sample].get('num_unique_chr', 0) 
        for sample in filtered_patient_ids
    }

    for patient_id in patient_counter.keys():
        value = nested_key_counts.get(patient_id, 0)
        cna_per_patient.append(value)

    k_vals = list(np.array(cna_types, dtype=int) - 2) if cna_types else []
    return {
        'cna_per_patient': cna_per_patient,
        'patient_idx': patient_idx,
        'num_sites': num_sites,
        'vals': vals,
        'patient_ages': patient_ages,
        'restrict_ages': restrict_ages,
        'epi_rates': epi_rates,
        'cna_types': cna_types,
        'k_vals': k_vals
    }
extracted_data = extract_data_from_dict_tri(dict_tri)
# print(extracted_data)
# print(extracted_data['patient_idx'], extracted_data['cna_per_patient'])
# print('asdfasdf',sum(extracted_data['cna_per_patient']))
# print(extracted_data['num_sites'],'fghjfhgjfg')
idx = extracted_data['patient_idx']
J = extracted_data['cna_per_patient']
P = max(idx)
cna_types_arr = np.array(extracted_data['cna_types'])
print(print(np.sum(cna_types_arr ==4)))
print(print(np.sum(cna_types_arr ==5)))
print(print(np.sum(cna_types_arr ==6)))
plt.hist(extracted_data['vals'])
plt.show()
unique_epi = list(dict.fromkeys(extracted_data['epi_rates']))
# model = cmdstanpy.CmdStanModel(stan_file='hier_real_mix_fix_sep_gen.stan')
model = cmdstanpy.CmdStanModel(stan_file='hier_real_mixed_fixed.stan')
counter = 0
print(extracted_data['cna_per_patient'])
for num in extracted_data['cna_per_patient']:   
    if num == 2:
        counter+=1
print(min(extracted_data['num_sites'])) 
print(max(extracted_data['num_sites'])) 
print(np.mean(extracted_data['num_sites']))

data = {
    'K': extracted_data['k_vals'], 
    'J': extracted_data['cna_per_patient'],
    'P': max(idx),
    'patient_idx': extracted_data['patient_idx'],
    'n': extracted_data['num_sites'], 
    'y': extracted_data['vals'],  
    'age': extracted_data['patient_ages'],
    'restrict_ages': extracted_data['restrict_ages'],
    'type': extracted_data['cna_types'],
    'mu': unique_epi

}

init_values = lambda: {
    # "mu_mean": np.random.normal(-4.3, 0.2),
    # "mu_sdd": np.random.uniform(0.001, 1),
    # "mu_raw": np.random.uniform(-1, 1, P).tolist(),
    # "gamma_raw": np.random.uniform(0.5, 1.5, P).tolist(),
    "t_raw": np.random.uniform(0.01, 0.99, sum(J)).tolist(),
    "kappa": np.random.uniform(10, 40, sum(J)).tolist(),
    "eta": np.random.uniform(0.8, 1.0, sum(J)).tolist(),
    "delta": np.random.uniform(0.0, 0.2, sum(J)).tolist(),
}

# init_values = [
#     {
#         't': [
#             np.random.uniform(0.01 * patient_ages[j], 0.99 * patient_ages[j])
#             for j in range(len(patient_ages))
#         ]
#     }
#     for _ in range(4)  
# ]

fit = model.sample(
    data=data,
    iter_sampling=1000,  
    iter_warmup=1000,  
    adapt_delta=0.8,
    max_treedepth=12,
    inits=init_values(),
    show_console=True,
    #save_warmup=True,
    # output_dir="output/"
)

summary_df = fit.summary()

import numpy as np

vals = np.array(extracted_data['vals'])
cna_types = np.array(extracted_data['cna_types'])  
num_sites = np.array(extracted_data['num_sites'])  

expanded_cna_types = np.repeat(cna_types, num_sites)  


vals_loh = vals[expanded_cna_types == 4]
vals_tri = vals[expanded_cna_types == 5]
vals_tet = vals[expanded_cna_types == 6]

az_fit = az.from_cmdstanpy(
    fit,
    posterior_predictive="y_rep",  
    observed_data={"y": vals} 
)




az_fit.to_netcdf(f"mixed_fixed2.nc")

with open(f'mixed_fixed2.pkl', 'wb') as f:
    pickle.dump(summary_df, f)
    




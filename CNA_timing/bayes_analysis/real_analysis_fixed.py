import cmdstanpy
import numpy as np
import pandas as pd
from scipy.special import softmax
import pickle
import arviz as az
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from plot import stat_creds, plot_event_time_uncertainty_hier_switched, chr_eval, mu_regplot, mu_inf_vs_evo, plot_event_time_uncertainty, plot_ppc, plot_pair

tri_az = 'mixed_fixed.nc'
tri_sum = 'mixed_fixed.pkl'
# tri_az = '/Volumes/T7/ICR1/real/tri_real(4.3)kap3.4,0.2_filt.nc'
# tri_sum = '/Volumes/T7/ICR1/real/tri_real(4.3)kap3.4,0.2_filt.pkl'
type = 2
# tri_az = 'loh_fixed.nc'
# tri_sum = 'loh_fixed.pkl'

if type == 2:
    az_fit = az.from_netcdf(tri_az)
   
    with open(tri_sum, 'rb') as f:
        summary_df = pickle.load(f)

    with open('dict_mixed.pkl', 'rb') as f:
        dict_tri = pickle.load(f)
    selected_dict = dict_tri 

all_beta_vals = [] 
counter=0

for patient_id, patient_data in selected_dict.items():
    chr_data = patient_data.get('chr_data', {})
    for chr_id, chr_info in chr_data.items():
        beta_vals = chr_info.get('beta', [])
        all_beta_vals.extend(beta_vals)
        print('sss',len(beta_vals))
        counter+=1


observed_data = {"y": np.array(all_beta_vals)}
posterior_predictive = {"y_rep": az_fit.posterior_predictive["y_rep"].values}
new_az_data = az.from_dict(observed_data=observed_data, posterior_predictive=posterior_predictive)

# fig, ax = plt.subplots(figsize=(8, 6))
# az.plot_ppc(new_az_data, data_pairs={"y": "y_rep"}, kind="kde", ax=ax)

# lines = ax.get_lines()  
# for line in lines:
#     if "Posterior predictive mean" in line.get_label():
#         x_data, y_data = line.get_data()
#         mask = x_data > 0.001
#         line.set_data(x_data[mask], y_data[mask])

# ax.set_xlabel("Fraction of alleles methylated", fontsize=14.5)
# ax.set_ylabel("Probability Density", fontsize=14.5)
# plt.legend(loc='upper center',fontsize=12)
# plt.savefig('ppc', dpi=300, bbox_inches="tight")
# plt.show()




# num_sites = []
# vals = []
# patient_ages = []
# patient_idx = []
# cna_per_patient = []

# def count_nested_keys(data):
#     return {sample: len(sample_data) for sample, sample_data in data.items()}

# nested_key_counts = count_nested_keys(selected_dict)

# for count, value in enumerate(nested_key_counts.values(), 1):
#     cna_per_patient.append(value)
#     for i in range(value):
#         patient_idx.append(count)

# for patient, chr_data in selected_dict.items():
#     for chr_values in chr_data.values():
#         beta_values = chr_values['beta']
#         num_sites.append(len(beta_values))
#         vals.extend(beta_values)
#         patient_ages.append(chr_values['age'])

# with open('/Users/finnkane/Desktop/ICR/dict_tri.pkl', 'rb') as f:
#         d = pickle.load(f)
# df_csv = pd.read_csv('CLL-fCpGs-by-CNA-v2.csv')
# df_csv.columns = df_csv.columns.str.strip()
# print(df_csv.columns)
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

ages_list = []

for patient_id, patient_data in dict_tri.items():
    for chr_id, chr_data in patient_data['chr_data'].items():
        ages_list.append(chr_data['age'])
print(ages_list)



unique_chr_counts = {patient_id: data['num_unique_chr'] for patient_id, data in dict_tri.items()}

for patient_id, num_unique_chr in unique_chr_counts.items():
    print(f"Patient {patient_id}: {num_unique_chr} unique chromosomes")



kappa_values = summary_df[summary_df.index.str.contains(r'^kappa\[\d+\]$')]['Mean'].values

sample_chr_counts = list(unique_chr_counts.items())   
split_kappa_values = np.split(
    kappa_values,
    np.cumsum([count for _, count in sample_chr_counts])[:-1]
)

epi_rates = [
    chr_data['epiRate']
    for sample_id, sample_data in dict_tri.items()
    if isinstance(sample_data, dict) and 'chr_data' in sample_data
    for chr_data in sample_data['chr_data'].values()
    if isinstance(chr_data, dict) and 'epiRate' in chr_data
]

epi_rates_unique = list(dict.fromkeys(epi_rates))

sample_ids = list(selected_dict.keys())
epi_rates_dict = {
    sample_id: epi_rates_unique[i]
    for i, sample_id in enumerate(sample_ids)
}

mu_summary = summary_df[summary_df.index.str.contains(r'^mu\[\d+\]$')]


mu_summary['Sample_id'] = sample_ids[:len(mu_summary)]
mu_dict = mu_summary.set_index('Sample_id')[['Mean', '5%', '95%']].to_dict(orient='index')






final_dict = {}

for (sample_id, nchr), kappa_vals in zip(sample_chr_counts, split_kappa_values):
    sample_data = dict_tri[sample_id]
    chr_data_dict = sample_data['chr_data'] 


    chr_entries = list(chr_data_dict.values()) 

    per_chr_info = []
    for kappa_val, chr_info in zip(kappa_vals, chr_entries):
        per_chr_info.append({
            'chr':       chr_info['chr'],
            'kappa':     kappa_val,
            'epiRate':   chr_info['epiRate'],
            'age':       chr_info['age'],
            'num_sites': chr_info['num_sites'],
            'cancerAge': chr_info['cancerAge'],
            'restrictAge': chr_info['restrictAge']
        })

    final_dict[sample_id] = {
        'epiRate':  epi_rates_dict.get(sample_id, None),  
        'chr_data': per_chr_info,
        'mu':       mu_dict.get(sample_id, {}).get('Mean', None),
        'mu_lower': mu_dict.get(sample_id, {}).get('5%', None),
        'mu_upper': mu_dict.get(sample_id, {}).get('95%', None),

    }

rows = []
sample_idx = 0

for sample_id, sample_info in final_dict.items():
    sample_idx += 1 
    for chr_info in sample_info['chr_data']:
        rows.append({
            'Sample_id':  sample_id,
            'index':      sample_idx,
            'chr':        chr_info['chr'],
            'num_kappa': chr_info['kappa'],
            'epiRate': chr_info['epiRate'], 
            'mu':       mu_dict.get(sample_id, {}).get('Mean', None),
            'mu_lower': mu_dict.get(sample_id, {}).get('5%', None),
            'mu_upper': mu_dict.get(sample_id, {}).get('95%', None),
        })

df_final = pd.DataFrame(rows)


patient_errors = {}

unique_true_mu = []
unique_inferred_mu = []



absolute_errors = [abs(inferred - true) for inferred, true in zip(unique_inferred_mu, unique_true_mu)]
mae = sum(absolute_errors) / len(absolute_errors) if absolute_errors else 0

print(f"Mean Absolute Error (MAE) Across Patients: {mae:.6f}")


t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]
eta = summary_df[summary_df.index.str.contains(r'^eta\[\d+\]$')]['Mean'].values
delta  = summary_df[summary_df.index.str.contains(r'^delta\[\d+\]$')]['Mean'].values
print(t_summary)
t = t_summary['Mean']

df_final['eta'] = eta
df_final['delta'] = delta
df_final['age'] = ages_list
df_final['t_mean'] = list(t)
df_final['bounds'] = list(t_summary['StdDev'])

# chr12_df = df_final[df_final['chr'] == 'chr12']
# print(chr12_df)
print(len(list(t)))
# chr_eval('chr_dist.pdf',df_final)

                                                                                                                                                                                                                                                                
t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]
kappa_summary = summary_df[summary_df.index.str.contains(r'^kappa\[\d+\]$')]




print(df_final)

event_times = df_final['t_mean']
patient_ages = df_final['age']
patient_idx = df_final['index']
chr = df_final['chr']
def calculate_bayesian_p_values(az_fit, patient_idx, event_types):

    event_time_samples = az_fit.posterior['t']  
    

    event_dim = 't_dim_0' if 't_dim_0' in event_time_samples.dims else 'event'

    patient_idx = np.array(patient_idx)
    event_types = np.array(event_types)

    unique_patients = np.unique(patient_idx)

    for patient in unique_patients:
        patient_mask = patient_idx == patient
        patient_event_types = event_types[patient_mask]
        patient_event_indices = np.where(patient_mask)[0]

        if len(patient_event_indices) > 1:
            print(f"\nPatient {patient}: Bayesian p-val")

            for i in range(len(patient_event_indices)):
                for j in range(i + 1, len(patient_event_indices)):
                    idx_i, idx_j = patient_event_indices[i], patient_event_indices[j]
                    event_i, event_j = patient_event_types[i], patient_event_types[j]

                    posterior_i = event_time_samples.isel({event_dim: idx_i}).values.flatten()
                    posterior_j = event_time_samples.isel({event_dim: idx_j}).values.flatten()

                    bayesian_p_value = np.mean(posterior_i > posterior_j)

                    print(f"  {event_i} vs {event_j}: p = {bayesian_p_value:.4f}")

calculate_bayesian_p_values(az_fit, patient_idx, chr)

# mu_regplot('reg_mixed',df_final)
a = stat_creds(az_fit, event_times, patient_ages, patient_idx, chr)
print(a)
event_types = np.array(extracted_data['cna_types'])

plot_event_time_uncertainty_hier_switched(
    't.png',
    az_fit,
    event_times,
    patient_ages,
    patient_idx,
    event_types,
    extracted_data,
    df_final
)

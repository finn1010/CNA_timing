import pickle
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plot import hist_plot
from scipy.special import softmax
# import logging
# logging.basicConfig(level=logging.DEBUG)

# with open('dict_loh.pkl', 'rb') as f:
#     dict_loh = pickle.load(f)
with open('dict_tri_cll.pkl', 'rb') as f:
    dict_tri = pickle.load(f)

# with open('dict_tri_mcl.pkl', 'rb') as f:
#     dict_tri = pickle.load(f)
type = 5
selected_dict = dict_tri
K=3


num_sites = []
vals = []
patient_ages = []
patient_idx = []
cna_per_patient = []
restrict_ages = []

def extract_data_from_dict_tri(dict_tri):
    cna_per_patient = []
    patient_idx = []
    num_sites = []
    vals = []
    patient_ages = []
    restrict_ages = []
    epi_rates = []

    nested_key_counts = {
        sample: data.get('num_unique_chr', 0) 
        for sample, data in dict_tri.items()
    }

    for count, value in enumerate(nested_key_counts.values(), start=1):
        cna_per_patient.append(value)
        patient_idx.extend([count] * value)

    for patient, data in dict_tri.items():
        chr_data = data.get('chr_data', {})
        for chr_values in chr_data.values():
            beta_values = chr_values.get('beta', [])
            num_sites.append(len(beta_values))
            vals.extend(beta_values)
            print(len(beta_values))
            patient_ages.append(chr_values.get('age', None))
            restrict_ages.append(chr_values.get('restrictAge', None))
            epi_rates.append(chr_values.get('epiRate', None))

    return {
        'cna_per_patient': cna_per_patient,
        'patient_idx': patient_idx,
        'num_sites': num_sites,
        'vals': vals,
        'patient_ages': patient_ages,
        'restrict_ages': restrict_ages,
        'epi_rates': epi_rates
    }

extracted_data = extract_data_from_dict_tri(dict_tri)
# print(extracted_data)
# print(extracted_data['patient_idx'], extracted_data['cna_per_patient'])
print('asdfasdf',sum(extracted_data['cna_per_patient']))
print(extracted_data['num_sites'],'fghjfhgjfg')
idx = extracted_data['patient_idx']
J = extracted_data['cna_per_patient']
P = max(idx)




def plot_cna_histogram(dict_tri, cna_number):

    cna_beta_values = []
    
    cna_counter = 0
    for patient, data in dict_tri.items():
        for chr_values in data['chr_data'].values():
            cna_counter += 1
            if cna_counter == cna_number:
                print(chr_values.get('sample_chr-id', 'N/A'))

                cna_beta_values.extend(chr_values['beta'])
                break 
    if not cna_beta_values:
        print(f"No beta values found for CNA {cna_number}. Ensure the CNA index is within range.")
        return

    plt.figure(figsize=(8, 5))
    plt.hist(cna_beta_values, bins=13, alpha=0.7, color='blue', edgecolor='black', density=True)
    plt.xlabel("Fraction of alleles methylated", fontsize=14.5)
    plt.ylabel("Probability Density", fontsize=14.5)
#     plt.title(f"Histogram of Beta Values for CNA {cna_number} (Sample Chr-ID: {sample_chr_id})")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim(0,1)
    plt.show()
for i in range(1, 170):
    plot_cna_histogram(dict_tri, i)
# plot_cna_histogram(dict_tri, 50)
# def plot_cna_histogram(dict_tri, cna_number, pdf_pages):
 
#     cna_beta_values = []
#     sample_chr_id = None
    
#     cna_counter = 0
#     for patient, data in dict_tri.items():
#         for chr_values in data['chr_data'].values():
#             cna_counter += 1
#             if cna_counter == cna_number:
#                 sample_chr_id = chr_values.get('sample_chr-id', 'N/A')
#                 # if sample_chr_id in target_chr_ids:
#                 #     cna_beta_values.extend(chr_values['beta'])
#                 # break  

#     # if not cna_beta_values or sample_chr_id not in target_chr_ids:
#         # return

#     plt.figure(figsize=(8, 5))
#     plt.hist(cna_beta_values, bins=40, alpha=0.7, color='blue', edgecolor='black', density=False)
#     plt.xlabel("Beta Value")
#     plt.ylabel("Density")
#     plt.xlim(0, 1)
#     plt.title(f"Histogram of Beta Values for CNA {cna_number} (Sample Chr-ID: {sample_chr_id})")
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     pdf_pages.savefig() 
#     plt.close()  

# with open('dict_tri_cll1.pkl', 'rb') as f:
#     dict_tri = pickle.load(f)

# with PdfPages('filtered_cna_histograms.pdf') as pdf_pages:


# print("PDF with filtered histograms has been created: 'filtered_cna_histograms.pdf'")


model = cmdstanpy.CmdStanModel(stan_file='hierarchical_real.stan')
data = {
    'K': K, 
    'J': extracted_data['cna_per_patient'],
    'P': max(idx),
    'patient_idx': extracted_data['patient_idx'],
    'n': extracted_data['num_sites'], 
    'y': extracted_data['vals'],  
    'age': extracted_data['patient_ages'],
    'restrict_ages': extracted_data['restrict_ages'],
    'type': type,
}

init_values = lambda: {
    "mu_mean": np.random.normal(-4.3, 0.2),
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

az_fit = az.from_cmdstanpy(
    fit,
    posterior_predictive="y_rep",  
    observed_data={"y": vals}      
)
az_fit.to_netcdf(f"cll_tri.nc")

with open(f'cll_tri.pkl', 'wb') as f:
    pickle.dump(summary_df, f)
    


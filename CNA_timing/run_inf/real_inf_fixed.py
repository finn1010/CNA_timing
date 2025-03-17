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

type = 4
selected_dict = dict_tri



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

unique_epi = list(dict.fromkeys(extracted_data['epi_rates']))
print(unique_epi)
model = cmdstanpy.CmdStanModel(stan_file='hier_real_mixed_fixed.stan')
K=[3]*13
type = [4]*13
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
az_fit.to_netcdf(f"loh_fixed.nc")

with open(f'loh_fixed.pkl', 'wb') as f:
    pickle.dump(summary_df, f)
    



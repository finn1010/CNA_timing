import pandas as pd
import pickle

target_chr_ids = {
    # "SCLL-030_chr2:224054-92326171",
    # "SCLL-148_chr12:145740-133779461",
    # "SCLL-182_chr2:10001-92285930",
    # "SCLL-287_chr12:145740-133779461",
    # "SCLL-330_chr12:145740-133779461",
    # "SCLL-626_chr12:145740-133779461",
    # "SCLL-642_chr12:145740-133779461",
    # "SCLL-754_chr12:145740-133779461",
    # "SCLL-057_chr20:60001-62965520",
    # "SCLL-1063_chr1:144145784-249233096",
    # "SCLL-1193_chr12:145740-133779461",
    # "SCLL-1295_chr12:145740-133779461",
    # "SCLL-1344_chr2:10001-91609170",
    # "SCLL-1371_chr2:10001-92326171",
    # "SCLL-297_chr12:145740-133779461",
    # "SCLL-356_chr2:16819-89129732",
    # "SCLL-1075_chr12:145740-133779461",
    # "SCLL-297_chr12:145740-133779461",
    # "SCLL-777_chr12:145740-133779461",
    # "SCLL-342_chr12:145740-133779461",
    # "SCLL-329_chr7:67357349-159128663",
    #mcl
    # "M248_chr12:190980-49142960",
    # "M014_chr12:51781754-133832256"
}

df = pd.read_csv('real_data/CLL-fCpGs-by-CNA-v3.csv')
# print(len(df['Name']))
# correlated_fCpGs = pd.read_csv('real_data/fcpgs_noncorrelated.csv')

# correlated_values = set(correlated_fCpGs.iloc[:, 0].astype(str))
# df['Name_1'] = df['Name_1'].astype(str)
# df_filtered = df[df['Name_1'].isin(correlated_values)]
# print(df_filtered)
df_filtered = df

chr_value_counts = df_filtered.groupby(['Sample_id', 'chr'])['beta'].size()

valid_chrs = chr_value_counts[chr_value_counts > 20].reset_index()['chr'].unique()
print(len(df['Name']))

filtered_df = df_filtered[df_filtered['chr'].isin(valid_chrs)]

gain_df = filtered_df[filtered_df['CNA.type'] == 'CNN-LOH']
print(len(gain_df['Name']))

dict_tri_summary = {}

cna_type_mapping = {
    'CNN-LOH': 4,
    'CN_Gain': 5,
    'High_Copy_Gain': 6
}

for sample_id, sample_group in gain_df.groupby('Sample_id'):
    chr_groups = {}
    
    for chr_, chr_group in sample_group.groupby('chr'):
        sample_chr_id = chr_group['sample_chr-id'].iloc[0]

        if sample_chr_id not in target_chr_ids and len(chr_group) >= 20:
            cna_type_str = chr_group['CNA.type'].iloc[0] 
            cna_type_numeric = cna_type_mapping.get(cna_type_str, None)

            chr_groups[chr_] = {
                'sample_chr-id': sample_chr_id,
                'beta': chr_group['beta'].tolist(),
                'chr': chr_group['chr'].iloc[0],
                'age': chr_group['Age'].iloc[0],
                'num_sites': len(chr_group),
                'cancerAge': chr_group['cancerAge'].iloc[0],
                'restrictAge': chr_group['Age'].iloc[0] - chr_group['cancerAge'].iloc[0],
                'epiRate': chr_group['epiRate'].iloc[0],
                'CNA.type': cna_type_numeric  
            }
    
    if chr_groups:
        dict_tri_summary[sample_id] = {
            'num_unique_chr': len(chr_groups),
            'chr_data': chr_groups
        }
# print(dict_tri_summary,'asdfasd')
with open('dict_loh_cll.pkl', 'wb') as f:
    pickle.dump(dict_tri_summary, f)

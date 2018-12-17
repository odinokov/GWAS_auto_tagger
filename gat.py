'''
To tag GWAS studies from http://ftp.ebi.ac.uk/pub/databases/gwas/
according to synonyms of the ranked top leading causes of death
'''

import pandas as pd
import re
from tqdm import tqdm # https://pypi.org/project/tqdm/

df = pd.read_csv('gwas-catalog-associations_ontology-annotated.zip',\
                 compression='zip', header=0, sep='\t', quotechar='"', low_memory=False)
df = df.copy()

GWAS_tags = {
    'Diseases of heart' : ['heart', 'coronary', 'cardiovascular', 'atherosclerosis'],
    'Malignant neoplasms' : ['malignant', 'neoplasms', 'cancer', 'glioma', 'glioblastoma', 'lymphoma', 'leukemia', 'sarcoma', 'melanoma', 'histiocytoma', 'osteosarcoma'],
    'Alzheimer\'s disease' : ['Alzheimer', 'AD', 'Alzheimer\s'],
    'Cerebrovascular diseases' : ['cerebrovascular', 'stroke', 'ischemic attack', 'aneurysms', 'aneurysm'],
    'Chronic lower respiratory diseases' : ['COPD', 'chronic obstructive pulmonary', 'lung function'],
    'Influenza and pneumonia' : ['pneumonia', 'pneumonitis', 'influenza'],
    'Nephritis, nephrotic syndrome and nephrosis' : ['nephritis', 'nephrotic', 'nephrosis', 'kidney disease', 'kidney function', 'renal function'],
    'Diabetes mellitus' : ['diabetes', 'insulin'],
    'Essential hypertension and hypertensive renal disease' : ['hypertension', 'hypertensive', 'blood pressure', 'atherosclerosis']
}

# update GWAS_causes' keys by REGEX template to look for a particular word or phrase in a string
for GWAS_keys, GWAS_values in GWAS_tags.items():
    new_GWAS_values_list = list()
    for GWAS_value in GWAS_values:
        new_GWAS_values_list.append(r'(?:^|\W)(' + GWAS_value.lower() + r')(?:$|\W)')
    GWAS_tags[GWAS_keys] = new_GWAS_values_list
    
# group all strings with trait for synonyms look up
df['TAG'] = df['STUDY'] + ' ' + df['DISEASE/TRAIT'] + ' ' + df['MAPPED_TRAIT']

df['RS_IDs'] = df['SNPS']

# df.columns.get_loc('SNPS')

# any rs...
RS_IDs_REGEX = r'(rs\d+)'    

# copy a header of modified df
df_clean = pd.DataFrame(data=None, columns=df.columns, index=None)

for df_index in tqdm(df.index):
     
    list_of_tags = list()
    
    # get all tags
    for GWAS_keys, GWAS_values in GWAS_tags.items():
        for GWAS_synonym in GWAS_values:
            if bool(re.search(GWAS_synonym, str(df.loc[df_index, 'TAG']).lower())):
                list_of_tags.append(GWAS_keys)
                
    # remove duplicates in the list of tags and store tags
    df.loc[df_index, 'TAG'] = '; '.join(list(set(list_of_tags)))
    
    # create a new record for each individual RS ID
    
    RS_IDs = re.findall(RS_IDs_REGEX, df.loc[df_index, 'RS_IDs'])
    
    if not bool(RS_IDs):
        RS_IDs = [None]
        
    for RS_ID in RS_IDs:
        df.loc[df_index, 'RS_IDs'] = RS_ID
        df_clean = df_clean.append(df.loc[[df_index,]], ignore_index=True, sort=False)
        
df_clean.to_csv('gwas_RS_IDs_tagged.csv', index=False)

'''
get a list of unique SNP/TAG pairs
'''

import pandas as pd
from tqdm import tqdm

df = pd.read_csv('gwas_RS_IDs_tagged_clean.csv', low_memory=False)
df_TGAG_RS_IDs = pd.DataFrame()

df_TGAG_RS_IDs = df[['TAG', 'RS_IDs']].dropna().reset_index(drop=True)
df_new = pd.DataFrame(data=None, columns=df_TGAG_RS_IDs.columns, index=None)
df_TGAG_RS_IDs = df_TGAG_RS_IDs.copy()
# df_TGAG_RS_IDs.sample(10)

for df_index in tqdm(df_TGAG_RS_IDs.index):
    TAGS_list = df_TGAG_RS_IDs.loc[df_index, 'TAG'].split(';')
    for TAG in TAGS_list:
        df_TGAG_RS_IDs['TAG'].loc[[df_index,]] = TAG.strip()
        df_new = df_new.append(df_TGAG_RS_IDs.loc[[df_index,]], ignore_index=True, sort=False)

df_new = df_new.drop_duplicates()

df_new.to_csv('RS_IDs_and_tags_only.csv', index=False)

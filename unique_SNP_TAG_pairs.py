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

'''
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 11})

df_new = pd.read_csv('RS_IDs_and_tags_only.csv', low_memory=False)

labels = df_new['TAG'].unique().tolist()

sizes = list()

for label in labels:
    sizes.append(df_new[df_new['TAG'] == label]['RS_IDs'].count())

fig1, ax1 = plt.subplots()

ax1.pie(sizes, explode=None, labels=labels, autopct='%1.f%%',
        pctdistance=0.8, labeldistance=1.1, radius = 1,
        shadow=False, startangle=30)

ax1.axis('equal') 

ax1.set_title('The distribution of 14K unique SNP/traits \n\
associated with the top causes of death\n', fontsize=16)

plt.savefig('test.png', dpi = (200), bbox_inches = 'tight')

plt.show()
'''

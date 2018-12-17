'''
get a list of unique SNP/TAG/Gene
'''

import pandas as pd
from tqdm import tqdm

df = pd.read_csv('gwas_RS_IDs_tagged.csv', low_memory=False)

df_TAG_RS_IDs_Gene = pd.DataFrame()

df_TAG_RS_IDs_Gene['TAG'], df_TAG_RS_IDs_Gene['RS_IDs'], df_TAG_RS_IDs_Gene['Gene'] = df['TAG'], df['RS_IDs'], df['MAPPED_GENE']

df_TAG_RS_IDs_Gene = df_TAG_RS_IDs_Gene.dropna(subset=['TAG']).reset_index(drop=True)

print('Total number of unique SNPS:', len(df_TAG_RS_IDs_Gene['RS_IDs'].value_counts()))

df_gene = df_tag = pd.DataFrame(data=None, columns=df_TAG_RS_IDs_Gene.columns, index=None)
df_TAG_RS_IDs_Gene.fillna('', inplace=True)
df_TAG_RS_IDs_Gene.sample(10)

# split the genes

df_ = df_TAG_RS_IDs_Gene.copy()

for df_index in tqdm(df_.index):
    genes_list = str(df_.loc[df_index, 'Gene']).replace(' - ', ', ').replace('; ', ', ').split(', ')
    for gene in genes_list:
        df_['Gene'].loc[[df_index,]] = gene
        df_gene = df_gene.append(df_.loc[[df_index,]], ignore_index=True, sort=False)
        
print('Total number of unique genes:', len(df_gene['Gene'].value_counts()))

# split the tags

df_ = df_gene.copy()

for df_index in tqdm(df_.index):
    tag_list = str(df_.loc[df_index, 'TAG']).split(';')
    for tag in tag_list:
        df_['TAG'].loc[[df_index,]] = tag.strip()
        df_tag = df_tag.append(df_.loc[[df_index,]], ignore_index=True, sort=False)
        
df_ = df_tag.drop_duplicates()

df_['TAG'].value_counts()

df_.to_csv('Tags_RS_IDs_and_genes.csv', index=False)

# plot a pie

import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 11})

# Sorting Pandas dataframe by frequency of TAG occurrence

df_dummy = df_.copy()
df_dummy['count'] = df_dummy.groupby('TAG')['TAG'].transform('count')
df_ = df_dummy.sort_values(by=['count']).reset_index(drop=True).drop(['count'], axis=1)

labels = df_['TAG'].unique().tolist()

sizes = list()

for label in labels:
    sizes.append(df_[df_['TAG'] == label]['RS_IDs'].count())

fig1, ax1 = plt.subplots()

ax1.pie(sizes, explode=None, labels=labels, autopct='%1.f%%',
        pctdistance=0.8, labeldistance=1.1, radius = 1,
        shadow=False, startangle=30)

ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.


ax1.set_title('The distribution of ' + str(len(df_['RS_IDs'])) + ' unique SNP/traits \n\
associated with the top causes of death\n', loc='right', fontsize=16)

# fig1.suptitle('\n © 2018', fontsize=8)

plt.savefig('test.png', dpi = (200), bbox_inches = 'tight')

plt.show()

df_new = df_.drop(['RS_IDs'], axis=1).drop_duplicates()
df_new = df_new.groupby('Gene')['TAG'].count().to_frame().sort_values(by=['TAG'], ascending = False)
df_new.reset_index().to_csv('traits_per_gene.csv', index=False)

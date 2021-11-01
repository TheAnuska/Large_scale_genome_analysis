import numpy as np
import pandas as pd
from collections import Counter

synonymous_codons = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTA', 'CTC', 'CTG', 'CTT'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'Y': ['TAC', 'TAT'],
}

#create list of codons on the base of dictionary, keeps both lists coupled.
list_of_codons = []

for key ,value  in synonymous_codons.items():
    for codon_name in value:
        list_of_codons.append(codon_name)

amino_acid_list = synonymous_codons.keys()


# number of best proteins selected from the database i every iteration
N = 25

MAX_RECORDS = 100000
print("max records", MAX_RECORDS)


def count_codons(data_frame):
    '''
    Separate nucleotide seq on codons, count codons
    :param data_frame
    :return: final_data_frame with the counted codons
    '''
    # split nucleotide seq on codons
    split_string = lambda x, n: [x[i:i + n] for i in range(0, len(x), n)]

    Nucleotide_loc = data_frame.columns.get_loc('CodonSeq')
    nucleotide_list = []
    for k in range(len(data_frame.index)):
        temp = split_string(data_frame.iloc[k, Nucleotide_loc], 3)
        nucleotide_list.append(temp)

    # count codons and add to new data frame
    temp_0 = Counter(nucleotide_list[0]).items()
    temp_df = pd.DataFrame(temp_0).set_index(0).T
    new_data_frame = temp_df

    for m in range(1, len(data_frame.index)):
        temp_m = Counter(nucleotide_list[m]).items()
        temp_m = pd.DataFrame(temp_m).set_index(0).T
        new_data_frame = new_data_frame.append(temp_m, ignore_index=True)

    codons_data_frame = pd.DataFrame(index=list_of_codons).T
    final_data_frame = pd.concat([codons_data_frame, new_data_frame],sort=True).fillna(0)
    final_data_frame = final_data_frame[final_data_frame.columns.intersection(list_of_codons)]
    return final_data_frame


def calculate_desc_CAI(df_withRA):
    '''
    Count CAI for each amino acid
    :param df_withRA: data frame with calculated relative adaptiveness
    :return: data frame with calculated CAI and sorted from the highest values
    '''

    CAI_list = []
    CAI_value = 1
    # n_values = df_counted.sum(axis=1)
    for rowIndex, row in df_withRA.iterrows():
        r = row.to_numpy()
        # number of nonzero elements
        n = np.count_nonzero(r)
        # Product of non zero elements
        p = np.prod(r[r.nonzero()])
        CAI_list.append(p**(1.0/n))

    new_data_frame = pd.DataFrame({'CAI': CAI_list})
    sort_data_frame = new_data_frame.sort_values('CAI', ascending=False)

    return sort_data_frame



def calculate_relative_adaptiveness(counted_codons_df):
    '''
    Count relative adaptiveness for synonymous codons
    :param counted_codons_df
    :return: data frame with calculated relative adaptiveness
    '''
    relative_adaptiveness_df = pd.DataFrame()

    sum_counted_codons_list = counted_codons_df.sum()
    sum_df = pd.DataFrame(sum_counted_codons_list).T
    sum_df = sum_df.replace(to_replace=0, value=1)

    for aa_name in amino_acid_list:
        aa_codon_list = synonymous_codons[aa_name]
        df_temp = sum_df[aa_codon_list]
        max_row = df_temp.max(axis=1)
        temp_div_df = df_temp.div(float(max_row))
        relative_adaptiveness_df = relative_adaptiveness_df.append(temp_div_df)

    relative_adaptiveness_df = pd.DataFrame(relative_adaptiveness_df.sum()).T

    return relative_adaptiveness_df



def selNbest(df_org,sorted_CAI_df,N):
    '''
    selects N proteins from the df_orgwih the highest CAI (taken from df_CAI)
    ensures that all codon are in this group,if not automatically increases N
    :param df_org: data frame with all proteins
    :param sorted_CAI_df: data frame with calculated CAI, sorted from highest values
    :param N: number of proteins
    :return: N proteins with highest CAI value
    '''

    top_CAI_index = (sorted_CAI_df.index.values.tolist())
    top_CAI_index = top_CAI_index[0:N]

    selected_proteins = df_org.iloc[top_CAI_index, :]

    return selected_proteins,top_CAI_index



def final_protein_index_list(counted_codons, init_RA, MAX_RECORDS):
    '''
    It is for search final proteins with optimal codons (top_CAI_index_new), and final weight table (selected_RA_df)
    :param counted_codons: data frame with counted codons for each protein
    :param init_RA:
    :param MAX_RECORDS:
    :return: new top CAI index, selected_RA_df, max(iteration_list)
    '''

    # computation is converged
    converged = False
    # maximum number of outer iterations
    max_iter = 20
    # iteration number
    it = 0
    iteration_list = []

    selected_RA_df = init_RA.copy()

    top_CAI_index_old = [0]*N
    while not converged and it < max_iter:
        it += 1
        iteration_list.append(it)
        # counted_selected_proteins = count_codons(selected_proteins)

        if it > 1:
            # since 2nd iteration we have recalculate RA
            #1 calculate RA for selected proteins
            selected_RA_df = calculate_relative_adaptiveness(counted_codons)

        # 2 Update database with RA values,
        df_wtihRA = replaceCodonWithRa(counted_codons,selected_RA_df)

        # 3 calculate CAI  in the updated database
        df_sorted_cai = calculate_desc_CAI(df_wtihRA[:MAX_RECORDS])

        # 4 select N best proteins
        selected_proteins, top_CAI_index_new = selNbest(df_bacterial_genome,df_sorted_cai,N)


        if top_CAI_index_old == top_CAI_index_new:
            converged = True
            # top_CAI_index_old.to_csv('/home/anna/PycharmProjects/Protein_domains/domains/TOP_25/top_CAI_index_' +organism_name+ '.csv', sep=',')
            print("CONVERGED!!!")
        else:
            top_CAI_index_old =top_CAI_index_new


    return top_CAI_index_new, selected_RA_df, max(iteration_list)


def replaceCodonWithRa(df_all_genes, df_calRa):
    '''
    Replace each codon with with calculated RA from top 25
    :param df_all_genes: data frame with all proteins genes
    :param df_calRa: data frame with calculated relative adativeness
    :return: new data frame with replaced codons
    '''
    df_out = df_all_genes.copy()
    for rowIndex, row in df_all_genes.iterrows():  # iterate over rows
        for columnName, elem in row.items():

            if elem>0:
                df_out.at[rowIndex,columnName] = df_calRa[columnName]
    return df_out


'''
import definitions from ribosomal_protein_genes
'''
import ribosomal_protein_genes

'''
select files to analyze: 
1. including initial ribosomal proteins genes
2. including all proteins genes 
'''
df_ribosomal_protein_genes = pd.read_csv(r'ribosomal_protein_genes.csv')
df_bacterial_genome = pd.read_csv(r'genes.csv')
count_codons_df = count_codons(df_bacterial_genome)

'''
It calculates number of codons in sequence.
'''
ribosomal_protein_count_codons_df, ribosomal_protein_RA = ribosomal_protein_genes.final_ribosomal_protein_genes(df_ribosomal_protein_genes)

'''
MAIN LOOP
'''
index_list, final_weight_table_df,iteration = final_protein_index_list(count_codons_df,ribosomal_protein_RA,MAX_RECORDS)

'''
save final table to file
'''
final_weight_table_df.to_csv(r'Final_weight_table.csv', index=False)

print('Program finished')

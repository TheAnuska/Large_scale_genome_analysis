from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import math

synonymous_codons = {
    'A': ['gct', 'gcc', 'gca', 'gcg'],
    'C': ['tgt', 'tgc'],
    'D': ['gat', 'gac'],
    'E': ['gaa', 'gag'],
    'F': ['ttc', 'ttt'],
    'G': ['gga', 'ggc', 'ggg', 'ggt'],
    'H': ['cac', 'cat'],
    'I': ['ata', 'atc', 'att'],
    'K': ['aaa', 'aag'],
    'L': ['tta', 'ttg', 'cta', 'ctc', 'ctg', 'ctt'],
    'N': ['aac', 'aat'],
    'P': ['cca', 'ccc', 'ccg', 'cct'],
    'Q': ['caa', 'cag'],
    'R': ['aga', 'agg', 'cga', 'cgc', 'cgg', 'cgt'],
    'S': ['agc', 'agt', 'tca', 'tcc', 'tcg', 'tct'],
    'T': ['aca', 'acc', 'acg', 'act'],
    'V': ['gta', 'gtc', 'gtg', 'gtt'],
    'Y': ['tac', 'tat'],
}

#create list of codons on the base of dictionary, keeps both lists coupled.
list_of_codons = []

for key ,value  in synonymous_codons.iteritems():
    for codon_name in value:
        list_of_codons.append(codon_name)

amino_acid_list = synonymous_codons.keys()


# number of best proteins selected from the database i every iteration
N = 25

MAX_RECORDS = 100000
print "max records",MAX_RECORDS

# -------------------------------------------------------------------
# select nucleotide seq without duplicates
# -------------------------------------------------------------------

def nucleotide_seq (data_frame):

    list_of_columns = ['protein', 'ProteinSeq', 'CodonSeq']

    new_data_frame = data_frame[list_of_columns].copy()
    new_data_frame = new_data_frame.drop_duplicates().reset_index(drop=True)

    return (new_data_frame)

# -------------------------------------------------------------------
# Separate nucleotide seq on codons, count codons
# -------------------------------------------------------------------

def count_codons(data_frame):
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
    final_data_frame = pd.concat([codons_data_frame, new_data_frame]).fillna(0)



    return (final_data_frame)

# -------------------------------------------------------------------
# Count CAI for each amino acid
# -------------------------------------------------------------------

def calculate_desc_CAI(df_counted, df_RA_values):

    CAI_list = []
    for i in range(len(df_counted.index)):
        CAI_value=1
        n = df_counted.loc[i].sum()
        for codon in list_of_codons:
            power = df_counted.ix[i, codon]
            root = df_RA_values[codon]

            CAI_value *= math.pow(root ** power, 1.0 / n)

        CAI_list.append(CAI_value)

    new_data_frame = pd.DataFrame({'CAI' : CAI_list})


    sort_data_frame = new_data_frame.sort_values('CAI', ascending=False)

    return sort_data_frame

# -------------------------------------------------------------------
# Count relative adaptivenes for synonymous codons
# -------------------------------------------------------------------


def calculate_relative_adaptiveness(counted_codons_df):
    relative_adaptivenes_df = pd.DataFrame()

    sum_counted_codons_list = counted_codons_df.sum()
    sum_df = pd.DataFrame(sum_counted_codons_list).T
    sum_df = sum_df.replace(to_replace=0, value=1) 

    for aa_name in amino_acid_list:
        aa_codon_list = synonymous_codons[aa_name]


        df_temp = sum_df[aa_codon_list]
        max_row = df_temp.max(axis=1)

        temp_div_df = df_temp.div(float(max_row))

        relative_adaptivenes_df = relative_adaptivenes_df.append(temp_div_df)

    relative_adaptivenes_df = pd.DataFrame(relative_adaptivenes_df.sum()).T

    return relative_adaptivenes_df

# -------------------------------------------------------------------
# selects N proteins from the df_orgwih the heigest CAI (taken from df_CAI)
# ensures that all codon are in this group,if not automaticaly increases N
# -------------------------------------------------------------------

def selNbest(df_org,sorted_CAI_df,N):

    top_CAI_index = (sorted_CAI_df.index.values.tolist())
    top_CAI_index = top_CAI_index[0:N]

    selected_proteins = df_org.iloc[top_CAI_index, :]

    return selected_proteins,top_CAI_index


# -------------------------------------------------------------------
# It is for search final proteins with optimal codons (top_CAI_index_new), and final weight table (selected_RA_df)
# -------------------------------------------------------------------
def final_protein_index_list(selected_priteins,top_CAI_index_old,MAX_RECORDS):

    # computation is converged
    converged = False
    # maximum number of outer iterations
    max_iter = 20
    # iteration number
    it = 0
    iteration_list = []

    while not converged and it < max_iter:
        it+=1
        iteration_list.append(it)
        counted_selected_proteins = count_codons(selected_priteins)

        #1 calculate RA for selected proteins
        selected_RA_df = calculate_relative_adaptiveness(counted_selected_proteins)

        #2 Update database with RA values,   #3 calculate CAI  in the updated database
        sorted_cai_df = calculate_desc_CAI(count_codons_df[:MAX_RECORDS], selected_RA_df)

        #4 select N best proteins
        selected_priteins, top_CAI_index_new = selNbest(nucleotide_seq_df,sorted_cai_df,N)


        if top_CAI_index_old == top_CAI_index_new:
            converged = True
            # top_CAI_index_old.to_csv('/home/anna/PycharmProjects/Protein_domains/domains/TOP_25/top_CAI_index_' +organism_name+ '.csv', sep=',')
            print("CONVERGED!!!")
        else:
            top_CAI_index_old =top_CAI_index_new


    return top_CAI_index_new, selected_RA_df, max(iteration_list)




# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

def organism_list(organism_file_path, row_number):

    organism_list = pd.read_csv(organism_file_path, sep=',')
    organism_id   = organism_list.ix[row_number, 3]
    organism_name = organism_list.ix[row_number, 2]
    organism_path = organism_list.ix[row_number, 1]

    return organism_id,organism_name,organism_path
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# -------------------------------------------------------------------
# Copy link from Two scripts
# -------------------------------------------------------------------

import bacteial_genomes
import ribosomal_protein_genes


begin = 3412
end = 3413


print range(begin,end)

#FINAL RESULT
all_final_weight_tabel = pd.read_csv('domains/final_weight_table.csv', sep=',')
iteration_df = pd.read_csv('domains/iteration_df.csv', sep=',')


for organism_number in range(begin,end):

    organism_id, organism_name, organism_path = organism_list('domains/organism_list.csv', organism_number)
    print('name: ', organism_id,organism_name)

    marge_genomesDf = pd.read_csv('domains/merge/' + organism_id + '.csv', sep=',')


    ribosomal_protein_RA = ribosomal_protein_genes.final_ribosomal_protein_genes(organism_id, organism_path,'http://10.117.11.77:7200/repositories/Prokaryotes')

    nucleotide_seq_df = nucleotide_seq(marge_genomesDf)

    # It calculates number of codons in sequence.
    count_codons_df = count_codons(nucleotide_seq_df)

    # #LOOP INITIALIZATION
    CAI_ribosomal_protein = calculate_desc_CAI(count_codons_df,ribosomal_protein_RA)
    selected_priteins, top_CAI_index_old = selNbest(nucleotide_seq_df,CAI_ribosomal_protein,N)
    index_list, final_weight_table_df,iteration = final_protein_index_list(selected_priteins,top_CAI_index_old,MAX_RECORDS)


    # add new row with finall weight for genome (organism)
    final_weight_table_df = final_weight_table_df[list_of_codons]
    print('final_weight_table_df: ', list(final_weight_table_df.loc[0,:]))
    all_final_weight_tabel.loc[begin,:] = list(final_weight_table_df.loc[0,:])

    # Count iteration
    print('max(iteration_list): ', iteration)
    max_iteration = pd.DataFrame([[begin, iteration]])

    all_25_best_protein = nucleotide_seq_df.iloc[index_list]

    begin+=1
    print 'organism number: ', (organism_number)

# write final weight table to output file
with open('domains/final_weight_table.csv', 'a') as file:
    (all_final_weight_tabel.to_csv(file, header=False, index=False))



print 'Program finished'








# Dziala :D
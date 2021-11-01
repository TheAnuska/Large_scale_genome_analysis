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
# -------------------------------------------------------------------
#create list of codons on the base of dictionary, keeps both lists coupled.
# -------------------------------------------------------------------

list_of_codons = []
for key, value in synonymous_codons.items():
    for codon_name in value:
        list_of_codons.append(codon_name)

amino_acid_list = synonymous_codons.keys()



def count_codons(data_frame):
    '''
    Separate nucleotide seq on codons, count codons, make sure that all of codons are used
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
    temp_0 = Counter(nucleotide_list[0]).items()
    temp_df = pd.DataFrame(temp_0).set_index(0).T
    new_data_frame = temp_df

    for m in range(1, len(data_frame.index)):
        temp_m = Counter(nucleotide_list[m]).items()
        temp_m = pd.DataFrame(temp_m).set_index(0).T
        new_data_frame = new_data_frame.append(temp_m, ignore_index=True)

    new_data_frame = pd.DataFrame(new_data_frame.sum()).T
    codons_data_frame = pd.DataFrame(index=list_of_codons).T
    final_data_frame = pd.concat([codons_data_frame,new_data_frame], sort=True).fillna(0)

    return final_data_frame


def calculate_relative_adaptiveness(counted_codons_df):
    '''
    Count relative adaptiveness for synonymous codons
    :param counted_codons_df
    :return: data frame with calculated relative adaptiveness
    :param counted_codons_df:
    :return:
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



def final_ribosomal_protein_genes(data_file):
    '''
    :param data_file: data frame with initial ribosomal proteins genes
    :return: initial conditions for the main analysis (df_count_codons, df_relative_adaptiveness)
    '''
    df_count_codons = count_codons(data_file)
    df_relative_adaptiveness = calculate_relative_adaptiveness(df_count_codons)

    return df_count_codons, df_relative_adaptiveness

print('End ribosomal protein genes file')

from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt


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
# amino_acid_list = []
for key ,value  in synonymous_codons.iteritems():
    for codon_name in value:
        list_of_codons.append(codon_name)

amino_acid_list = synonymous_codons.keys()


def get_sparql_dataframe(service, query):
    """
    Helper function to convert SPARQL results into a Pandas data frame.
    """
    sparql = SPARQLWrapper(service)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    result = sparql.query()

    processed_results = json.load(result.response)
    cols = processed_results['head']['vars']

    out = []
    for row in processed_results['results']['bindings']:
        item = []
        for c in cols:
            item.append(row.get(c, {}).get('value'))
        out.append(item)

    return pd.DataFrame(out, columns=cols)


# last_csv = 'CAI_test.csv'
last_csv = 'domains/CAI_2.csv'
# domain_csv = 'domains/domain_4.csv'
middle_domain_csv = 'domains/CAI_middle_domain_2.csv'


def getGenomes(endpoint,organism):
    query = \
       '''
PREFIX gbol: <http://gbol.life/0.1/>
SELECT DISTINCT ?protein  
WHERE {
VALUES ?genome{'''+organism+'''}
?genome a gbol:Sample .
?contig gbol:sample ?genome .
?contig gbol:feature ?gene .
?gene gbol:transcript ?transcript .
?transcript gbol:sequence ?CDS .
?transcript gbol:feature ?cds .
?cds gbol:protein ?protein .
?protein gbol:feature ?feature .
?feature gbol:provenance ?provenance .
?feature gbol:signatureDesc ?domain_description .
?protein gbol:sequence ?ProteinSeq .
?protein gbol:feature ?domain .
?protein gbol:sha384 ?sha .
?domain gbol:xref ?go .
?domain gbol:location ?loc .
?domain gbol:xref ?pfam .
?pfam gbol:accession ?pfamAccession .
?pfam gbol:db ?database1 .
?database1 gbol:id 'pfam' .

FILTER(regex(?domain_description, 'ribosomal protein', 'i'))
} 
ORDER BY ?protein
# LIMIT 25

    '''
    resultsDf = get_sparql_dataframe(endpoint, query)
    return (resultsDf)


# -------------------------------------------------------------------
# select nucleotide seq without duplicates
# -------------------------------------------------------------------

def nucleotide_seq (data_frame):

    list_of_columns = ['protein', 'ProteinSeq', 'CodonSeq']

    new_data_frame = data_frame[list_of_columns].copy()
    new_data_frame = new_data_frame.drop_duplicates().reset_index(drop=True)

    return (new_data_frame)

# -------------------------------------------------------------------
# Separate nucleotide seq on codons, count codons, make sure that all of codons are used
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

    # DF with final codons marge with empty DF (with all codons as a column names)
    new_data_frame = pd.DataFrame(new_data_frame.sum()).T
    codons_data_frame = pd.DataFrame(index=list_of_codons).T
    final_data_frame = pd.concat([codons_data_frame,new_data_frame]).fillna(0)
    # final_data_frame = final_data_frame[list_of_codons].replace(to_replace=0., value=1.) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return (final_data_frame)


def calculate_relative_adaptiveness(counted_codons_df):
    relative_adaptivenes_df = pd.DataFrame()

    sum_counted_codons_list = counted_codons_df.sum()
    sum_df = pd.DataFrame(sum_counted_codons_list).T
    # print sum_df
    sum_df = sum_df.replace(to_replace=0, value=1) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for aa_name in amino_acid_list:
        aa_codon_list = synonymous_codons[aa_name]


        df_temp = sum_df[aa_codon_list]
        max_row = df_temp.max(axis=1)

        temp_div_df = df_temp.div(float(max_row))

        relative_adaptivenes_df = relative_adaptivenes_df.append(temp_div_df)

    relative_adaptivenes_df = pd.DataFrame(relative_adaptivenes_df.sum()).T

    return relative_adaptivenes_df

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

def organism_list(organism_file_path, row_number):

    organism_list = pd.read_csv(organism_file_path, sep=',')
    organism_id   = organism_list.ix[row_number, 3]
    organism_name = organism_list.ix[row_number, 2]
    organism_path = organism_list.ix[row_number, 1]


    return organism_id,organism_name,organism_path

# ------------------------------------------------------------------------------------
#  Final def ribosomal protein genes;
# return data frame with list of ribosomal genes, which initiates computing weight table
# ------------------------------------------------------------------------------------


def final_ribosomal_protein_genes(organism_id,organism_path, endpoint):


    ribosomal_protein_genes = getGenomes(endpoint, '<' + organism_path + '>')
    genomesList = ribosomal_protein_genes["protein"]

    genomesDf = pd.read_csv('domains/merge/' + organism_id + '.csv', sep=',')
    genomesDf = genomesDf[['protein', 'ProteinSeq', 'CodonSeq']]

    ribosomal_protein_df = pd.merge(ribosomal_protein_genes, genomesDf, how='left', on=['protein'])
    ribosomal_protein_df = ribosomal_protein_df.dropna()
    ribosomal_protein_df_nucleotide = nucleotide_seq(ribosomal_protein_df)
    ribosomal_protein_df_nucleotide.to_csv('domains/RA/RA_' + organism_id + '.csv',sep=',', index=False)

    ribosomal_protein_count_codons = count_codons(ribosomal_protein_df_nucleotide)
    ribosomal_protein_RA = calculate_relative_adaptiveness(ribosomal_protein_count_codons)
    ribosomal_protein_RA = ribosomal_protein_RA[list_of_codons]

    return ribosomal_protein_RA


print('end ribosomal_protein_genes')

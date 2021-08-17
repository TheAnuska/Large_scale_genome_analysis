from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt


# from tabulate import tabulate


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


def getGenomes(endpoint,organism):
    query = \
        '''
PREFIX gbol: <http://gbol.life/0.1/>
    SELECT DISTINCT ?protein ?ProteinSeq ?CodonSeq ?beginpos ?endpos ?pfamAccession
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
?loc gbol:begin ?begin .
?loc gbol:end ?end .
?begin gbol:position ?beginpos.
?end gbol:position ?endpos .
?domain gbol:xref ?pfam .
?pfam gbol:accession ?pfamAccession .
?pfam gbol:db ?database1 .
?database1 gbol:id 'pfam' .
?mrna gbol:feature ?cds ;
      gbol:sequence ?CodonSeq .
        }    

    ORDER BY ?organism ?protein ?pfamAccession ?beginpos

    '''
    resultsDf = get_sparql_dataframe(endpoint, query)
    return (resultsDf)


# ------------------------------------------------------------------------------------
# 1. Select the rows with the same length of protein ('ProteinSeq') and nucleotide seq ('CodonSeq')
# ------------------------------------------------------------------------------------

def the_same_length(data_frame):
    data_frame['ProtSeq_len'] = data_frame['ProteinSeq'].apply(len)
    data_frame['CodonSeq_len'] = data_frame['CodonSeq'].apply(len)
    data_frame['D'] = data_frame['CodonSeq_len'] / 3 - 1
    data_frame = data_frame.loc[data_frame['D'] == data_frame['ProtSeq_len']]

    # delete columns needed for calculation
    new_data_frame = data_frame.drop(['ProtSeq_len', 'CodonSeq_len', 'D'], axis=1)

    return (new_data_frame)


# ------------------------------------------------------------------------------------
# 2. MERGE proteins with tis same pf and sequence
# Join domains with this same pfamAccession and sequence of nucleotides and amino acids
# ------------------------------------------------------------------------------------

# merge lists overlaping - def helping in marge ovelaping domains
def merge_lists(l):
    s = map(set, l)
    i, n = 0, len(s)
    while i < n - 1:
        for j in xrange(i + 1, n):
            if s[i].intersection(s[j]):
                s[i].update(s[j])
                del s[j]
                n -= 1
                break
        else:
            i += 1
    return [sorted(i) for i in s]


def marge_domains(data_frame):

    # 1. Marge rows for this same 'pfamAccession'
    data_frame.loc[:, 'beginpos'] = pd.to_numeric(data_frame['beginpos'], errors='coerce')
    data_frame.loc[:, 'endpos'] = pd.to_numeric(data_frame['endpos'], errors='coerce')

    temp_data_frame = data_frame.groupby(
        ['protein', 'ProteinSeq', 'CodonSeq', 'pfamAccession']).agg(
        {'beginpos': np.min, 'endpos': np.max}).reset_index()

    # 2. Marge rows for ovelaping domains
    col_to_remove = [r'Unnamed: 0', r'pfamAccession', r'Codon_domain']

    for col_name in col_to_remove:
        if col_name in list(temp_data_frame.columns):
            temp_data_frame.drop(col_name, axis=1, inplace=True)

    grouped_data_frame = temp_data_frame.groupby(['protein', 'ProteinSeq', 'CodonSeq'])

    new_data_frame = pd.DataFrame(columns=list(temp_data_frame.columns.values))

    added_groups = 0
    for name, group in grouped_data_frame:
        sorted = group.sort_values('beginpos')

        protein_loc = sorted.columns.get_loc('protein')
        ProteinSeq_loc = sorted.columns.get_loc('ProteinSeq')
        CodonSeq_loc = sorted.columns.get_loc('CodonSeq')
        beginpos_loc = sorted.columns.get_loc('beginpos')
        endpos_loc = sorted.columns.get_loc('endpos')

        protein_1 = sorted.iloc[0, protein_loc]
        ProteinSeq_1 = sorted.iloc[0, ProteinSeq_loc]
        CodonSeq_1 = sorted.iloc[0, CodonSeq_loc]

        # merge domain which are intersecting,(beginpo,endpos)
        domain_parts = []
        for row_idx in range(len(sorted)):
            beginpos_1 = sorted.iloc[row_idx, beginpos_loc]
            endpos_1 = sorted.iloc[row_idx, endpos_loc]

            domain_parts.append(range(beginpos_1, endpos_1 + 1))

        merged_domain_parts = merge_lists(domain_parts)

        for part in merged_domain_parts:
            values = [protein_1, ProteinSeq_1, CodonSeq_1, part[0], part[-1]]
            new_data_frame.loc[added_groups] = values
            added_groups += 1

    return (temp_data_frame)

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


# ------------------------------------------------------------------------------------
# 3. Final results from this script give a data frame with selected proteins for a genome
# ------------------------------------------------------------------------------------

def final_result_organism(organism_id,organism_path,endpoint):

    # organism_id,organism_name = organism_list(organism__list_file_path)

    genomesDf = getGenomes(endpoint, '<' + organism_path + '>')
    genomesList = genomesDf["protein"]

    genomesDf['CodonSeq'] = genomesDf['CodonSeq'].str.lower()
    genomesDf_temp = genomesDf.drop_duplicates().reset_index(drop=True)

    genomesDf_the_same_length = the_same_length(genomesDf_temp)

    marge_genomesDf = marge_domains(genomesDf_the_same_length)
    marge_genomesDf = marge_genomesDf.sort_values(['protein', 'ProteinSeq', 'CodonSeq', 'beginpos'], inplace=False)
    marge_genomesDf.to_csv('domains/merge/merge_' + organism_id + '.csv', sep=',', index=False)

    return marge_genomesDf



print('end_Organism.py')


































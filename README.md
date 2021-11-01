Analysis was based on the conceptual framework proposed by Puigbo et al. (2007). The OPTIMIZER framework is a non-line application where codon usage was used to optimize a DNA sequence to predict a group of highly expressed genes. The modified algorithm starts from randomly chosen bacterial genes to create the initial weight table. It then goes through a series of iterations to find a group of highly expressed genes and a table with relative adaptiveness values.

Before using the algorithm, each genome was reduced by the codons whose role is not significant regarding the codon bias (start or stop codons) or codon is related with one amino acid.

The algorithm describes how the final relative adaptiveness table, the collection of Relative Adaptiveness (RA) of the most highly expressed genes from each bacterial genome, was performed.  The algorithm for a  chosen genome was initiated by the group of selected genes corresponding to this genome and the RA calculation.  After calculation,  each codon in the whole genome was replaced with received RA. The next step was to calculate Codon Adaptation  Index  (CAI)  of each gene and select 25 genes with the highest CAI value. The following step was to compare the given 25 genes with the initiating genes. If they were different, then the new genes become initiating genes, and the process was repeated until the new group of genes were equal to genes in the previous iteration. Once the final group of genes was established, RA's value was saved into the final relative adaptiveness table.

[1] Pere  Puigb`o,  Eduard  Guzm ́an,  Antoni  Romeu,  and  Santiago  Garcia-Vallv ́e.   OPTIMIZER: A web server for optimizing the codon usage ofDNA sequences.Nucleic Acids Research, 35(SUPPL.2), jul 2007.


ribosomal_protein_genes.csv - initial file with selected ribosomal proteins genes

genes.csv - file with all proteins genes

Final_weight_table.csv - file with the results


1. ribosomal_protein_genes.py - script to analyze initial file
2. relative_adaptiveness.py - main script


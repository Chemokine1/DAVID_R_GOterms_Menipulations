# DAVID_R_GOterms_Menipulations
This repository provides different pipelines, written in the R programming language, for visualizing Gene Ontology (GO) terms associated with sets of genes analyzed by DAVID.

The repository contains the following scripts:

"GOterms_DAVID_Chart": This script returns a columns chart that displays the significant GO terms (Benjamini < 0.05), their fold enrichment, and the number of genes found for each term in the analyzed gene set.

"GOterms_DAVID_Shared2Files_Chart": This script returns a columns chart that shows the shared significant GO terms (Benjamini < 0.05) and their fold enrichment for two different files. It also displays the number of genes associated with each term.

"GOterms_In_File1_NOT_In_FIle2_DAVID_Chart": This script returns a columns chart that displays the significant GO terms (Benjamini < 0.05), their fold enrichment, and the number of genes found for each term in file1 but not in file2. (Only the GO terms that are unique for gene set number 1)

"GOterms_Shared2Files_afterAntiJoin2Files": This script combines the pipelines from scripts 3 and 2 and requires four files. It returns a columns chart that displays the significant GO terms (Benjamini < 0.05) and their fold enrichment for the shared GO terms between two files after performing an antijoin operation.

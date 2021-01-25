DEPARTMENT OF BIOINFORMATICS 
THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS

All the R scripts were developed in R studio R version 3.6.1 (2019-07-05).

The scripts were used for the analysis and visualisation of results in 

Tomazou, M., Bourdakou, M., Minadakis, G., Zachariou, M., Oulas, A., Karatzas, E., ... & Spyrou, G. (2020). Multi-omics data integration and network-based analysis drives a multiplex drug repurposing approach to a shortlist of candidate drugs against COVID-19.
Corresponding Author George Spyrou email: georges@cing.ac.cy

Scripts description:

integration_github_part1.R
By Margarita Zachariou
email: margaritaz@cing.ac.cy
The script takes as an input lists of genes from multiple sources which include gene symbols, 
entrez ids, original scores and ranks (if available) and normalised scores (weights)
All gene lists have been manually curated to annotate which genes can not be found in the Genemania database
	 *  ranked genes lists in terms of absolute logFC from three serum transcriptomic (T) datasets (Series 15, BALF, PBMC)[6],
	 *  ranked genes lists in terms of absolute log fold from one proteomics (P) dataset[5]
	 *  ranked genes lists in terms of p-value from one metabolomics (M) dataset[5]
	 *  unranked list of host proteins (PPI) which interact with SARS-CoV-2 from Gordon et al.[8]
	 *  unranked unique gene list from HPA, excluding the genes identified in Gordon et al.
The output of this script is used as an input to Cytoscape to generate the genemania network 
and as an input to integration_part2 to complete the integration process.

integration_github_part2.R
By Margarita Zachariou
email: margaritaz@cing.ac.cy 
Takes as input the integration based on the node score table (df_int_filtered) and the 
Genemania network files
The edge score is callculated and combined to the node score to obtain the overall integration score.
The output of this script is used as an input to IntegrationNet.R to plot the integration network 
and as an input to CODRES and Pathwalks tools.

IntegrationNet.R:
By Marios Tomazou
email: mariost@cing.ac.cy
This script takes as an input the resulting files from the integration step comprising lists of ranked genes, their MIGe, MIGn and MIG scores and creates a network with binned node sizes and coloured based on their originating dataset. In addition and calculates the degree, Betweeness and Closeness and plots their distributions.

DrugsSankey_Human.R:
By Marios Tomazou
email: mariost@cing.ac.cy
This script takes as an input the drug list indicating the shortlisted drugs, the mapping of KEGG pathways with their IDs, the mapping of genes identified as drug targets with their KEGG pathways, a table of drugs and their originating list and the output of the PWalks2net script.           The output is a sankey plot visualising the mapping of lists to shortlisted drugs to genes to KEGG pathways.    

DrugsSankey_Pathogen.R:
By Marios Tomazou
email: mariost@cing.ac.cy
This script takes as an input the drug list indicating the shortlisted drugs, a table of drugs indicating their originating list and their viral target proteins. The output is a sankey plot visualising the mapping of lists to shortlisted drugs to genes to KEGG pathways.     

mainAnalysis_HPAV.R:
By George Minadakis
email: georgem@cing.ac.cy
This script takes as an input a drug information file from DrugBank, a list of pathogens and their NCBI taxonomy IDs, the taxonomy distance matrix of these pathogens, a network of host pathogen interactions file and a list of host proteins interacting with SARS-CoV-2. The output is the HPAV list which include drugs targeting SARS-COV-2 proteins scored based on the commonality of other pathogens with SARS-COV-2 in terms of their taxonomy distance and host-pathogen interactions.


mainAnalysis_HPH.R:
By George Minadakis 
email: georgem@cing.ac.cy
This script takes as an input a drug information file from DrugBank, a list of pathogens and their NCBI taxonomy IDs, the taxonomy distance matrix of these pathogens, a network of host pathogen interactions file and a list of host proteins interacting with SARS-CoV-2. The output is the HPH list which include drugs targeting the host proteins commonly affected by SARS-CoV-2 and other pathogens

taxonomyDrugR.R:
By Marios Tomazou
email: mariost@cing.ac.cy
This script takes as an input a list of pathogens and the targets, polypeptides and drugs files from Drugbank. The output is a data frame of anti-pathogen drugs scored with respect to the taxonomic distance of the target pathogens to the pathogen of interest (herein SARS-CoV-2).

taxaplot.R:
By Marios Tomazou
email: mariost@cing.ac.cy
The script takes as an input a list of pathogens with their NCBI taxonomy IDs and the DrugBank drugs targeting their proteins. Loose and exact match refer to drugs that target the species in which the taxonomy ID belogns or exactly the taxonomy ID, respectively. In addition, a second with the pathogens and their taxonomy IDs and the number of the corresponding human-pathogen protein interactions is also required. The output is a circular dendrogram based on the taxonomy distance highlighting pathogens of interest and showing as a heatmaps the number of pathogen-human protein interactions and number of drugs targeting these pathogens.

PWalks2net.R:
By Marios Tomazou
email: mariost@cing.ac.cy
This script takes as an input the result of the Pathwalks analysis, a fully random walk results, and optionally a layout matrix. Then it performs an Odds ratio analysis to identify the most enriched pathways which are highlighted in the resulting pathway community network and     
bar chart. Finally the script saves the results in table format.

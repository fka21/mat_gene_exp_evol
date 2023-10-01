# mat_gene_exp_evol
Scripts used for the "Comparative Analysis of Maternal Gene Expression Patterns: Unraveling Evolutionary Signatures Across Reproductive Modes" publication.

Scripts used for de novo assemblies can be found in the Assemblies directory. These scripts contain the commands used for single-end read data and paired-end data, also a script used iteratively for each species for running EviGene. Finally an R script can be found for generating a BUSCO analysis figure.

A separate R script (functions.R) is included, which contains custom functions used throughout the analysis. 

The previous directory can be used separately with the appropriate data. 

The rest of the scripts rely for an initial step of differential gene expression (DGE) analysis. The scripts used for this can be found in the Quantification directory. Within log subdirectory the salmon mapping rates are found. Within scripts subdirectory the scripts used for DGE are found. Command for pre-processing of the raw reads (fastp.sh), the indexing of transcriptomes (read_lengths.sh, salmon_indexing.sh) and the alignmnet (quantification.sh) are found here.  The salmon quantification files are imported with tximport.R script. As it would require alot of memory only tximport objects are saved in the enviroment. Maternally expressed genes are defined here aswell for each gene. After running this script, DGE itself is done be all_species_deseq.R script. Within the script the Drosophila melanogaster maternally expressed genes based from RNAseq is compared to in situ hybridization data from FLY-FISH database.  A misc script is also included for checking the maternally expressed gene categorizations. 

To amalgamate the gene expression data in a single gene expression matrix containing the scripts in Orthology_inference directory were used. OrthoFinder run can be replicated with orthofinder.sh script (prerequisite to this is the longest.py ancilliary script provided by the OrthoFinder developers to extract the longest proteins from the provided sequences). The species tree calibration was done using species_tree_calibration.R script. To combine the gene expression data with orthology data the OG_tbl_generator.R was used.

The gene ontological (GO) analysis was performed using scripts found in the GO_analysis directory. To retrieve GO annotation the uniprot_downloader.R script was used. It requires the annotated orthogroups as input. The GO analysis itself is performed by the GO.R script. Intermediate data are also attached used by GO.R 

Gene feature length comparisons were analysed using the script found in Feature_length_analysis directory. 

Modeling was done in the Phylogenetic_logistic_regression and Model_fittings directories. Initially the Model_fittings directory should be run. Here the phylosignal.R script is used to assess phylogenetic signal present in the dataset.  Two separate subdirectories (multiregime and multiregime_fc) are used to model gene expression data and fold change data separately. To run model fitting use the model_fitting.R script in each subdirectory. To build the model parameter matrix use the parameter_estimator.R scripts. Parameters are analysed using the model_param_analysis.R script.

Intermediate datasets can be accessed on Zenodo ( DOI: 10.5281/zenodo.8374018 ).

Note: please adjust the paths towards all data/scripts according to your own setup. 

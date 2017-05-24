################################################################################
################################################################################
#
# How do the type of QTL effect and the form of the residual term influence QTL
# detection in multi-parent population? A case study in the maize EU-NAM
# population
#
# Vincent Garin, Valentin Wimmer, Sofiane Mezmouk, Marcos Malosetti
# Fred van Eeuwijk
#
################################################################################
################################################################################


*** The repository MPP_EUNAM contains all the data (except the raw genotype array data)
and the script (~\MPP_EUNAM\scripts\Complete_analysis.R) to reproduce all results
presented in the study.


*** The original data come from: 


-phenotype (trait) data : http://www.genetics.org/content/198/1/3/suppl/DC1

-map : http://maizegdb.org/cgi-bin/displayrefrecord.cgi?id=9024747

-genotype matrix : http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50558

!!! This dataset is to big to be loaded in a Github repository. Please load it and
put it in the following folder: ~MPP_EUNAM/data/geno using the following
name: geno_array_EUNAM.csv

-> ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50558/matrix/



*** Except these data all origninal and intermediate datasets necessary
to run the full analyses are present in the folowing folders of the repository

~MPP_EUNAM/data/clustering

~MPP_EUNAM/data/geno

~MPP_EUNAM/data/map

~MPP_EUNAM/data/mpp_data

~MPP_EUNAM/data/pheno


*** to run the analysis you will need a list of package mentioned in the script.
Some of these package are not available on CRAN. You can find them at the
following locations

mppR: -> ~\MPP_EUNAM\software\mppR_1.0.tar.gz

clustaplo: -> ~\MPP_EUNAM\software\clusthaplo_1.2.tar.gz

or https://cran.r-project.org/src/contrib/Archive/clusthaplo/

! clustering results have been done using Rhmm package. To be able
to run these two packages simulataneously We had to use version
R 2.14.



*** The script ~\MPP_EUNAM\scripts\Complete_analysis.R allows to reproduce
all the results of the research. The different steps are listed bellow.
Each step can be run independently of the others.

*** The script ~\MPP_EUNAM\scripts\Table_Figure_Article.R allow to reproduce
the table and figures of the article

*** The script ~\MPP_EUNAM\scripts\Table_Figure_Supp_Mat.R allow to reproduce
the table and figures of the appendix

*** concerning the production of results (QTL analysis) (step 8 - 12)
the original results of the study are already saved in the folder
~\MPP_EUNAM\results.


1. Adjusted means computation


2. Map processing


3. Genotype matrix sorting


4. Determination of the three subsets


5. short subset data processing


6. heterogeneous subset data processing


7. Long subset data processing


8. Significance threshold determination full datasets


9. QTL analysis full datasets


10. Significance threshold determination reduced datasets


11. Cross-validation (reduced datasets)


12. Multi QTL effect model full dataset


13. Multi QTL effect model CV (reduced datasets)
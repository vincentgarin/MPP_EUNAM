################################################################################
################################################################################
#
# The effect of genetic distance between parents in MPP QTL analysis
#
# Vincent Garin, Marcos Malosetti, Fred van Eeuwijk, 2016
#
################################################################################
################################################################################

# Library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppR build for this research


library(mppR)
library(asreml)
library(outliers)


#####################
#  Raw data sources #
#####################

# the raw data used come from the following sources:

# genotype matrix : http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50558

# phenotype (trait) data : http://www.genetics.org/content/198/1/3/suppl/DC1

# map : http://maizegdb.org/cgi-bin/displayrefrecord.cgi?id=9024747

# You can find the necessary data in the folowing folders of the repository

# ~MPP_EUNAM/data/geno

# ~MPP_EUNAM/data/pheno

# ~MPP_EUNAM/data/map


##################################
# 1. Adjusted means computation  #
##################################

# 1.1 load data 
##############

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

pheno <- read.csv("./data/pheno/PhenotypicDataDent.csv")

# change the column names
col.names <- colnames(pheno)
col.names[c(2, 7)] <- c("Fam", "Block")
colnames(pheno) <- col.names
rm(col.names)

par_names <- c("B73", "D06", "D09", "EC169", "F252", "F353", "F618", "F98902",
               "Mo17", "UH250", "UH304", "W117")

# change MO17 for Mo17
pheno$Genotype <- as.character(pheno$Genotype)
pheno$Fam <- as.character(pheno$Fam)
pheno$Pedigree <- as.character(pheno$Pedigree)
pheno[which(pheno[, 1] == "MO17"), 1:3] <- "Mo17"
pheno$Genotype <- as.factor(pheno$Genotype)
pheno$Fam <- as.factor(pheno$Fam)
pheno$Pedigree <- as.factor(pheno$Pedigree)


# 1.2 Put plot information as missing for DMC and DMY where number of plant
# is bellow 0.7 the median score per location
#############################################


# compute the median number of plant per location

loc.id <- levels(pheno$LOC)

median.NBPL <- c()

for(i in 1:4){
  
  median.i <- median(pheno$NBPL[pheno$LOC == loc.id[i]])
  median.NBPL <- c(median.NBPL, median.i)
  
}

# determine per location which plot contain 70% less plant than the number of
# data

missing.plot <- c()

for(i in 1:4){
  
  missing.plot.i <- pheno$NBPL[pheno$LOC == loc.id[i]] < (median.NBPL[i] * 0.7)
  missing.plot <- c(missing.plot, missing.plot.i)
  
}


# Put missins values for DMY and DMC

pheno[missing.plot, 9:10] <- NA

# 1.3 Perform the Grubb test to detect outlying plot
###################################################

# Grubb test applied on the residual of a mixed model. If the residual of the
# mixed model contained an outlying observation (Grubb test pval< 0.05), the 
# plot value for this specific trait with the largest residual is put as missing.
# This procedure is applied iteratively until no extreme value is detected.
# The test is done for the 5 traits.

# iteration on all the five traits

trait.id <- c("DMY", "DMC", "PH", "DtTAS", "DtSILK")
pos.id <- c(9, 10, 11, 12, 13)

library(outliers)

for (i in 1:5){
  
  fix_formula <- paste(trait.id[i],"~1",sep="")
  
  # set initial value for p_val
  
  p_val <- 0.01
  count <- 0
  
  while(p_val < 0.05){
    
    # compute the mixed model
    
    model <- asreml(fixed = as.formula(fix_formula), random = ~ Genotype + LOC +
                      Genotype:LOC + LOC:Rep +  LOC:Rep:Block,
                    rcov= ~ idv(units), data = pheno, trace = FALSE)
    
    # compute the residuals
    res <- model$residuals
    
    # do the Grubbs test on the residuals
    test <- grubbs.test(res, type = 10, two.sided = TRUE)
    
    # get the p-val of the test
    p_val <- test$p.value
    
    # test if p_val is lower than 0.05 (presence of an outlier)
    if(p_val < 0.05){
      
      # remove the outlying value from the data
      
      # check if there are missing values within the residuals
      
      if(sum(is.na(res))>0){
        
        pheno[which(abs(res) == max(abs(res[-which(is.na(res))]))),
              pos.id[i]] <- NA
        
        count <- count + 1
        
        
      } else {
        
        
        pheno[which(abs(res) == max(abs(res))), pos.id[i]] <- NA
        count <- count + 1
        
      }
      
    }
    
  }
  
}


# 1.4 subset the data that were used in Lehermeier et al. research
#################################################################

# Several data have not been used for the final analysis. These lines comprise
# the cross CFD08. Other lines that were not genotyped were also excluded.
# lines with a call rate <0.9 were also removed. Finally lines were discarded
# due to problem due to the map construction. These information have been
# provided by Christina Lehermeier.

# To separate the line that will be further used and to the other, I use
# a list transmitted by Christina Lehermeier which correspond to the list
# of lines used in Lehermeier et al. (2014)

lines_used <- read.csv("./data/pheno/List_lines_Dent_Lehermeier.csv")

# assign specific family indicator for the non used lines

ind <- !(pheno[,1] %in% lines_used[,1])

pheno$Fam <- as.character(pheno$Fam)
pheno$Fam[ind] <- "Extra"
pheno$Fam <- as.factor(pheno$Fam)

# save data

write.csv(pheno,"./data/pheno/pheno_red.csv")


# 1.5 Adjusted means computation
#################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

pheno <- read.csv("./data/pheno/pheno_red.csv", row.names = 1)

lines_used <- read.csv("./data/pheno/List_lines_Dent_Lehermeier.csv")

# suppress the factor levels with 0 occurences

pheno$Genotype <- as.factor(as.character(pheno$Genotype))
pheno$Fam <- as.factor(as.character(pheno$Fam))

# order data by family

pheno <- pheno[order(pheno$Fam), ]

par_names <- c("B73", "D06", "D09", "EC169", "F252", "F353", "F618", "F98902",
               "Mo17", "UH250", "UH304", "W117")


trait.id <- c("DMY", "DMC", "PH", "DtTAS", "DtSILK")
Adj_means <- matrix(0, dim(lines_used)[1], 5)


for (i in 1:5){
  
  fix_formula <- paste(trait.id[i], " ~ 1 + Genotype", sep = "")
  model <- asreml(fixed = as.formula(fix_formula),
                  random = ~  LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                  rcov = ~ idv(units), data = pheno,
                  na.method.X = "omit", na.method.Y="omit")
  
  pred <- predict(model, classify = "Genotype")
  index <- pred$predictions$pvals$Genotype %in% lines_used[, 1]
  Adj_means[, i] <- pred$predictions$pvals$predicted.value[index]
  
} 

colnames(Adj_means) <- c("DMY","DMC","PH","DtTAS","DtSILK")
rownames(Adj_means) <- lines_used[, 1]

# save the adjusted means

write.csv(Adj_means,"./data/pheno/Adj_means.csv")

# remove all useless variables

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

################################################################################
######################## End computation adjusted means ########################
################################################################################

#####################
# 2. Map processing #
#####################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# original complete map

map_dent <- read.table("./data/map/Complete_map_Dent.txt", h = TRUE)

# information concerning 50k marker of the array

array <- read.csv("./data/geno/Plateform_Data.csv")
array <- array[, c(1, 3, 8)]

# 2.1 Remove positions in the array or in the map that have the same identifier
###############################################################################

# check if there are some markers identifier that are used twice or more

length(table(array[, 2]))
length(table(map_dent[, 2]))

# In both case, there are less identifier than positions.
# This means that some marker identifier are repeated on the map and on the array

# in order to make a correspondance between the map position and the marker
# in the data matrix we need to only have one possible identifier
# the position with two or more identifier will be removed from the map
# and from the marker matrix

# list of markers that are positioned twice or more in the map

rep_mk <- which(table(map_dent[, 2]) > 1)
ID_mk <- attr(table(map_dent[, 2])[rep_mk], "dimnames")

# remove the list of misplaced markers from the full map

map_dent_red <- map_dent[!(map_dent[, 2] %in% ID_mk[[1]]), ]


# list of markers that are positioned twice or more in the map

rep_mk <- which(table(array[, 2]) > 1)
ID_mk <- attr(table(array[, 2])[rep_mk], "dimnames")

# remove the list of misplaced markers from the full map

array_red <- array[!(array[, 2] %in% ID_mk[[1]]), ]

# 2.2 match SNP_id and rsid
###########################

inter.mk.map <- intersect(array_red[, 2], map_dent_red[, 2])

new_map_dent <- map_dent_red[map_dent_red[, 2] %in% inter.mk.map, ]
array_sel <- array_red[array_red[, 2] %in% inter.mk.map, ]

# combine the rs ID to have the same SNP id in the map and the array

# put array_sel in the same order as map data

array_sel_sort <- array_sel[match(new_map_dent$Locus, array_sel$RSID), ]

new_map_dent2 <- cbind(array_sel_sort, new_map_dent)

map <- new_map_dent2[, c(1, 2, 3, 4, 6)]

# save map with rsid and SNPID

write.table(map, "./data/map/map_Dent_SNPID.txt")

# 2.3 select only the Panzea markers like in Giraud et al. (2014)
#################################################################

mk.origin <- substring(text = map[, 1], 1, 2)
map <- map[which(mk.origin == "PZ"), ]
map <- map[, c(1, 4, 5)]

mk.names <- as.character(map[, 1])
mk.names <- gsub("-", ".", mk.names)
map[, 1] <- mk.names

# save map with all Panzea markers

write.table(map, "./data/map/map_panzea_full.txt")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list = variable)
rm(variable)

################################################################################
############################## End map processing ##############################
################################################################################


##############################
# 3. Genotype matrix sorting #
##############################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# 3.1 General sorting
#####################

# genotype matrix (! full array data, heavy file)

geno <- read.csv("./data/geno/geno_array_EUNAM.csv")

# parents of the Dent panel

parents <- c("F353", "B73", "D06", "D09", "EC169", "F252", "F618", "Mo17",
             "UH250", "UH304", "W117")

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)

# map

map <- read.table("./data/map/map_panzea_full.txt", h = TRUE)

# get the genotype identifier

geno.id <- geno[36, ]

geno.id2 <- c()

for (i in 2:2291) {
  
  name <- strsplit(as.character(geno.id[, i]), split = ",", fixed = TRUE)[[1]][2]
  
  geno.id2 <- c(geno.id2, name)
  
}

geno.id3 <- substring(geno.id2, 2, nchar(geno.id2))

# Remove the row that are not the genotypic information

geno <- geno[-(1:73), ]
colnames(geno) <- c("ID_ref", geno.id3)

### 3.1.1 select only the marker position that are present in the map

geno <- geno[geno[, 1] %in% map[, 1], ]

# put the marker in the same order as in the map

geno <- geno[match(map[, 1], geno[, 1]), ]

# reorder genotypes as row and markers as column and modify scores

geno <- t(geno)
geno <- as.matrix(geno)
colnames(geno) <- geno[1, ]
geno <- geno[-1, ]

geno[geno == "AB"] <- NA
geno[geno == "NC"] <- NA

# save the genotype matrix of the parents for Clusthaplo

geno.pz.par <- geno[rownames(geno) %in% parents, ]

write.csv(geno.pz.par, "./data/geno/geno_panzea_par.csv")

### 3.1.2 select only the genotypes that have been phenotyped

# select the genotype that have been phenotyped and the parents

entries.sel <- c(parents, rownames(pheno))

geno <- geno[rownames(geno) %in% entries.sel, ]

write.csv(geno, "./data/geno/geno_panzea_full.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list = variable)
rm(variable)

################################################################################
######################## End genotype matrix sorting ###########################
################################################################################

#############################################
# 4. Determine the three subset populations #
#############################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)
geno.par <- as.matrix(geno.par)


# SM coefficient distance between the central parent (F353)

kin.mat.par <- SM_comp(mk.mat = geno.par)

dist.cent.par <- kin.mat.par[, 6]
dist.cent.par <- dist.cent.par[-6]

# relate information distance between the two parents to the cross
# information

crosses <- c("CFD02","CFD03","CFD04","CFD05","CFD06","CFD07",
             "CFD09","CFD10","CFD11","CFD12")

parent2 <- c("B73","D06","D09","EC169","F252","F618","Mo17",
             "UH250","UH304","W117")

dist.cent.par <- data.frame(crosses, parent2, dist.cent.par,
                            stringsAsFactors = FALSE)

colnames(dist.cent.par) <- c("cross","parent2","SM")

# order the crosses according to genetic distance between parents

dist.ord <- dist.cent.par[order(dist.cent.par$SM, decreasing = TRUE), ]
dist.ord

#  form the subsets

short.cross <- dist.ord$cross[1:5]
long.cross <- dist.ord$cross[6:10]
hetero.cross <- dist.ord$cross[c(1, 3, 5, 8, 10)]

par.short <- dist.ord$parent2[1:5]
par.long <- dist.ord$parent2[6:10]
par.hetero <- dist.ord$parent2[c(1, 3, 5, 8, 10)]

cross.partition <- data.frame(short.cross, long.cross, hetero.cross,par.short,
                              par.long,par.hetero)

# save

write.table(cross.partition,"./data/geno/cross.partition.txt")

# removed unused variable

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

################################################################################
######################### End subset determination #############################
################################################################################


###################################
# 5. short subset data processing #
###################################

library(mppR)

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv", row.names = 1)
geno <- as.matrix(geno)

# select the subset of 361 genotypes

geno_used <- read.table("./data/geno/geno_sel_list.txt")
parent_used <- read.table("./data/geno/cross.partition.txt",
                          stringsAsFactors = FALSE)

geno.off <- geno[rownames(geno) %in% geno_used[, "short"], ]
geno.par <- geno[rownames(geno) %in% c("F353", parent_used[, "par.short"]), ]

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)

trait <- pheno[rownames(pheno) %in% geno_used[, "short"], 1, drop = FALSE]
trait <- data.frame(rownames(trait), trait[, 1], stringsAsFactors = FALSE)

cross.ind <- substr(trait[, 1], 1, 5)

# map

map <- read.table("./data/map/map_panzea_full.txt", h = TRUE,
                  stringsAsFactors = FALSE)

# par.per.cross

par.per.cross <- cbind(unique(cross.ind), rep("F353", 5),
                       parent_used[, "par.short"])

colnames(par.per.cross) <- c("cr", "par1", "par2")


# 5.1 IBD mppData object
########################

data_IBD <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
                    trait = trait, cross.ind = cross.ind,
                    par.per.cross = par.per.cross, MAF.pop.lim = 0.01,
                    MAF.cr.lim =  c(0, 0, 0, 0, 0), MAF.cr.miss = TRUE,
                    ABH = TRUE)

# save data for parental clustering later

write.table(x = data_IBD$map.par.clu, file = "./data/map/hap_map_sh.txt")
write.table(x = data_IBD$geno.par.clu, file = "./data/geno/hap_geno_sh.txt")


data_short_IBD <- mppData_form(geno.off = data_IBD$geno.off,
                               geno.par = data_IBD$geno.par, IBS = FALSE,
                               type = "dh", map = data_IBD$map,
                               trait = data_IBD$trait,
                               cross.ind = data_IBD$cross.ind,
                               par.per.cross = data_IBD$par.per.cross,
                               step = 50, dir = "./data/geno")

saveRDS(data_short_IBD, file = "./data/mpp_data/data_short_IBD.rds")
# previously data.short.ABH

# 5.2 IBS mppData object
########################

data_IBS <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
                    trait = trait, cross.ind = cross.ind,
                    par.per.cross = par.per.cross, MAF.pop.lim = 0.01,
                    MAF.cr.lim =  c(0, 0, 0, 0, 0), MAF.cr.miss = TRUE,
                    ABH = FALSE)


# Imputation with Beagle

library(synbreed)

# form a gpData object

family <- data.frame(data_IBS$cross.ind)
rownames(family) <- rownames(data_IBS$geno.off)

map2 <- data_IBS$map
map2[, 1] <- as.character(map2[, 1])
rownames(map2) <- map2[, 1]
map2 <- map2[, 2:3]
colnames(map2) <- c("chr", "pos")

trait <- data_IBS$trait[, 2, drop = FALSE]

gp <- create.gpData(pheno = trait, geno = data_IBS$geno.off, map = map2,
                    family = family, map.unit = "Mb")
# impute

# For the moment as Beagle is not working. I will keep family when it is working
# I will change family for BeagleAfterFamily.

gp.imp <- codeGeno(gpData = gp, impute = TRUE, impute.type = "family",
                   maf = .001, nmiss = .2, label.heter = NULL, verbose = TRUE)

data_short_IBS <- mppData_form(geno.off = gp.imp$geno,
                               geno.par = data_IBS$geno.par, IBS = TRUE,
                               IBS.format = "012", type = "dh",
                               map = data_IBS$map, trait = data_IBS$trait,
                               cross.ind = data_IBS$cross.ind,
                               par.per.cross = data_IBS$par.per.cross)

# here I can fill the geno.par and the allele.ref

allele.ref <- geno_012(mk.mat = data_IBS$geno.off)[[2]]
data_short_IBS$allele.ref <- allele.ref
data_short_IBS$geno.par <- data_IBS$geno.par

# then save the IBS mppData object

saveRDS(data_short_IBS,file = "./data/mpp_data/data_short_IBS.rds")
#previously data.short.biall

# 5.3 parental lines clustering
###############################

### This part need to be executed on R 2.14 to be able to run clusthaplo
# and RHmm at the same time !!!

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

library("clusthaplo")
library(RHmm)
library(mppR)

# load data

hap_map <- read.table("./data/map/hap_map_sh.txt", stringsAsFactors = FALSE)
mk_data <- read.table("./data/geno/hap_geno_sh.txt", stringsAsFactors = FALSE)

data <- readRDS("./data/mpp_data/data_short_IBD.rds")  
cons_map <- data$map[, -3]

set.seed(26589)

cluster <- parent_cluster(haplo.map = hap_map, consensus.map = cons_map,
                          marker.data = t(mk_data), na.strings = NA,
                          step.size = 50, clustering.method = "hmm", K = 10,
                          window = 2, plot = FALSE)

cluster[[2]] #average number of ancestors

par.clu.short <- cluster[[1]] # clustering results

# save the parent clustering results

write.table(par.clu.short,"./data/clustering/par_clu_short2.txt")


################################################################################
############################## End short subset ################################
################################################################################

###########################################
# 6. heterogeneous subset data processing #
###########################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

library(mppR)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv", row.names = 1)
geno <- as.matrix(geno)

# select the subset of 361 genotypes

geno_used <- read.table("./data/geno/geno_sel_list.txt")
parent_used <- read.table("./data/geno/cross.partition.txt",
                          stringsAsFactors = FALSE)

geno.off <- geno[rownames(geno) %in% geno_used[, "het"], ]
geno.par <- geno[rownames(geno) %in% c("F353", parent_used[, "par.hetero"]), ]

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)

trait <- pheno[rownames(pheno) %in% geno_used[, "het"], 1, drop = FALSE]
trait <- data.frame(rownames(trait), trait[, 1], stringsAsFactors = FALSE)

cross.ind <- substr(trait[, 1], 1, 5)

# map

map <- read.table("./data/map/map_panzea_full.txt", h = TRUE,
                  stringsAsFactors = FALSE)

# par.per.cross

par.per.cross <- cbind(unique(cross.ind), rep("F353", 5),
                       parent_used[, "par.hetero"])

colnames(par.per.cross) <- c("cr", "par1", "par2")


# 6.1 IBD mppData object
########################

data_IBD <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
                    trait = trait, cross.ind = cross.ind,
                    par.per.cross = par.per.cross, MAF.pop.lim = 0.01,
                    MAF.cr.lim =  c(0, 0, 0, 0, 0), MAF.cr.miss = TRUE,
                    ABH = TRUE)

# save data for parental clustering later

write.table(x = data_IBD$map.par.clu, file = "./data/map/hap_map_het.txt")
write.table(x = data_IBD$geno.par.clu, file = "./data/geno/hap_geno_het.txt")


data_het_IBD <- mppData_form(geno.off = data_IBD$geno.off,
                             geno.par = data_IBD$geno.par, IBS = FALSE,
                             type = "dh", map = data_IBD$map,
                             trait = data_IBD$trait,
                             cross.ind = data_IBD$cross.ind,
                             par.per.cross = data_IBD$par.per.cross,
                             step = 50, dir = "./data/geno")

saveRDS(data_het_IBD, file = "./data/mpp_data/data_hetero_IBD.rds")
# previously data.hetero.ABH

# 6.2 IBS mppData object
########################

data_IBS <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
                    trait = trait, cross.ind = cross.ind,
                    par.per.cross = par.per.cross, MAF.pop.lim = 0.01,
                    MAF.cr.lim =  c(0, 0, 0, 0, 0), MAF.cr.miss = TRUE,
                    ABH = FALSE)


# Imputation with Beagle

library(synbreed)

# form a gpData object

family <- data.frame(data_IBS$cross.ind)
rownames(family) <- rownames(data_IBS$geno.off)

map2 <- data_IBS$map
map2[, 1] <- as.character(map2[, 1])
rownames(map2) <- map2[, 1]
map2 <- map2[, 2:3]
colnames(map2) <- c("chr", "pos")

trait <- data_IBS$trait[, 2, drop = FALSE]

gp <- create.gpData(pheno = trait, geno = data_IBS$geno.off, map = map2,
                    family = family, map.unit = "Mb")
# impute

# For the moment as Beagle is not working. I will keep family when it is working
# I will change family for BeagleAfterFamily.

gp.imp <- codeGeno(gpData = gp, impute = TRUE, impute.type = "family",
                   maf = .001, nmiss = .2, label.heter = NULL, verbose = TRUE)

data_het_IBS <- mppData_form(geno.off = gp.imp$geno,
                             geno.par = data_IBS$geno.par, IBS = TRUE,
                             IBS.format = "012", type = "dh",
                             map = data_IBS$map, trait = data_IBS$trait,
                             cross.ind = data_IBS$cross.ind,
                             par.per.cross = data_IBS$par.per.cross)

# here I can fill the geno.par and the allele.ref

allele.ref <- geno_012(mk.mat = data_IBS$geno.off)[[2]]
data_het_IBS$allele.ref <- allele.ref
data_het_IBS$geno.par <- data_IBS$geno.par

# then save the IBS mppData object

saveRDS(data_het_IBS,file = "./data/mpp_data/data_hetero_IBS.rds")
#previously data.hetero.biall

# 6.3 parental lines clustering
###############################

### This part need to be executed on R 2.14 to be able to run clusthaplo
# and RHmm at the same time !!!

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

library("clusthaplo")
library(RHmm)
library(mppR)

# load data

hap_map <- read.table("./data/map/hap_map_het.txt", stringsAsFactors = FALSE)
mk_data <- read.table("./data/geno/hap_geno_het.txt", stringsAsFactors = FALSE)

data <- readRDS("./data/mpp_data/data_het_IBD.rds")  
cons_map <- data$map[, -3]

set.seed(3687918)

cluster <- parent_cluster(haplo.map = hap_map, consensus.map = cons_map,
                          marker.data = t(mk_data), na.strings = NA,
                          step.size = 50, clustering.method = "hmm", K = 10,
                          window = 2, plot = FALSE)

cluster[[2]] #average number of ancestors

par.clu.het <- cluster[[1]] # clustering results

# save the parent clustering results

write.table(par.clu.het,"./data/clustering/par_clu_hetero2.txt")


################################################################################
########################### End heterogeneous subset ###########################
################################################################################

##################################
# 7. Long subset data processing #
##################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

library(mppR)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv", row.names = 1)
geno <- as.matrix(geno)

# select the subset of 361 genotypes

geno_used <- read.table("./data/geno/geno_sel_list.txt")
parent_used <- read.table("./data/geno/cross.partition.txt",
                          stringsAsFactors = FALSE)

geno.off <- geno[rownames(geno) %in% geno_used[, "long"], ]
geno.par <- geno[rownames(geno) %in% c("F353", parent_used[, "par.long"]), ]

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)

trait <- pheno[rownames(pheno) %in% geno_used[, "long"], 1, drop = FALSE]
trait <- data.frame(rownames(trait), trait[, 1], stringsAsFactors = FALSE)

cross.ind <- substr(trait[, 1], 1, 5)

# map

map <- read.table("./data/map/map_panzea_full.txt", h = TRUE,
                  stringsAsFactors = FALSE)

# par.per.cross

par.per.cross <- cbind(unique(cross.ind), rep("F353", 5),
                       parent_used[, "par.long"])

colnames(par.per.cross) <- c("cr", "par1", "par2")


# 7.1 IBD mppData object
########################

data_IBD <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
                    trait = trait, cross.ind = cross.ind,
                    par.per.cross = par.per.cross, MAF.pop.lim = 0.01,
                    MAF.cr.lim =  c(0, 0, 0, 0, 0), MAF.cr.miss = TRUE,
                    ABH = TRUE)

# save data for parental clustering later

write.table(x = data_IBD$map.par.clu, file = "./data/map/hap_map_lg.txt")
write.table(x = data_IBD$geno.par.clu, file = "./data/geno/hap_geno_lg.txt")


data_long_IBD <- mppData_form(geno.off = data_IBD$geno.off,
                              geno.par = data_IBD$geno.par, IBS = FALSE,
                              type = "dh", map = data_IBD$map,
                              trait = data_IBD$trait,
                              cross.ind = data_IBD$cross.ind,
                              par.per.cross = data_IBD$par.per.cross,
                              step = 50, dir = "./data/geno")

saveRDS(data_long_IBD, file = "./data/mpp_data/data_long_IBD.rds")
# previously data.hetero.ABH

# 7.2 IBS mppData object
########################

data_IBS <- QC_proc(geno.off = geno.off, geno.par = geno.par, map = map,
                    trait = trait, cross.ind = cross.ind,
                    par.per.cross = par.per.cross, MAF.pop.lim = 0.01,
                    MAF.cr.lim =  c(0, 0, 0, 0, 0), MAF.cr.miss = TRUE,
                    ABH = FALSE)


# Imputation with Beagle

library(synbreed)

# form a gpData object

family <- data.frame(data_IBS$cross.ind)
rownames(family) <- rownames(data_IBS$geno.off)

map2 <- data_IBS$map
map2[, 1] <- as.character(map2[, 1])
rownames(map2) <- map2[, 1]
map2 <- map2[, 2:3]
colnames(map2) <- c("chr", "pos")

trait <- data_IBS$trait[, 2, drop = FALSE]

gp <- create.gpData(pheno = trait, geno = data_IBS$geno.off, map = map2,
                    family = family, map.unit = "Mb")
# impute

# For the moment as Beagle is not working. I will keep family when it is working
# I will change family for BeagleAfterFamily.

gp.imp <- codeGeno(gpData = gp, impute = TRUE, impute.type = "family",
                   maf = .001, nmiss = .2, label.heter = NULL, verbose = TRUE)

data_long_IBS <- mppData_form(geno.off = gp.imp$geno,
                              geno.par = data_IBS$geno.par, IBS = TRUE,
                              IBS.format = "012", type = "dh",
                              map = data_IBS$map, trait = data_IBS$trait,
                              cross.ind = data_IBS$cross.ind,
                              par.per.cross = data_IBS$par.per.cross)

# here I can fill the geno.par and the allele.ref

allele.ref <- geno_012(mk.mat = data_IBS$geno.off)[[2]]
data_long_IBS$allele.ref <- allele.ref
data_long_IBS$geno.par <- data_IBS$geno.par

# then save the IBS mppData object

saveRDS(data_long_IBS,file = "./data/mpp_data/data_long_IBS.rds")
#previously data.hetero.biall

# 7.3 parental lines clustering
###############################

### This part need to be executed on R 2.14 to be able to run clusthaplo
# and RHmm at the same time !!!

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

library("clusthaplo")
library(RHmm)
library(mppR)

# load data

hap_map <- read.table("./data/map/hap_map_lg.txt", stringsAsFactors = FALSE)
mk_data <- read.table("./data/geno/hap_geno_lg.txt", stringsAsFactors = FALSE)

data <- readRDS("./data/mpp_data/data_long_IBD.rds")  
cons_map <- data$map[, -3]

set.seed(6873)

cluster <- parent_cluster(haplo.map = hap_map, consensus.map = cons_map,
                          marker.data = t(mk_data), na.strings = NA,
                          step.size = 50, clustering.method = "hmm", K = 10,
                          window = 2, plot = FALSE)

cluster[[2]] #average number of ancestors

par.clu.lg <- cluster[[1]] # clustering results

# save the parent clustering results

write.table(par.clu.lg,"./data/clustering/par_clu_long2.txt")


################################################################################
################################ End long subset ###############################
################################################################################

############# stop there...

#########################################################
# 8. Significance threshold determination (permutation) #
#########################################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppR)


# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)


partition.perm <- expand.grid(c("short","het","long"),
                              c("par","anc","biall"),c("DMY","PH"))

partition.perm  <- partition.perm[order(partition.perm[,1]),]
partition.perm <- data.frame(partition.perm, rep(c("IBD","IBD","IBS"),6))
colnames(partition.perm) <- c("subset", "Q.eff","trait","ABH.bi")

N <- 1000
q.val <- c(0.9,0.95)

for(i in 1:18){
  
  par.exc <- as.character(unlist(par.exc[i, ]))
  
  file.name <- paste0("./data/mpp_data/data_",par.exc[1],"_", par.exc[4],".rds")
  
  # load the mppData object
  
  data <- readRDS(file = file.name)
  
  # put the right trait
  
  pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
  pheno.red <- pheno[rownames(pheno) %in% data$geno.id, par.exc[3], drop = FALSE]
  
  trait <- data.frame(rownames(pheno.red), pheno.red[, 1],
                      stringsAsFactors = FALSE)
  
  data <- mppData_chgPheno(trait = trait, mppData = data)
  
  # load the par.clu object
  
  if(par.exc[2] == "anc"){
    
    file.name <- paste0("./data/clustering/par_clu_", par.exc[1],"2.txt")
    
    par.clu <- as.matrix(read.table(file.name))
    
  }
  
  
  
  library(parallel)
  n.cores <- detectCores()
  cluster <- makeCluster((n.cores - 2))
  
  # compute the permutation test
  
  perm.thre <- mpp_perm(mppData = data, Q.eff = par.exc[2], VCOV = "h.err",
                        N = N, q.val = q.val, parallel = TRUE,
                        par.clu = par.clu, cluster = cluster)

  stopCluster(cl = cluster)
  
# save the results

# Origninal threshold result of the study are already located in the following
# folder ~/MPP_EUNAM/results/Threshold_permutation. The permutation results
# are therefore not saved otherwise they will overwrite the original results
# from the study. If you want nevertheless to save these results uncomment the
# next lines
  
  file.name <- paste0(paste(par.exc[1:3], collapse = "_"), "2.txt")
  folder <- paste0("./results/Threshold_permutation/", file.name)

  write.table(cbind(perm.thre$max.pval, perm.thre$seed), folder)
  
}

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

###############################################################################
############################ End threshold permutation ########################
###############################################################################

##################################
# 10. QTL analysis full datasets #
##################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppR)
library(asreml)

# # need to load version 0.7.1 of igraph library
# # https://cran.r-project.org/web/packages/igraph/index.html
# 
# library(igraph)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# partition

subset <- rep(c("short", "hetero", "long"), each = 12)
Q.eff <- rep(c("par", "anc", "biall"), times = 6, each = 2)
trait <- rep(c("DMY", "PH"), times = 3, each = 6)
VCOV <- rep(c("h.err", "cr.err"), times = 18)
data.type <- rep(c("IBD", "IBD", "IBS"), times = 6, each = 2)

part.QTL <- data.frame(subset, Q.eff, trait, VCOV, data.type)

folder <- "./results/Threshold_permutation/"
my.loc <- "./results/QTL_analyses"


for(i in 1:36){
  
  # extract the ith partition
  part <- as.character(unlist(part.QTL[i,]))
  
  # get the threshold value
  
  max.val <- read.table(paste0(folder, paste0(part[1],"_",part[2],"_",
                                              part[3],".txt")))
  
  thre <- quantile(max.val[, 1], 0.95)
  
  data <- readRDS(file = paste0("./data/mpp_data/data_", part[1],"_",
                                part[5],".rds"))
  
  # put the right trait
  
  pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
  pheno.red <- pheno[rownames(pheno) %in% data$geno.id, part[3], drop = FALSE]
  
  trait <- data.frame(rownames(pheno.red), pheno.red[, 1],
                      stringsAsFactors = FALSE)
  
  data <- mppData_chgPheno(trait = trait, mppData = data)
  
  
  # load the par.clu object
  
  if(part[2] == "anc"){
    
    file.name <- paste0("./data/clustering/par_clu_", part[1],"2.txt")
    
    par.clu <- as.matrix(read.table(file.name))
    
  }
  
  
  
  if(part[4] == "h.err") {
    
    library(parallel)
    parallel = TRUE
    n.cores <- detectCores()
    cluster <- makeCluster((n.cores-1))

  } else {

    parallel = FALSE

  }

  
  proc <- mpp_proc(pop.name = "Test", trait.name = part[3], mppData = data,
                   Q.eff = part[2], par.clu = par.clu, VCOV = part[4],
                   est.gen.eff = TRUE, thre.cof = thre,
                   N.cim = 2, thre.QTL = thre, alpha.bk = 0.01,
                   parallel = parallel, cluster = cluster,
                   silence.print = TRUE, output.loc = my.loc)

  if(parallel==TRUE){stopCluster(cl = cluster)}
  
}

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

###############################################################################
################################ End QTL analysis #############################
###############################################################################

############ stop there

################################
# 11. Cross-validation process #
################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppRDraft)

library(qtl)
library(asreml)
library(MASS)
library(stringr)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)


# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)


# partition

subset <- rep(c("short","hetero","long"),each=6)
Q.eff <- rep(c("par","anc","biall"),times = 6,each = 1)
trait <- rep(c("DMY","PH"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)
data.type <- rep(c("ABH","ABH","biall"), times=6, each=1)

part.CV <- data.frame(subset,Q.eff,trait,VCOV,data.type)

folder <- "./results/Threshold_permutation/"
my.loc <- "./results/CV"

# start the loop here

for(i in 1:18){
  
  # extract the ith partition
  part <- as.character(unlist(part.CV[i,]))
  
  # get the threshold value
  max.val <- read.table(paste0(folder, paste0(part[1],"_",part[2],"_",
                                              part[3],".txt")))
  
  thre <- quantile(max.val[,1],0.95)
  
  data <- readRDS(file = paste0("./data/mpp_data/data.",part[1],".",part[5],".rds"))
  
  # put the right trait
  
  pheno <- read.csv(paste0("./data/pheno/pheno.",part[1],".sorted.csv"))
  
  if(part[3]=="DMY"){
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$DMY,stringsAsFactors = F)
    heritability <- 0.57
    
  } else {
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$PH,stringsAsFactors = F)
    heritability <- 0.81
  }
  
  data <- sub.trait(data,trait)
  
  # load the par.clu object
  
  if(part[2]=="anc"){
    
    par.clu <- read.table(paste0("./data/clustering/par_clu_", part[1],".txt"))
    
  }
  
  # Cross-validation results are already saved in the following folder
  # ~/MPP_EUNAM/results/CV. If however you want to reproduce the
  # results, uncomment the part bellow. Be carefull that the code is programmed
  # in parallel! If you do not want to run it in parallel modify option parllel in
  # mpp.proc. Be also carefull with the construction of cluster
  # (cluster <- makeCluster((n.cores-1)); stopCluster(cl=cluser)).
  # The cross-validation results are product of randomisation. Therefore
  # reproduction will not give exactly the same results as the one of the 
  # study.
  
  # library(parallel)
  # n.cores <- detectCores()
  # cluster <- makeCluster((n.cores-1))
  # 
  # CV <- mpp.CV(pop.name = part[1], trait.name = part[3],data = data,
  #              parallel = TRUE,cluster = cluster, her = heritability,
  #              Rep = 10, k = 5,Q.eff = part[2],VCOV = part[4],
  #              par.clu = par.clu,cof.sel = "SIM",thre.cof = thre,
  #              cim = "CIM",N.cim = 2,thre.QTL = thre,output.loc = my.loc)
  # 
  # stopCluster(cl = cluster)
  
}

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

###############################################################################
################################ End CV procedure #############################
###############################################################################


##################################
# 12. Multi QTL effect model fit #
##################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppRDraft)

library(qtl)
library(asreml)
library(MASS)
library(stringr)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(xtable)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)


# partition

subset <- rep(c("short","hetero","long"),times=2)
trait <- rep(c("DMY","PH"),each=3)
VCOV <- rep("u.err", 6)

part.step <- data.frame(subset,trait,VCOV)

folder <- "./results/Threshold_permutation/"
my.loc <- "./results/MQeff"

# options(warn=2)


for (i in 1:6){
  
  part <- as.character(unlist(part.step[i,]))
  
  # get the threshold value
  
  # we do the average of the three type of QTL incidence matrices
  
  max.val.par <- read.table(paste0(folder, paste0(part[1],"_","par","_",
                                                  part[2],".txt")))
  
  max.val.anc <- read.table(paste0(folder, paste0(part[1],"_","anc","_",
                                                  part[2],".txt")))
  
  max.val.biall <- read.table(paste0(folder, paste0(part[1],"_","biall","_",
                                                    part[2],".txt")))
  
  thre1 <- quantile(max.val.par[,1],0.95)
  thre2 <- quantile(max.val.anc[,1],0.95)
  thre3 <- quantile(max.val.biall[,1],0.95)
  
  # take the average value of the three subsets thresholds
  
  thre <- mean(c(thre1,thre2,thre3))
  
  data.ABH <- readRDS(file = paste0("./data/mpp_data/data.",part[1],".","ABH",
                                    ".rds"))
  
  data.bi <- readRDS(file = paste0("./data/mpp_data/data.",part[1],".",
                                   "biall",".rds"))
  
  # put the right trait
  
  pheno <- read.csv(paste0("./data/pheno/pheno.",part[1],".sorted.csv"))
  
  if(part[2]=="DMY"){
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$DMY,stringsAsFactors = F)
    
  } else {
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$PH,stringsAsFactors = F)
    
  }
  
  data.ABH <- sub.trait(data.ABH,trait)
  data.bi <- sub.trait(data.bi,trait)
  
  
  # par.clu object
  
  par.clu <- read.table(paste0("./data/clustering/par_clu_", part[1],".txt"))
  
  # datasets
  
  Q.eff <- c("par","anc","biall")
  
  # Multi QTL effect model results are already saved in the following folder
  # ~/MPP_EUNAM/results/MQeff. If however you want to reproduce the
  # results, uncomment the part bellow. Be carefull that the code is programmed
  # in parallel! If you do not want to run it in parallel modify option parllel in
  # mpp.multi.Qeff.proc. Be also carefull with the construction of cluster
  # (cluster <- makeCluster((n.cores-1)); stopCluster(cl=cluser)).
  
  
  # library(parallel)
  # n.cores <- detectCores()
  # cluster <- makeCluster((n.cores-1))
  # 
  # QTL <- mpp.multi.Qeff.proc(pop.name = part[1], trait.name = part[2],
  #                            parallel = TRUE, cluster = cluster,
  #                            data.ABH = data.ABH,data.bi = data.bi,
  #                            Q.eff = Q.eff,par.clu = par.clu,VCOV = part[3],
  #                            threshold = thre, output.loc = my.loc)
  # 
  # stopCluster(cl = cluster)
  
  
}


variable <- ls()
rm(list=variable)
rm(variable)


################################################################################
###################### End stepwise multi-QTL effect model #####################
################################################################################

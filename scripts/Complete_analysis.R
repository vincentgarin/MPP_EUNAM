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

path <- "~/MPP_EUNAM/"

# path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

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

path <- "~/MPP_EUNAM/"

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

path <- "~/MPP_EUNAM/"

setwd(path)

# original complete map

map <- read.table("./data/map/Complete_map_Dent.txt", h = TRUE)

# information concerning 50k marker of the array

array <- read.csv("./data/geno/Plateform_Data.csv")
array <- array[, c(1, 3, 5, 6)]


# 2.1 Remove positions in the array or in the map that have the same identifier
###############################################################################

# check if there are some markers identifier that are used twice or more

length(table(array[, 2]))
length(table(map[, 2]))

# In both case, there are less identifier than positions.
# This means that some marker identifier are repeated on the map and on the array

# in order to make a correspondance between the map position and the marker
# in the data matrix we need to only have one possible identifier
# the position with two or more identifier will be removed from the map
# and from the marker matrix

# list of markers that are positioned twice or more in the map

rep_mk <- which(table(map[, 2]) > 1)
ID_mk <- attr(table(map[, 2])[rep_mk], "dimnames")

# remove the list of misplaced markers from the full map

map <- map[!(map[, 2] %in% ID_mk[[1]]), ]


# list of markers that are positioned twice or more in the map

rep_mk <- which(table(array[, 2]) > 1)
ID_mk <- attr(table(array[, 2])[rep_mk], "dimnames")

# remove the list of misplaced markers from the full map

array <- array[!(array[, 2] %in% ID_mk[[1]]), ]

# 2.2 match SNP_id and rsid
###########################

inter.mk.map <- intersect(array[, 2], map[, 2])

map <- map[map[, 2] %in% inter.mk.map, ]
array <- array[array[, 2] %in% inter.mk.map, ]

# combine the rs ID to have the same SNP id in the map and the array

# put array_sel in the same order as map data

array <- array[match(map$Locus, array$RSID), ]

map <- cbind(array, map)

map <- map[, c(1, 2, 5, 7, 3, 4)]

# 2.3 check the concordance between genetic and physical positioning
####################################################################

# remove the 108 markers that are physically misplaced

map <- map[map[, 6] != 99, ]

identical(map[, 3], map[, 5])

sum((map[, 3] - map[, 5]) != 0)

map[(map[, 3] - map[, 5]) != 0, ]

# remove also 3 markers that are placed on no chromosome

map <- map[map[, 5] !=0, ]

identical(map[, 3], map[, 5])

colnames(map) <- c("SNP_id", "RSID", "chr_gen", "cM", "chr_ph", "bp")


# 2.3 select only the Panzea markers
####################################

mk.origin <- substring(text = map[, 1], 1, 2)
map <- map[which(mk.origin == "PZ"), ]


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

path <- "~/MPP_EUNAM/"

setwd(path)

# 3.1 General sorting
#####################

# genotype matrix (! full array data, heavy file)

geno <- read.csv("./data/geno/geno_array_EUNAM.csv")

# geno <- read.csv("F:/EU_NAM/Raw_geno/geno_raw_EUNAM.csv")

# parents of the Dent panel

parents <- c("F353", "B73", "D06", "D09", "EC169", "F252", "F618", "Mo17",
             "UH250", "UH304", "W117")

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)

# map

map <- read.table("./data/map/map_panzea_full.txt", h = TRUE,
                  stringsAsFactors = FALSE)

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

# select only the used lines

entries.sel <- c(parents, rownames(pheno))

geno <- geno[rownames(geno) %in% entries.sel, ]

# 3.2 genotype marker matrix imputation
#######################################

library(synbreed)
options(scipen = 999)

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)

map <- read.table("./data/map/map_panzea_full.txt", row.names = 1,
                  stringsAsFactors = FALSE)

geno.imp <- geno[12:852, ]

# form a gpData object

family <- data.frame(substr(rownames(geno.imp), 1, 5))
rownames(family) <- rownames(geno.imp)

# order the markers within chromosome according to the physical map order

mk.list.ph <- c()

for(i in 1:10){

  # subset chr i

  map_i <- map[map$chr_ph == i,]
  map_i_ord <- map_i[order(map_i$bp), ]

  mk.list.ph <- rbind(mk.list.ph, map_i_ord[, c(1, 5, 6)])

}

# remove the marker that are at the same position

diff_mk <- c(diff(mk.list.ph[, 3]))
mk.same.pos <- mk.list.ph[diff_mk == 0, 1]
mk.list.ph <- mk.list.ph[!(mk.list.ph[, 1] %in% mk.same.pos), ]

# reorder the genotypes marker matrix according to the physical order

geno.imp <- geno.imp[, mk.list.ph[, 1]]

map_bp <- mk.list.ph[, 2:3]
colnames(map_bp) <- c("chr", "pos")
rownames(map_bp) <- mk.list.ph[, 1]

# imputation

gp <- create.gpData(geno = geno.imp, map = map_bp, family = family,
                    map.unit = "bp")

gp.imp <- codeGeno(gpData = gp, impute = TRUE, impute.type = "beagleAfterFamily",
                   maf = 0, nmiss = .5, label.heter = NULL, verbose = TRUE)

geno.imp <- gp.imp$geno

# 3.3 save data
###############

# select the marker used for imputation

geno <- geno[, colnames(geno.imp)]
map <- map[map[, 1] %in% colnames(geno.imp), ]

# order the markers of the genotype matrix according to genetic distances

geno <- geno[, map[, 1]]
geno.imp <- geno.imp[, map[, 1]]

# save the row parent genotypes

parents <- c("F353", "B73", "D06", "D09", "EC169", "F252", "F618", "Mo17",
             "UH250", "UH304", "W117")

geno.par <- geno[rownames(geno) %in% parents, ]

write.csv(geno.par, "./data/geno/geno_panzea_par.csv")

# save the raw genotype data

write.csv(geno, "./data/geno/geno_panzea_full.csv")

# save the imputed genotype data

write.csv(geno.imp, "./data/geno/geno_panzea_full_imp.csv")

# save modified map data

write.table(map, "./data/map/map_panzea_full.txt")


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

path <- "~/MPP_EUNAM/"

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

path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv", row.names = 1,
                 check.names = FALSE)
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
map <- map[, c(1, 3, 4)]
colnames(map) <- c("mk.id", "chr", "cM")

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


map <- data_short_IBD$map

# divide the map in bin of 1.05 cM

chr.ind <- unique(map[, 2])

bin.ind <- c()

for(i in chr.ind){

  map_i <- map[map[, 2] == i, ]

  int.mk.dist_i <- c(0, diff(map_i[, 4]))
  cum.dist <- 0
  bin.ind_i <- c()
  bin.nb <- 1

  for(j in 1:dim(map_i)[1]){

    cum.dist <- cum.dist + int.mk.dist_i[j]

    if(cum.dist < 1.05){

      bin.ind_i <- c(bin.ind_i, paste0("chr", i, "_", bin.nb))

    } else{

      bin.nb <- bin.nb + 1
      bin.ind_i <- c(bin.ind_i, paste0("chr",i, "_", bin.nb))
      cum.dist <- 0

    }

  }

  bin.ind <- c(bin.ind, bin.ind_i)

}


geno.off.red <- geno.off[, map[, 1]]
maf <- QC_MAF(mk.mat = geno.off.red) # marker allele frequency
maf.aug <- cbind(maf, 1:length(maf))

mk.sel <- c()

bin.ind.id <- unique(bin.ind)

for (i in 1:length(bin.ind.id)){

  maf_i <- maf.aug[(bin.ind == bin.ind.id[i]), , drop = FALSE]
  mk.sel <- c(mk.sel, maf_i[which.max(maf_i[, 1]), 2])

}

data_short_IBD_red <- mppData_subset(mk.list = mk.sel, mppData = data_short_IBD)

saveRDS(data_short_IBD_red, file = "./data/mpp_data/data_short_IBD_red.rds")


# 5.2 IBS mppData object
########################


# load imputed data

geno.imp <- read.csv("./data/geno/geno_panzea_full_imp.csv", row.names = 1,
                             check.names = FALSE)
geno.imp <- as.matrix(geno.imp)

geno.off <- geno.imp[data_short_IBD$geno.id, data_IBD$map[, 1]]

data_short_IBS <- mppData_form(geno.off = geno.off,
                               geno.par = data_IBD$geno.par, IBS = TRUE,
                               IBS.format = "012", type = "dh",
                               map = data_IBD$map, trait = data_IBD$trait,
                               cross.ind = data_IBD$cross.ind,
                               par.per.cross = data_IBD$par.per.cross)

# here I can fill the geno.par and the allele.ref

geno.raw <- geno[rownames(geno.off), colnames(geno.off)]
geno.par <- geno.par[, colnames(geno.off)]

allele.ref <- geno_012(mk.mat = geno.raw)[[2]]
data_short_IBS$allele.ref <- allele.ref
data_short_IBS$geno.par <- data.frame(data_short_IBS$map, t(geno.par),
                                      stringsAsFactors = FALSE)

# then save the IBS mppData object

saveRDS(data_short_IBS,file = "./data/mpp_data/data_short_IBS.rds")


data_short_IBS_red <- mppData_subset(mk.list = mk.sel, mppData = data_short_IBS)

saveRDS(data_short_IBS_red, file = "./data/mpp_data/data_short_IBS_red.rds")

# 5.3 parental lines clustering
###############################

### This part need to be executed on R 2.14 to be able to run clusthaplo
# and RHmm at the same time !!!

path <- "~/MPP_EUNAM/"

setwd(path)

library("clusthaplo")
library(RHmm)
library(mppR)

# load data

hap_map <- read.table("./data/map/hap_map_sh.txt", stringsAsFactors = FALSE)
mk_data <- read.table("./data/geno/hap_geno_sh.txt", stringsAsFactors = FALSE,
                      check.names = FALSE)

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

write.table(par.clu.short,"./data/clustering/par_clu_short.txt")


################################################################################
############################## End short subset ################################
################################################################################

###########################################
# 6. heterogeneous subset data processing #
###########################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

library(mppR)

path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv", row.names = 1,
                 check.names = FALSE)
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

map <- map[, c(1, 3, 4)]
colnames(map) <- c("mk.id", "chr", "cM")

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


map <- data_het_IBD$map

# divide the map in bin of 1.05 cM

# do it per chromosome

chr.ind <- unique(map[, 2])

bin.ind <- c()

for(i in chr.ind){

  map_i <- map[map[, 2] == i, ]

  int.mk.dist_i <- c(0, diff(map_i[, 4]))
  cum.dist <- 0
  bin.ind_i <- c()
  bin.nb <- 1

  for(j in 1:dim(map_i)[1]){

    cum.dist <- cum.dist + int.mk.dist_i[j]

    if(cum.dist < 1.05){

      bin.ind_i <- c(bin.ind_i, paste0("chr", i, "_", bin.nb))

    } else{

      bin.nb <- bin.nb + 1
      bin.ind_i <- c(bin.ind_i, paste0("chr",i, "_", bin.nb))
      cum.dist <- 0

    }

  }

  bin.ind <- c(bin.ind, bin.ind_i)

}

geno.off.red <- geno.off[, map[, 1]]
maf <- QC_MAF(mk.mat = geno.off.red) # marker allele frequency
maf.aug <- cbind(maf, 1:length(maf))

mk.sel <- c()

bin.ind.id <- unique(bin.ind)

for (i in 1:length(bin.ind.id)){

  maf_i <- maf.aug[(bin.ind == bin.ind.id[i]), , drop = FALSE]
  mk.sel <- c(mk.sel, maf_i[which.max(maf_i[, 1]), 2])

}

data_het_IBD_red <- mppData_subset(mk.list = mk.sel, mppData = data_het_IBD)

saveRDS(data_het_IBD_red, file = "./data/mpp_data/data_hetero_IBD_red.rds")


# 6.2 IBS mppData object
########################

# load imputed data

geno.imp <- read.csv("./data/geno/geno_panzea_full_imp.csv", row.names = 1,
                     check.names = FALSE)
geno.imp <- as.matrix(geno.imp)

geno.off <- geno.imp[data_het_IBD$geno.id, data_IBD$map[, 1]]

data_het_IBS <- mppData_form(geno.off = geno.off,
                               geno.par = data_IBD$geno.par, IBS = TRUE,
                               IBS.format = "012", type = "dh",
                               map = data_IBD$map, trait = data_IBD$trait,
                               cross.ind = data_IBD$cross.ind,
                               par.per.cross = data_IBD$par.per.cross)

# here I can fill the geno.par and the allele.ref

geno.raw <- geno[rownames(geno.off), colnames(geno.off)]
geno.par <- geno.par[, colnames(geno.off)]

allele.ref <- geno_012(mk.mat = geno.raw)[[2]]
data_het_IBS$allele.ref <- allele.ref
data_het_IBS$geno.par <- data.frame(data_het_IBS$map, t(geno.par),
                                      stringsAsFactors = FALSE)

# then save the IBS mppData object

saveRDS(data_het_IBS,file = "./data/mpp_data/data_hetero_IBS.rds")


data_het_IBS_red <- mppData_subset(mk.list = mk.sel, mppData = data_het_IBS)

saveRDS(data_het_IBS_red, file = "./data/mpp_data/data_hetero_IBS_red.rds")

# 6.3 parental lines clustering
###############################

### This part need to be executed on R 2.14 to be able to run clusthaplo
# and RHmm at the same time !!!

path <- "~/MPP_EUNAM/"

setwd(path)

library("clusthaplo")
library(RHmm)
library(mppR)

# load data

hap_map <- read.table("./data/map/hap_map_het.txt", stringsAsFactors = FALSE)
mk_data <- read.table("./data/geno/hap_geno_het.txt", stringsAsFactors = FALSE,
                      check.names = FALSE)

data <- readRDS("./data/mpp_data/data_hetero_IBD.rds")
cons_map <- data$map[, -3]

set.seed(3687918)

cluster <- parent_cluster(haplo.map = hap_map, consensus.map = cons_map,
                          marker.data = t(mk_data), na.strings = NA,
                          step.size = 50, clustering.method = "hmm", K = 10,
                          window = 2, plot = FALSE)

cluster[[2]] #average number of ancestors

par.clu.het <- cluster[[1]] # clustering results

# save the parent clustering results

write.table(par.clu.het,"./data/clustering/par_clu_hetero.txt")


################################################################################
########################### End heterogeneous subset ###########################
################################################################################

##################################
# 7. Long subset data processing #
##################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

library(mppR)

path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv", row.names = 1,
                 check.names = FALSE)
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

map <- map[, c(1, 3, 4)]
colnames(map) <- c("mk.id", "chr", "cM")

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


map <- data_long_IBD$map

# divide the map in bin of 1 cM

# do it per chromosome

chr.ind <- unique(map[, 2])

bin.ind <- c()

for(i in chr.ind){

  map_i <- map[map[, 2] == i, ]

  int.mk.dist_i <- c(0, diff(map_i[, 4]))
  cum.dist <- 0
  bin.ind_i <- c()
  bin.nb <- 1

  for(j in 1:dim(map_i)[1]){

    cum.dist <- cum.dist + int.mk.dist_i[j]

    if(cum.dist < 1){

      bin.ind_i <- c(bin.ind_i, paste0("chr", i, "_", bin.nb))

    } else{

      bin.nb <- bin.nb + 1
      bin.ind_i <- c(bin.ind_i, paste0("chr",i, "_", bin.nb))
      cum.dist <- 0

    }

  }

  bin.ind <- c(bin.ind, bin.ind_i)

}

geno.off.red <- geno.off[, map[, 1]]
maf <- QC_MAF(mk.mat = geno.off.red) # marker allele frequency
maf.aug <- cbind(maf, 1:length(maf))

mk.sel <- c()

bin.ind.id <- unique(bin.ind)

for (i in 1:length(bin.ind.id)){

  maf_i <- maf.aug[(bin.ind == bin.ind.id[i]), , drop = FALSE]
  mk.sel <- c(mk.sel, maf_i[which.max(maf_i[, 1]), 2])

}

data_long_IBD_red <- mppData_subset(mk.list = mk.sel, mppData = data_long_IBD)

saveRDS(data_long_IBD_red, file = "./data/mpp_data/data_long_IBD_red.rds")


# 7.2 IBS mppData object
########################

# load imputed data

geno.imp <- read.csv("./data/geno/geno_panzea_full_imp.csv", row.names = 1,
                     check.names = FALSE)
geno.imp <- as.matrix(geno.imp)

geno.off <- geno.imp[data_long_IBD$geno.id, data_IBD$map[, 1]]

data_long_IBS <- mppData_form(geno.off = geno.off,
                               geno.par = data_IBD$geno.par, IBS = TRUE,
                               IBS.format = "012", type = "dh",
                               map = data_IBD$map, trait = data_IBD$trait,
                               cross.ind = data_IBD$cross.ind,
                               par.per.cross = data_IBD$par.per.cross)

# here I can fill the geno.par and the allele.ref

geno.raw <- geno[rownames(geno.off), colnames(geno.off)]
geno.par <- geno.par[, colnames(geno.off)]

allele.ref <- geno_012(mk.mat = geno.raw)[[2]]
data_long_IBS$allele.ref <- allele.ref
data_long_IBS$geno.par <- data.frame(data_long_IBS$map, t(geno.par),
                                      stringsAsFactors = FALSE)

# then save the IBS mppData object

saveRDS(data_long_IBS,file = "./data/mpp_data/data_long_IBS.rds")

# data_long_IBS <- readRDS("./data/mpp_data/data_long_IBS.rds")

data_long_IBS_red <- mppData_subset(mk.list = mk.sel, mppData = data_long_IBS)

saveRDS(data_long_IBS_red, file = "./data/mpp_data/data_long_IBS_red.rds")

# 7.3 parental lines clustering
###############################

### This part need to be executed on R 2.14 to be able to run clusthaplo
# and RHmm at the same time !!!

path <- "~/MPP_EUNAM/"

setwd(path)

library("clusthaplo")
library(RHmm)
library(mppR)

# load data

hap_map <- read.table("./data/map/hap_map_lg.txt", stringsAsFactors = FALSE)
mk_data <- read.table("./data/geno/hap_geno_lg.txt", stringsAsFactors = FALSE,
                      check.names = FALSE)

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

write.table(par.clu.lg,"./data/clustering/par_clu_long.txt")


################################################################################
################################ End long subset ###############################
################################################################################


#########################################################
# 8. Significance threshold determination (permutation) #
#########################################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppR)


# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)

# partition

subset <- rep(c("short", "hetero", "long"), each = 12)
Q.eff <- rep(c("par", "anc", "biall"), times = 6, each = 2)
trait <- rep(c("DMY", "PH"), times = 3, each = 6)
VCOV <- rep(c("h.err", "cr.err_fast"), times = 18)
data.type <- rep(c("IBD", "IBD", "IBS"), times = 6, each = 2)

part.perm <- data.frame(subset, Q.eff, trait, VCOV, data.type)

N <- 1000
q.val <- c(0.9, 0.95)

library(parallel)
cluster <- makeCluster(6)


for(i in 1:36){

  par.exc <- as.character(unlist(part.perm[i, ]))

  file.name <- paste0("./data/mpp_data/data_", par.exc[1], "_",
                      par.exc[5], ".rds")

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

    file.name <- paste0("./data/clustering/par_clu_", par.exc[1],".txt")

    par.clu <- as.matrix(read.table(file.name))

  }

  # compute the permutation test

  print(par.exc[1:4])

  perm.thre <- mpp_perm(mppData = data, Q.eff = par.exc[2], VCOV = par.exc[4],
                        N = N, q.val = q.val, parallel = TRUE,
                        par.clu = par.clu, cluster = cluster)


  # save the results

  file.name <- paste0(paste(par.exc[1:4], collapse = "_"), ".txt")
  folder <- paste0("./results/Threshold_permutation/", file.name)

  write.table(cbind(perm.thre$max.pval, perm.thre$seed), folder)

}


###############################################################################
############################ End threshold permutation ########################
###############################################################################

##################################
# 9. QTL analysis full datasets #
##################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppR)
library(asreml)
library(parallel)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)

# partition

subset <- rep(c("short", "hetero", "long"), each = 12)
Q.eff <- rep(c("par", "anc", "biall"), times = 6, each = 2)
trait <- rep(c("DMY", "PH"), times = 3, each = 6)
VCOV <- rep(c("h.err", "cr.err"), times = 18)
data.type <- rep(c("IBD", "IBD", "IBS"), times = 6, each = 2)

part.QTL <- data.frame(subset, Q.eff, trait, VCOV, data.type)

folder <- "./results/Threshold_permutation/"
res.loc <- "./results/QTL_analyses"


for(i in 1:36){

  # extract the ith partition
  part <- as.character(unlist(part.QTL[i,]))

  print(part[1:4])

  # get the threshold value

  max.val <- read.table(paste0(folder, paste0(paste(part[1:3], collapse = "_"),
                                              "_h.err.txt")))

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

    file.name <- paste0("./data/clustering/par_clu_", part[1], ".txt")
    par.clu <- as.matrix(read.table(file.name))

  }

  if(part[4] == "h.err") {

    parallel = TRUE
    cluster <- makeCluster(4)

  } else {

    parallel = FALSE

  }

  if(part[2] == "biall") est.gen.eff <- FALSE else est.gen.eff <- TRUE

  proc <- mpp_proc(pop.name = part[1], trait.name = part[3], mppData = data,
                   Q.eff = part[2], par.clu = par.clu, VCOV = part[4],
                   est.gen.eff = est.gen.eff, thre.cof = thre,
                   N.cim = 2, thre.QTL = thre, alpha.bk = 0.01,
                   parallel = parallel, cluster = cluster,
                   silence.print = TRUE, output.loc = res.loc)

  if(parallel==TRUE){stopCluster(cl = cluster)}

}



###############################################################################
################################ End QTL analysis #############################
###############################################################################

##############################################################
# 10. Significance threshold determination on reduced subset #
##############################################################

# library

# you can find in the ~MPP_EUNAM/software a tar.gz version of the package
# mppRDraft build for this research

library(mppR)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)

# partition

subset <- rep(c("short", "hetero", "long"), each = 12)
Q.eff <- rep(c("par", "anc", "biall"), times = 6, each = 2)
trait <- rep(c("DMY", "PH"), times = 3, each = 6)
VCOV <- rep(c("h.err", "cr.err_fast"), times = 18)
data.type <- rep(c("IBD", "IBD", "IBS"), times = 6, each = 2)

part.perm <- data.frame(subset, Q.eff, trait, VCOV, data.type)

N <- 1000
q.val <- c(0.9, 0.95)

library(parallel)
cluster <- makeCluster(6)


for(i in 1:36){

  par.exc <- as.character(unlist(part.perm[i, ]))

  file.name <- paste0("./data/mpp_data/data_", par.exc[1], "_",
                      par.exc[5], "_red.rds")

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

    file.name <- paste0("./data/clustering/par_clu_", par.exc[1],".txt")

    par.clu <- as.matrix(read.table(file.name))
    par.clu <- par.clu[data$map[, 1], ]

  }

  # compute the permutation test

  print(par.exc[1:4])

  perm.thre <- mpp_perm(mppData = data, Q.eff = par.exc[2], VCOV = par.exc[4],
                        N = N, q.val = q.val, parallel = TRUE,
                        par.clu = par.clu, cluster = cluster)


  # save the results

  file.name <- paste0(paste(par.exc[1:4], collapse = "_"), ".txt")
  folder <- paste0("./results/Threshold_CV/", file.name)

  write.table(cbind(perm.thre$max.pval, perm.thre$seed), folder)

}

###############################################################################
########################## End threshold reduced subset #######################
###############################################################################

#########################
# 11. Cross-validation  #
#########################

library(mppR)
library(asreml)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)

# partition

subset <- rep(c("short", "hetero", "long"), each = 12)
Q.eff <- rep(c("par", "anc", "biall"), times = 6, each = 2)
trait <- rep(c("DMY", "PH"), times = 3, each = 6)
VCOV <- rep(c("h.err", "cr.err_fast"), times = 18)
data.type <- rep(c("IBD", "IBD", "IBS"), times = 6, each = 2)

part.QTL <- data.frame(subset, Q.eff, trait, VCOV, data.type)

folder <- "./results/Threshold_CV/"
res.loc <- "./results/CV"

library(parallel)
parallel <-  TRUE
cluster <- makeCluster(6)

# heritability and cross size information

heritability.DMY <- c(0.59, 0.78, 0.42, 0.31, 0.71, 0.61, 0.52, 0.6, 0.56, 0.64)
heritability.PH <- c(0.69, 0.8, 0.87, 0.8, 0.74, 0.85, 0.76, 0.85, 0.85, 0.81)

cr.names <- c("CFD11", "CFD06", "CFD04", "CFD07", "CFD03", "CFD10", "CFD09",
              "CFD12", "CFD05","CFD02")

heritability <- cbind(heritability.DMY, heritability.PH)
rownames(heritability) <- cr.names
colnames(heritability) <- c("DMY", "PH")


for(i in 1:36){

  # extract the ith partition
  part <- as.character(unlist(part.QTL[i,]))

  print(i)
  print(part[1:4])

  # get the threshold value

  max.val <- read.table(paste0(folder, paste0(part[1:4], collapse = "_"), ".txt"))

  thre <- quantile(max.val[, 1], 0.95)

  data <- readRDS(file = paste0("./data/mpp_data/data_", part[1],"_",
                                part[5],"_red.rds"))

  # put the right trait

  pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
  pheno.red <- pheno[rownames(pheno) %in% data$geno.id, part[3], drop = FALSE]

  trait <- data.frame(rownames(pheno.red), pheno.red[, 1],
                      stringsAsFactors = FALSE)

  data <- mppData_chgPheno(trait = trait, mppData = data)

  her <- heritability[unique(data$cross.ind), part[3]]


  # load the par.clu object

  if(part[2] == "anc"){

    file.name <- paste0("./data/clustering/par_clu_", part[1],".txt")

    par.clu <- as.matrix(read.table(file.name))
    par.clu <- par.clu[data$map[, 1], ]


  } else {par.clu <- NULL}


  CV <- mpp_CV(pop.name = part[1], trait.name = part[3],
               mppData = data, her = her, Rep = 20, k = 5, Q.eff = part[2],
               par.clu = par.clu, VCOV = part[4], thre.cof = thre, N.cim = 2,
               thre.QTL = thre, alpha.bk = 0.01, parallel = TRUE,
               cluster = cluster, silence.print = TRUE,
               output.loc = res.loc)

}

###############################################################################
################################ End CV procedure #############################
###############################################################################

###########################################
# 13. Multi QTL effect model full dataset #
###########################################

library(mppR)
library(nlme)
library(asreml)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)


# partition

part.MQE <- expand.grid(c("short", "hetero", "long"), c("DMY", "PH"),
                        c("h.err", "cr.err"), stringsAsFactors = FALSE)


folder <- "./results/Threshold_permutation/"
res.loc <- "./results/MQE"

library(parallel)
cluster <- makeCluster(6)


for (i in 1:dim(part.MQE)[1]){

  part <- as.character(unlist(part.MQE[i,]))

  print(part)

  # get the threshold value

  # we do the average of the three type of QTL incidence matrices

  max.val.par <- read.table(paste0(folder, paste0(part[1],"_","par","_",
                                                  part[2],"_","h.err.txt")))

  max.val.anc <- read.table(paste0(folder, paste0(part[1],"_","anc","_",
                                                  part[2],"_","h.err.txt")))

  max.val.biall <- read.table(paste0(folder, paste0(part[1],"_","biall","_",
                                                    part[2],"_","h.err.txt")))

  thre1 <- quantile(max.val.par[, 1],0.95)
  thre2 <- quantile(max.val.anc[, 1],0.95)
  thre3 <- quantile(max.val.biall[, 1],0.95)

  # take the average value of the three subsets thresholds

  thre <- mean(c(thre1, thre2, thre3))

  data <- readRDS(file = paste0("./data/mpp_data/data_", part[1],
                                "_IBD.rds"))

  data.bi <- readRDS(file = paste0("./data/mpp_data/data_", part[1],
                                   "_IBS.rds"))

  # put the right trait

  pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
  pheno.red <- pheno[rownames(pheno) %in% data$geno.id, part[2], drop = FALSE]

  trait <- data.frame(rownames(pheno.red), pheno.red[, 1],
                      stringsAsFactors = FALSE)

  data <- mppData_chgPheno(trait = trait, mppData = data)
  data.bi <- mppData_chgPheno(trait = trait, mppData = data.bi)


  # par.clu object

  par.clu <- read.table(paste0("./data/clustering/par_clu_", part[1],".txt"))
  par.clu <- as.matrix(par.clu)

  Q.eff <- c("par","anc","biall")

  if(part[3] == "h.err") parallel <- TRUE else parallel <- FALSE


  MQE <- MQE_proc(pop.name = part[1], trait.name = part[2], mppData = data,
                  mppData_bi = data.bi, Q.eff = Q.eff, par.clu = par.clu,
                  VCOV = part[3], threshold = thre, parallel = parallel,
                  cluster = cluster, output.loc = res.loc)

}


###############################################################################
################################ End MQE model ################################
###############################################################################


#################################
# 13. Multi QTL effect model CV #
#################################

library(mppR)
library(nlme)
library(asreml)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "~/MPP_EUNAM/"

setwd(path)


# partition

part.MQE <- expand.grid(c("short", "hetero", "long"), c("DMY", "PH"),
                        c("h.err", "cr.err_fast"), stringsAsFactors = FALSE)


folder <- "./results/Threshold_CV/"
res.loc <- "./results/MQE_CV"

library(parallel)
cluster <- makeCluster(6)

# heritability information

heritability.DMY <- c(0.59, 0.78, 0.42, 0.31, 0.71, 0.61, 0.52, 0.6, 0.56, 0.64)
heritability.PH <- c(0.69, 0.8, 0.87, 0.8, 0.74, 0.85, 0.76, 0.85, 0.85, 0.81)

cr.names <- c("CFD11", "CFD06", "CFD04", "CFD07", "CFD03", "CFD10", "CFD09",
              "CFD12", "CFD05","CFD02")

heritability <- cbind(heritability.DMY, heritability.PH)
rownames(heritability) <- cr.names
colnames(heritability) <- c("DMY", "PH")


for (i in 1:dim(part.MQE)[1]){

  part <- as.character(unlist(part.MQE[i,]))

  print(part)

  # get the threshold value

  # we do the average of the three type of QTL incidence matrices

  max.val.par <- read.table(paste0(folder, paste0(part[1],"_","par","_",
                                                  part[2],"_", part[3],".txt")))

  max.val.anc <- read.table(paste0(folder, paste0(part[1],"_","anc","_",
                                                  part[2],"_", part[3],".txt")))

  max.val.biall <- read.table(paste0(folder, paste0(part[1],"_","biall","_",
                                                    part[2],"_",part[3],".txt")))

  thre1 <- quantile(max.val.par[,1],0.95)
  thre2 <- quantile(max.val.anc[,1],0.95)
  thre3 <- quantile(max.val.biall[,1],0.95)

  # take the average value of the three subsets thresholds

  thre <- mean(c(thre1, thre2, thre3))

  data <- readRDS(file = paste0("./data/mpp_data/data_", part[1],
                                "_IBD_red.rds"))

  data.bi <- readRDS(file = paste0("./data/mpp_data/data_", part[1],
                                   "_IBS_red.rds"))

  # put the right trait

  pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
  pheno.red <- pheno[rownames(pheno) %in% data$geno.id, part[2], drop = FALSE]

  trait <- data.frame(rownames(pheno.red), pheno.red[, 1],
                      stringsAsFactors = FALSE)

  data <- mppData_chgPheno(trait = trait, mppData = data)
  data.bi <- mppData_chgPheno(trait = trait, mppData = data.bi)


  # par.clu object

  par.clu <- read.table(paste0("./data/clustering/par_clu_", part[1], ".txt"))
  par.clu <- as.matrix(par.clu)
  par.clu <- par.clu[data$map[, 1], ]

  Q.eff <- c("par","anc","biall")

  her <- heritability[unique(data$cross.ind), part[2]]


  MQE <- MQE_CV(pop.name = part[1], trait.name = part[2], mppData = data,
                mppData_bi = data.bi, her = her, Rep = 20, k = 5, Q.eff = Q.eff,
                par.clu = par.clu, VCOV = part[3], threshold = thre,
                parallel = TRUE, cluster = cluster, silence.print = TRUE,
                output.loc = res.loc)

}

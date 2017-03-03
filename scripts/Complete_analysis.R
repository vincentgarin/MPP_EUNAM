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
# mppRDraft build for this research

library(mppRDraft)

library(qtl)
library(asreml)
library(MASS)
library(outliers)
library(stringr)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)

# need to load version 0.7.1 of igraph library
# https://cran.r-project.org/web/packages/igraph/index.html

library(igraph) 

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
col.names[c(2,7)] <- c("Fam","Block")
colnames(pheno) <- col.names
rm(col.names)

par_names <- c("B73","D06", "D09", "EC169", "F252", "F353", "F618", "F98902",
               "Mo17", "UH250", "UH304","W117")

# change MO17 for Mo17
pheno$Genotype <- as.character(pheno$Genotype)
pheno$Fam <- as.character(pheno$Fam)
pheno$Pedigree <- as.character(pheno$Pedigree)
pheno[which(pheno[,1]=="MO17"),1:3] <- "Mo17"
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
  
  median.i <- median(pheno$NBPL[pheno$LOC==loc.id[i]])
  median.NBPL <- c(median.NBPL,median.i)
  
}

# determine per location which plot contain 70% less plant than the number of
# data

missing.plot <- c()

for(i in 1:4){
  
  missing.plot.i <- pheno$NBPL[pheno$LOC==loc.id[i]]<(median.NBPL[i]*0.7)
  missing.plot <- c(missing.plot, missing.plot.i)
  
}


# Put missins values for DMY and DMC

pheno[missing.plot,9:10] <- NA

# 1.3 Perform the Grubb test to detect outlying plot
###################################################

# Grubb test applied on the residual of a mixed model. If the residual of the
# mixed model contained an outlying observation (Grubb test pval< 0.05), the 
# plot value for this specific trait with the largest residual is put as missing.
# This procedure is applied iteratively until no extreme value is detected.
# The test is done for the 5 traits.


# iteration on all the five traits

trait.id <- c("DMY","DMC","PH","DtTAS","DtSILK")
pos.id <- c(9,10,11,12,13)

for (i in 1:5){
  
  fix_formula <- paste(trait.id[i],"~1",sep="")
  
  # set initial value for p_val
  
  p_val <- 0.01
  count <- 0
  
  while(p_val<0.05){
    
    # compute the mixed model
    
    model <- asreml(fixed = as.formula(fix_formula), random= ~ Genotype + LOC+
                      Genotype:LOC + LOC:Rep +  LOC:Rep:Block,
                    rcov= ~ idv(units),
                    data = pheno, trace=F)
    
    # compute the residuals
    res <- model$residuals
    
    # do the Grubbs test on the residuals
    test <- grubbs.test(res,type=10,two.sided = TRUE)
    
    # get the p-val of the test
    p_val <- test$p.value
    
    # test if p_val is lower than 0.05 (presence of an outlier)
    if(p_val<0.05){
      
      # remove the outlying value from the data
      
      # check if there are missing values within the residuals
      
      if(sum(is.na(res))>0){
        
        pheno[which(abs(res)==max(abs(res[-which(is.na(res))]))),pos.id[i]] <- NA
        count <- count+1
        
        
      } else {
        
        
        pheno[which(abs(res)==max(abs(res))),pos.id[i]] <- NA
        count <- count+1
        
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

pheno <- read.csv("./data/pheno/pheno_red.csv")
pheno <- pheno[,-1]

lines_used <- read.csv("./data/pheno/List_lines_Dent_Lehermeier.csv")

# suppress the factor levels with 0 occurences

pheno$Genotype <- as.factor(as.character(pheno$Genotype))
pheno$Fam <- as.factor(as.character(pheno$Fam))

# order data by family

pheno <- pheno[order(pheno$Fam),]

par_names <- c("B73","D06", "D09", "EC169", "F252", "F353", "F618", "F98902",
               "Mo17", "UH250", "UH304","W117")

# 1.5.1 DMY
###########

# the model correspond to model 1 in Lehermeier et al. (2014)

model <- asreml(fixed = DMY~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

var.DMY <- summary(model)$varcomp
var.DMY <-round(cbind(var.DMY[c(2:11),c(2:3)], var.DMY[c(13:22),c(2:3)],
                      var.DMY[c(26:35),2]),2)
her <- round(100*(var.DMY[,1]/(var.DMY[,1] + (var.DMY[,3]/4) + (var.DMY[,5]/5))),2)
var.DMY <- cbind(var.DMY[,-5], her)
rownames(var.DMY) <- par_names[-c(6,8)]
colnames(var.DMY) <- c("sigma.g", "std.err", "sigma.ge", "std.err","heritability")
var.DMY

# adjusted means

model <- asreml(fixed = DMY~1 + Genotype, random= ~  LOC:Genotype + LOC:Rep
                + LOC:Rep:Block, rcov=~idv(units), data = pheno,
                na.method.X = "omit", na.method.Y="omit")

pred <- predict(model,classify="Genotype")
pred.DMY <- pred$predictions$pvals$predicted.value

# select same lines as Lehermeier et al.

pred.DMY <- data.frame(pred$predictions$pvals$Genotype,pred.DMY)
pred.DMY <- pred.DMY[pred.DMY[,1] %in% lines_used[,1],]

# 1.5.2 DMC
###########

model <- asreml(fixed = DMC~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

var.DMC <- summary(model)$varcomp
var.DMC <-round(cbind(var.DMC[c(2:11),c(2:3)], var.DMC[c(13:22),c(2:3)],
                      var.DMC[c(26:35),2]),2)
her <- round(100*(var.DMC[,1]/(var.DMC[,1] + (var.DMC[,3]/4) + (var.DMC[,5]/5))),2)
var.DMC <- cbind(var.DMC[,-5], her)
rownames(var.DMC) <- par_names[-c(6,8)]
colnames(var.DMC) <- c("sigma.g", "std.err", "sigma.ge", "std.err","heritability")
var.DMC

# adjusted means

model <- asreml(fixed = DMC~1 + Genotype, random= ~  LOC:Genotype + LOC:Rep
                + LOC:Rep:Block, rcov=~idv(units), data = pheno,
                na.method.X = "omit", na.method.Y="omit")

pred <- predict(model,classify="Genotype")
pred.DMC <- pred$predictions$pvals$predicted.value

# select same lines as Lehermeier et al.

pred.DMC <- data.frame(pred$predictions$pvals$Genotype,pred.DMC)
pred.DMC <- pred.DMC[pred.DMC[,1] %in% lines_used[,1],]

# 1.5.3 PH
###########

model <- asreml(fixed = PH~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

var.PH <- summary(model)$varcomp
var.PH <-round(cbind(var.PH[c(2:11),c(2:3)], var.PH[c(13:22),c(2:3)],
                     var.PH[c(26:35),2]),2)
her <- round(100*(var.PH[,1]/(var.PH[,1] + (var.PH[,3]/4) + (var.PH[,5]/5))),2)
var.PH <- cbind(var.PH[,-5], her)
rownames(var.PH) <- par_names[-c(6,8)]
colnames(var.PH) <- c("sigma.g", "std.err", "sigma.ge", "std.err","heritability")
var.PH

# adjusted means

model <- asreml(fixed = PH~1 + Genotype, random= ~  LOC:Genotype + LOC:Rep
                + LOC:Rep:Block, rcov=~idv(units), data = pheno,
                na.method.X = "omit", na.method.Y="omit")

pred <- predict(model,classify="Genotype")
pred.PH <- pred$predictions$pvals$predicted.value

# select same lines as Lehermeier et al.

pred.PH <- data.frame(pred$predictions$pvals$Genotype,pred.PH)
pred.PH <- pred.PH[pred.PH[,1] %in% lines_used[,1],]

# 1.5.4 DtTAS
##############

model <- asreml(fixed = DtTAS~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

var.DtTAS <- summary(model)$varcomp
var.DtTAS <-round(cbind(var.DtTAS[c(2:11),c(2:3)], var.DtTAS[c(13:22),c(2:3)],
                        var.DtTAS[c(26:35),2]),2)
her <- round(100*(var.DtTAS[,1]/(var.DtTAS[,1] + (var.DtTAS[,3]/4) + (var.DtTAS[,5]/5))),2)
var.DtTAS <- cbind(var.DtTAS[,-5], her)
rownames(var.DtTAS) <- par_names[-c(6,8)]
colnames(var.DtTAS) <- c("sigma.g", "std.err", "sigma.ge", "std.err","heritability")
var.DtTAS

# adjusted means

model <- asreml(fixed = DtTAS~1 + Genotype, random= ~  LOC:Genotype + LOC:Rep
                + LOC:Rep:Block, rcov=~idv(units), data = pheno,
                na.method.X = "omit", na.method.Y="omit")

pred <- predict(model,classify="Genotype")
pred.DtTAS <- pred$predictions$pvals$predicted.value

# select same lines as Lehermeier et al.

pred.DtTAS <- data.frame(pred$predictions$pvals$Genotype,pred.DtTAS)
pred.DtTAS <- pred.DtTAS[pred.DtTAS[,1] %in% lines_used[,1],]

# 1.5.5 DtSILK
##############

model <- asreml(fixed = DtSILK~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

var.DtSILK <- summary(model)$varcomp
var.DtSILK <-round(cbind(var.DtSILK[c(2:11),c(2:3)], var.DtSILK[c(13:22),c(2:3)],
                         var.DtSILK[c(26:35),2]),2)
her <- round(100*(var.DtSILK[,1]/(var.DtSILK[,1] + (var.DtSILK[,3]/4) + (var.DtSILK[,5]/5))),2)
var.DtSILK <- cbind(var.DtSILK[,-5], her)
rownames(var.DtSILK) <- par_names[-c(6,8)]
colnames(var.DtSILK) <- c("sigma.g", "std.err", "sigma.ge", "std.err","heritability")
var.DtSILK

# adjusted means

model <- asreml(fixed = DtSILK~1 + Genotype, random= ~  LOC:Genotype + LOC:Rep
                + LOC:Rep:Block, rcov=~idv(units), data = pheno,
                na.method.X = "omit", na.method.Y="omit")

pred <- predict(model,classify="Genotype")
pred.DtSILK <- pred$predictions$pvals$predicted.value

# select same lines as Lehermeier et al.

pred.DtSILK <- data.frame(pred$predictions$pvals$Genotype,pred.DtSILK)
pred.DtSILK <- pred.DtSILK[pred.DtSILK[,1] %in% lines_used[,1],]


# adjusted means matrix

Adj_means <- data.frame(pred.DMY[,2],pred.DMC[,2],pred.PH[,2],pred.DtTAS[,2],
                        pred.DtSILK[,2])
colnames(Adj_means) <- c("DMY","DMC","PH","DtTAS","DtSILK")
rownames(Adj_means) <- pred.DMY[,1]

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

map_dent <- read.table("./data/map/Complete_map_Dent.txt",h=T)

# information concerning 50k marker of the array

array <- read.csv("./data/geno/Plateform_Data.csv")
array <- array[,c(1,3,8)]

# 2.1 Remove positions in the array or in the map that have the same identifier
###############################################################################

# check if there are some markers identifier that are used twice or more

length(table(array[,2]))
length(table(map_dent[,2]))

# In both case, there are less identifier than positions.
# This means that some marker identifier are repeated on the map and on the array

# in order to make a correspondance between the map position and the marker
# in the data matrix we need to only have one possible identifier
# the position with two or more identifier will be removed from the map
# and from the marker matrix

# list of markers that are positioned twice or more in the map

rep_mk <- which(table(map_dent[,2])>1)
ID_mk <- attr(table(map_dent[,2])[rep_mk],"dimnames")

# remove the list of misplaced markers from the full map

map_dent_red <- map_dent[!(map_dent[,2]%in%ID_mk[[1]]),]


# list of markers that are positioned twice or more in the map

rep_mk <- which(table(array[,2])>1)
ID_mk <- attr(table(array[,2])[rep_mk],"dimnames")

# remove the list of misplaced markers from the full map

array_red <- array[!(array[,2]%in%ID_mk[[1]]),]

# 2.2 match SNP_id and rsid
###########################

inter.mk.map <- intersect(array_red[,2], map_dent_red[,2])

new_map_dent <- map_dent_red[map_dent_red[,2]%in%inter.mk.map,]
array_sel <- array_red[array_red[,2]%in%inter.mk.map,]

# combine the rs ID to have the same SNP id in the map and the array

# put array_sel in the same order as map data

array_sel_sort <- array_sel[match(new_map_dent$Locus,array_sel$RSID),]

new_map_dent2 <- cbind(array_sel_sort,new_map_dent)

map <- new_map_dent2[,c(1,2,3,4,6)]

# save map with rsid and SNPID

write.table(map,"./data/map/map_Dent_SNPID.txt")

# 2.3 select only the Panzea markers like in Giraud et al. (2014)
#################################################################

mk.origin <- substring(text = map[,1],1,2)
map <- map[which(mk.origin=="PZ"),]
map <- map[,c(1,4,5)]

# save map with all Panzea markers

write.table(map,"./data/map/map_panzea_full.txt")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
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

parents <- c("F353","B73","D06","D09","EC169","F252","F618","Mo17","UH250",
             "UH304","W117")

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv")
rownames(pheno) <- as.character(pheno[,1])
pheno <- pheno[,-1]

# map

map <- read.table("./data/map/map_panzea_full.txt",h=T)

# get the genotype identifier

geno.id <- geno[36,]

geno.id2 <- c()

for (i in 2:2291){
  
  name <- strsplit(as.character(geno.id[,i]), split= ",",fixed=TRUE)[[1]][2]
  
  geno.id2 <- c(geno.id2, name)
  
}

geno.id3 <- substring(geno.id2,2,nchar(geno.id2))

# Remove the row that are not the genotypic information

geno <- geno[-(1:73),]
colnames(geno) <- c("ID_ref",geno.id3)

### 3.1.1 select only the marker position that are present in the map

geno <- geno[geno[,1]%in%map[,1],]

# put the marker in the same order as in the map

geno <- geno[match(map[,1],geno[,1]),]

# reorder genotypes as row and markers as column and modify scores

geno <- t(geno)
geno <- as.matrix(geno)
colnames(geno) <- geno[1,]
geno <-geno[-1,]

geno <- gsub("AA","A",geno)
geno <- gsub("BB","B",geno)
geno <- gsub("AB","-",geno)
geno <- gsub("NC","-",geno)

# save the genotype matrix of the parents for Clusthaplo

geno.pz.par <- geno[rownames(geno) %in% parents,]

write.csv(geno.pz.par,"./data/geno/geno_panzea_par.csv")

### 3.1.2 select only the genotypes that have been phenotyped

# select the genotype that have been phenotyped and the parents

entries.sel <- c(parents, rownames(pheno))

geno <- geno[rownames(geno)%in%entries.sel,]

write.csv(geno,"./data/geno/geno_panzea_full.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

################################################################################
######################## End genotype matrix sorting ###########################
################################################################################

############################################
# 4. Determine the three subset populations #
############################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)

# SM coefficient distance between the central parent (F353)

kin.mat.par <- kinship.matrix(mk.mat = geno.par,method = "SM")

dist.cent.par <- kin.mat.par[,6]
dist.cent.par <- dist.cent.par[-6]

# relate information distance between the two parents to the cross
# information

crosses <- c("CFD02","CFD03","CFD04","CFD05","CFD06","CFD07",
             "CFD09","CFD10","CFD11","CFD12")

# parent per cross

parent1 <- as.character(rep("F353",10))
parent2 <- c("B73","D06","D09","EC169","F252","F618","Mo17",
             "UH250","UH304","W117")
par.per.cross <- cbind(crosses, parent1, parent2)

dist.cent.par <- data.frame(par.per.cross[,1],par.per.cross[,3],dist.cent.par)
colnames(dist.cent.par) <- c("cross","parent2","cor")
dist.cent.par

# order the crosses according to genetic distance between parents

dist.ord <- dist.cent.par[order(dist.cent.par$cor,decreasing = T),]
dist.ord

#  form the subsets

short.cross <- dist.ord$cross[1:5]
long.cross <- dist.ord$cross[6:10]
hetero.cross <- dist.ord$cross[c(1,3,5,8,10)]

par.short <- dist.ord$parent2[1:5]
par.long <- dist.ord$parent2[6:10]
par.hetero <- dist.ord$parent2[c(1,3,5,8,10)]

cross.partition <- data.frame(short.cross, long.cross, hetero.cross,par.short,
                              par.long,par.hetero)

# save

write.table(cross.partition,"./data/geno/cross.partition.txt")

# 4.1 PC bi-plot different subsets parents
##########################################

mpp.plot.PC(mk.mat = geno.par[rownames(geno.par)%in%c("F353",as.character(par.short)),],
            cross.ind = c("F353",as.character(par.short)), main = "Short subset")

mpp.plot.PC(mk.mat = geno.par[rownames(geno.par)%in%c("F353",as.character(par.long)),],
            cross.ind = c("F353",as.character(par.long)), main = "Long subset")

mpp.plot.PC(mk.mat = geno.par[rownames(geno.par)%in%c("F353",as.character(par.hetero)),],
            cross.ind = c("F353",as.character(par.hetero)), main = "Heterogeneous subset")

# 4.2 average distance between the parents of the different subsets
####################################################################

kin.mat.short <- kin.mat.par[rownames(kin.mat.par) %in% c("F353",as.character(par.short)),]
kin.mat.short <- kin.mat.short[,colnames(kin.mat.short) %in% c("F353",as.character(par.short))]

kin.mat.hetero <- kin.mat.par[rownames(kin.mat.par) %in% c("F353",as.character(par.hetero)),]
kin.mat.hetero <- kin.mat.hetero[,colnames(kin.mat.hetero) %in% c("F353",as.character(par.hetero))]

kin.mat.long <- kin.mat.par[rownames(kin.mat.par) %in% c("F353",as.character(par.long)),]
kin.mat.long <- kin.mat.long[,colnames(kin.mat.long) %in% c("F353",as.character(par.long))]

mean(kin.mat.short[lower.tri(kin.mat.short)])
mean(kin.mat.hetero[lower.tri(kin.mat.hetero)])
mean(kin.mat.long[lower.tri(kin.mat.long)])

sort(kin.mat.short[lower.tri(kin.mat.short)], decreasing = T)
sort(kin.mat.hetero[lower.tri(kin.mat.hetero)], decreasing = T)

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

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv")
rownames(pheno) <- as.character(pheno[,1])
pheno <- pheno[,-1]

cross.ind <- substr(rownames(pheno),1,5)

# cross partition

cross.partition <- read.table("./data/geno/cross.partition.txt")


# index for selection with parent and genotype belonging to the selected cross

index.sel <- c(rownames(geno)[1:11] %in% c("F353",as.character(cross.partition[,4])),
               cross.ind %in% cross.partition[,1])

geno.short <- geno[index.sel,]

# 5.1 subset line to equalize with long subset
##########################################

# take 361 lines to have same population size with
# the long subset

sel.ind <- sort(sample(7:486,361))

geno.short <- rbind(geno.short[1:6,], geno.short[sel.ind,])


# 5.2 Remove monomorphic and missing value positions
####################################################

# detect monomorphic position in the parents and remove

par.scores <- allele.scores(mk.mat = geno.short[1:6,])

par.mono <- which(par.scores[1,]=="mono")

geno.short <- geno.short[,-par.mono]

# detect monomorphic position in the offspring and remove

off.scores <- allele.scores(mk.mat = geno.short[7:dim(geno.short)[1],])

off.mono <- which(off.scores[1,]=="mono")

geno.short <- geno.short[,-off.mono]

# rare alleles with MAF <0.01

rare.all <- rare.allele.ind(mk.mat = geno.short[7:dim(geno.short)[1],],
                            threshold = 0.01)

geno.short <- geno.short[,-rare.all[,2]]

# missing markers with more than 10%

miss.ind.mk <- missing.indicator(mk.mat = geno.short[7:dim(geno.short)[1],],
                                 threshold=0.1,col=2)

geno.short  <- geno.short[,-miss.ind.mk[,2]]

# genotypes with more than 25% missing values

miss.ind.gen <- missing.indicator(mk.mat = geno.short[7:dim(geno.short)[1],],
                                  threshold=0.25,col=1)

# no line with more than 25% missing values

# save the variable

write.csv(geno.short,"./data/geno/geno.short.sorted.csv")

# remove all useless variables

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 5.3 match phenotype, genotype and map, subset 1 marker per position 
######################################################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# load data

# phenotype

pheno <- read.csv("./data/pheno/Adj_means.csv")
row.names <- as.character(pheno[,1])
rownames(pheno) <- row.names
pheno <- pheno[,-1]

# genotype

geno <- read.csv("./data/geno/geno.short.sorted.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]

# map

map <- read.table("./data/map/map_panzea_full.txt",h=T)
mk.names <- as.character(map[,1])
mk.names <- gsub("-",".",mk.names)
map[,1] <- mk.names


# match genotype and map

match <- match.marker(mk.mat = geno, map = map)
geno <- match$new.mk.mat
map <- match$new.map
rm(match)


# match genotype and phenotype

geno.par <- geno[1:6,] # separate parents from offspring

match <- match.genotype(mk.mat = geno, pheno = pheno)
pheno <- match$new.pheno
geno <-  match$new.mk.mat
rm(match)

# subset single unique position

difference <- diff(map[,3])

# add a 1 for the first position

difference <- c(1,difference)

map.short.sg.pos <- map[-which(difference==0),]

match <- match.marker(mk.mat = geno.par, map = map.short.sg.pos)
geno.short.par.sg.pos <- match$new.mk.mat

match <- match.marker(mk.mat = geno, map = map.short.sg.pos)
geno.short.off.sg.pos <- match$new.mk.mat
rm(match)


# save data

write.table(map.short.sg.pos,"./data/map/map.short.sg.pos.txt")

write.csv(geno.short.par.sg.pos,"./data/geno/geno.short.par.sg.pos.csv")
write.csv(geno.short.off.sg.pos,"./data/geno/geno.short.off.sg.pos.csv")

write.csv(pheno,"./data/pheno/pheno.short.sorted.csv")


variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 5.4 form genotype matrix for ABH  models 
##########################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# load data

# genotypes

geno.off <- read.csv("./data/geno/geno.short.off.sg.pos.csv",
                     stringsAsFactors = FALSE,row.names = 1,)
geno.off <- as.matrix(geno.off)


geno.par <- read.csv("./data/geno/geno.short.par.sg.pos.csv",
                     stringsAsFactors = FALSE, row.names=1)
geno.par <- as.matrix(geno.par)


# make a cross indicator

cross.ind <- substring(rownames(geno.off),1,5)
cross.id <- levels(as.factor(cross.ind))

# make a parent per cross indicator

parent1 <- as.character(rep("F353",5))
parent2 <- rownames(geno.par)
parent2 <- parent2[-which(parent2=="F353")]

par.per.cross <- cbind(as.character(cross.id), parent1, parent2)

# ABH assignement

geno.ABH <- cross.ABH(par.sc = geno.par,off.sc = geno.off,
                      cross.ind = cross.ind, par.per.cross=par.per.cross)

# save data

write.csv(geno.ABH,"./data/geno/geno.short.ABH.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

# 5.5 formation of the mpp.data object (ABH, biall)
###################################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Load data

# Genotypes

geno.ABH <- read.csv("./data/geno/geno.short.ABH.csv", row.names = 1,
                     stringsAsFactors = FALSE)

geno.biall <- read.csv("./data/geno/geno.short.off.sg.pos.csv", row.names = 1,
                       stringsAsFactors = FALSE)

geno.par <- read.csv("./data/geno/geno.short.par.sg.pos.csv", row.names = 1,
                     stringsAsFactors = FALSE)

# cross.ind

cross.ind <- substring(rownames(geno.ABH),1,5)
cross.id <- levels(as.factor(cross.ind))

# make a parent per cross indicator

parent1 <- as.character(rep("F353",5))
parent2 <- rownames(geno.par)
parent2 <- parent2[-which(parent2=="F353")]


par.per.cross <- cbind(as.character(cross.id), parent1, parent2)

# phenotype

pheno <- read.csv("./data/pheno/pheno.short.sorted.csv",  row.names = 1)

trait.DMY <- data.frame(rownames(pheno), pheno[,1],stringsAsFactors = F)

# map

map <- read.table("./data/map/map.short.sg.pos.txt")

# 5.5.1 mpp.data (par, anc)
#########################

data.short.ABH <- mpp.data(geno = geno.ABH,geno.par = geno.par,biall = FALSE,
                           type = "dh",map = map,trait = trait.DMY,
                           cross.ind = cross.ind,par.per.cross = par.per.cross,
                           step = 50,dir = "./data/geno")

saveRDS(data.short.ABH,file="./data/mpp_data/data.short.ABH.rds")


# 5.5.2 mpp.data (biallelic)
##########################

data.short.biall <- mpp.data(geno = geno.biall,geno.par = geno.par,biall = TRUE,
                             type = "dh",map = map,trait = trait.DMY,
                             cross.ind = cross.ind,par.per.cross = par.per.cross)

saveRDS(data.short.biall,file="./data/mpp_data/data.short.biall.rds")


variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 5.6 Descriptive statistics : mean, variance
##############################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# phenotype

pheno <- read.csv("./data/pheno/pheno.short.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

cross.ind <- substring(rownames(pheno),1,5)
cross.id <- levels(as.factor(cross.ind))

# boxplot

boxplot(pheno$DMY~cross.ind, main="General boxplot (DMY)",
        names= cross.id)

boxplot(pheno$PH~cross.ind, main="General boxplot (PH)",
        names= cross.id)

# mean and variance

variances.DMY <- c()
mean.DMY <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno$DMY[cross.ind==cross.id[i]])
  
  variances.DMY <- c(variances.DMY,var.i)
  
  mean.i <- mean(pheno$DMY[cross.ind==cross.id[i]])
  mean.DMY <- c(mean.DMY,mean.i)
  
}

cbind(mean.DMY,sqrt(variances.DMY))

variances.PH <- c()
mean.PH <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno$PH[cross.ind==cross.id[i]])
  variances.PH <- c(variances.PH,var.i)
  
  mean.i <- mean(pheno$PH[cross.ind==cross.id[i]])
  mean.PH <- c(mean.PH,mean.i)
  
}

cbind(mean.PH,sqrt(variances.PH))


variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

################################################################################
############################## End short subset ################################
################################################################################

##################
# 6. Long subset #
##################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv")
rownames(pheno) <- as.character(pheno[,1])
pheno <- pheno[,-1]

cross.ind <- substr(rownames(pheno),1,5)

# cross partition

cross.partition <- read.table("./data/geno/cross.partition.txt")


# index for selection with parent and genotype belonging to the selected cross

index.sel <- c(rownames(geno)[1:11] %in% c("F353",as.character(cross.partition[,5])),
               cross.ind %in% cross.partition[,2])

geno.long <- geno[index.sel,]


# 6.1 Remove monomorphic and missing value positions
####################################################

par.scores <- allele.scores(mk.mat = geno.long[1:6,])

par.mono <- which(par.scores[1,]=="mono")
par.miss <- which(par.scores[1,]=="miss")

prob.ind <- c(par.mono,par.miss)

geno.long <- geno.long[,-prob.ind]

# detect monomorphic position in the offspring and remove

off.scores <- allele.scores(mk.mat = geno.long[7:dim(geno.long)[1],])

off.mono <- which(off.scores[1,]=="mono")

geno.long <- geno.long[,-off.mono]

# rare alleles with MAF <0.01

rare.all <- rare.allele.ind(mk.mat = geno.long[7:dim(geno.long)[1],],
                            threshold = 0.01)

# no marker with MAF<0.01

# missing markers with more than 10%

miss.ind.mk <- missing.indicator(mk.mat = geno.long[7:dim(geno.long)[1],],
                                 threshold=0.1,col=2)

geno.long  <- geno.long[,-miss.ind.mk[,2]]

# genotypes with more than 25% missing values

miss.ind.gen <- missing.indicator(mk.mat = geno.long[7:dim(geno.long)[1],],
                                  threshold=0.25,col=1)

# no line with more than 25% missing values

# save the variable

write.csv(geno.long,"./data/geno/geno.long.sorted.csv")

# remove all useless variables

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)



# 6.2. match phenotype, genotype and map, subset 1 marker per position
######################################################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# load data

# phenotype

pheno <- read.csv("./data/pheno/Adj_means.csv")
row.names <- as.character(pheno[,1])
rownames(pheno) <- row.names
pheno <- pheno[,-1]

# genotype

geno <- read.csv("./data/geno/geno.long.sorted.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]

# map

map <- read.table("./data/map/map_panzea_full.txt",h=T)
mk.names <- as.character(map[,1])
mk.names <- gsub("-",".",mk.names)
map[,1] <- mk.names


# match genotype and map

match <- match.marker(mk.mat = geno, map = map)
geno <- match$new.mk.mat
map <- match$new.map
rm(match)


# match genotype and phenotype

# separate parents from offspring

geno.par <- geno[1:6,]

match <- match.genotype(mk.mat = geno, pheno = pheno)
pheno <- match$new.pheno
geno <-  match$new.mk.mat
rm(match)

# subset single unique position

difference <- diff(map[,3])

# add a 1 for the first position

difference <- c(1,difference)

map.long.sg.pos <- map[-which(difference==0),]

match <- match.marker(mk.mat = geno.par, map = map.long.sg.pos)
geno.long.par.sg.pos <- match$new.mk.mat

match <- match.marker(mk.mat = geno, map = map.long.sg.pos)
geno.long.off.sg.pos <- match$new.mk.mat
rm(match)


# save data

write.table(map.long.sg.pos,"./data/map/map.long.sg.pos.txt")

write.csv(geno.long.par.sg.pos,"./data/geno/geno.long.par.sg.pos.csv")
write.csv(geno.long.off.sg.pos,"./data/geno/geno.long.off.sg.pos.csv")

write.csv(pheno,"./data/pheno/pheno.long.sorted.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 6.3 form genotype matrix for ABH  models 
###########################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)


# load data

# genotypes

geno.off <- read.csv("./data/geno/geno.long.off.sg.pos.csv")
rownames(geno.off) <- geno.off[,1]
geno.off <- geno.off[,-1]

geno.par <- read.csv("./data/geno/geno.long.par.sg.pos.csv")
rownames(geno.par) <- geno.par[,1]
geno.par <- geno.par[,-1]

# make a cross indicator

cross.ind <- substring(rownames(geno.off),1,5)
cross.id <- levels(as.factor(cross.ind))

# make a parent per cross indicator

parent1 <- as.character(rep("F353",5))
parent2 <- rownames(geno.par)
parent2 <- parent2[-which(parent2=="F353")]

par.per.cross <- cbind(as.character(cross.id), parent1, parent2)

# ABH assignement

geno.ABH <- cross.ABH(par.sc = geno.par,off.sc = geno.off,
                      cross.ind = cross.ind, par.per.cross=par.per.cross)

# save data

write.csv(geno.ABH,"./data/geno/geno.long.ABH.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 6.4 formation of the mpp.data object (ABH, biall)
###################################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Load data

# Genotypes

geno.ABH <- read.csv("./data/geno/geno.long.ABH.csv")
rownames(geno.ABH) <- geno.ABH[,1]
geno.ABH <- geno.ABH[,-1]

geno.biall <- read.csv("./data/geno/geno.long.off.sg.pos.csv")
rownames(geno.biall) <- geno.biall[,1]
geno.biall <- geno.biall[,-1]

geno.par <- read.csv("./data/geno/geno.long.par.sg.pos.csv")
rownames(geno.par) <- geno.par[,1]
geno.par <- geno.par[,-1]

# cross.ind

cross.ind <- substring(rownames(geno.ABH),1,5)
cross.id <- levels(as.factor(cross.ind))

# make a parent per cross indicator

parent1 <- as.character(rep("F353",5))
parent2 <- rownames(geno.par)
parent2 <- parent2[-which(parent2=="F353")]


par.per.cross <- cbind(as.character(cross.id), parent1, parent2)

# phenotype

pheno <- read.csv("./data/pheno/pheno.long.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

trait.DMY <- data.frame(rownames(pheno), pheno[,1],stringsAsFactors = F)

# map

map <- read.table("./data/map/map.long.sg.pos.txt")

# 6.4.1 mpp.data (par, anc)
#########################

data.long.ABH <- mpp.data(geno = geno.ABH,geno.par = geno.par,biall = FALSE,
                          type = "dh",map = map,trait = trait.DMY,
                          cross.ind = cross.ind,par.per.cross = par.per.cross,
                          step = 50,dir = "./data/geno")

saveRDS(data.long.ABH,file="./data/mpp_data/data.long.ABH.rds")

# 6.4.2 mpp.data (biallelic)
##########################

data.long.biall <- mpp.data(geno = geno.biall,geno.par = geno.par,biall = TRUE,
                            type = "dh",map = map,trait = trait.DMY,
                            cross.ind = cross.ind,par.per.cross = par.per.cross)

saveRDS(data.long.biall,file="./data/mpp_data/data.long.biall.rds")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

# 6.5 Descriptive statistics : mean, variance
#############################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# phenotype

pheno <- read.csv("./data/pheno/pheno.long.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

cross.ind <- substring(rownames(pheno),1,5)
cross.id <- levels(as.factor(cross.ind))

# boxplot

boxplot(pheno$DMY~cross.ind, main="General boxplot (DMY)",
        names= cross.id)

boxplot(pheno$PH~cross.ind, main="General boxplot (PH)",
        names= cross.id)

# mean and variance

variances.DMY <- c()
mean.DMY <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno$DMY[cross.ind==cross.id[i]])
  
  variances.DMY <- c(variances.DMY,var.i)
  
  mean.i <- mean(pheno$DMY[cross.ind==cross.id[i]])
  mean.DMY <- c(mean.DMY,mean.i)
  
}

cbind(mean.DMY,sqrt(variances.DMY))

variances.PH <- c()
mean.PH <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno$PH[cross.ind==cross.id[i]])
  variances.PH <- c(variances.PH,var.i)
  
  mean.i <- mean(pheno$PH[cross.ind==cross.id[i]])
  mean.PH <- c(mean.PH,mean.i)
  
}

cbind(mean.PH,sqrt(variances.PH))

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

################################################################################
############################## End long subset ################################
################################################################################

####################
# 7. hetero subset #
####################

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno <- read.csv("./data/geno/geno_panzea_full.csv",row.names = 1)

# phenotypic data

pheno <- read.csv("./data/pheno/Adj_means.csv",row.names = 1)

cross.ind <- substr(rownames(pheno),1,5)

# cross partition

cross.partition <- read.table("./data/geno/cross.partition.txt")

# index for selection with parent and genotype belonging to the selected cross

index.sel <- c(rownames(geno)[1:11] %in% c("F353",as.character(cross.partition[,6])),
               cross.ind %in% cross.partition[,3])

geno.hetero <- geno[index.sel,]

# 7.1 subset line to equalize with long subset
##############################################

# take 361 lines to have same population size with
# the long subset

sel.ind <- sort(sample(7:434,361))

geno.hetero <- rbind(geno.hetero[1:6,], geno.hetero[sel.ind,])

# 7.2 Remove monomorphic and missing value positions
####################################################

par.scores <- allele.scores(mk.mat = geno.hetero[1:6,])

par.mono <- which(par.scores[1,]=="mono")
par.miss <- which(par.scores[1,]=="miss")

prob.ind <- c(par.mono,par.miss)

geno.hetero <- geno.hetero[,-prob.ind]

off.scores <- allele.scores(mk.mat = geno.hetero[7:dim(geno.hetero)[1],])

off.mono <- which(off.scores[1,]=="mono")

geno.hetero <- geno.hetero[,-off.mono]

# rare alleles with MAF <0.01

rare.all <- rare.allele.ind(mk.mat = geno.hetero[7:dim(geno.hetero)[1],],
                            threshold = 0.01)

geno.hetero <- geno.hetero[,-rare.all[,2]]

# missing markers with more than 10%

miss.ind.mk <- missing.indicator(mk.mat = geno.hetero[7:dim(geno.hetero)[1],],
                                 threshold=0.1,col=2)

geno.hetero  <- geno.hetero[,-miss.ind.mk[,2]]

# genotypes with more than 25% missing values

miss.ind.gen <- missing.indicator(mk.mat = geno.hetero[7:dim(geno.hetero)[1],],
                                  threshold=0.25,col=1)

# no line with more than 25% missing values

# save the variable

write.csv(geno.hetero,"./data/geno/geno.hetero.sorted.csv")

# remove all useless variables

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 7.3 match phenotype, genotype and map, subset 1 marker per position
#####################################################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# load data

# phenotype

pheno <- read.csv("./data/pheno/Adj_means.csv",row.names = 1)

# genotype

geno <- read.csv("./data/geno/geno.hetero.sorted.csv",row.names = 1)

# map

map <- read.table("./data/map/map_panzea_full.txt",h=T)
mk.names <- as.character(map[,1])
mk.names <- gsub("-",".",mk.names)
map[,1] <- mk.names

# match genotype and map

match <- match.marker(mk.mat = geno, map = map)
geno <- match$new.mk.mat
map <- match$new.map
rm(match)


# match genotype and phenotype

# separate parents from offspring

geno.par <- geno[1:6,]

match <- match.genotype(mk.mat = geno, pheno = pheno)
pheno <- match$new.pheno
geno <-  match$new.mk.mat
rm(match)

# subset single unique position

difference <- diff(map[,3])

# add a 1 for the first position

difference <- c(1,difference)

map.hetero.sg.pos <- map[-which(difference==0),]

match <- match.marker(mk.mat = geno.par, map = map.hetero.sg.pos)
geno.hetero.par.sg.pos <- match$new.mk.mat

match <- match.marker(mk.mat = geno, map = map.hetero.sg.pos)
geno.hetero.off.sg.pos <- match$new.mk.mat
rm(match)


# save data

write.table(map.hetero.sg.pos,"./data/map/map.hetero.sg.pos.txt")

write.csv(geno.hetero.par.sg.pos,"./data/geno/geno.hetero.par.sg.pos.csv")
write.csv(geno.hetero.off.sg.pos,"./data/geno/geno.hetero.off.sg.pos.csv")

write.csv(pheno,"./data/pheno/pheno.hetero.sorted.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 7.4 form genotype matrix for ABH  models
##########################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)


# load data

# genotypes

geno.off <- read.csv("./data/geno/geno.hetero.off.sg.pos.csv", row.names = 1,
                     stringsAsFactors = FALSE)
geno.off <- as.matrix(geno.off)


geno.par <- read.csv("./data/geno/geno.hetero.par.sg.pos.csv", row.names = 1,
                     stringsAsFactors = FALSE)
geno.par <- as.matrix(geno.par)

# make a cross indicator

cross.ind <- substring(rownames(geno.off),1,5)
cross.id <- levels(as.factor(cross.ind))

# make a parent per cross indicator

parent1 <- as.character(rep("F353",5))
parent2 <- rownames(geno.par)
parent2 <- parent2[-which(parent2=="F353")]

par.per.cross <- cbind(as.character(cross.id), parent1, parent2)

#  ABH assignement

geno.ABH <- cross.ABH(par.sc = geno.par,off.sc = geno.off,
                      cross.ind = cross.ind, par.per.cross=par.per.cross)

# save data

write.csv(geno.ABH,"./data/geno/geno.hetero.ABH.csv")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


# 7.5 formation of the mpp.data object (ABH, biall)
######################################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Load data

# Genotypes

geno.ABH <- read.csv("./data/geno/geno.hetero.ABH.csv",row.names = 1)

geno.biall <- read.csv("./data/geno/geno.hetero.off.sg.pos.csv",row.names = 1)

geno.par <- read.csv("./data/geno/geno.hetero.par.sg.pos.csv",row.names = 1)


# cross.ind

cross.ind <- substring(rownames(geno.ABH),1,5)
cross.id <- levels(as.factor(cross.ind))

# make a parent per cross indicator

parent1 <- as.character(rep("F353",5))
parent2 <- rownames(geno.par)
parent2 <- parent2[-which(parent2=="F353")]

par.per.cross <- cbind(as.character(cross.id), parent1, parent2)

# phenotype

pheno <- read.csv("./data/pheno/pheno.hetero.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

trait.DMY <- data.frame(rownames(pheno), pheno[,1],stringsAsFactors = F)

# map

map <- read.table("./data/map/map.hetero.sg.pos.txt")

# 7.4.1 mpp.data (par, anc)
#########################

data.hetero.ABH <- mpp.data(geno = geno.ABH,geno.par = geno.par,biall = FALSE,
                            type = "dh",map = map,trait = trait.DMY,
                            cross.ind = cross.ind,par.per.cross = par.per.cross,
                            step = 50,dir = "./data/geno")

saveRDS(data.hetero.ABH,file="./data/mpp_data/data.hetero.ABH.rds")


# 7.4.2 mpp.data (biallelic)
##########################

data.hetero.biall <- mpp.data(geno = geno.biall,geno.par = geno.par,biall = TRUE,
                              type = "dh",map = map,trait = trait.DMY,
                              cross.ind = cross.ind,par.per.cross = par.per.cross)

saveRDS(data.hetero.biall,file="./data/mpp_data/data.hetero.biall.rds")

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

# 7.5 Descriptive statistics : mean, variance
#############################################

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Load data

# phenotype

pheno <- read.csv("./data/pheno/pheno.hetero.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

cross.ind <- substring(rownames(pheno),1,5)
cross.id <- levels(as.factor(cross.ind))


# boxplot

boxplot(pheno$DMY~cross.ind, main="General boxplot (DMY)",
        names= cross.id)

boxplot(pheno$PH~cross.ind, main="General boxplot (PH)",
        names= cross.id)

# mean and variance

variances.DMY <- c()
mean.DMY <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno$DMY[cross.ind==cross.id[i]])
  
  variances.DMY <- c(variances.DMY,var.i)
  
  mean.i <- mean(pheno$DMY[cross.ind==cross.id[i]])
  mean.DMY <- c(mean.DMY,mean.i)
  
}

cbind(mean.DMY,sqrt(variances.DMY))

variances.PH <- c()
mean.PH <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno$PH[cross.ind==cross.id[i]])
  variances.PH <- c(variances.PH,var.i)
  
  mean.i <- mean(pheno$PH[cross.ind==cross.id[i]])
  mean.PH <- c(mean.PH,mean.i)
  
}

cbind(mean.PH,sqrt(variances.PH))

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

# determination of the number of marker in common between the three different
# subsets

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

map.short <- read.table("./data/map/map.short.sg.pos.txt")
map.hetero <- read.table("./data/map/map.hetero.sg.pos.txt")
map.long <- read.table("./data/map/map.long.sg.pos.txt")

short_N_hetero <- intersect(map.short[, 1], map.hetero[, 1])
short_N_hetero_N_long <- intersect(short_N_hetero, map.long[, 1])

################################################################################
############################## End hetero subset ###############################
################################################################################

#########################
# 8. Clusthaplo subsets #
#########################

# This part of the script must be run on a version of R where function
# clusthaplo and RHmm in order to use HMM for the clustering process.
# R 2.14 is an option

library(clusthaplo)
library(RHmm)
library(mppRDraft)

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# cross partition

cross.partition <- read.table("./data/geno/cross.partition.txt")

par.short <- c("F353", as.character(cross.partition[,4]))
par.long <- c("F353", as.character(cross.partition[,5]))
par.hetero <- c("F353", as.character(cross.partition[,6]))

# map with all markers (22K)

map.full <- read.table("./data/map/map_panzea_full.txt")

# genotype data corresponding to these markers for the parental lines

geno.par <- read.csv("./data/geno/geno_panzea_par.csv",row.names = 1)
geno.par <- t(geno.par)

# subset geno par (short, long, hetero)

geno.par.short <- geno.par[,colnames(geno.par) %in% par.short]
geno.par.long <- geno.par[,colnames(geno.par) %in% par.long]
geno.par.hetero <- geno.par[,colnames(geno.par) %in% par.hetero]

# map from the IBD computation with only single position

map.short.sg.pos <- read.table("./data/map/map.short.sg.pos.txt")
map.hetero.sg.pos <- read.table("./data/map/map.hetero.sg.pos.txt")
map.long.sg.pos <- read.table("./data/map/map.long.sg.pos.txt")


# 8.1 clustering process (take around 5 min)
#########################


par.clu.short <- parent.cluster(haplo.map=map.full, consensus.map=map.short.sg.pos,
                                marker.data=geno.par.short, na.strings="-",
                                step.size=50, window=2, K=10,
                                clustering.method="hmm")

par.clu.long <- parent.cluster(haplo.map=map.full, consensus.map=map.long.sg.pos,
                               marker.data=geno.par.long, na.strings="-",
                               step.size=50, window=2, K=10,
                               clustering.method="hmm")

par.clu.hetero <- parent.cluster(haplo.map=map.full, consensus.map=map.hetero.sg.pos,
                                 marker.data=geno.par.hetero, na.strings="-",
                                 step.size=50, window=2, K=10,
                                 clustering.method="hmm")


# average number of cluster

par.clu.short$av.cl

par.clu.hetero$av.cl

par.clu.long$av.cl


# save clustering object

par.clu.short <- par.clu.short[[1]]
colnames(par.clu.short) <- par.short

par.clu.long <- par.clu.long[[1]]
colnames(par.clu.long) <- par.long

par.clu.hetero <- par.clu.hetero[[1]]
colnames(par.clu.hetero) <- par.hetero

# the results of the clustering process using clusthaplo were not
# the same when we repeated the analysis. If you want to get the 
# same results use the clustering results already saved in
# ~MPP_EUNAM/data/clustering

# write.table(par.clu.short,"./data/clustering/par_clu_short.txt")
# write.table(par.clu.long,"./data/clustering/par_clu_long.txt")
# write.table(par.clu.hetero,"./data/clustering/par_clu_hetero.txt")

# end of the part that must be done in R version 2.14

# removed unused variable

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)


###############################################################################
############################ End parental line clustering #####################
###############################################################################

#########################################################
# 9. Significance threshold determination (permutation) #
#########################################################

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


partition.perm <- expand.grid(c("short","hetero","long"),
                              c("par","anc","biall"),c("DMY","PH"))

partition.perm  <- partition.perm[order(partition.perm[,1]),]
partition.perm <- data.frame(partition.perm, rep(c("ABH","ABH","biall"),6))
colnames(partition.perm) <- c("subset", "Q.eff","trait","ABH.bi")

N <- 2000
q.val <- c(0.9,0.95)

for(i in 1:18){
  
  par.exc <- partition.perm[i,]
  par.exc <- as.character(unlist(par.exc))
  
  file.name <- paste("./data/mpp_data/data.",par.exc[1],".",
                     par.exc[4],".rds",sep="")
  
  # load the mpp.data object
  
  data <- readRDS(file = file.name)
  
  # put the right trait
  
  file.name <- paste0("./data/pheno/pheno.",par.exc[1],".sorted.csv")
  
  pheno <- read.csv(file.name)
  
  if(par.exc[3]=="DMY"){
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$DMY,stringsAsFactors = F)
    
  } else {
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$PH,stringsAsFactors = F)
    
  }
  
  
  data <- sub.trait(data,trait)
  
  # load the par.clu object
  
  if(par.exc[2]=="anc"){
    
    file.name <- paste("./data/clustering/par_clu_", par.exc[1],".txt",sep="")
    
    par.clu <- read.table(file.name)
    
  }
  
  # If you want/can run the code on parallel in your computer
  # you can uncomment an run the following code. Do not forget
  # to remove the part of the code that compute the permutation
  # test without using parallel functionalities.
  
  # perm.test <- library(parallel)
  # n.cores <- detectCores()
  # cluster <- makeCluster((n.cores-1))
  # 
  # # compute the permutation test
  # 
  # perm.thre <- mpp.perm(data = data, Q.eff = par.exc[2],VCOV = "u.err",
  #                       N = N, q.val = q.val, parallel = TRUE,
  #                       par.clu=par.clu, cluster = cluster)
  # 
  # 
  # stopCluster(cl = cluster)
  

  # save the results

# Origninal threshold result of the study are already located in the following
# folder ~/MPP_EUNAM/results/Threshold_permutation. The permutation results
# are therefore not saved otherwise they will overwrite the original results
# from the study. If you want nevertheless to save these results uncomment the
# next lines
  
  file.name <- paste(par.exc[1],"_",par.exc[2],"_",par.exc[3],sep="")
  file.name <- paste(file.name,".txt",sep="")

  folder <- "./results/Threshold_permutation/"
  folder <- paste(folder,file.name,sep="")

  write.table(cbind(perm.thre$max.pval,perm.thre$seed),folder)
  
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

library(mppRDraft)

library(qtl)
library(asreml)
library(MASS)
library(stringr)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)

# need to load version 0.7.1 of igraph library
# https://cran.r-project.org/web/packages/igraph/index.html

library(igraph)

# load data

# SPECIFY YOUR PATH HERE (the location of ~/MPP_EUNAM)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# partition

subset <- rep(c("short","hetero","long"),each=12)
Q.eff <- rep(c("par","anc","biall"),times = 6,each = 2)
trait <- rep(c("DMY","PH"),times = 3,each = 6)
VCOV <- rep(c("u.err","cr.err"),times = 18)
data.type <- rep(c("ABH","ABH","biall"), times=6, each=2)

part.QTL <- data.frame(subset,Q.eff,trait,VCOV,data.type)

folder <- "./results/Threshold_permutation/"
my.loc <- "./results/QTL_analyses"


for(i in 1:36){
  
  # extract the ith partition
  part <- as.character(unlist(part.QTL[i,]))
  
  # get the threshold value
  
  max.val <- read.table(paste0(folder, paste0(part[1],"_",part[2],"_",
                                              part[3],".txt")))
  
  thre <- quantile(max.val[,1],0.95)
  
  data <- readRDS(file = paste0("./data/mpp_data/data.",part[1],".",part[5],".rds"))
  
  # put the right trait
  
  pheno <- read.csv(paste0("./data/pheno/pheno.",part[1],".sorted.csv"))
  
  if(part[3]=="DMY"){
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$DMY,stringsAsFactors = F)
    
  } else {
    
    trait <-  data.frame(as.character(pheno[,1]), pheno$PH,stringsAsFactors = F)
    
  }
  
  data <- sub.trait(data,trait)
  
  # load the par.clu object
  
  if(part[2]=="anc"){
    
    par.clu <- read.table(paste0("./data/clustering/par_clu_", part[1],".txt"))
    
  }
  
  # QTL analysis results are already saved in the following folder
  # ~/MPP_EUNAM/results/QTL_analyses. If however you want to reproduce the
  # results, uncomment the part bellow. Be carefull that the code is programmed
  # in parallel! If you do not want to run it in parallel modify option parllel in
  # mpp.proc. Be also carefull with the construction of cluster
  # (cluster <- makeCluster((n.cores-1)); stopCluster(cl=cluser)).
  
  
  # if(part[4]=="u.err") {
  # 
  #   parallel=TRUE
  #   library(parallel)
  #   n.cores <- detectCores()
  #   cluster <- makeCluster((n.cores-1))
  # 
  # } else {
  # 
  #   parallel=FALSE
  # 
  # }
  # 
  # proc <- mpp.proc(pop.name = part[1], trait.name = part[3],
  #                   parallel = parallel, cluster = cluster ,
  #                   data = data, Q.eff = part[2],par.clu = par.clu,
  #                   VCOV = part[4], est.gen.eff = TRUE,thre.cof = thre,
  #                   N.cim = 2, effect.comp = TRUE, thre.QTL = thre,
  #                   CI = TRUE, output.loc = my.loc)
  # 
  # if(parallel==TRUE){stopCluster(cl = cluster)}
  
}

variable <- ls()
variable <- variable[-which(variable == "path")]

rm(list=variable)
rm(variable)

###############################################################################
################################ End QTL analysis #############################
###############################################################################

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

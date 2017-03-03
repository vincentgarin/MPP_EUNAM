##############################
# table and figures Appendix #
##############################

# S2 Kinship matrix GS
######################


library(mppRDraft)
library(xtable)

path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv")
rownames(geno.par) <- geno.par[,1]
geno.par <- geno.par[,-1]

# SM coefficient distance between the central parent (F353)

kin.mat.par <- kinship.matrix(mk.mat = geno.par, method = "SM")

kin.mat.par <- round(kin.mat.par,3)
kin.mat.par2 <- kin.mat.par

kin.mat.par2[lower.tri(kin.mat.par2)] <- ""

kin.mat.par2

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

print(xtable(kin.mat.par2,digits = c(0,3,3,3,3,3,3,3,3,3,3,3),
             caption = "SM table", table.placement="!ht",
             caption.placement="top",
             align = c("c","|","c","c","c","c","c","c","c","c","c","c","c"),),
      include.rownames=TRUE)

write(x = print(xtable(kin.mat.par2,digits = c(0,3,3,3,3,3,3,3,3,3,3,3),
                       caption = "SM table", table.placement="!ht",
                       caption.placement="top",
                       align = c("c","|","c","c","c","c","c","c","c","c","c","c","c"),),
                include.rownames=TRUE),
      file=paste(folder,"SM_table.txt",sep = ""))

################################################################################
######################### End Kinship matrix GS ################################
################################################################################

# S4 Genetic map subsets
########################

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/"

path <- "~/MPP_EUNAM/data"

setwd(path)

# Load data

# Genotypes

geno <- read.csv("./geno/geno.short.ABH.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]


# phenotype

pheno <- read.csv("./pheno/pheno.short.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

trait <- data.frame(rownames(pheno), pheno[,1],stringsAsFactors = F)

# map

map <- read.table("./map/map.short.sg.pos.txt")

dir <-  "E:/PhD/Test" 


# genotype matrix

# form the chromosome and the position (in cM) information

chr.info <- t(map[,2:3])
colnames(chr.info) <- map[,1]

# add the map information to the marker information

geno <- rbind(chr.info,geno)

# addition of the phenotype column
# first extend the pheno column with two empty

trait <- c(c("",""),trait[,2])

geno <- cbind(trait, geno)

# Export the data in a .csv file in the specified directory

file.name <- paste(dir,"/","Cross_object.csv",sep="")

write.csv(geno,file=file.name,row.names=F)

# read the data to form after the cross object

# select what type of cross we want to read

# Form the cross according to the type of cross specified by the user

cross.object  <-  read.cross("csv", file = file.name, genotypes=c("A","B"),
                             alleles=c("A","B"))

# map projection

map.sum.short <- summary.map(cross.object)
# map.sum.short <- map.sum.short[,-1]
map.sum.short <- data.frame(c(1:10,"Overal"),map.sum.short)
colnames(map.sum.short) <- c("Chromosome","N","Length(cM)","Average spacing",
                             "maximum spacing")

map.sum.short.Ltx <- print(xtable(map.sum.short,digits = c(0,0,0,1,1,1)),
                           include.rownames=FALSE)

write(x = print(xtable(map.sum.short,digits = c(0,0,0,1,1,1)),
                include.rownames=FALSE),
      file=paste(folder,"table/","map_sum_short.txt",sep = ""))

jpeg(filename = paste(folder,"figures/","map_short.jpeg",sep = ""),
     width = 1000, height = 1000, units = "px", pointsize = 36,
     quality = 100)
plot.map(cross.object,main="EU-NAM Dent short subset map",)
dev.off()


########
# long #
########

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/"

path <- "~/MPP_EUNAM/data"

setwd(path)

# Load data

# Genotypes

geno <- read.csv("./geno/geno.long.ABH.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]


# phenotype

pheno <- read.csv("./pheno/pheno.long.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

trait <- data.frame(rownames(pheno), pheno[,1],stringsAsFactors = F)

# map

map <- read.table("./map/map.long.sg.pos.txt")

dir <-  "E:/PhD/Test" 


# genotype matrix

# form the chromosome and the position (in cM) information

chr.info <- t(map[,2:3])
colnames(chr.info) <- map[,1]

# add the map information to the marker information

geno <- rbind(chr.info,geno)

# addition of the phenotype column
# first extend the pheno column with two empty

trait <- c(c("",""),trait[,2])

geno <- cbind(trait, geno)

# Export the data in a .csv file in the specified directory

file.name <- paste(dir,"/","Cross_object.csv",sep="")

write.csv(geno,file=file.name,row.names=F)

# read the data to form after the cross object

# select what type of cross we want to read

# Form the cross according to the type of cross specified by the user

cross.object  <-  read.cross("csv", file = file.name, genotypes=c("A","B"),
                             alleles=c("A","B"))

# map projection

map.sum.long <- summary.map(cross.object)
map.sum.long <- data.frame(c(1:10,"Overal"),map.sum.long)
colnames(map.sum.long) <- c("Chromosome","N","Length(cM)","Average spacing",
                            "maximum spacing")

map.sum.short.Ltx <- print(xtable(map.sum.long,digits = c(0,0,0,1,1,1)),
                           include.rownames=FALSE)

write(x = print(xtable(map.sum.long,digits = c(0,0,0,1,1,1)),
                include.rownames=FALSE),
      file=paste(folder,"table/","map_sum_long.txt",sep = ""))

jpeg(filename = paste(folder,"figures/","map_long.jpeg",sep = ""),
     width = 1000, height = 1000, units = "px", pointsize = 36,
     quality = 100)
plot.map(cross.object,main="EU-NAM Dent long subset map",)

dev.off()


##########
# hetero #
##########

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/"

path <- "~/MPP_EUNAM/data"

setwd(path)

# Load data

# Genotypes

geno <- read.csv("./geno/geno.hetero.ABH.csv")
rownames(geno) <- geno[,1]
geno <- geno[,-1]


# phenotype

pheno <- read.csv("./pheno/pheno.hetero.sorted.csv")
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

trait <- data.frame(rownames(pheno), pheno[,1],stringsAsFactors = F)

# map

map <- read.table("./map/map.hetero.sg.pos.txt")

dir <-  "E:/PhD/Test" 


# genotype matrix

# form the chromosome and the position (in cM) information

chr.info <- t(map[,2:3])
colnames(chr.info) <- map[,1]

# add the map information to the marker information

geno <- rbind(chr.info,geno)

# addition of the phenotype column
# first extend the pheno column with two empty

trait <- c(c("",""),trait[,2])

geno <- cbind(trait, geno)

# Export the data in a .csv file in the specified directory

file.name <- paste(dir,"/","Cross_object.csv",sep="")

write.csv(geno,file=file.name,row.names=F)

# read the data to form after the cross object

# select what type of cross we want to read

# Form the cross according to the type of cross specified by the user

cross.object  <-  read.cross("csv", file = file.name, genotypes=c("A","B"),
                             alleles=c("A","B"))

# map projection

map.sum.hetero <- summary.map(cross.object)
map.sum.hetero <- data.frame(c(1:10,"Overal"),map.sum.hetero)
colnames(map.sum.hetero) <- c("Chromosome","N","Length(cM)","Average spacing",
                              "maximum spacing")

map.sum.short.Ltx <- print(xtable(map.sum.hetero,digits = c(0,0,0,1,1,1)),
                           include.rownames=FALSE)

write(x = print(xtable(map.sum.hetero,digits = c(0,0,0,1,1,1)),
                include.rownames=FALSE),
      file=paste(folder,"table/","map_sum_hetero.txt",sep = ""))

jpeg(filename = paste(folder,"figures/","map_hetero.jpeg",sep = ""),
     width = 1000, height = 1000, units = "px", pointsize = 36,
     quality = 100)
plot.map(cross.object,main="EU-NAM Dent heterogeneous subset map",)

dev.off()

################################################################################
########################### End map subsets ####################################
################################################################################


# S5 variance components table 
##############################

library(asreml)

path <- "~/MPP_EUNAM/data"

setwd(path)

pheno <- read.csv("./pheno/pheno_red.csv")
pheno <- pheno[,-1]

lines_used <- read.csv("./pheno/List_lines_Dent_Lehermeier.csv")

# suppress the factor levels with 0 occurences

pheno$Genotype <- as.factor(as.character(pheno$Genotype))
pheno$Fam <- as.factor(as.character(pheno$Fam))

# order data by family

pheno <- pheno[order(pheno$Fam),]

par_names <- c("B73","D06", "D09", "EC169", "F252", "F353", "F618", "F98902",
               "Mo17", "UH250", "UH304","W117")

# DMY
#####

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


var.DMY.Ltx <- print(xtable(var.DMY, align = c("l","c","c","c","c","c"),
                            digits = c(2,2,2,2,2,2)))

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

write(x = print(xtable(var.DMY, align = c("l","c","c","c","c","c"),
                       digits = c(2,2,2,2,2,2))),
      file=paste0(folder,"var.comp.DMY.txt"))

# PH
####

model <- asreml(fixed = PH~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

var.PH <- summary(model)$varcomp
var.PH <-round(cbind(var.PH[c(2:11), c(2:3)], var.PH[c(13:22), c(2:3)],
                     var.PH[c(26:35), 2]), 2)
her <- round(100*(var.PH[, 1]/(var.PH[, 1] + (var.PH[, 3]/4) + (var.PH[, 5]/5))), 2)
var.PH <- cbind(var.PH[, -5], her)
rownames(var.PH) <- par_names[-c(6, 8)]
colnames(var.PH) <- c("sigma.g", "std.err", "sigma.ge", "std.err","heritability")
var.PH

var.PH.Ltx <- print(xtable(var.PH, align = c("l","c","c","c","c","c"),
                           digits = c(2,2,2,2,2,2)))

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

write(x = print(xtable(var.PH, align = c("l","c","c","c","c","c"),
                       digits = c(2,2,2,2,2,2))),
      file=paste0(folder,"var.comp.PH.txt"))

################################################################################
########################### End var comp table #################################
################################################################################

# S6 Adjusted means descriptive statistics
##########################################

# phenotype

path <- "~/MPP_EUNAM/data"

setwd(path)

# short

pheno.sh <- read.csv("./pheno/pheno.short.sorted.csv", row.names = 1)

cross.ind.sh <- substring(rownames(pheno.sh),1,5)
cross.id.sh <- levels(as.factor(cross.ind.sh))

variances.DMY.sh <- c()
mean.DMY.sh <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno.sh$DMY[cross.ind.sh==cross.id.sh[i]])
  
  variances.DMY.sh <- c(variances.DMY.sh,var.i)
  
  mean.i <- mean(pheno.sh$DMY[cross.ind.sh==cross.id.sh[i]])
  mean.DMY.sh <- c(mean.DMY.sh,mean.i)
  
}

variances.PH.sh <- c()
mean.PH.sh <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno.sh$PH[cross.ind.sh==cross.id.sh[i]])
  
  variances.PH.sh <- c(variances.PH.sh,var.i)
  
  mean.i <- mean(pheno.sh$PH[cross.ind.sh==cross.id.sh[i]])
  mean.PH.sh <- c(mean.PH.sh,mean.i)
  
}


# heterogeneous

pheno.het <- read.csv("./pheno/pheno.hetero.sorted.csv", row.names = 1)

cross.ind.het <- substring(rownames(pheno.het),1,5)
cross.id.het <- levels(as.factor(cross.ind.het))

variances.DMY.het <- c()
mean.DMY.het <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno.het$DMY[cross.ind.het==cross.id.het[i]])
  
  variances.DMY.het <- c(variances.DMY.het,var.i)
  
  mean.i <- mean(pheno.het$DMY[cross.ind.het==cross.id.het[i]])
  mean.DMY.het <- c(mean.DMY.het,mean.i)
  
}

variances.PH.het <- c()
mean.PH.het <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno.het$PH[cross.ind.het==cross.id.het[i]])
  
  variances.PH.het <- c(variances.PH.het,var.i)
  
  mean.i <- mean(pheno.het$PH[cross.ind.het==cross.id.het[i]])
  mean.PH.het <- c(mean.PH.het,mean.i)
  
}


# long

pheno.lg <- read.csv("./pheno/pheno.long.sorted.csv", row.names = 1)

cross.ind.lg <- substring(rownames(pheno.lg),1,5)
cross.id.lg <- levels(as.factor(cross.ind.lg))

variances.DMY.lg <- c()
mean.DMY.lg <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno.lg$DMY[cross.ind.lg==cross.id.lg[i]])
  
  variances.DMY.lg <- c(variances.DMY.lg,var.i)
  
  mean.i <- mean(pheno.lg$DMY[cross.ind.lg==cross.id.lg[i]])
  mean.DMY.lg <- c(mean.DMY.lg,mean.i)
  
}

variances.PH.lg <- c()
mean.PH.lg <- c()

for (i in 1:5) {
  
  # subset cross
  
  var.i<- var(pheno.lg$PH[cross.ind.lg==cross.id.lg[i]])
  
  variances.PH.lg <- c(variances.PH.lg,var.i)
  
  mean.i <- mean(pheno.lg$PH[cross.ind.lg==cross.id.lg[i]])
  mean.PH.lg <- c(mean.PH.lg,mean.i)
  
}


# table descriptive statistics
##############################

mean.DMY <- c(mean.DMY.sh, mean.DMY.het, mean.DMY.lg)
mean.PH <- c(mean.PH.sh, mean.PH.het, mean.PH.lg)

variances.DMY <- c(variances.DMY.sh, variances.DMY.het, variances.DMY.lg)
variances.PH <- c(variances.PH.sh, variances.PH.het, variances.PH.lg)


crosses <- c(cross.id.sh, cross.id.het, cross.id.lg)
parents <- c("D06","D09","F252","F618","UH304", "B73","D06","D09", "UH304",
             "W117", "B73","EC169","Mo17","UH250","W117")

table <- data.frame(crosses, parents, mean.DMY,
                    variances.DMY, mean.PH, variances.PH)

# export xtable

colnames(table) <- c("Cross","Parent2","mean DMY","var. DMY",
                     "mean PH","var. PH")

table.Ltx <- print(xtable(table,digits = c(1,1,1,1,1,1,1),
                          align = c("l","l","l","c","c","c","c")),
                   include.rownames=FALSE)

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

write(x = print(xtable(table,digits = c(1,1,1,1,1,1,1),
                       align = c("l","l","l","c","c","c","c")),
                include.rownames=FALSE),
      file=paste(folder,"des_stat_trait.txt",sep = ""))

################################################################################
########################### End descriptive stat ###############################
################################################################################

# S7 Genetic variance versus parental relatedness
#################################################

path <- "~/MPP_EUNAM/"

setwd(path)

library(mppRDraft)
library(xtable)
library(ggplot2)

crosses <- c("CFD02","CFD03","CFD04","CFD05","CFD06","CFD07",
             "CFD09","CFD10","CFD11","CFD12")

parent2 <- c("B73","D06","D09","EC169","F252","F618","Mo17",
             "UH250","UH304","W117")

# compute SM coefficients

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)

# SM coefficient distance between the central parent (F353)

kin.mat.par <- kinship.matrix(mk.mat = geno.par, method = "SM")

SM <- kin.mat.par[6, ]
SM <- SM[-6]


# compute genetic variance
##########################

library(asreml)

path <- "~/MPP_EUNAM/data"

setwd(path)

pheno <- read.csv("./pheno/pheno_red.csv")
pheno <- pheno[,-1]

lines_used <- read.csv("./pheno/List_lines_Dent_Lehermeier.csv")

# suppress the factor levels with 0 occurences

pheno$Genotype <- as.factor(as.character(pheno$Genotype))
pheno$Fam <- as.factor(as.character(pheno$Fam))

# order data by family

pheno <- pheno[order(pheno$Fam),]

par_names <- c("B73","D06", "D09", "EC169", "F252", "F353", "F618", "F98902",
               "Mo17", "UH250", "UH304","W117")

# DMY

model <- asreml(fixed = DMY~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")


DMY.var <- summary(model)$varcomp
DMY.var <- DMY.var[c(3:12), 2]


# PH

model <- asreml(fixed = PH~1, random= ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                rcov=~at(Fam):units, data = pheno,
                na.method.X = "omit", na.method.Y="omit")

# variance component table

PH.var <- summary(model)$varcomp
PH.var <- PH.var[c(3:12), 2]

data <- data.frame(crosses, parent2,SM,DMY.var,PH.var)

plot(SM, DMY.var)

plot(SM, PH.var)


# DMY

p0 <- ggplot(data, aes(x=SM, y=DMY.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "DMY full population")

# save plot

jpeg(filename = "DMY_pop_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p0

dev.off()


# DMY short

par.short <- c("UH304", "D06","D09","F252","F618")

data.short <- data[data$parent2 %in% par.short,]


p1 <- ggplot(data.short, aes(x=SM, y=DMY.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "DMY short population")

# save plot

jpeg(filename = "DMY_short_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p1

dev.off()

# DMY long

par.long <- c("B73","EC169","Mo17", "UH250","W117")

data.long <- data[data$parent2 %in% par.long,]

p2 <- ggplot(data.long, aes(x=SM, y=DMY.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "DMY long population")

# save plot

jpeg(filename = "DMY_long_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p2

dev.off()

# DMY hetero

par.hetero <- c("UH304", "D06","D09","W117","B73")

data.hetero <- data[data$parent2 %in% par.hetero,]

p3 <- ggplot(data.hetero, aes(x=SM, y=DMY.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "DMY heterogeneous population")

# save plot

jpeg(filename = "DMY_hetero_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p3

dev.off()

################################################################################
################## Plant height ################################################
################################################################################


p0 <- ggplot(data, aes(x=SM, y=PH.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "PH full population")

# save plot

jpeg(filename = "PH_pop_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p0

dev.off()


# PH short

par.short <- c("UH304", "D06","D09","F252","F618")

data.short <- data[data$parent2 %in% par.short,]


p1 <- ggplot(data.short, aes(x=SM, y=PH.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "PH short population")

# save plot

jpeg(filename = "PH_short_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p1

dev.off()

# PH long

par.long <- c("B73","EC169","Mo17", "UH250","W117")

data.long <- data[data$parent2 %in% par.long,]

p2 <- ggplot(data.long, aes(x=SM, y=PH.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "PH long population")

# save plot

jpeg(filename = "PH_long_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p2

dev.off()

# PH hetero

par.hetero <- c("UH304", "D06","D09","W117","B73")

data.hetero <- data[data$parent2 %in% par.hetero,]

p3 <- ggplot(data.hetero, aes(x=SM, y=PH.var, label=parent2)) +
  
  geom_point(shape=1, size=4) +    # Use hollow circles
  
  geom_text(hjust = 0.5, vjust=-1) + 
  
  stat_smooth(method = "lm", aes(colour="linear"), formula = y ~ poly(x, 1),
              size = 1, se = FALSE) +
  stat_smooth(method = "lm", aes(colour="quadratic"), formula = y ~ poly(x, 2),
              size = 1,se = FALSE) +
  
  theme_bw() +
  
  scale_colour_brewer(name = 'Trend', palette = 'Set2') +
  
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=16),
    axis.title.y = element_text(face="bold", color="black", size=16),
    plot.title = element_text(face="bold", color = "black", size=20),
    legend.position=c(1,1),
    legend.justification=c(1,1),
    legend.title = element_text(colour=1, size=14, face="bold"),
    legend.text = element_text(colour=1, size=14)) +
  
  
  labs(x="genetic relatedness (SM)", 
       y = "genetic variance", 
       title= "PH heterogeneous population")

# save plot

jpeg(filename = "PH_hetero_gen_vs_rel.jpg",width = 1200,height = 1200,
     quality = 100, res = 100)

p3

dev.off()

################################################################################
#################### End Genetic variance vs relatedness  ######################
################################################################################

# S11 threshold table
######################

library(xtable)

path <- "~/MPP_EUNAM/data"

setwd(path)


# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("par","anc","biall"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 9)

partition <- data.frame(subset,Q.eff,trait)

folder <- "~/MPP_EUNAM/results/Threshold_permutation/"


# need to split for each elements


# first the average percentage of detection of each QTL

thre.tab <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  max.val <- read.table(paste0(folder, paste0(part[1],"_",part[2],"_",
                                              part[3],".txt")))
  
  thre.tab[i.ind[i],j.ind[i]] <-round(quantile(max.val[,1],0.95),2)
  
  
}

av.Q.inc <- round(colMeans(thre.tab), 2)
av.subset <- round(rowMeans(thre.tab), 2)

thre.tab <- rbind(thre.tab, rep("",3))


thre.tab <- rbind(thre.tab, av.Q.inc)
thre.tab <- cbind(thre.tab, rep("",5), c(av.subset,"",""))


rownames(thre.tab) <- c("short", "het.", "long","\n", "Average")
colnames(thre.tab) <- c("parental", "ancestal", "bi-allelic","\n", "Average")

thre.tab

# export xtable
##############

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"
file=paste0(folder,"threshold_DMY.txt")

print.xtable(thre.tab)

write(x = print(xtable(thre.tab, align = c("l","c","c","c","c","c"),),
                table.placement = "H", include.rownames=TRUE,
                caption.placement = "top"), file=file)


# Plant height
##############

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("par","anc","biall"),times = 3,each = 1)
trait <- rep(c("PH"),times = 9)

partition <- data.frame(subset,Q.eff,trait)

folder <- "~/MPP_EUNAM/results/Threshold_permutation/"


# need to split for each elements


# first the average percentage of detection of each QTL

thre.tab <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  max.val <- read.table(paste0(folder, paste0(part[1],"_",part[2],"_",
                                              part[3],".txt")))
  
  thre.tab[i.ind[i],j.ind[i]] <-round(quantile(max.val[,1],0.95),2)
  
  
}

av.Q.inc <- round(colMeans(thre.tab), 2)
av.subset <- round(rowMeans(thre.tab), 2)

thre.tab <- rbind(thre.tab, rep("",3))


thre.tab <- rbind(thre.tab, av.Q.inc)
thre.tab <- cbind(thre.tab, rep("",5), c(av.subset,"",""))


rownames(thre.tab) <- c("short", "het.", "long","\n", "Average")
colnames(thre.tab) <- c("parental", "ancestal", "bi-allelic","\n", "Average")

thre.tab

# export xtable
##############

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"
file=paste0(folder,"threshold_PH.txt")

print.xtable(thre.tab)

write(x = print(xtable(thre.tab, align = c("l","c","c","c","c","c"),),
                table.placement = "H", include.rownames=TRUE,
                caption.placement = "top"), file=file)

################################################################################
########################### End threshold computation ##########################
################################################################################

# S12 Cross-validation general results
######################################

# DMY
#####

# grand average table
#####################

library(xtable)


# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"


# N.QTL DMY with parenthesis

N.QTL.DMY.txt <- matrix(0,3,3)
N.QTL.DMY.num <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # NQTL
  
  value <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                            part[2],"_",part[4],"_CV/N_QTL.txt")))
  
  N.QTL.DMY.num[i.ind[i],j.ind[i]] <- round(mean(unlist(value)),1)
  
  N.QTL.DMY.txt[i.ind[i],j.ind[i]] <- paste0("(",round(mean(unlist(value)),1),")")
  
  
}


col.m <- paste0("(",round(colMeans(N.QTL.DMY.num), 1),")")
row.n <- paste0("(",round(rowMeans(N.QTL.DMY.num), 1),")")

N.QTL.DMY.txt <- rbind(N.QTL.DMY.txt, rep("",3))

N.QTL.DMY.txt <- rbind(N.QTL.DMY.txt, col.m)
N.QTL.DMY.txt <- cbind(N.QTL.DMY.txt, rep("",5), c(row.n,"",""))


rownames(N.QTL.DMY.txt) <- c("short", "het.", "long","\n", "av.")
colnames(N.QTL.DMY.txt) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

N.QTL.DMY.txt

# p.TS

p.TS.DMY <- matrix(0,3,3)
p.VS.DMY <- matrix(0,3,3)
bias.DMY <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  pes <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  p.TS.DMY[i.ind[i],j.ind[i]] <- round(mean(unlist(pes),na.rm = T),2)
  
  pts <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  p.VS.DMY[i.ind[i],j.ind[i]] <- round(mean(unlist(pts),na.rm = T),2)
  
  bias <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/bias.txt")))
  
  bias.DMY[i.ind[i],j.ind[i]] <- round(mean(unlist(bias),na.rm = T),2)
  
}

# p.TS

col.m <- round(colMeans(p.TS.DMY), 2)
row.n <- round(rowMeans(p.TS.DMY), 2)

p.TS.DMY <- rbind(p.TS.DMY, rep("",3))

p.TS.DMY <- rbind(p.TS.DMY, col.m)
p.TS.DMY <- cbind(p.TS.DMY, rep("",5), c(row.n,"",""))


rownames(p.TS.DMY) <- c("short", "het.", "long","\n", "av.")
colnames(p.TS.DMY) <- c("parental", "ancestal", "bi-allelic","\n", "av.")


p.TS.DMY
N.QTL.DMY.txt

# paste together all the elements of these table

i.ind <- rep(1:5,each=5)
j.ind <- rep(1:5,time=5)

# start the loop here

p.TS.DMY.aug <- matrix(0,5,5)

for(i in 1:25){
  
  p.TS.DMY.aug[i.ind[i], j.ind[i]] <- paste(p.TS.DMY[i.ind[i], j.ind[i]],
                                            N.QTL.DMY.txt[i.ind[i], j.ind[i]])
  
}


# p.VS

col.m <- round(colMeans(p.VS.DMY), 2)
row.n <- round(rowMeans(p.VS.DMY), 2)

p.VS.DMY <- rbind(p.VS.DMY, rep("",3))

p.VS.DMY <- rbind(p.VS.DMY, col.m)
p.VS.DMY <- cbind(p.VS.DMY, rep("",5), c(row.n,"",""))


rownames(p.VS.DMY) <- c("short", "het.", "long","\n", "av.")
colnames(p.VS.DMY) <- c("parental", "ancestal", "bi-allelic","\n", "av.")


p.VS.DMY

# bias

col.m <- round(colMeans(bias.DMY), 2)
row.n <- round(rowMeans(bias.DMY), 2)

bias.DMY <- rbind(bias.DMY, rep("",3))

bias.DMY <- rbind(bias.DMY, col.m)
bias.DMY <- cbind(bias.DMY, rep("",5), c(row.n,"",""))


rownames(bias.DMY) <- c("short", "het.", "long","\n", "av.")
colnames(bias.DMY) <- c("parental", "ancestal", "bi-allelic","\n", "av.")


bias.DMY

# make the complete table

DMY.CV.table <- rbind(p.TS.DMY.aug[1,],p.VS.DMY[1,],bias.DMY[1,],
                      rep("",5),
                      p.TS.DMY.aug[2,],p.VS.DMY[2,],bias.DMY[2,],
                      rep("",5),
                      p.TS.DMY.aug[3,],p.VS.DMY[3,],bias.DMY[3,],
                      rep("",5),
                      rbind(p.TS.DMY.aug[5,], p.VS.DMY[5,], bias.DMY[5,]))

row.names1 <- c("p.TS","p.VS","bias","\n",
                "p.TS","p.VS","bias","\n",
                "p.TS","p.VS","bias","\n",
                "p.TS","p.VS","bias")

row.names2 <- c("\n","short","\n","\n",
                "\n","het.","\n","\n",
                "\n","long","\n","\n",
                "\n","av.","\n")


DMY.CV.table <- cbind(row.names2, row.names1, DMY.CV.table)

colnames(DMY.CV.table)[1:2] <- c("\n", "\n")

###### export the result table

exp.fold <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

write(x = print.xtable(xtable(DMY.CV.table, 
                              caption = "CV results DMY",
                              digits = c(2,2,2,2,2,2,0,2), 
                              align = c("l","l","c","c","c","c","c","c")),
                       table.placement="!ht", caption.placement="top",
                       include.rownames = FALSE),
      file=paste0(exp.fold,"CV_res_DMY.txt"))

############################## end average CV results DMY ######################

# (ind) average table
#####################


library(xtable)

path <- "~/MPP_EUNAM/"

setwd(path)


# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 9)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# first the average percentage of detection of each QTL

av.det.num <- matrix(0,3,3)
av.det.txt <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # select positions detected 10 times or more (20 %)
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  if(dim(QTLs)[1]==0) {
    
    av.det.num[i.ind[i],j.ind[i]] <- 0
    
    av.det.txt[i.ind[i],j.ind[i]] <- paste0("(", 0," \\%",")")    
    
  } else {
    
    av.det.num[i.ind[i],j.ind[i]] <- round(100*mean(QTLs$N)/50)
    
    av.det.txt[i.ind[i],j.ind[i]] <- paste0("(",round(100*mean(QTLs$N)/50)," \\%",")")  
    
  }
  
  
}


col.m <- paste0("(",round(colMeans(av.det.num))," \\%",")")
row.n <- paste0("(",round(rowMeans(av.det.num))," \\%",")")

av.det.txt <- rbind(av.det.txt, rep("",3))

av.det.txt <- rbind(av.det.txt, col.m)
av.det.txt <- cbind(av.det.txt, rep("",5), c(row.n,"",""))


rownames(av.det.txt) <- c("short", "het.", "long","\n", "av.")
colnames(av.det.txt) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

av.det.txt


# p.TS, p.VS, bias

N.20 <- matrix(0,3,3)
p.TS.d <- matrix(0,3,3)
p.TS.sg <- matrix(0,3,3)
p.VS.d <- matrix(0,3,3)
p.VS.sg <- matrix(0,3,3)
bias.d <- matrix(0,3,3)
bias.sg <- matrix(0,3,3)
N.tot <- matrix(0,3,3)


i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  N.tot[i.ind[i],j.ind[i]] <- dim(QTLs)[1]
  
  # select positions detected 10 times or more
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  if(dim(QTLs)[1]==0) {
    
    N.20[i.ind[i],j.ind[i]] <- 0
    p.TS.d[i.ind[i],j.ind[i]] <- 0
    p.TS.sg[i.ind[i],j.ind[i]] <- 0
    p.VS.d[i.ind[i],j.ind[i]] <- 0
    p.VS.sg[i.ind[i],j.ind[i]] <- 0
    bias.d[i.ind[i],j.ind[i]] <- 0
    bias.sg[i.ind[i],j.ind[i]] <- 0
    
  } else {
    
    N.20[i.ind[i],j.ind[i]] <- dim(QTLs)[1]
    p.TS.d[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pes.d),1)
    p.TS.sg[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pes.s),1)
    p.VS.d[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pts.d),1)
    p.VS.sg[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pts.s),1)
    bias.d[i.ind[i],j.ind[i]] <- round(mean(QTLs$bias.d),1)
    bias.sg[i.ind[i],j.ind[i]] <- round(mean(QTLs$bias.s),1)  
    
  }
  
}

# N.20

col.m <- round(colMeans(N.20), 1)
row.n <- round(rowMeans(N.20), 1)

N.20 <- rbind(N.20, rep("",3))

N.20 <- rbind(N.20, col.m)
N.20 <- cbind(N.20, rep("",5), c(row.n,"",""))


rownames(N.20) <- c("short", "het.", "long","\n", "av.")
colnames(N.20) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

N.20

# fuse the number of QTL detected 10 times and the percentage of detection

# paste together all the elements of these table

i.ind <- rep(1:5,each=5)
j.ind <- rep(1:5,time=5)

# start the loop here

N.20.aug <- matrix(0,5,5)

for(i in 1:25){
  
  N.20.aug[i.ind[i], j.ind[i]] <- paste(N.20[i.ind[i], j.ind[i]],
                                        av.det.txt[i.ind[i], j.ind[i]])
  
}


N.20.aug

# N.tot

col.m <- round(colMeans(N.tot), 1)
row.n <- round(rowMeans(N.tot), 1)

N.tot <- rbind(N.tot, rep("",3))

N.tot <- rbind(N.tot, col.m)
N.tot <- cbind(N.tot, rep("",5), c(row.n,"",""))


rownames(N.tot) <- c("short", "het.", "long","\n", "av.")
colnames(N.tot) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

N.tot

# p.TS.d

col.m <- round(colMeans(p.TS.d), 2)
row.n <- round(rowMeans(p.TS.d), 2)

p.TS.d <- rbind(p.TS.d, rep("",3))

p.TS.d <- rbind(p.TS.d, col.m)
p.TS.d <- cbind(p.TS.d, rep("",5), c(row.n,"",""))


rownames(p.TS.d) <- c("short", "het.", "long","\n", "av.")
colnames(p.TS.d) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.TS.d

# p.VS.d

col.m <- round(colMeans(p.VS.d), 2)
row.n <- round(rowMeans(p.VS.d), 2)

p.VS.d <- rbind(p.VS.d, rep("",3))

p.VS.d <- rbind(p.VS.d, col.m)
p.VS.d <- cbind(p.VS.d, rep("",5), c(row.n,"",""))


rownames(p.VS.d) <- c("short", "het.", "long","\n", "av.")
colnames(p.VS.d) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.VS.d

# bias.d

col.m <- round(colMeans(bias.d), 2)
row.n <- round(rowMeans(bias.d), 2)

bias.d <- rbind(bias.d, rep("",3))

bias.d <- rbind(bias.d, col.m)
bias.d <- cbind(bias.d, rep("",5), c(row.n,"",""))


rownames(bias.d) <- c("short", "het.", "long","\n", "av.")
colnames(bias.d) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

bias.d


# p.TS.sg

col.m <- round(colMeans(p.TS.sg), 2)
row.n <- round(rowMeans(p.TS.sg), 2)

p.TS.sg <- rbind(p.TS.sg, rep("",3))

p.TS.sg <- rbind(p.TS.sg, col.m)
p.TS.sg <- cbind(p.TS.sg, rep("",5), c(row.n,"",""))


rownames(p.TS.sg) <- c("short", "het.", "long","\n", "av.")
colnames(p.TS.sg) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.TS.sg

# p.VS.sg

col.m <- round(colMeans(p.VS.sg), 2)
row.n <- round(rowMeans(p.VS.sg), 2)

p.VS.sg <- rbind(p.VS.sg, rep("",3))

p.VS.sg <- rbind(p.VS.sg, col.m)
p.VS.sg <- cbind(p.VS.sg, rep("",5), c(row.n,"",""))


rownames(p.VS.sg) <- c("short", "het.", "long","\n", "av.")
colnames(p.VS.sg) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.VS.sg

# bias.sg

col.m <- round(colMeans(bias.sg), 2)
row.n <- round(rowMeans(bias.sg), 2)

bias.sg <- rbind(bias.sg, rep("",3))

bias.sg <- rbind(bias.sg, col.m)
bias.sg <- cbind(bias.sg, rep("",5), c(row.n,"",""))


rownames(bias.sg) <- c("short", "het.", "long","\n", "av.")
colnames(bias.sg) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

bias.sg

# then combine everything together

# make the complete table

DMY.CV.ind.tab <- rbind(N.tot[1,],N.20.aug[1,],p.TS.d[1,],p.VS.d[1,], bias.d[1,],
                        rep("",5),p.TS.sg[1,],p.VS.sg[1,], bias.sg[1,],rep("",5),
                        N.tot[2,],N.20.aug[2,],p.TS.d[2,],p.VS.d[2,], bias.d[2,],
                        rep("",5),p.TS.sg[2,],p.VS.sg[2,], bias.sg[2,],rep("",5),
                        N.tot[3,],N.20.aug[3,],p.TS.d[3,],p.VS.d[3,], bias.d[3,],
                        rep("",5),p.TS.sg[3,],p.VS.sg[3,], bias.sg[3,],rep("",5),
                        rbind(N.tot[5,],N.20.aug[5,],p.TS.d[5,],p.VS.d[5,],
                              bias.d[5,],rep("",5),p.TS.sg[5,],p.VS.sg[5,],
                              bias.sg[5,]))

row.names1 <- c("$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$","\n",
                "$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$","\n",
                "$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$","\n",
                "$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$")

row.names2 <- c("\n","\n","\n","\n","short","\n","\n","\n","\n","\n",
                "\n","\n","\n","\n","het.","\n","\n","\n","\n","\n",
                "\n","\n","\n","\n","long","\n","\n","\n","\n","\n",
                "\n","\n","\n","\n","av.","\n","\n","\n","\n")

DMY.CV.ind.tab <- cbind(row.names2,row.names1, DMY.CV.ind.tab)

colnames(DMY.CV.ind.tab)[1:2] <- c("\n", "\n")

###### export the result table

exp.fold <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

print.xtable(xtable(DMY.CV.ind.tab,
                    caption = "CV results individual positions DMY",
                    digits = c(2,2,2,2,2,2,0,2),
                    align = c("l","l","c","c","c","c","c","c")),
             table.placement="!ht", caption.placement="top",
             include.rownames = FALSE,sanitize.text.function=function(x){x},
             file=paste0(exp.fold,"CV_ind_res_DMY.txt"))


################################################################################
################################################################################
################################################################################

# plant height
##############

# grand average table
#####################

library(xtable)


# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"


# N.QTL PH with parenthesis

N.QTL.PH.txt <- matrix(0,3,3)
N.QTL.PH.num <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # NQTL
  
  value <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                            part[2],"_",part[4],"_CV/N_QTL.txt")))
  
  N.QTL.PH.num[i.ind[i],j.ind[i]] <- round(mean(unlist(value)),1)
  
  N.QTL.PH.txt[i.ind[i],j.ind[i]] <- paste0("(",round(mean(unlist(value)),1),")")
  
  
}


col.m <- paste0("(",round(colMeans(N.QTL.PH.num), 1),")")
row.n <- paste0("(",round(rowMeans(N.QTL.PH.num), 1),")")

N.QTL.PH.txt <- rbind(N.QTL.PH.txt, rep("",3))

N.QTL.PH.txt <- rbind(N.QTL.PH.txt, col.m)
N.QTL.PH.txt <- cbind(N.QTL.PH.txt, rep("",5), c(row.n,"",""))


rownames(N.QTL.PH.txt) <- c("short", "het.", "long","\n", "av.")
colnames(N.QTL.PH.txt) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

N.QTL.PH.txt

# p.TS

p.TS.PH <- matrix(0,3,3)
p.VS.PH <- matrix(0,3,3)
bias.PH <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  pes <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  p.TS.PH[i.ind[i],j.ind[i]] <- round(mean(unlist(pes),na.rm = T),2)
  
  pts <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  p.VS.PH[i.ind[i],j.ind[i]] <- round(mean(unlist(pts),na.rm = T),2)
  
  bias <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/bias.txt")))
  
  bias.PH[i.ind[i],j.ind[i]] <- round(mean(unlist(bias),na.rm = T),2)
  
}

# p.TS

col.m <- round(colMeans(p.TS.PH), 2)
row.n <- round(rowMeans(p.TS.PH), 2)

p.TS.PH <- rbind(p.TS.PH, rep("",3))

p.TS.PH <- rbind(p.TS.PH, col.m)
p.TS.PH <- cbind(p.TS.PH, rep("",5), c(row.n,"",""))


rownames(p.TS.PH) <- c("short", "het.", "long","\n", "av.")
colnames(p.TS.PH) <- c("parental", "ancestal", "bi-allelic","\n", "av.")


p.TS.PH
N.QTL.PH.txt

# paste together all the elements of these table

i.ind <- rep(1:5,each=5)
j.ind <- rep(1:5,time=5)

# start the loop here

p.TS.PH.aug <- matrix(0,5,5)

for(i in 1:25){
  
  p.TS.PH.aug[i.ind[i], j.ind[i]] <- paste(p.TS.PH[i.ind[i], j.ind[i]],
                                           N.QTL.PH.txt[i.ind[i], j.ind[i]])
  
}


# p.VS

col.m <- round(colMeans(p.VS.PH), 2)
row.n <- round(rowMeans(p.VS.PH), 2)

p.VS.PH <- rbind(p.VS.PH, rep("",3))

p.VS.PH <- rbind(p.VS.PH, col.m)
p.VS.PH <- cbind(p.VS.PH, rep("",5), c(row.n,"",""))


rownames(p.VS.PH) <- c("short", "het.", "long","\n", "av.")
colnames(p.VS.PH) <- c("parental", "ancestal", "bi-allelic","\n", "av.")


p.VS.PH

# bias

col.m <- round(colMeans(bias.PH), 2)
row.n <- round(rowMeans(bias.PH), 2)

bias.PH <- rbind(bias.PH, rep("",3))

bias.PH <- rbind(bias.PH, col.m)
bias.PH <- cbind(bias.PH, rep("",5), c(row.n,"",""))


rownames(bias.PH) <- c("short", "het.", "long","\n", "av.")
colnames(bias.PH) <- c("parental", "ancestal", "bi-allelic","\n", "av.")


bias.PH

# make the complete table

PH.CV.table <- rbind(p.TS.PH.aug[1,],p.VS.PH[1,],bias.PH[1,],
                     rep("",5),
                     p.TS.PH.aug[2,],p.VS.PH[2,],bias.PH[2,],
                     rep("",5),
                     p.TS.PH.aug[3,],p.VS.PH[3,],bias.PH[3,],
                     rep("",5),
                     rbind(p.TS.PH.aug[5,], p.VS.PH[5,], bias.PH[5,]))

row.names1 <- c("p.TS","p.VS","bias","\n",
                "p.TS","p.VS","bias","\n",
                "p.TS","p.VS","bias","\n",
                "p.TS","p.VS","bias")

row.names2 <- c("\n","short","\n","\n",
                "\n","het.","\n","\n",
                "\n","long","\n","\n",
                "\n","av.","\n")


PH.CV.table <- cbind(row.names2, row.names1, PH.CV.table)

colnames(PH.CV.table)[1:2] <- c("\n", "\n")

###### export the result table

exp.fold <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

write(x = print.xtable(xtable(PH.CV.table, 
                              caption = "CV results PH",
                              digits = c(2,2,2,2,2,2,0,2), 
                              align = c("l","l","c","c","c","c","c","c")),
                       table.placement="!ht", caption.placement="top",
                       include.rownames = FALSE),
      file=paste0(exp.fold,"CV_res_PH.txt"))

# (ind) average table
#####################

library(xtable)

path <- "~/MPP_EUNAM/"

setwd(path)


# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 9)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# first the average percentage of detection of each QTL

av.det.num <- matrix(0,3,3)
av.det.txt <- matrix(0,3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # select positions detected 10 times or more (20 %)
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  if(dim(QTLs)[1]==0) {
    
    av.det.num[i.ind[i],j.ind[i]] <- 0
    
    av.det.txt[i.ind[i],j.ind[i]] <- paste0("(", 0," \\%",")")    
    
  } else {
    
    av.det.num[i.ind[i],j.ind[i]] <- round(100*mean(QTLs$N)/50)
    
    av.det.txt[i.ind[i],j.ind[i]] <- paste0("(",round(100*mean(QTLs$N)/50)," \\%",")")  
    
  }
  
  
}


col.m <- paste0("(",round(colMeans(av.det.num))," \\%",")")
row.n <- paste0("(",round(rowMeans(av.det.num))," \\%",")")

av.det.txt <- rbind(av.det.txt, rep("",3))

av.det.txt <- rbind(av.det.txt, col.m)
av.det.txt <- cbind(av.det.txt, rep("",5), c(row.n,"",""))


rownames(av.det.txt) <- c("short", "het.", "long","\n", "av.")
colnames(av.det.txt) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

av.det.txt


# p.TS, p.VS, bias

N.20 <- matrix(0,3,3)
p.TS.d <- matrix(0,3,3)
p.TS.sg <- matrix(0,3,3)
p.VS.d <- matrix(0,3,3)
p.VS.sg <- matrix(0,3,3)
bias.d <- matrix(0,3,3)
bias.sg <- matrix(0,3,3)
N.tot <- matrix(0,3,3)


i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  N.tot[i.ind[i],j.ind[i]] <- dim(QTLs)[1]
  
  # select positions detected 10 times or more
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  if(dim(QTLs)[1]==0) {
    
    N.20[i.ind[i],j.ind[i]] <- 0
    p.TS.d[i.ind[i],j.ind[i]] <- 0
    p.TS.sg[i.ind[i],j.ind[i]] <- 0
    p.VS.d[i.ind[i],j.ind[i]] <- 0
    p.VS.sg[i.ind[i],j.ind[i]] <- 0
    bias.d[i.ind[i],j.ind[i]] <- 0
    bias.sg[i.ind[i],j.ind[i]] <- 0
    
  } else {
    
    N.20[i.ind[i],j.ind[i]] <- dim(QTLs)[1]
    p.TS.d[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pes.d),1)
    p.TS.sg[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pes.s),1)
    p.VS.d[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pts.d),1)
    p.VS.sg[i.ind[i],j.ind[i]] <- round(mean(QTLs$av.pts.s),1)
    bias.d[i.ind[i],j.ind[i]] <- round(mean(QTLs$bias.d),1)
    bias.sg[i.ind[i],j.ind[i]] <- round(mean(QTLs$bias.s),1)  
    
  }
  
}

# N.20

col.m <- round(colMeans(N.20), 1)
row.n <- round(rowMeans(N.20), 1)

N.20 <- rbind(N.20, rep("",3))

N.20 <- rbind(N.20, col.m)
N.20 <- cbind(N.20, rep("",5), c(row.n,"",""))


rownames(N.20) <- c("short", "het.", "long","\n", "av.")
colnames(N.20) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

N.20

# fuse the number of QTL detected 10 times and the percentage of detection

# paste together all the elements of these table

i.ind <- rep(1:5,each=5)
j.ind <- rep(1:5,time=5)

# start the loop here

N.20.aug <- matrix(0,5,5)

for(i in 1:25){
  
  N.20.aug[i.ind[i], j.ind[i]] <- paste(N.20[i.ind[i], j.ind[i]],
                                        av.det.txt[i.ind[i], j.ind[i]])
  
}


N.20.aug

# N.tot

col.m <- round(colMeans(N.tot), 1)
row.n <- round(rowMeans(N.tot), 1)

N.tot <- rbind(N.tot, rep("",3))

N.tot <- rbind(N.tot, col.m)
N.tot <- cbind(N.tot, rep("",5), c(row.n,"",""))


rownames(N.tot) <- c("short", "het.", "long","\n", "av.")
colnames(N.tot) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

N.tot

# p.TS.d

col.m <- round(colMeans(p.TS.d), 2)
row.n <- round(rowMeans(p.TS.d), 2)

p.TS.d <- rbind(p.TS.d, rep("",3))

p.TS.d <- rbind(p.TS.d, col.m)
p.TS.d <- cbind(p.TS.d, rep("",5), c(row.n,"",""))


rownames(p.TS.d) <- c("short", "het.", "long","\n", "av.")
colnames(p.TS.d) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.TS.d

# p.VS.d

col.m <- round(colMeans(p.VS.d), 2)
row.n <- round(rowMeans(p.VS.d), 2)

p.VS.d <- rbind(p.VS.d, rep("",3))

p.VS.d <- rbind(p.VS.d, col.m)
p.VS.d <- cbind(p.VS.d, rep("",5), c(row.n,"",""))


rownames(p.VS.d) <- c("short", "het.", "long","\n", "av.")
colnames(p.VS.d) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.VS.d

# bias.d

col.m <- round(colMeans(bias.d), 2)
row.n <- round(rowMeans(bias.d), 2)

bias.d <- rbind(bias.d, rep("",3))

bias.d <- rbind(bias.d, col.m)
bias.d <- cbind(bias.d, rep("",5), c(row.n,"",""))


rownames(bias.d) <- c("short", "het.", "long","\n", "av.")
colnames(bias.d) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

bias.d


# p.TS.sg

col.m <- round(colMeans(p.TS.sg), 2)
row.n <- round(rowMeans(p.TS.sg), 2)

p.TS.sg <- rbind(p.TS.sg, rep("",3))

p.TS.sg <- rbind(p.TS.sg, col.m)
p.TS.sg <- cbind(p.TS.sg, rep("",5), c(row.n,"",""))


rownames(p.TS.sg) <- c("short", "het.", "long","\n", "av.")
colnames(p.TS.sg) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.TS.sg

# p.VS.sg

col.m <- round(colMeans(p.VS.sg), 2)
row.n <- round(rowMeans(p.VS.sg), 2)

p.VS.sg <- rbind(p.VS.sg, rep("",3))

p.VS.sg <- rbind(p.VS.sg, col.m)
p.VS.sg <- cbind(p.VS.sg, rep("",5), c(row.n,"",""))


rownames(p.VS.sg) <- c("short", "het.", "long","\n", "av.")
colnames(p.VS.sg) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

p.VS.sg

# bias.sg

col.m <- round(colMeans(bias.sg), 2)
row.n <- round(rowMeans(bias.sg), 2)

bias.sg <- rbind(bias.sg, rep("",3))

bias.sg <- rbind(bias.sg, col.m)
bias.sg <- cbind(bias.sg, rep("",5), c(row.n,"",""))


rownames(bias.sg) <- c("short", "het.", "long","\n", "av.")
colnames(bias.sg) <- c("parental", "ancestal", "bi-allelic","\n", "av.")

bias.sg

# then combine everything together

# make the complete table

PH.CV.ind.tab <- rbind(N.tot[1,],N.20.aug[1,],p.TS.d[1,],p.VS.d[1,], bias.d[1,],
                       rep("",5),p.TS.sg[1,],p.VS.sg[1,], bias.sg[1,],rep("",5),
                       N.tot[2,],N.20.aug[2,],p.TS.d[2,],p.VS.d[2,], bias.d[2,],
                       rep("",5),p.TS.sg[2,],p.VS.sg[2,], bias.sg[2,],rep("",5),
                       N.tot[3,],N.20.aug[3,],p.TS.d[3,],p.VS.d[3,], bias.d[3,],
                       rep("",5),p.TS.sg[3,],p.VS.sg[3,], bias.sg[3,],rep("",5),
                       rbind(N.tot[5,],N.20.aug[5,],p.TS.d[5,],p.VS.d[5,],
                             bias.d[5,],rep("",5),p.TS.sg[5,],p.VS.sg[5,],
                             bias.sg[5,]))

row.names1 <- c("$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$","\n",
                "$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$","\n",
                "$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$","\n",
                "$N_{tot}$","$N(\\geq 20\\%)$","$\\bar{p.TS.d}$",
                "$\\bar{p.VS.d}$","$\\bar{bias.d}$","\n","$\\bar{p.TS.sg}$",
                "$\\bar{p.VS.sg}$","$\\bar{bias.sg}$")

row.names2 <- c("\n","\n","\n","\n","short","\n","\n","\n","\n","\n",
                "\n","\n","\n","\n","het.","\n","\n","\n","\n","\n",
                "\n","\n","\n","\n","long","\n","\n","\n","\n","\n",
                "\n","\n","\n","\n","av.","\n","\n","\n","\n")

PH.CV.ind.tab <- cbind(row.names2,row.names1, PH.CV.ind.tab)

colnames(PH.CV.ind.tab)[1:2] <- c("\n", "\n")

###### export the result table

exp.fold <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

print.xtable(xtable(PH.CV.ind.tab,
                    caption = "CV results individual positions PH",
                    digits = c(2,2,2,2,2,2,0,2),
                    align = c("l","l","c","c","c","c","c","c")),
             table.placement="!ht", caption.placement="top",
             include.rownames = FALSE,sanitize.text.function=function(x){x},
             file=paste0(exp.fold,"CV_ind_res_PH.txt"))


################################################################################
############################# End CV general results ###########################
################################################################################

# S13 Cross-validation detailed results
#######################################

# Library

library(mppRDraft)
library(qtl)
library(asreml)
library(MASS)
library(outliers)
library(stringr)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(igraph)
library(xtable)

# cross-validation results


setwd(path)


# partition

subset <- rep(c("short","hetero","long"),each=6)
title <- rep(c("Short","Heterogeneous","Long"),each=6)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 6,each = 1)
title2 <- rep(c("parental","ancestral","bi-allelic"),times = 6,each = 1)
trait <- rep(c("DMY","PH"),times = 3,each = 3)

partition <- data.frame(subset,trait,Q.eff,title,title2)

my.loc <- "~/MPP_EUNAM/results/CV/"
fig.loc <- "E:/PhD/Manuscript/1st_chapter/Latex/figures/"
CV.sup.file <- "E:/PhD/Manuscript/1st_chapter/Latex/CV_sup_mat.txt"
fig.cmd.line <- "\\includegraphics[width=\\textwidth, natwidth=610, natheight=642]{figures/"
tab.loc <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

# start the loop here

for(i in 1:18){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  # load the CV results
  
  file <- paste0(my.loc,paste(part[1:3],collapse = "_"),"_u.err_CV/QTL_profiles.txt")
  
  CV.res <- read.table(file)
  
  # put the results in a list format
  
  res <- list(QTL.profiles=CV.res)
  
  # produce the CV plot
  
  # save in the figures location
  
  pdf(paste0(fig.loc,paste(part[1:3],collapse = "_"),"_CV.pdf"),
      height = 10, width = 16)
  
  mpp.plot.CV(res = res, main = paste(part[1:3],collapse = " "))
  
  dev.off()
  
  # table CV QTL results
  
  file <- paste0(my.loc,paste(part[1:3],collapse = "_"),"_u.err_CV/QTL.txt")
  
  QTL.res <- read.table(file)
  
  # select only the positins detected 10% or more
  
  QTL.res <- QTL.res[QTL.res$N >= 5, ]
  QTL.res <- QTL.res[,-c(3,9,10,11,15,16,17)]
  QTL.res <- QTL.res[, c(1:4, 8:10, 5:7)]
  colnames(QTL.res)[c(1,3,5,6,8,9)] <- c("Mk. names","pos [cM]","av.pts.d",
                                         "av.pvs.d","av.pts.s","av.pvs.s")
  
  # make the latex table
  
  print.xtable(xtable(QTL.res, digits = c(0,0,0,1,0,1,1,1,1,1,1),
                      align = c("l","l","c","c","c","c","c","c","c","c","c")),
               table.placement="H", include.rownames = FALSE,
               sanitize.text.function=function(x){x},
               file=paste0(tab.loc,paste(part[1:3],collapse = "_"),
                           "_tab.txt"))
  
  # write the text for the figure supplemental material
  
  el <- paste0("\\subsection{",paste(part[c(4,2,5)],collapse = " "),"}")
  
  write(el,file = CV.sup.file,append = TRUE)
  write("\n",file = CV.sup.file,append = TRUE)
  write("\\iffalse",file = CV.sup.file,append = TRUE)
  write("\n",file = CV.sup.file,append = TRUE)
  write("\\begin{figure}[H]",file = CV.sup.file,append = TRUE)
  write(paste0(fig.cmd.line,paste(part[1:3],collapse = "_"),"_CV.pdf","}"),
        file = CV.sup.file,append = TRUE)
  write("\\end{figure}",file = CV.sup.file,append = TRUE)
  write("\n",file = CV.sup.file,append = TRUE)
  write("\\fi",file = CV.sup.file,append = TRUE)
  write("\n",file = CV.sup.file,append = TRUE)
  
  el <- paste0("\\input{table/", paste(part[1:3],collapse = "_"), "_tab.txt", "}")
  
  write(el, file = CV.sup.file,append = TRUE)
  write("\n",file = CV.sup.file,append = TRUE)
  write("\\newpage",file = CV.sup.file,append = TRUE)
  write("\n",file = CV.sup.file,append = TRUE)
  
  
}

################################################################################
############################# End CV detailed results ##########################
################################################################################

# S14 Cross-validation detailed results
#######################################

# Library

library(mppRDraft)
library(qtl)
library(asreml)
library(MASS)
library(outliers)
library(stringr)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(igraph)
library(xtable)


# SPECIFY YOUR PATH HERE (the location of EU_NAM_DENT folder)

path <- "~/MPP_EUNAM/"

setwd(path)


# partition

sub <- rep(x = c("short","hetero","long"),each=6,times=2)
title <- rep(x = c("Short","Heterogeneous","Long"),each=6,times=2)
trait <- rep(c("DMY","PH"),each=18)
model <- rep(c("parental", "ancestral", "biallelic"),each=2,times=6)
effect <- c(rep("CIM",12), rep("effect_CIM",6), rep("CIM",12),
            rep("effect_CIM",6)) 
VCOV <- rep(c("u.err","cr.err"),times=18)

folder <- "./results/Threshold_permutation/"

part.fold <- cbind(sub,trait,model,effect,VCOV)

sub <- rep(x = c("short","hetero","long"),each=6,times=2)
title <- rep(x = c("Short","Heterogeneous","Long"),each=6,times=2)
Q.eff <- rep(c("par","anc","biall"),times = 6,each = 2)
title2 <- rep(c("parental","ancestral","bi-allelic"),times = 6,each = 2)
trait <- rep(c("DMY","PH"),each=18)
# VCOV <- rep(c("u.err","cr.err"),times=18)
VCOV <- rep(c("uerr","crerr"),times=18)
title3 <- rep(c("HRT","CSRT"),times=18)
data.type <- rep(c("ABH","ABH","biall"), each=2, times=6)

part.data <- data.frame(sub,Q.eff,trait,VCOV,data.type, title, title2, title3)


my.loc <- "~/MPP_EUNAM/results/QTL_analyses/"
folder <- "./results/Threshold_permutation/"

fig.loc <- "E:/PhD/Manuscript/1st_chapter/Latex/figures/"
QTL.sup.file <- "E:/PhD/Manuscript/1st_chapter/Latex/QTL_sup_mat2.txt"
fig.cmd.line <- "\\includegraphics[width=\\textwidth, natwidth=610, natheight=642]{figures/"
tab.loc <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"

# start the loop here


for(i in 1:36){
  
  # extract the ith partition
  part <- as.character(unlist(part.data[i,]))
  
  # get threshold
  
  max.val <- read.table(paste0(folder, paste0(part[1],"_",part[2],"_",
                                              part[3],".txt")))
  
  thre <- quantile(max.val[,1],0.95)
  
  # open the data set
  data <- readRDS(file = paste0("./data/mpp_data/data.",part[1],".",part[5],".rds"))
  
  # load cim results
  
  file <- paste0(my.loc,paste0(part.fold[i,],collapse = "_"),"/CIM.txt")
  
  CIM <- read.table(file)
  
  # plot the CIM profile
  
  if(part[2]=="biall"){
    
    # # CIM profile  
    # 
    # pdf(paste0(fig.loc,paste(part[1:4],collapse = "_"),"_QTL_prof.pdf"),
    #     height = 10, width = 16)  
    # 
    # print(mpp.plot(mpp.res = CIM,type = "h", threshold = thre,
    #                main = paste(part[1:4],collapse=" ")))
    # 
    # dev.off()
    
    # write the text for the figure supplemental material
    
    el <- paste0("\\subsection{",paste(part[c(6,3,7,8)],collapse = " "),"}")
    
    write(el,file = QTL.sup.file,append = TRUE)
    
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\iffalse",file = QTL.sup.file, append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    
    write("\\begin{figure}[H]",file = QTL.sup.file,append = TRUE)
    write(paste0(fig.cmd.line,paste(part[1:4],collapse = "_"),"_QTL_prof.pdf","}"),
          file = QTL.sup.file,append = TRUE)
    write("\\end{figure}",file = QTL.sup.file,append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\newpage",file = QTL.sup.file,append = TRUE)
    
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\fi",file = QTL.sup.file, append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    
  } else {
    
    # CIM profile
    
    # pdf(paste0(fig.loc,paste(part[1:4],collapse = "_"),"_QTL_prof.pdf"),
    #     height = 10, width = 16)
    # 
    # print(mpp.plot(mpp.res = CIM,type = "l", threshold = thre,
    #                main = paste(part[1:4],collapse=" ")))
    # 
    # dev.off()
    
    # write text for sup material
    
    el <- paste0("\\subsection{",paste(part[c(6,3,7,8)],collapse = " "),"}")
    
    write(el,file = QTL.sup.file,append = TRUE)
    
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\iffalse",file = QTL.sup.file, append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    
    write("\\begin{figure}[H]",file = QTL.sup.file,append = TRUE)
    write(paste0(fig.cmd.line,paste(part[1:4],collapse = "_"),"_QTL_prof.pdf","}"),
          file = QTL.sup.file,append = TRUE)
    write("\\end{figure}",file = QTL.sup.file,append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    
    
    # genetic effect profile
    
    # pdf(paste0(fig.loc,paste(part[1:4],collapse = "_"),"_gen_eff.pdf"),
    #     height = 10, width = 16)
    # 
    # gen.eff.plot(data = data, mpp.res = CIM, Q.eff = part[2], threshold = thre,
    #              main = paste(part[1:4],collapse=" "))
    # 
    # dev.off()
    
    # write text for sup material
    
    write("\\begin{figure}[H]",file = QTL.sup.file,append = TRUE)
    write(paste0(fig.cmd.line,paste(part[1:4],collapse = "_"),"_gen_eff.pdf","}"),
          file = QTL.sup.file,append = TRUE)
    write("\\end{figure}",file = QTL.sup.file,append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\newpage",file = QTL.sup.file, append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\fi",file = QTL.sup.file, append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    
    
  }
  
  # produce the table
  
  if((i != 7) & (i != 9)){
    
    file <- paste0(my.loc,paste0(part.fold[i,],collapse = "_"),"/table.QTL.effects.txt")
    
    QTL.table <- read.table(file, sep = "\t",strip.white = TRUE,h=T)
    
    colnames(QTL.table)[c(2,4,7,8)] <- c("Mk names","pos [cM]","-log10pval","r2")
    
    print.xtable(scalebox = 0.75,xtable(QTL.table, digits = c(0,0,0,0,1,1,1,1,1,0,2,2,2,2,2),
                                        align = c("l","l","l",rep("c",12))),
                 table.placement="H", include.rownames = FALSE,
                 sanitize.text.function=function(x){x},
                 file=paste0(tab.loc,paste(part[1:4],collapse = "_"),
                             "_Qeff_tab.txt"))
    
    # write commande for the table
    
    
    el <- paste0("\\input{table/", paste(part[1:4],collapse = "_"),
                 "_Qeff_tab.txt", "}")
    
    write(el, file = QTL.sup.file,append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    write("\\newpage",file = QTL.sup.file,append = TRUE)
    write("\n",file = QTL.sup.file,append = TRUE)
    
  }
  
  
  
}


################################################################################
############################# End CV detailed results ##########################
################################################################################
###########################################
# table and figures Supplemental material #
###########################################

# S2 Kinship matrix GS
######################


library(mppR)
library(xtable)

path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)
geno.par <- as.matrix(geno.par)

# SM coefficient distance between parents

kin.mat.par <- SM_comp(mk.mat = geno.par)

kin.mat.par <- round(kin.mat.par, 3)
kin.mat.par2 <- kin.mat.par

kin.mat.par2[lower.tri(kin.mat.par2)] <- ""

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

par.short <- c("UH304", "F252", "D09", "F618", "D06")
par.hetero <- c("UH304", "D09", "D06", "W117", "B73")
par.long <-  c("UH250", "Mo17",  "W117", "EC169", "B73")

kin.mat.short <- kin.mat.par[rownames(kin.mat.par) %in%
                               c("F353",as.character(par.short)), ]
kin.mat.short <- kin.mat.short[, colnames(kin.mat.short) %in%
                                 c("F353",as.character(par.short))]

kin.mat.hetero <- kin.mat.par[rownames(kin.mat.par) %in%
                                c("F353",as.character(par.hetero)), ]
kin.mat.hetero <- kin.mat.hetero[, colnames(kin.mat.hetero) %in%
                                   c("F353",as.character(par.hetero))]

kin.mat.long <- kin.mat.par[rownames(kin.mat.par) %in%
                              c("F353",as.character(par.long)), ]
kin.mat.long <- kin.mat.long[, colnames(kin.mat.long) %in%
                               c("F353",as.character(par.long))]

mean(kin.mat.short[lower.tri(kin.mat.short)])
mean(kin.mat.hetero[lower.tri(kin.mat.hetero)])
mean(kin.mat.long[lower.tri(kin.mat.long)])

################################################################################
######################### End Kinship matrix GS ################################
################################################################################

# S3 PC biplot of the EU-NAM Dent parents
#########################################

path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)
geno.par <- as.matrix(geno.par)

mk.mat <-  geno.par
cross.ind <- rownames(geno.par)

col.vect <- c("red", "blue", "green", "cyan", "magenta", "black", "grey",
              "black", "brown", "orange", "violet")

shape.vect <- c(3, 3, 3, 3, 3, 15, 3, 3, 3, 3, 3)

# remove monomorphic or missing columns

MAF <- QC_MAF(mk.mat = mk.mat)
index <- which((MAF == 0)| is.na(MAF))
mk.mat <- mk.mat[, -index]

# impute to have complete data

Sample  <-  function(x, ...) if (length(x) == 1) { x
} else {sample(x, replace = TRUE, ...)}
count  <-  function(v) sum(v)
impute  <-  function(d) {
  r  <-  apply(d, 2, function(col) {
    col0  <-  col
    col[is.na(col)] = Sample(na.omit(col), count(is.na(col)))
    col
  })
  r
}

mk.mat <- impute(mk.mat)

mk.mat <- geno_012(mk.mat)[[1]]

# compute PCs

res <- prcomp(mk.mat, center = TRUE, scale. = TRUE)

var.PC <- round((res$sdev[1:2]^2)/sum(res$sdev^2)*100, 2)

# Plot PCs

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)

plot(res$x[,1], res$x[,2], col = col.vect, main = "EU-NAM Dent parents",
     xlab = paste("PC1", paste("(", var.PC[1], "%", ")", sep = "")),
     ylab = paste("PC2", paste("(",var.PC[2], "%",")", sep = "")),
     pch = shape.vect)

# Add legend to top right, outside plot region
legend("topright", inset = c(-0.2, 0), legend = levels(as.factor(cross.ind)),
       title = "Parents", pch = shape.vect, col = col.vect)


# save the plot

setwd("E:/PhD/Manuscript/1st_chapter/Latex/figures")

jpeg(filename = "Parents_PC_biplot.jpg", width = 758, height = 585,
     quality = 100, res = 100)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)

plot(res$x[,1], res$x[,2], col = col.vect, main = "EU-NAM Dent parents",
     xlab = paste("PC1", paste("(", var.PC[1], "%", ")", sep = "")),
     ylab = paste("PC2", paste("(",var.PC[2], "%",")", sep = "")),
     pch = shape.vect)

# Add legend to top right, outside plot region
legend("topright", inset = c(-0.2, 0), legend = levels(as.factor(cross.ind)),
       title = "Parents", pch = shape.vect, col = col.vect)

dev.off()

################################################################################
############################## End PC biplot ###################################
################################################################################

# S4 Genetic map subsets
########################

library(qtl)
library(xtable)

path <- "~/MPP_EUNAM/"

setwd(path)

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/"

# short
#######

# load the mppData objects

mppData <- readRDS("./data/mpp_data/data_short_IBD.rds")

map <- summaryMap(mppData$geno)

map <- data.frame(c(1:10, "Overal"), map)
colnames(map) <- c("Chromosome","N","Length(cM)","Average spacing",
                   "maximum spacing")

map.Ltx <- print(xtable(map, digits = c(0,0,0,1,1,1)), include.rownames = FALSE)

write(x = print(xtable(map, digits = c(0,0,0,1,1,1)),
                include.rownames = FALSE),
      file = paste(folder, "table/", "map_sum_short.txt", sep = ""))

jpeg(filename = paste(folder, "figures/", "map_short.jpeg", sep = ""),
     width = 1000, height = 1000, units = "px", pointsize = 36, quality = 100)

plot.map(mppData$geno, main = "EU-NAM Dent short subset map")

dev.off()

# heterogeneous
###############

# load the mppData objects

mppData <- readRDS("./data/mpp_data/data_hetero_IBD.rds")

map <- summaryMap(mppData$geno)

map <- data.frame(c(1:10, "Overal"), map)
colnames(map) <- c("Chromosome","N","Length(cM)","Average spacing",
                   "maximum spacing")

map.Ltx <- print(xtable(map, digits = c(0,0,0,1,1,1)), include.rownames = FALSE)

write(x = print(xtable(map, digits = c(0,0,0,1,1,1)),
                include.rownames = FALSE),
      file = paste(folder, "table/", "map_sum_hetero.txt", sep = ""))

jpeg(filename = paste(folder, "figures/", "map_hetero.jpeg", sep = ""),
     width = 1000, height = 1000, units = "px", pointsize = 36, quality = 100)

plot.map(mppData$geno, main = "EU-NAM Dent hetero subset map")

dev.off()

# long
######

# load the mppData objects

mppData <- readRDS("./data/mpp_data/data_long_IBD.rds")

map <- summaryMap(mppData$geno)

map <- data.frame(c(1:10, "Overal"), map)
colnames(map) <- c("Chromosome","N","Length(cM)","Average spacing",
                   "maximum spacing")

map.Ltx <- print(xtable(map, digits = c(0,0,0,1,1,1)), include.rownames = FALSE)

write(x = print(xtable(map, digits = c(0,0,0,1,1,1)),
                include.rownames = FALSE),
      file = paste(folder, "table/", "map_sum_long.txt", sep = ""))

jpeg(filename = paste(folder, "figures/", "map_long.jpeg", sep = ""),
     width = 1000, height = 1000, units = "px", pointsize = 36, quality = 100)

plot.map(mppData$geno, main = "EU-NAM Dent long subset map")

dev.off()

################################################################################
########################### End map subsets ####################################
################################################################################

# S7 Genetic variance versus parental relatedness
#################################################

library(asreml)
library(mppR)
library(xtable)
library(ggplot2)


path <- "~/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)
geno.par <- as.matrix(geno.par)

# SM coefficient distance between parents

kin.mat.par <- SM_comp(mk.mat = geno.par)

SM <- kin.mat.par[6, ]
SM <- SM[-6]

crosses <- c("CFD02","CFD03","CFD04","CFD05","CFD06","CFD07",
             "CFD09","CFD10","CFD11","CFD12")

parent2 <- c("B73","D06","D09","EC169","F252","F618","Mo17",
             "UH250","UH304","W117")


# compute genetic variance
##########################

pheno <- read.csv("./data/pheno/pheno_red.csv", row.names = 1)

lines_used <- read.csv("./data/pheno/List_lines_Dent_Lehermeier.csv")

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

data <- data.frame(crosses, parent2, SM, DMY.var, PH.var)

plot(SM, DMY.var)

plot(SM, PH.var)


# DMY

p0 <- ggplot(data, aes(x = SM, y = DMY.var, label = parent2)) +
  
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

setwd("E:/PhD/Manuscript/1st_chapter/Latex/figures")

jpeg(filename = "DMY_pop_gen_vs_rel.jpg", width = 1200, height = 1200,
     quality = 100, res = 100)

p0

dev.off()


# PH

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


################################################################################
#################### End Genetic variance vs relatedness  ######################
################################################################################


# S9 threshold table
######################

library(xtable)

path <- "~/MPP_EUNAM/"

setwd(path)


# partition

subset <- rep(c("short", "hetero", "long"),each=3)
Q.eff <- rep(c("par", "anc", "biall"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 9)

partition <- data.frame(subset,Q.eff,trait)

folder <- "./results/Threshold_permutation/"


# need to split for each elements


# first the average percentage of detection of each QTL

thre.tab <- matrix(0, 3, 3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# start the loop here

for(i in 1:9){
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  max.val <- read.table(paste0(folder, paste0(part, collapse = "_"),
                               "_h.err.txt"))
  
  thre.tab[i.ind[i],j.ind[i]] <-round(quantile(max.val[, 1], 0.95), 2)
  
  
}

av.Q.inc <- round(colMeans(thre.tab), 2)
av.subset <- round(rowMeans(thre.tab), 2)

thre.tab <- rbind(thre.tab, rep("",3))


thre.tab <- rbind(thre.tab, av.Q.inc)
thre.tab <- cbind(thre.tab, rep("",5), c(av.subset,"",""))


rownames(thre.tab) <- c("short", "het.", "long","\n", "Average")
colnames(thre.tab) <- c("parental", "ancestal", "bi-allelic","\n", "Average - MQE")

thre.tab

# export xtable
##############

folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"
file = paste0(folder, "threshold_DMY.txt")


write(x = print(xtable(thre.tab, align = c("l","c","c","c","c","c"),),
                table.placement = "H", include.rownames = TRUE,
                caption.placement = "top"), file = file)


# Plant height
##############

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("par","anc","biall"),times = 3,each = 1)
trait <- rep(c("PH"),times = 9)

partition <- data.frame(subset,Q.eff,trait)

folder <- "./results/Threshold_permutation/"


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
  max.val <- read.table(paste0(folder, paste0(part, collapse = "_"),
                               "_h.err.txt"))
  
  thre.tab[i.ind[i],j.ind[i]] <-round(quantile(max.val[,1],0.95),2)
  
  
}

av.Q.inc <- round(colMeans(thre.tab), 2)
av.subset <- round(rowMeans(thre.tab), 2)

thre.tab <- rbind(thre.tab, rep("",3))


thre.tab <- rbind(thre.tab, av.Q.inc)
thre.tab <- cbind(thre.tab, rep("",5), c(av.subset,"",""))


rownames(thre.tab) <- c("short", "het.", "long","\n", "Average")
colnames(thre.tab) <- c("parental", "ancestal", "bi-allelic","\n", "Average - MQE")

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
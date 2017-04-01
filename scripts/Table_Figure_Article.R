#############################
# table and figures article #
#############################

# Table 1
#########

library(mppR)
library(xtable)

path <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/"

setwd(path)

# Genotype

geno.par <- read.csv("./data/geno/geno_panzea_par.csv", row.names = 1)
geno.par <- as.matrix(geno.par)

# SM coefficient distance between the central parent (F353)

kin.mat.par <- SM_comp(mk.mat = geno.par)

dist.cent.par <- kin.mat.par[, 6]
dist.cent.par <- dist.cent.par[ -6]

# relate information distance between the two parents to the cross
# information

crosses <- c("CFD02","CFD03","CFD04","CFD05","CFD06","CFD07",
             "CFD09","CFD10","CFD11","CFD12")

parent2 <- c("B73","D06","D09","EC169","F252","F618","Mo17",
             "UH250","UH304","W117")

# sample sizes

geno_used <- read.table("./data/geno/geno_sel_list.txt")


short.cr <- substr(geno_used[, 1], 1, 5)
hetero.cr <- substr(geno_used[, 2], 1, 5)
long.cr <- substr(geno_used[, 3], 1, 5)

short.cr.sz <- table(short.cr)
hetero.cr.sz <- table(hetero.cr)
long.cr.sz <- table(long.cr)

N.sh <- rep(0, 10)
N.sh[crosses %in% names(short.cr.sz)] <- short.cr.sz
N.het <- rep(0, 10)
N.het[crosses %in% names(hetero.cr.sz)] <- hetero.cr.sz
N.lg <- rep(0, 10)
N.lg[crosses %in% names(long.cr.sz)] <- long.cr.sz

# genetic variance components

library(asreml)

pheno <- read.csv("./data/pheno/pheno_red.csv", row.names = 1)
lines_used <- read.csv("./data/pheno/List_lines_Dent_Lehermeier.csv")
adj_means <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
cross.ind <- substr(rownames(adj_means), 1, 5)

# suppress the factor levels with 0 occurences

pheno$Genotype <- as.factor(as.character(pheno$Genotype))
pheno$Fam <- as.factor(as.character(pheno$Fam))

# order data by family

pheno <- pheno[order(pheno$Fam), ]


trait.id <- c("DMY", "PH")
trait.stats <- c()


for (i in 1:2){
  
  trait.av <- tapply(X = adj_means[, trait.id[i]], INDEX = as.factor(cross.ind),
                     FUN = function(x) mean(x, na.rm = TRUE))
  
  trait.av <- round(trait.av[-11], 1)
  
  fix_formula <- paste(trait.id[i], " ~ 1", sep = "")
  
  model <- asreml(fixed = as.formula(fix_formula),
                  random = ~ Fam + at(Fam):Genotype + 
                  at(Fam):LOC:Genotype + LOC:Rep + LOC:Rep:Block,
                  rcov=~at(Fam):units, data = pheno,
                  na.method.X = "omit", na.method.Y = "omit")
  
  var.comp <- summary(model)$varcomp
  var.comp <-round(cbind(var.comp[c(3:12), 2], var.comp[c(14:23), 2],
                         var.comp[c(27:36), 2]), 2)
  
  her <- round(100*(var.comp[, 1]/(var.comp[, 1] + (var.comp[, 2]/4) +
                                     (var.comp[, 3]/5))), 1)
  
  trait.stats <- cbind(trait.stats, cbind(trait.av, var.comp[, 1], her))
  
} 

blank <- rep("", 10)

subset.part <- data.frame(crosses, parent2, dist.cent.par, trait.stats[, 1:3],
                          blank, trait.stats[, 4:6], N.sh, N.het, N.lg)

subset.part <- subset.part[order(subset.part[, 3], decreasing = TRUE), ]


colnames(subset.part) <- c("Cross","Parent","SM","X.DMY", "var.g.DMY",
                           "her.DMY", "", "X.PH", "var.g.PH", "her.PH" ,"short",
                           "het.", "long")


# export to Latex


folder <- "E:/PhD/Manuscript/1st_chapter/Latex/table/"
file <- paste0(folder,"subset_partition.txt")

write(x = print(xtable(subset.part,digits = c(0,0,0,3,1,1,1,0,1,1,1,0,0,0),
                       caption = "Subsets of the population",
                       align = c("l","l","l","c","c","c","c","c","c","c","c","c",
                                 "c","c")),
                table.placement = NULL,
                include.rownames=FALSE, caption.placement = "top"), file=file)

################################################################################
############################# End Table 1 ######################################
################################################################################

# Table 2
#########

library(xtable)
library(ggplot2)


folder <- "F:/EU_NAM/EU_NAM_DENT/MPP_EUNAM/results2/QTL_analyses/"

# DMY

subset <- c(rep("short", 8), rep("hetero", 8), rep("long", 8))
VCOV <- rep(c(rep("h.err", 4), rep("cr.err", 4)), 3)
Q.eff <- rep(c("par", "anc", "biall", "MQE"), 6)

partition <- cbind(subset, Q.eff, VCOV)

res_DMY <- matrix("-", 6, 4)

i.ind <- rep(1:6, each = 4)
j.ind <- rep(1:4, time = 6)


for (i in 1:dim(partition)[1]){
  
  part <- partition[i, ]
  
  # load the data
  
  if(part[2] == "MQE"){
    
    res_DMY[i.ind[i], j.ind[i]] <- "-"
    
  } else {
    
    file <- paste0(folder, paste0("QTLan_", part[1], "_DMY_", part[2],
                                  "_", part[3]), "/QTL_genResults.txt")
    
    res <- tryCatch(read.table(file, row.names = 1), error = function(e) NULL)
    
    if(!is.null(res)){
      
      N.QTL <- round(res[1, ], 0)
      
      R2 <- round(res[3, ], 1)
      
      res_DMY[i.ind[i], j.ind[i]] <- paste(N.QTL, paste0("(", R2, ")"))
      
    }
    
  }
  
}

res_DMY

res_PH <- matrix("-", 6, 4)

i.ind <- rep(1:6, each = 4)
j.ind <- rep(1:4, time = 6)


for (i in 1:dim(partition)[1]){
  
  part <- partition[i, ]
  
  # load the data
  
  if(part[2] == "MQE"){
    
    res_PH[i.ind[i], j.ind[i]] <- "-"
    
  } else {
    
    file <- paste0(folder, paste0("QTLan_", part[1], "_PH_", part[2],
                                  "_", part[3]), "/QTL_genResults.txt")
    
    res <- tryCatch(read.table(file, row.names = 1), error = function(e) NULL)
    
    if(!is.null(res)){
      
      N.QTL <- round(res[1, ], 0)
      
      R2 <- round(res[3, ], 1)
      
      res_PH[i.ind[i], j.ind[i]] <- paste(N.QTL, paste0("(", R2, ")"))
      
    }
    
  }
  
}

res_PH

# for the whole table

Tab2 <- cbind(res_DMY, rep("", 6), res_PH)

################################################################################
############################# End Table 2 ######################################
################################################################################

# Figure 3
##########

library(xtable)
library(ggplot2)

# multiplot function
####################

# From Cookbook for R: "Graphs" Multiple graphs on one pages (ggplot2)

# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# plot A
########

# Average pTS and pVS over the three QTL effects per subset

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# short subset

pTS.short <- c()
pVS.short <- c()

for(i in 1:3){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.short <- c(pTS.short,c(t(pTS)))
  
  # pVS vector
  
  pVS.short <- c(pVS.short,c(t(pVS)))
  
}

# het subset

pTS.het <- c()
pVS.het <- c()

for(i in 4:6){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.het <- c(pTS.het,c(t(pTS)))
  
  # pVS vector
  
  pVS.het <- c(pVS.het,c(t(pVS)))
  
}

# long subset

pTS.lg <- c()
pVS.lg <- c()

for(i in 7:9){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.lg <- c(pTS.lg,c(t(pTS)))
  
  # pVS vector
  
  pVS.lg <- c(pVS.lg,c(t(pVS)))
  
}

av.pTS.sh_D <- mean(pTS.short)
av.pVS.sh_D <- mean(pVS.short)
se.pTS.sh_D <- sd(pTS.short)
se.pVS.sh_D <- sd(pVS.short)

av.pTS.het_D <- mean(pTS.het)
av.pVS.het_D <- mean(pVS.het)
se.pTS.het_D <- sd(pTS.het)
se.pVS.het_D <- sd(pVS.het)

av.pTS.lg_D <- mean(pTS.lg)
av.pVS.lg_D <- mean(pVS.lg)
se.pTS.lg_D <- sd(pTS.lg)
se.pVS.lg_D <- sd(pVS.lg)

#### PH

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"


# table A: Average pTS and pVS over the three QTL effects per subset

# short subset

pTS.short <- c()
pVS.short <- c()

for(i in 1:3){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.short <- c(pTS.short,c(t(pTS)))
  
  # pVS vector
  
  pVS.short <- c(pVS.short,c(t(pVS)))
  
}

# het subset

pTS.het <- c()
pVS.het <- c()

for(i in 4:6){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.het <- c(pTS.het,c(t(pTS)))
  
  # pVS vector
  
  pVS.het <- c(pVS.het,c(t(pVS)))
  
}

# long subset

pTS.lg <- c()
pVS.lg <- c()

for(i in 7:9){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.lg <- c(pTS.lg,c(t(pTS)))
  
  # pVS vector
  
  pVS.lg <- c(pVS.lg,c(t(pVS)))
  
}

av.pTS.sh_P <- mean(pTS.short)
av.pVS.sh_P <- mean(pVS.short)
se.pTS.sh_P <- sd(pTS.short)
se.pVS.sh_P <- sd(pVS.short)

av.pTS.het_P <- mean(pTS.het)
av.pVS.het_P <- mean(pVS.het)
se.pTS.het_P <- sd(pTS.het)
se.pVS.het_P <- sd(pVS.het)

av.pTS.lg_P <- mean(pTS.lg)
av.pVS.lg_P <- mean(pVS.lg)
se.pTS.lg_P <- sd(pTS.lg)
se.pVS.lg_P <- sd(pVS.lg)

# data.frame

pTS <- c(av.pTS.sh_P, av.pTS.het_P, av.pTS.lg_P,av.pTS.sh_D,
         av.pTS.het_D, av.pTS.lg_D)

se.pTS <- c(se.pTS.sh_P, se.pTS.het_P, se.pTS.lg_P, se.pTS.sh_D,
            se.pTS.het_D, se.pTS.lg_D)

pVS <- c(av.pVS.sh_P, av.pVS.het_P, av.pVS.lg_P,av.pVS.sh_D,
         av.pVS.het_D, av.pVS.lg_D)

se.pVS <- c(se.pVS.sh_P, se.pVS.het_P, se.pVS.lg_P, se.pVS.sh_D,
            se.pVS.het_D, se.pVS.lg_D)

# arange data in data.frame

data <- data.frame(rep(c("PH","DMY"), each=6), rep(c("short","het.","long"),4),
                   rep(c("pTS","pVS"), 2,each=3), c(pTS,pVS), c(se.pTS,se.pVS))

colnames(data) <- c("Trait","subset","set","var.ex","se")

data$Trait <- as.factor(data$Trait)
data$subset <- factor(data$subset, levels = c("short", "het.", "long"))
data$set <- as.factor(data$set)

# graph A
#########

pd <- position_dodge(width = 0.4)

pA <- ggplot(data, aes(x = subset, y = var.ex, group = set, colour = set,
                       shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group=set), position = pd) + # draw a line between point
  
  
  
  # modification of the legend elements
  
  ylab("Explained genetic variance (%)") +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        axis.text.x  = element_text(size=16),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size=16)) + 
  theme(legend.text = element_text(size=16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ Trait,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size=16, face="bold"))


#################### End graph A

# plot D
########

# Average pTS and pVS over the three different subset per QTL effects

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# parental effect

pTS.par <- c()
pVS.par <- c()

for(i in c(1,4,7)){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.par <- c(pTS.par,c(t(pTS)))
  
  # pVS vector
  
  pVS.par <- c(pVS.par,c(t(pVS)))
  
}

# ancestral effect

pTS.anc <- c()
pVS.anc <- c()

for(i in c(2,5,8)){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.anc <- c(pTS.anc,c(t(pTS)))
  
  # pVS vector
  
  pVS.anc <- c(pVS.anc,c(t(pVS)))
  
}

# bi-allelic effect

pTS.bi <- c()
pVS.bi <- c()

for(i in c(3,6,9)){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.bi <- c(pTS.bi,c(t(pTS)))
  
  # pVS vector
  
  pVS.bi <- c(pVS.bi,c(t(pVS)))
  
}

av.pTS.par_D <- mean(pTS.par)
av.pVS.par_D <- mean(pVS.par)
se.pTS.par_D <- sd(pTS.par)
se.pVS.par_D <- sd(pVS.par)

av.pTS.anc_D <- mean(pTS.anc)
av.pVS.anc_D <- mean(pVS.anc)
se.pTS.anc_D <- sd(pTS.anc)
se.pVS.anc_D <- sd(pVS.anc)

av.pTS.bi_D <- mean(pTS.bi)
av.pVS.bi_D <- mean(pVS.bi)
se.pTS.bi_D <- sd(pTS.bi)
se.pVS.bi_D <- sd(pVS.bi)

#### PH

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# parental effect

pTS.par <- c()
pVS.par <- c()

for(i in c(1,4,7)){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.par <- c(pTS.par,c(t(pTS)))
  
  # pVS vector
  
  pVS.par <- c(pVS.par,c(t(pVS)))
  
}

# ancestral effect

pTS.anc <- c()
pVS.anc <- c()

for(i in c(2,5,8)){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.anc <- c(pTS.anc,c(t(pTS)))
  
  # pVS vector
  
  pVS.anc <- c(pVS.anc,c(t(pVS)))
  
}

# bi-allelic effect

pTS.bi <- c()
pVS.bi <- c()

for(i in c(3,6,9)){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  
  pTS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_es.txt")))
  
  
  pVS <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                          part[2],"_",part[4],"_CV/p_ts.txt")))
  
  # pTS vector
  
  pTS.bi <- c(pTS.bi,c(t(pTS)))
  
  # pVS vector
  
  pVS.bi <- c(pVS.bi,c(t(pVS)))
  
}

av.pTS.par_P <- mean(pTS.par)
av.pVS.par_P <- mean(pVS.par)
se.pTS.par_P <- sd(pTS.par)
se.pVS.par_P <- sd(pVS.par)

av.pTS.anc_P <- mean(pTS.anc)
av.pVS.anc_P <- mean(pVS.anc)
se.pTS.anc_P <- sd(pTS.anc)
se.pVS.anc_P <- sd(pVS.anc)

av.pTS.bi_P <- mean(pTS.bi)
av.pVS.bi_P <- mean(pVS.bi)
se.pTS.bi_P <- sd(pTS.bi)
se.pVS.bi_P <- sd(pVS.bi)

# data.frame

var.ex <- c(av.pTS.par_P, av.pTS.anc_P, av.pTS.bi_P, av.pVS.par_P, av.pVS.anc_P,
            av.pVS.bi_P, av.pTS.par_D, av.pTS.anc_D, av.pTS.bi_D, av.pVS.par_D,
            av.pVS.anc_D, av.pVS.bi_D)

se <- c(se.pTS.par_P, se.pTS.anc_P, se.pTS.bi_P, se.pVS.par_P, se.pVS.anc_P,
        se.pVS.bi_P, se.pTS.par_D, se.pTS.anc_D, se.pTS.bi_D, se.pVS.par_D,
        se.pVS.anc_D, se.pVS.bi_D)


# arange data in data.frame

data <- data.frame(rep(c("PH","DMY"), each=6), rep(c("par","anc","bi"),4),
                   rep(c("pTS","pVS"), 2,each=3), var.ex, se)

colnames(data) <- c("Trait","Qeff","set","var.ex","se")

data$Trait <- as.factor(data$Trait)
data$Qeff <- factor(data$Qeff, levels = c("par", "anc", "bi"))
data$set <- as.factor(data$set)

# graph D
#########

pd <- position_dodge(width = 0.4)

pD <- ggplot(data, aes(x=Qeff, y=var.ex,group=set,colour=set,shape=set)) +
  
  # plot expl.var (y) by Qeff, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 2* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group=set), position = pd) + # draw a line between point
  
  # modification of the legend elements
  
  ylab("Explained genetic variance (%)") +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        axis.text.x  = element_text(size=16),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size=16)) + 
  theme(legend.text = element_text(size=16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ Trait,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size=16, face="bold"))

# graph B
#########

# the average number of QTL that are detected at least 20% of the time over the
# three type of QTL effect per subset.

# DMY
#####

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# short subset

Q20 <- c()
Qtot <- c()

for(i in 1:3){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # total number of QTLs
  
  Qtot <- c(Qtot, as.character(QTLs$Mk_names))
  
  # Number of QTL detected at least 20% of the time
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  Q20 <- c(Q20, as.character(QTLs$Mk_names))
  
}

Qtot.short_D <- length(unique(Qtot))
Q20.short_D <- length(unique(Q20))

# heterogeneous subset

Q20 <- c()
Qtot <- c()

for(i in 4:6){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # total number of QTLs
  
  Qtot <- c(Qtot, as.character(QTLs$Mk_names))
  
  # Number of QTL detected at least 20% of the time
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  Q20 <- c(Q20, as.character(QTLs$Mk_names))
  
}

Qtot.hetero_D <- length(unique(Qtot))
Q20.hetero_D <- length(unique(Q20))

# Long subset

Q20 <- c()
Qtot <- c()

for(i in 7:9){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # total number of QTLs
  
  Qtot <- c(Qtot, as.character(QTLs$Mk_names))
  
  # Number of QTL detected at least 20% of the time
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  Q20 <- c(Q20, as.character(QTLs$Mk_names))
  
}

Qtot.long_D <- length(unique(Qtot))
Q20.long_D <- length(unique(Q20))

# PH
####

# partition

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- data.frame(subset,Q.eff,trait,VCOV)

folder <- "~/MPP_EUNAM/results/CV/"

# short subset

Q20 <- c()
Qtot <- c()

for(i in 1:3){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # total number of QTLs
  
  Qtot <- c(Qtot, as.character(QTLs$Mk_names))
  
  # Number of QTL detected at least 20% of the time
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  Q20 <- c(Q20, as.character(QTLs$Mk_names))
  
}

Qtot.short_P <- length(unique(Qtot))
Q20.short_P <- length(unique(Q20))

# heterogeneous subset

Q20 <- c()
Qtot <- c()

for(i in 4:6){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # total number of QTLs
  
  Qtot <- c(Qtot, as.character(QTLs$Mk_names))
  
  # Number of QTL detected at least 20% of the time
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  Q20 <- c(Q20, as.character(QTLs$Mk_names))
  
}

Qtot.hetero_P <- length(unique(Qtot))
Q20.hetero_P <- length(unique(Q20))

# Long subset

Q20 <- c()
Qtot <- c()

for(i in 7:9){
  
  # load data
  
  # extract the ith partition
  part <- as.character(unlist(partition[i,]))
  
  # get the threshold value
  
  QTLs <- read.table(paste0(folder, paste0(part[1],"_",part[3],"_",
                                           part[2],"_",part[4],"_CV/QTL.txt")))
  
  # total number of QTLs
  
  Qtot <- c(Qtot, as.character(QTLs$Mk_names))
  
  # Number of QTL detected at least 20% of the time
  
  QTLs <- QTLs[QTLs$N>=10,]
  
  Q20 <- c(Q20, as.character(QTLs$Mk_names))
  
}

Qtot.long_P <- length(unique(Qtot))
Q20.long_P <- length(unique(Q20))

# compose the dataset

Qtot <- c(Qtot.short_P,Qtot.hetero_P,Qtot.long_P,
          Qtot.short_D,Qtot.hetero_D,Qtot.long_D)

Q20 <- c(Q20.short_P,Q20.hetero_P,Q20.long_P,
         Q20.short_D,Q20.hetero_D,Q20.long_D)

data <- data.frame(rep(c("PH","DMY"), each=3), rep(c("short","het.","long"),2),
                   Qtot, Q20)

colnames(data)[1:2] <- c("Trait","subset")

data$Trait <- as.factor(data$Trait)
data$subset <- factor(data$subset, levels = c("short", "het.", "long"))

# plot

# N.QTL tot

pC <- ggplot(data = data, aes(x=subset, y=Qtot, group=Trait,
                              colour=Trait, shape=Trait)) +
  
  geom_line() +
  
  geom_point(size=7) +
  
  theme(legend.title = element_text(size=14)) + 
  
  theme(legend.text = element_text(size=14)) +
  
  ylab("N QTL") +
  
  theme(axis.title.x = element_text(face="bold",size=16),
        axis.title.y = element_text(face="bold",size=18),
        axis.text.x  = element_text(size=16),
        axis.text = element_text(size = 18))

# graph B
#########

# plot

# Q 20

pB <- ggplot(data = data, aes(x=subset, y=Q20, group=Trait,
                              colour=Trait, shape=Trait)) +
  
  geom_line() +
  
  geom_point(size=7) +
  
  theme(legend.title = element_text(size=14)) + 
  
  theme(legend.text = element_text(size=14)) +
  
  ylab("N QTL >= 20%") +
  
  theme(axis.title.x = element_text(face="bold",size=16),
        axis.title.y = element_text(face="bold",size=18),
        axis.text.x  = element_text(size=16),
        axis.text = element_text(size = 18))

# combination all graphs


multiplot(pA, pC, pB, pD, cols=2)

# export the result

setwd("E:/PhD/Manuscript/1st_chapter/Latex/figures")

jpeg(filename = "CV_sum_res2.jpg",width = 1200,height = 900, quality = 100,
     res = 100)

multiplot(pA, pB, pC, pD, cols=2)

dev.off()

################################################################################
############################# End Figure 3 #####################################
################################################################################

# Table 2
#########

# DMY HRT
#########

folder <- "~/MPP_EUNAM/results/QTL_analyses/"

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- cbind(subset, Q.eff, trait, VCOV)

res_DMY_uerr <- matrix("-",3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# iteration does not goes through 4 and 5 because no QTL was detected for
# these configuration (DMY, heterogeneous, parental and ancestral)

for (i in c(1:3,6:9)){
  
  part <- partition[i,]
  
  # number of detected QTLs
  
  if(part[1]=="long") {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2],
                                  "_effect", "_CIM_", part[4],
                                  "/summary_results.csv")) 
    
  } else {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2], "_CIM_",
                                  part[4],"/summary_results.csv")) 
    
  }
  
  
  
  res <- read.csv(file, row.names = 1)
  
  N.QTL <- round(res[1, ], 0)
  
  R2 <- round(res[3, ], 1)
  
  res_DMY_uerr[i.ind[i],j.ind[i]] <- paste(N.QTL, paste0("(", R2, ")"))
  
  
}

res_DMY_uerr

# DMY CSRT
###########

folder <- "~/MPP_EUNAM/results/QTL_analyses/"

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("DMY"),times = 3,each = 3)
VCOV <- rep("cr.err",times = 9)

partition <- cbind(subset, Q.eff, trait, VCOV)

res_DMY_crerr <- matrix("-",3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)


for (i in 1:9){
  
  part <- partition[i,]
  
  # number of detected QTLs
  
  if(part[1]=="long") {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2],
                                  "_effect", "_CIM_", part[4],
                                  "/summary_results.csv")) 
    
  } else {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2], "_CIM_",
                                  part[4],"/summary_results.csv")) 
    
  }
  
  
  
  res <- read.csv(file, row.names = 1)
  
  N.QTL <- round(res[1, ], 0)
  
  R2 <- round(res[3, ], 1)
  
  res_DMY_crerr[i.ind[i],j.ind[i]] <- paste(N.QTL, paste0("(", R2, ")"))
  
  
}


# results of the MQE model
##########################

subset <- c("short", "hetero", "long")

folder <- "~/MPP_EUNAM/results/MQeff/"

MQE_DMY <- matrix("",3,3)


for(i in 1:3){
  
  if(subset[i] == "long"){
    
    
    QTL.file <- paste0(folder, subset[i],"_DMY_mQeff_proc_u.err","/table.QTL.txt")
    
    QTL <- read.table(QTL.file, h=T)
    n.QTL <- dim(QTL)[1]
    detail <- table(QTL$Q.eff)[c(3,1,2)]
    detail <- paste0("(", paste(detail, collapse = "/"), ")")
    
    # R squared
    
    QTL_report <- read.table(paste0(folder, subset[i],"_DMY_mQeff_proc_u.err",
                                    "/QTL_REPORT.txt"), fill = T)
    
    R2 <- round(as.numeric(as.character(QTL_report[4, 4])), 1)
    
    R2 <- paste0("(", R2, ")")
    
    MQE_DMY[i,2] <- paste(n.QTL, detail, R2)
    
  } else {
    
    QTL.file <- paste0(folder, subset[i],"_DMY_mQeff_proc_u.err","/QTL.txt")
    
    QTL <- read.table(QTL.file, h=T)
    n.QTL <- dim(QTL)[1]
    detail <- table(QTL$QTL.inc)[c(3,1,2)]
    detail <- paste0("(", paste(detail, collapse = "/"), ")")
    
    # R squared
    
    QTL_report <- read.table(paste0(folder, subset[i],"_DMY_mQeff_proc_u.err",
                                    "/QTL_REPORT.txt"), fill = T)
    
    R2 <- round(as.numeric(as.character(QTL_report[4, 4])), 1)
    
    R2 <- paste0("(", R2, ")")
    
    MQE_DMY[i,2] <- paste(n.QTL, detail, R2)
    
  }
  
  
}

MQE_DMY

# combine all results together

DMY_res <- rbind(res_DMY_uerr[1, ], res_DMY_crerr[1, ], MQE_DMY[1, ],
                 res_DMY_uerr[2, ], res_DMY_crerr[2, ], MQE_DMY[2, ],
                 res_DMY_uerr[3, ], res_DMY_crerr[3, ], MQE_DMY[3, ])

# PH
####

# PH HRT
#########

folder <- "~/MPP_EUNAM/results/QTL_analyses/"

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 3,each = 3)
VCOV <- rep("u.err",times = 9)

partition <- cbind(subset, Q.eff, trait, VCOV)

res_PH_uerr <- matrix("-",3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)

# iteration does not goes through 4 and 5 because no QTL was detected for
# these configuration (PH, heterogeneous, parental and ancestral)

for (i in 1:9){
  
  part <- partition[i,]
  
  # number of detected QTLs
  
  if(part[1]=="long") {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2],
                                  "_effect", "_CIM_", part[4],
                                  "/summary_results.csv")) 
    
  } else {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2], "_CIM_",
                                  part[4],"/summary_results.csv")) 
    
  }
  
  
  
  res <- read.csv(file, row.names = 1)
  
  N.QTL <- round(res[1, ], 0)
  
  R2 <- round(res[3, ], 1)
  
  res_PH_uerr[i.ind[i],j.ind[i]] <- paste(N.QTL, paste0("(", R2, ")"))
  
  
}

res_PH_uerr

subset <- rep(c("short","hetero","long"),each=3)
Q.eff <- rep(c("parental","ancestral","biallelic"),times = 3,each = 1)
trait <- rep(c("PH"),times = 3,each = 3)
VCOV <- rep("cr.err",times = 9)

partition <- cbind(subset, Q.eff, trait, VCOV)

res_PH_crerr <- matrix("-",3,3)

i.ind <- rep(1:3,each=3)
j.ind <- rep(1:3,time=3)


for (i in 1:9){
  
  part <- partition[i,]
  
  # number of detected QTLs
  
  if(part[1]=="long") {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2],
                                  "_effect", "_CIM_", part[4],
                                  "/summary_results.csv")) 
    
  } else {
    
    file <- paste0(folder, paste0(part[1], "_", part[3], "_", part[2], "_CIM_",
                                  part[4],"/summary_results.csv")) 
    
  }
  
  
  
  res <- read.csv(file, row.names = 1)
  
  N.QTL <- round(res[1, ], 0)
  
  R2 <- round(res[3, ], 1)
  
  res_PH_crerr[i.ind[i],j.ind[i]] <- paste(N.QTL, paste0("(", R2, ")"))
  
  
}


# results of the MQE model
##########################

subset <- c("short", "hetero", "long")

folder <- "~/MPP_EUNAM/results/MQeff/"

MQE_PH <- matrix("",3,3)


for(i in 1:3){
  
  if(subset[i] == "long"){
    
    
    QTL.file <- paste0(folder, subset[i],"_PH_mQeff_proc_u.err","/table.QTL.txt")
    
    QTL <- read.table(QTL.file, h=T)
    n.QTL <- dim(QTL)[1]
    detail <- table(QTL$Q.eff)[c(3,1,2)]
    detail <- paste0("(", paste(detail, collapse = "/"), ")")
    
    # R squared
    
    QTL_report <- read.table(paste0(folder, subset[i],"_PH_mQeff_proc_u.err",
                                    "/QTL_REPORT.txt"), fill = T)
    
    R2 <- round(as.numeric(as.character(QTL_report[4, 4])), 1)
    
    R2 <- paste0("(", R2, ")")
    
    MQE_PH[i,2] <- paste(n.QTL, detail, R2)
    
  } else {
    
    QTL.file <- paste0(folder, subset[i],"_PH_mQeff_proc_u.err","/QTL.txt")
    
    QTL <- read.table(QTL.file, h=T)
    n.QTL <- dim(QTL)[1]
    detail <- table(QTL$QTL.inc)[c(3,1,2)]
    detail <- paste0("(", paste(detail, collapse = "/"), ")")
    
    # R squared
    
    QTL_report <- read.table(paste0(folder, subset[i],"_PH_mQeff_proc_u.err",
                                    "/QTL_REPORT.txt"), fill = T)
    
    R2 <- round(as.numeric(as.character(QTL_report[4, 4])), 1)
    
    R2 <- paste0("(", R2, ")")
    
    MQE_PH[i,2] <- paste(n.QTL, detail, R2)
    
  }
  
  
}

MQE_PH

# combine all results together

PH_res <- rbind(res_PH_uerr[1, ], res_PH_crerr[1, ], MQE_PH[1, ],
                res_PH_uerr[2, ], res_PH_crerr[2, ], MQE_PH[2, ],
                res_PH_uerr[3, ], res_PH_crerr[3, ], MQE_PH[3, ])

# combine both traits togheter

res_QTL_analysis <- cbind(DMY_res, rep("", 3), PH_res)

# add row and columns names

model <- rep(c("HRT", "CSRT", "MQE"), times=3)
subset <- c("", "short", "", "", "hetero", "", "", "long", "")

res_QTL_analysis <- cbind(subset, model, res_QTL_analysis)

colnames(res_QTL_analysis) <- c("\n", "\n", c("parental","ancestral","bi-allelic"),
                                "\n", c("parental","ancestral","bi-allelic"))

# export xtable format
######################

library(xtable)

folder <- "F:/PhD/Manuscript/1st_chapter/Latex/table/"
file=paste0(folder,"QTL_analyses.txt")

write(x = print(xtable(res_QTL_analysis,digits = c(0, 0, 0, 1,  1, 1, 1, 1, 1, 1),
                       caption = "QTL detection results full subsets",
                       align = c("l","l","l", "c", "c", "c","c","c","c","c")), 
                table.placement = NULL, hline.after = c(0, 3, 6, 9), 
                include.rownames=FALSE, include.colnames = TRUE, caption.placement = "top"),
      file=file)

################################################################################
############################# End Table 2 #####################################
################################################################################

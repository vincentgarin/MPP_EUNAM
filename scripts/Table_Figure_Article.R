#############################
# table and figures article #
#############################

# Table 1
#########

library(mppR)
library(xtable)

path <- "~/MPP_EUNAM/"

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


folder <- "~/MPP_EUNAM/results/QTL_analyses/"
folder.MQE <- "~/MPP_EUNAM/results/MQE/"

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
    
    file <- paste0(folder.MQE, paste0("MQE_", part[1], "_DMY_", part[3]),
                   "/QTL_genResults.txt")
    
    file2 <- paste0(folder.MQE, paste0("MQE_", part[1], "_DMY_", part[3]),
                    "/QTL.txt")
    
    res <- tryCatch(read.table(file, row.names = 1), error = function(e) NULL)
    
    QTL <- tryCatch(read.table(file2, row.names = 1), error = function(e) NULL)
    
    if(!is.null(res)){
      
      N.QTL <- round(res[1, ], 0)
      
      R2 <- round(res[3, ], 1)
      
      n.par <- sum(QTL[, 4] == "par"); if(n.par == 0) {n.par <- "-"}
      n.anc <- sum(QTL[, 4] == "anc"); if(n.anc == 0) {n.anc <- "-"}
      n.biall <- sum(QTL[, 4] == "biall"); if(n.biall == 0) {n.biall <- "-"}
      
      QTL.type <- paste0("(", n.par, "/", n.anc, "/", n.biall,")")
      
      res_DMY[i.ind[i], j.ind[i]] <- paste(N.QTL, QTL.type, paste0("(", R2, ")"))
      
    }
    
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


for (i in 1:dim(partition)[1]){
  
  part <- partition[i, ]
  
  # load the data
  
  if(part[2] == "MQE"){
    
    file <- paste0(folder.MQE, paste0("MQE_", part[1], "_PH_", part[3]),
                   "/QTL_genResults.txt")
    
    file2 <- paste0(folder.MQE, paste0("MQE_", part[1], "_PH_", part[3]),
                    "/QTL.txt")
    
    res <- tryCatch(read.table(file, row.names = 1), error = function(e) NULL)
    
    QTL <- tryCatch(read.table(file2, row.names = 1), error = function(e) NULL)
    
    if(!is.null(res)){
      
      N.QTL <- round(res[1, ], 0)
      
      R2 <- round(res[3, ], 1)
      
      n.par <- sum(QTL[, 4] == "par"); if(n.par == 0) {n.par <- "-"}
      n.anc <- sum(QTL[, 4] == "anc"); if(n.anc == 0) {n.anc <- "-"}
      n.biall <- sum(QTL[, 4] == "biall"); if(n.biall == 0) {n.biall <- "-"}
      
      QTL.type <- paste0("(", n.par, "/", n.anc, "/", n.biall,")")
      
      res_PH[i.ind[i], j.ind[i]] <- paste(N.QTL, QTL.type, paste0("(", R2, ")"))
      
    }
    
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
colnames(Tab2) <- c("parental", "ancestral", "bi-allelic", "MQE", " ",
                    "parental", "ancestral", "bi-allelic", "MQE")

file.name <- paste0("E:/PhD/Manuscript/1st_chapter/Latex/table/Tab2.txt")

# export table 2

# write(x = print(xtable(Tab2, digits = c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1),
#                        caption = "Tab2",
#                        align = c("l", rep("c", 9))), 
#                 table.placement = NULL, hline.after = c(0, 2, 4, 6), 
#                 include.rownames = TRUE, include.colnames = TRUE,
#                 caption.placement = "top"),
#       file = file.name)

################################################################################
############################# End Table 2 ######################################
################################################################################

# Figure 2
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

# short DMY
###########

setwd("~/MPP_EUNAM/results/")

res_info <- data.frame(rep("DMY", 16), rep("short", 16),
                       rep(c("h.err", "cr.err_fast"), each = 8), 
                       rep(c("par", "anc", "biall", "MQE"), time = 4),
                       rep(c("pTS", "pVS"), each = 4), stringsAsFactors = FALSE)

res_info <- res_info[c(1:4, 9:12, 5:8, 13:16), ]

res_tab <- matrix(0, 16, 2)

colnames(res_info) <- c("trait", "subset", "VCOV", "Q.eff", "set")

res <- data.frame(res_info, res_tab, stringsAsFactors = FALSE)

for(i in 1:(dim(res_info)[1]/2)){
  
  part <- res_info[i, ]
  part <- part[-5]
  part <- part[c(2, 1, 4, 3)]
  
  if (part[3] == "MQE"){
    
    file.name <- paste0("./MQE_CV/MQE_CV_", paste(part[-3], collapse = "_"),
                        "/CV_res.txt")
    
  } else {
    
    file.name <- paste0("./CV/CV_", paste(part, collapse = "_"), "/CV_res.txt")
    
  }
  
  res.i <- read.table(file.name, h = TRUE)
  
  # pTS results
  
  res[i, 6] <- mean(res.i[, 2], na.rm = TRUE)
  res[i, 7] <- sd(res.i[, 2], na.rm = TRUE)
    
  # pVS results
    
  res[(i + 8), 6] <- mean(res.i[, 3], na.rm = TRUE)
  res[(i + 8), 7] <- sd(res.i[, 3], na.rm = TRUE)
  
}

# format data.frame

data <- res
colnames(data)[6:7] <- c("var.ex", "se")
data$VCOV <- factor(rep(c("HRT", "HRT", "HRT", "HRT", "CSRT",
                          "CSRT", "CSRT", "CSRT"), 2), levels = c("HRT", "CSRT"))
data$Q.eff[data$Q.eff == "biall"] <- "bi"
data$Q.eff <- factor(data$Q.eff, levels = c("par", "anc", "bi", "MQE"))
data$set <- as.factor(data$set)

# plot

pd <- position_dodge(width = 0.4)

short_DMY <- ggplot(data, aes(x = Q.eff, y = var.ex, group = set, colour = set,
                       shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group = set), position = pd) + # draw a line between point
  
  
  # modification of the legend elements
  
  xlab("QTL effect") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ VCOV,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size = 16, face = "bold")) +
  
  theme(legend.position = "none")

short_DMY


# short PH
##########

setwd("~/MPP_EUNAM/results/")

res_info <- data.frame(rep("PH", 16), rep("short", 16),
                       rep(c("h.err", "cr.err_fast"), each = 8), 
                       rep(c("par", "anc", "biall", "MQE"), time = 4),
                       rep(c("pTS", "pVS"), each = 4), stringsAsFactors = FALSE)

res_info <- res_info[c(1:4, 9:12, 5:8, 13:16), ]

res_tab <- matrix(0, 16, 2)

colnames(res_info) <- c("trait", "subset", "VCOV", "Q.eff", "set")

res <- data.frame(res_info, res_tab, stringsAsFactors = FALSE)

for(i in 1:(dim(res_info)[1]/2)){
  
  part <- res_info[i, ]
  part <- part[-5]
  part <- part[c(2, 1, 4, 3)]
  
  if (part[3] == "MQE"){
    
    file.name <- paste0("./MQE_CV/MQE_CV_", paste(part[-3], collapse = "_"),
                        "/CV_res.txt")
    
  } else {
    
    file.name <- paste0("./CV/CV_", paste(part, collapse = "_"), "/CV_res.txt")
    
  }
  
  res.i <- read.table(file.name, h = TRUE)
  
  # pTS results
  
  res[i, 6] <- mean(res.i[, 2], na.rm = TRUE)
  res[i, 7] <- sd(res.i[, 2], na.rm = TRUE)
  
  # pVS results
  
  res[(i + 8), 6] <- mean(res.i[, 3], na.rm = TRUE)
  res[(i + 8), 7] <- sd(res.i[, 3], na.rm = TRUE)
  
}

# format data.frame

data <- res
colnames(data)[6:7] <- c("var.ex", "se")
data$VCOV <- factor(rep(c("HRT", "HRT", "HRT", "HRT", "CSRT",
                          "CSRT", "CSRT", "CSRT"), 2), levels = c("HRT", "CSRT"))
data$Q.eff[data$Q.eff == "biall"] <- "bi"
data$Q.eff <- factor(data$Q.eff, levels = c("par", "anc", "bi", "MQE"))
data$set <- as.factor(data$set)

# plot

pd <- position_dodge(width = 0.4)

short_PH <- ggplot(data, aes(x = Q.eff, y = var.ex, group = set, colour = set,
                              shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group = set), position = pd) + # draw a line between point
  
  
  # modification of the legend elements
  
  xlab("QTL effect") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ VCOV,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size = 16, face = "bold")) +
  
  theme(legend.position = "none")

short_PH

# hetero DMY
############

setwd("~/MPP_EUNAM/results/")

res_info <- data.frame(rep("DMY", 16), rep("hetero", 16),
                       rep(c("h.err", "cr.err_fast"), each = 8), 
                       rep(c("par", "anc", "biall", "MQE"), time = 4),
                       rep(c("pTS", "pVS"), each = 4), stringsAsFactors = FALSE)

res_info <- res_info[c(1:4, 9:12, 5:8, 13:16), ]

res_tab <- matrix(0, 16, 2)

colnames(res_info) <- c("trait", "subset", "VCOV", "Q.eff", "set")

res <- data.frame(res_info, res_tab, stringsAsFactors = FALSE)

for(i in 1:(dim(res_info)[1]/2)){
  
  part <- res_info[i, ]
  part <- part[-5]
  part <- part[c(2, 1, 4, 3)]
  
  if (part[3] == "MQE"){
    
    file.name <- paste0("./MQE_CV/MQE_CV_", paste(part[-3], collapse = "_"),
                        "/CV_res.txt")
    
  } else {
    
    file.name <- paste0("./CV/CV_", paste(part, collapse = "_"), "/CV_res.txt")
    
  }
  
  res.i <- read.table(file.name, h = TRUE)
  
  # pTS results
  
  res[i, 6] <- mean(res.i[, 2], na.rm = TRUE)
  res[i, 7] <- sd(res.i[, 2], na.rm = TRUE)
  
  # pVS results
  
  res[(i + 8), 6] <- mean(res.i[, 3], na.rm = TRUE)
  res[(i + 8), 7] <- sd(res.i[, 3], na.rm = TRUE)
  
}

# format data.frame

data <- res
colnames(data)[6:7] <- c("var.ex", "se")
data$VCOV <- factor(rep(c("HRT", "HRT", "HRT", "HRT", "CSRT",
                          "CSRT", "CSRT", "CSRT"), 2), levels = c("HRT", "CSRT"))
data$Q.eff[data$Q.eff == "biall"] <- "bi"
data$Q.eff <- factor(data$Q.eff, levels = c("par", "anc", "bi", "MQE"))
data$set <- as.factor(data$set)

# plot

pd <- position_dodge(width = 0.4)

hetero_DMY <- ggplot(data, aes(x = Q.eff, y = var.ex, group = set, colour = set,
                             shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group = set), position = pd) + # draw a line between point
  
  
  # modification of the legend elements
  
  xlab("QTL effect") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ VCOV,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size = 16, face = "bold")) +
  
  theme(legend.position = "none")

hetero_DMY

# hetero PH
############

setwd("~/MPP_EUNAM/results/")

res_info <- data.frame(rep("PH", 16), rep("hetero", 16),
                       rep(c("h.err", "cr.err_fast"), each = 8), 
                       rep(c("par", "anc", "biall", "MQE"), time = 4),
                       rep(c("pTS", "pVS"), each = 4), stringsAsFactors = FALSE)

res_info <- res_info[c(1:4, 9:12, 5:8, 13:16), ]

colnames(res_info) <- c("trait", "subset", "VCOV", "Q.eff", "set")

res <- data.frame(res_info, res_tab, stringsAsFactors = FALSE)

for(i in 1:(dim(res_info)[1]/2)){
  
  part <- res_info[i, ]
  part <- part[-5]
  part <- part[c(2, 1, 4, 3)]
  
  if (part[3] == "MQE"){
    
    file.name <- paste0("./MQE_CV/MQE_CV_", paste(part[-3], collapse = "_"),
                        "/CV_res.txt")
    
  } else {
    
    file.name <- paste0("./CV/CV_", paste(part, collapse = "_"), "/CV_res.txt")
    
  }
  
  res.i <- read.table(file.name, h = TRUE)
  
  # pTS results
  
  res[i, 6] <- mean(res.i[, 2], na.rm = TRUE)
  res[i, 7] <- sd(res.i[, 2], na.rm = TRUE)
  
  # pVS results
  
  res[(i + 8), 6] <- mean(res.i[, 3], na.rm = TRUE)
  res[(i + 8), 7] <- sd(res.i[, 3], na.rm = TRUE)
  
}

# format data.frame

data <- res
colnames(data)[6:7] <- c("var.ex", "se")
data$VCOV <- factor(rep(c("HRT", "HRT", "HRT", "HRT", "CSRT",
                          "CSRT", "CSRT", "CSRT"), 2), levels = c("HRT", "CSRT"))
data$Q.eff[data$Q.eff == "biall"] <- "bi"
data$Q.eff <- factor(data$Q.eff, levels = c("par", "anc", "bi", "MQE"))
data$set <- as.factor(data$set)

# plot

pd <- position_dodge(width = 0.4)

hetero_PH <- ggplot(data, aes(x = Q.eff, y = var.ex, group = set, colour = set,
                               shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group = set), position = pd) + # draw a line between point
  
  
  # modification of the legend elements
  
  xlab("QTL effect") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ VCOV,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size = 16, face = "bold")) +
  
  theme(legend.position = "none")


hetero_PH

# long DMY
##########

setwd("~/MPP_EUNAM/results/")

res_info <- data.frame(rep("DMY", 16), rep("long", 16),
                       rep(c("h.err", "cr.err_fast"), each = 8), 
                       rep(c("par", "anc", "biall", "MQE"), time = 4),
                       rep(c("pTS", "pVS"), each = 4), stringsAsFactors = FALSE)

res_info <- res_info[c(1:4, 9:12, 5:8, 13:16), ]

colnames(res_info) <- c("trait", "subset", "VCOV", "Q.eff", "set")

res <- data.frame(res_info, res_tab, stringsAsFactors = FALSE)

for(i in 1:(dim(res_info)[1]/2)){
  
  part <- res_info[i, ]
  part <- part[-5]
  part <- part[c(2, 1, 4, 3)]
  
  if (part[3] == "MQE"){
    
    file.name <- paste0("./MQE_CV/MQE_CV_", paste(part[-3], collapse = "_"),
                        "/CV_res.txt")
    
  } else {
    
    file.name <- paste0("./CV/CV_", paste(part, collapse = "_"), "/CV_res.txt")
    
  }
  
  res.i <- read.table(file.name, h = TRUE)
  
  # pTS results
  
  res[i, 6] <- mean(res.i[, 2], na.rm = TRUE)
  res[i, 7] <- sd(res.i[, 2], na.rm = TRUE)
  
  # pVS results
  
  res[(i + 8), 6] <- mean(res.i[, 3], na.rm = TRUE)
  res[(i + 8), 7] <- sd(res.i[, 3], na.rm = TRUE)
  
}

# format data.frame

data <- res
colnames(data)[6:7] <- c("var.ex", "se")
data$VCOV <- factor(rep(c("HRT", "HRT", "HRT", "HRT", "CSRT",
                          "CSRT", "CSRT", "CSRT"), 2), levels = c("HRT", "CSRT"))
data$Q.eff[data$Q.eff == "biall"] <- "bi"
data$Q.eff <- factor(data$Q.eff, levels = c("par", "anc", "bi", "MQE"))
data$set <- as.factor(data$set)

# plot

pd <- position_dodge(width = 0.4)

long_DMY <- ggplot(data, aes(x = Q.eff, y = var.ex, group = set, colour = set,
                              shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group = set), position = pd) + # draw a line between point
  
  
  # modification of the legend elements
  
  xlab("QTL effect") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ VCOV,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size = 16, face = "bold")) +
  
  theme(legend.position = "none")

long_DMY

# long PH
##########

setwd("~/MPP_EUNAM/results/")

res_info <- data.frame(rep("PH", 16), rep("long", 16),
                       rep(c("h.err", "cr.err_fast"), each = 8), 
                       rep(c("par", "anc", "biall", "MQE"), time = 4),
                       rep(c("pTS", "pVS"), each = 4), stringsAsFactors = FALSE)

res_info <- res_info[c(1:4, 9:12, 5:8, 13:16), ]

colnames(res_info) <- c("trait", "subset", "VCOV", "Q.eff", "set")

res <- data.frame(res_info, res_tab, stringsAsFactors = FALSE)

for(i in 1:(dim(res_info)[1]/2)){
  
  part <- res_info[i, ]
  part <- part[-5]
  part <- part[c(2, 1, 4, 3)]
  
  if (part[3] == "MQE"){
    
    file.name <- paste0("./MQE_CV/MQE_CV_", paste(part[-3], collapse = "_"),
                        "/CV_res.txt")
    
  } else {
    
    file.name <- paste0("./CV/CV_", paste(part, collapse = "_"), "/CV_res.txt")
    
  }
  
  res.i <- read.table(file.name, h = TRUE)
  
  # pTS results
  
  res[i, 6] <- mean(res.i[, 2], na.rm = TRUE)
  res[i, 7] <- sd(res.i[, 2], na.rm = TRUE)
  
  # pVS results
  
  res[(i + 8), 6] <- mean(res.i[, 3], na.rm = TRUE)
  res[(i + 8), 7] <- sd(res.i[, 3], na.rm = TRUE)
  
}

# format data.frame

data <- res
colnames(data)[6:7] <- c("var.ex", "se")
data$VCOV <- factor(rep(c("HRT", "HRT", "HRT", "HRT", "CSRT",
                          "CSRT", "CSRT", "CSRT"), 2), levels = c("HRT", "CSRT"))
data$Q.eff[data$Q.eff == "biall"] <- "bi"
data$Q.eff <- factor(data$Q.eff, levels = c("par", "anc", "bi", "MQE"))
data$set <- as.factor(data$set)

# plot

pd <- position_dodge(width = 0.4)

long_PH <- ggplot(data, aes(x = Q.eff, y = var.ex, group = set, colour = set,
                             shape = set)) +
  
  # plot expl.var (y) by subset, differentiate two lines by colour and shape
  # (pTS and pVS)
  
  geom_errorbar(aes(ymin = var.ex-(2*se), ymax = var.ex+(2*se)), width = .2,
                position = pd) +
  
  # add error bar : value +- 1* standard error
  
  geom_point(size = 6, position = pd) + # draw point and specify size of the points
  geom_line(aes(group = set), position = pd) + # draw a line between point
  
  
  # modification of the legend elements
  
  xlab("QTL effect") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 16)) +
  
  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) +
  
  # divide the graph in several panels
  
  facet_grid(. ~ VCOV,) +
  
  # legend of the panels
  
  theme(strip.text.x = element_text(size = 16, face = "bold")) +

  theme(legend.position = "none")
  
long_PH


multiplot(short_DMY, hetero_DMY, long_DMY, short_PH, hetero_PH, long_PH, cols = 2)

# export the result

setwd("E:/PhD/Manuscript/1st_chapter/Latex/figures")

jpeg(filename = "CV_res.jpg",width = 1200, height = 900, quality = 100,
     res = 100)

multiplot(short_DMY, hetero_DMY, long_DMY, short_PH, hetero_PH, long_PH, cols = 2)

dev.off()

################################################################################
############################# End Table 2 ######################################
################################################################################
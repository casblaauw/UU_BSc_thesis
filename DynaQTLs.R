### Import, reshape and explore starting data --------------------------------------------------------------------------------------------
## Libraries ----
library(openxlsx)
library(ggplot2)
library(gplots)
library(cowplot)
library(tidyverse)
library(forcats)
library(reshape2)
library(randomForest)
library(randomForestExplainer)

## Needed data --------------------------------------------------------------------------------
# Random forest results
load("RF_results2/data_expvar_perm.out")
load("RF_results2/data_expvar_dev.out")
load("RF_results2/data_expvar_hs.out")
load("RF_results2/data_expvar_rec.out")

load("RF_results2/data_incmse_perm.out")
load("RF_results2/data_incmse_dev.out")
load("RF_results2/data_incmse_hs.out")
load("RF_results2/data_incmse_rec.out")

load("RF_results2/data_pur_perm.out")
load("RF_results2/data_pur_dev.out")
load("RF_results2/data_pur_hs.out")
load("RF_results2/data_pur_rec.out")

# Data on the expressed genes
WormGenes <- read.csv("WormGenes_all.csv", header=FALSE, row.names=1)
colnames(WormGenes) <- c("name", "start", "end", "chr")

# marker info
mrk.info <- read.delim("marker.txt",sep="\t")
head(mrk.info)

# genetic map
gm.dev <- read.delim("map_CT.txt",sep="\t")
head(gm.dev)
gm.hs <- read.delim("map_HS.txt",sep="\t")
head(gm.hs)
gm.rec <- read.delim("map_REC.txt", sep="\t")
head(gm.rec)

# gene expression
gexp.dev <- read.delim("Snoek_Sterken_CT_gene_expression_RILs.txt",sep="\t")
head(gexp.dev)
gexp.hs <- read.delim("Snoek_Sterken_HS_gene_expression_RILs.txt",sep="\t")
head(gexp.hs)
gexp.rec <- read.delim("Snoek_Sterken_REC_gene_expression_RILs.txt",sep="\t")
head(gexp.rec)


# projection axis 
load("projection_alldata.rdata")
head(projection_all)
projection_all <- data.frame(projection_all)

## Reshape data
    # select the specific RILs only
    proj.dev <- projection_all[projection_all$experimentype == "Dev" &
                                 projection_all$line == "RIL" & 
                                 projection_all$genotype %in% colnames(gm.dev),]
    proj.hs <- projection_all[projection_all$experimentype == "HS" &
                                projection_all$line == "RIL" & 
                                projection_all$genotype %in% colnames(gm.hs),]
    proj.rec <- projection_all[projection_all$experimentype == "Rec" &
                                 projection_all$line == "RIL" & 
                                 projection_all$genotype %in% colnames(gm.rec),]
    dim(proj.dev); dim (proj.hs); dim(proj.rec)
    # remove duplicated rils (select first entry....)
    proj.dev <- proj.dev[!duplicated(proj.dev$genotype),]
    proj.hs <- proj.hs[!duplicated(proj.hs$genotype),]
    proj.rec <- proj.rec[!duplicated(proj.rec$genotype),]
    
    # set ril names as rownames for easier selection and ordering
    rownames(proj.dev) <- proj.dev$genotype
    head(proj.dev)
    rownames(proj.hs) <- proj.hs$genotype
    head(proj.hs)
    rownames(proj.rec) <- proj.rec$genotype
    head(proj.rec)
    
    # get the genetic map and proj.axis in the same order 
        # first check if all RILs are in gm and proj
    colnames(gm.dev) %in% rownames(proj.dev)
    rownames(proj.dev) %in% colnames(gm.dev)
    
    colnames(gm.hs) %in% rownames(proj.hs)
    rownames(proj.hs) %in% colnames(gm.hs)
    
    colnames(gm.rec) %in% rownames(proj.rec)
    rownames(proj.rec) %in% colnames(gm.rec)
    
    # adjust gm to only have the RILS in proj
    gm.dev.ril <- gm.dev[, colnames(gm.dev) %in% rownames(proj.dev)]
    proj.dev <- proj.dev[colnames(gm.dev.ril),]
    colnames(gm.dev.ril) == rownames(proj.dev)
    head(as.matrix(proj.dev))
    
    gm.hs.ril <- gm.hs[, colnames(gm.hs) %in% rownames(proj.hs)]
    proj.hs <- proj.hs[colnames(gm.hs.ril),]
    colnames(gm.hs.ril) == rownames(proj.hs)
    head(as.matrix(proj.hs))
    
    gm.rec.ril <- gm.rec[, colnames(gm.rec) %in% rownames(proj.rec)]
    proj.rec <- proj.rec[colnames(gm.rec.ril),]
    colnames(gm.rec.ril) == rownames(proj.rec)
    head(as.matrix(proj.rec))

# remove the redundant data from the genetic map
    # many markers have no differences due to the limited size of the population
    gm.dev.nd <- as.matrix(gm.dev.ril[!duplicated(gm.dev.ril),] )
    head(gm.dev.nd)
    
    gm.hs.nd <- as.matrix(gm.hs.ril[!duplicated(gm.hs.ril),] )
    head(gm.hs.nd)
    
    gm.rec.nd <- as.matrix(gm.rec.ril[!duplicated(gm.rec.ril),] )
    head(gm.rec.nd)
    
# prepare gene expression data
    gexp.dev.nd <- gexp.dev[,colnames(gm.dev.nd)]
    colnames(gexp.dev.nd) == colnames(gm.dev.nd)
    
    gexp.hs.nd <- gexp.hs[,colnames(gm.hs.nd)]
    colnames(gexp.hs.nd) == colnames(gm.hs.nd)
    
    gexp.rec.nd <- gexp.rec[,colnames(gm.rec.nd)]
    colnames(gexp.rec.nd) == colnames(gm.rec.nd)

# Prepare data for in random forest
    use.gm.dev <- (t(gm.dev.nd)) ; dim(use.gm.dev)
    prax.dev <- proj.dev$devprojection ; length(prax.dev)

  
    use.gm.hs <- (t(gm.hs.nd)) ; dim(use.gm.hs)
    prax.hs <- proj.hs$devprojection ; length(prax.hs)

    
    use.gm.rec <- (t(gm.rec.nd)) ; dim(use.gm.rec)
    prax.rec <- proj.rec$devprojection ; length(prax.rec)

    
    
## Exploratory plots -----------------------------------------------------------
    #Distributions of total developing worms and developing RILs among development axis
    # Prepare the data
    ril.all <- projection_all
    ril.all <- ril.all[!duplicated(ril.all$genotype),]
    ril.all[ril.all == "NIL" | ril.all == "PL"] <- "Other"
    
    ril.dev <- projection_all[projection_all$experimentype == "Dev",]
    ril.dev <- ril.dev[!duplicated(ril.dev$genotype),]
    ril.dev[ril.dev == "NIL" | ril.dev == "PL"] <- "Other"
    
    ril.hs <- projection_all[projection_all$experimentype == "HS",]
    ril.hs <- ril.hs[!duplicated(ril.hs$genotype),]
    ril.hs[ril.hs == "NIL" | ril.hs == "PL"] <- "Other"
    
    ril.rec <- projection_all[projection_all$experimentype == "Rec",]
    ril.rec <- ril.rec[!duplicated(ril.rec$genotype),]
    ril.rec[ril.rec == "NIL" | ril.rec == "PL"] <- "Other"
    
    
    ril.devplot = ggplot(data=ril.all, aes(x=line, y=devprojection, fill = line)) +
      ggtitle("Distribution of developing worms on development axis by line type") +
      ylab("Projection on development axis [arbitrary units]") + 
      xlab("Line type") +
      scale_fill_manual(values = c("slateblue2", "tomato2")) +
      geom_violin() +
      coord_flip()
    ril.devplot   
    
    ril.allplot = ggplot(data=ril.dev, aes(x=line, y=devprojection, fill = line)) +
      ggtitle("Distribution of all worms on development axis by line type") +
      ylab("Projection on development axis [arbitrary units]") + 
      xlab("Line type") +
      scale_fill_manual(values = c("slateblue2", "tomato2")) +
      geom_violin() +
      coord_flip()
    ril.allplot
    
## test random forest run ---------------------------------------------------
    
    ## example gene from Francesconi & Lehner 
    ## fpn-1.3 .. WBGene00019969 
    ## Y75B7B.2 .. WBGene00022289 
    ## clec-147 .. WBGene00007313
    ## hsp-16.41 .. WBGene00002018

    use.gexp <- as.numeric(gexp.dev.nd["WBGene00019969",]) ; length(use.gexp)
    # use.gexp <- log2( (2^use.gexp)/mean(2^use.gexp))
    prax <- proj.dev$devprojection ; length(prax)
    use.gm <- (t(gm.dev.nd)) ; dim(use.gm)
    # use.gm <- apply(use.gm,2,as.factor)
    rf.res <- randomForest(use.gexp~prax + .,data = use.gm,proximity = T, importance =T,ntree=5000,nodesize=2,maxnodes=8)
    
head(rf.res$importance)
        
plot(rf.res$importance[,1])

top.mrk <- which(rf.res$importance[-1,1] == max(rf.res$importance[-1,1]))
qtl.allele <- (use.gm[,top.mrk])
qtl.allele[qtl.allele == -1 ] <- "CB"
qtl.allele[qtl.allele == 1 ] <- "N2"

to.pl <- data.frame(prax,qtl.allele,use.gexp)

this <- ggplot(to.pl,aes(prax,use.gexp,col=qtl.allele))+
      geom_smooth(se=F) +
      geom_point()+
      ggtitle(colnames(use.gm)[top.mrk])+ 
      labs(x = "Projection on development axis [arbitrary units]", y = "Expression level") +
      scale_color_manual(values = c("blue","red"))

this

head(to.pl)
    

# projection_all[projection_all$line == "NIL" &
#                  projection_all$experimentype == "Dev",c("genotype","devprojection")]

## Loop over all genes in dev worms ------------------------------------------------------------------------

#Create empty matrix
rf.all.dev <- matrix(,nrow = nrow(gexp.dev.nd), ncol = ncol(use.gm.dev), dimnames = list(rownames(gexp.dev.nd), colnames(use.gm.dev)))

    #Optional: slim down datasets for testing
    gexp.dev.nd <- gexp.dev[,colnames(gm.dev.nd)] #Remake for safety
    gexp.dev.nd <- gexp.dev.nd[seq(1, nrow(gexp.dev.nd), 100),] #Select every 100th gene for sample dataset
    rf.all.dev <- matrix(,nrow = nrow(gexp.dev.nd), ncol = ncol(use.gm.dev), dimnames = list(rownames(gexp.dev.nd), colnames(use.gm.dev))) #Resize matrix for smaller dataset



for (row in 1:nrow(gexp.dev.nd)) {
  use.gexp <- as.numeric(gexp.dev.nd[row,])
  rf.res <- randomForest(use.gexp~prax.dev + .,data = use.gm.dev,proximity = T, importance =T,ntree=5000,nodesize=2,maxnodes=8)
  rf.all.dev[row,] <- t(rf.res$importance[-1,1])
  print(paste(row, "/", nrow(gexp.dev.nd)))
}
write.csv(rf.all.dev, file = "RF_devgenes.csv")



## Loop over all genes in HS worms ------------------------------------------------------------------------

#Create empty matrix
rf.all.hs <- matrix(,nrow = nrow(gexp.hs.nd), ncol = ncol(use.gm.hs), dimnames = list(rownames(gexp.hs.nd), colnames(use.gm.hs)))

    #Optional: slim down datasets for testing
    gexp.hs.nd <- gexp.hs[,colnames(gm.hs.nd)] #Remake for safety
    gexp.hs.nd <- gexp.hs.nd[seq(1, nrow(gexp.hs.nd), 100),] #Select every 100th gene for sample dataset
    rf.all.hs <- matrix(,nrow = nrow(gexp.hs.nd), ncol = ncol(use.gm.hs), dimnames = list(rownames(gexp.hs.nd), colnames(use.gm.hs))) #Resize matrix for smaller dataset

for (row in 1:nrow(gexp.hs.nd)) {
  use.gexp <- as.numeric(gexp.hs.nd[row,])
  rf.res <- randomForest(use.gexp~prax.hs + .,data = use.gm.hs,proximity = T, importance =T,ntree=5000,nodesize=2,maxnodes=8)
  rf.all.hs[row,] <- t(rf.res$importance[-1,1])
  print(paste(row, "/", nrow(gexp.hs.nd)))
}
write.csv(rf.all.hs, file = "RF_hsgenes.csv")




## Loop over all genes in rec worms ------------------------------------------------------------------------

#Create empty matrix
rf.all.rec <- matrix(,nrow = nrow(gexp.rec.nd), ncol = ncol(use.gm.rec), dimnames = list(rownames(gexp.rec.nd), colnames(use.gm.rec)))

    #Optional: slim down datasets to test
    gexp.rec.nd <- gexp.rec[,colnames(gm.rec.nd)] #Remake for safety
    gexp.rec.nd <- gexp.rec.nd[seq(1, nrow(gexp.rec.nd), 100),] #Select every 100th gene for sample dataset
    rf.all.rec <- matrix(,nrow = nrow(gexp.rec.nd), ncol = ncol(use.gm.rec), dimnames = list(rownames(gexp.rec.nd), colnames(use.gm.rec))) #Resize matrix for smaller dataset


for (row in 1:nrow(gexp.rec.nd)) {
  use.gexp <- as.numeric(gexp.rec.nd[row,])
  rf.res <- randomForest(use.gexp~prax.rec + .,data = use.gm.rec,proximity = T, importance =T,ntree=5000,nodesize=2,maxnodes=8)
  rf.all.rec[row,] <- t(rf.res$importance[-1,1])
  print(paste(row, "/", nrow(gexp.rec.nd)))
}
write.csv(rf.all.rec, file = "RF_recgenes.csv")

### Analyse all eQTLs genome-wide --------------------------------------------------------------------------------------
## Prepare genome-wide data/functions -----------------------

#Function to create a full dataset with gene/marker locations from eQTL pairs
qtlGenome <- function(df, cutoff = qtl_cutoff, genedata = WormGenes, mrkdata = mrk.info){
  qtl_GxM <- df[,-1] #Separate the gene x marker matrix
  qtl_GxPrax <- df[,1]# from the gene x developmental axis data
  qtl_GxPrax <- rownames_to_column(data.frame(qtl_GxPrax))
  colnames(qtl_GxPrax) <- c("geneid", "devImportance")
  
  qtl_pairs <- melt(qtl_GxM, #Reshape data from wide to long
                    varnames = c("geneid", "mrkid"), 
                    value.name = "importance") %>%
    filter(importance >= cutoff) %>% #Select the values above the cutoff
    left_join(y = qtl_GxPrax, by = "geneid") %>% #Attach dev importance per gene
    mutate(totalImportance = importance/mean(importance) + devImportance/mean(devImportance))
  
  #Attach gene data
  qtl_merge1 <- merge(qtl_pairs, genedata, by.x = "geneid", by.y = "row.names") %>%
    rename(genename = name ,
           genestart = start, 
           geneend = end, 
           genechr = chr)
  
  #Attach marker data
  qtl_merge2 <- merge(qtl_merge1, mrkdata, by.x = "mrkid", by.y = "marker") %>%
    rename(mrkchr = chr,
           mrkstart = start,
           mrkend = end) 
  
  #Classify and clean up
  qtl_merged <- qtl_merge2 %>%
    filter(genechr != "MtDNA" & mrkchr != "MtDNA") %>%
    mutate(distance = abs(genestart-mrkstart), #calculate distance between gene and marker
           type = ifelse(mrkchr == genechr, #If the two chromosomes are the same, check distance
                         yes = ifelse(distance < 1000000, #Check for 1 Mb apart if same chromosome
                                      "local", 
                                      "distant"),
                         no = "distant") #If two different chromosomes, always distant
    ) %>%
    select(geneid, genename, genestart, geneend, genechr, mrkid, mrkstart, mrkend, mrkchr, distance, type, importance, devImportance, totalImportance) %>% #Reorder columns
    mutate(genechr = fct_relevel(genechr, "X", "V", "IV", "III", "II", "I"), #Needed for cis-trans plot
           mrkchr = fct_relevel(mrkchr, "I", "II", "III", "IV", "V", "X")) %>%
    arrange(desc(importance))
  
  return(qtl_merged)
}

### Plotting functions

## Cis-trans plots
#Function for a standard cis-trans QTL plot    
qtlCTplot <- function(df){
  ggplot(df, aes(mrkstart, genestart, shape = type)) + 
    geom_point() + 
    facet_grid(genechr~mrkchr) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = "Marker position", y = "Gene position")
}

#Function for cis-trans plots of different datasets, splits on 'compare'
qtlCompare <- function(df, compare){
  ggplot(df, aes_string("mrkstart", "genestart", col = compare)) + 
    geom_point(alpha = 0.3) + 
    facet_grid(genechr~mrkchr) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = "Marker position", y = "Gene position")
}

## Hotspot plots
#Function to plot a histogram of distant eQTLs to find hotspots
qtlHotspot <- function(df){
  df <- filter(df, type == "distant")
  ggplot(data = df, aes(mrkstart)) + 
    geom_histogram(binwidth = 2e6) + #Not sure about best binwidth
    facet_wrap(~mrkchr,nrow = 1) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"), #Mimic c-tplot style
          axis.text = element_blank(), 
          axis.ticks = element_blank()) +
    labs(x = "Marker position", y = "#eQTLs") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) #Actually start from 0
}

#Function to plot a distant eQTLs histogram with multiple datasets
qtlCompHotspot <- function(df, compare){
  df <- filter(df, type == "distant")
  ggplot(data = df, aes_string("mrkstart", fill = compare)) + 
    geom_histogram(binwidth = 2e6) + #Not sure about best binwidth
    facet_wrap(~mrkchr,nrow = 1) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"), #Mimic c-tplot style
          axis.text = element_blank(), 
          axis.ticks = element_blank()) +
    labs(x = "Marker position", y = "#eQTLs") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) #Actually start from 0
}


## Combined cis-trans and hotspot plots
qtlBoth <- function(df){
  plot_grid(qtlCTplot(df) + theme(panel.grid.minor = element_blank(),
                                  plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm",),
                                  axis.title.x = element_blank()
  ), 
  qtlHotspot(df) + theme(panel.grid.minor.x = element_blank(),
                         plot.margin = unit(c(0,0,0.2,0), "cm"),
                         strip.background = element_blank(),
                         strip.text.x = element_blank()
  ), 
  ncol = 1, 
  align = 'v', 
  axis = 'lr',
  rel_heights = c(3,1))
}

qtlCompBoth <- function(df, compare){
  plot_grid(qtlCompare(df, compare) + theme(panel.grid.minor = element_blank(),
                                            plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm",),
                                            axis.title.x = element_blank(),
                                            legend.title = element_text(size = 11)
  ), 
  qtlCompHotspot(df, compare) + theme(panel.grid.minor.x = element_blank(),
                                      plot.margin = unit(c(0,0,0.2,0), "cm"),
                                      strip.background = element_blank(),
                                      strip.text.x = element_blank()
  ), 
  ncol = 1, 
  align = 'v', 
  axis = 'lr',
  rel_heights = c(3,1))
}


## Start analysis ----------------  
#Get the genome-wide data for all three conditions
qtl_dev <- qtlGenome(data.incmse.dev, cutoff = 0.2)
qtl_hs <- qtlGenome(data.incmse.hs, cutoff = 0.2)
qtl_rec <- qtlGenome(data.incmse.rec, cutoff = 0.2)
head(qtl_dev)

#Count total eQTLs
matrix(c(
  sum(qtl_dev$type == "local"), sum(qtl_hs$type == "local"), sum(qtl_rec$type == "local"),
  sum(qtl_dev$type == "distant"), sum(qtl_hs$type == "distant"), sum(qtl_rec$type == "distant"),
  nrow(qtl_dev), nrow(qtl_hs), nrow(qtl_rec)),
  nrow = 3, dimnames = list(c("Dev", "HS", "Rec"), c("Local", "Distant", "Total"))
)
## Start cis-trans plot ----
qtl_devplot <- qtlCTplot(qtl_dev) + ggtitle("eQTLs found during development")
qtl_devplot

qtlBoth(qtl_dev)

print(c("Total number of eQTLs:", nrow(qtl_dev)))
print(c("Number of local eQTLs:", sum(qtl_dev$type == "local")))
print(c("Number of distant eQTLs:", sum(qtl_dev$type == "distant")))
ggsave(plot = qtl_ctplot,
       filename = paste0("cistransplot_",
                         "t", nrow(qtl_dev), #Total eQTLs
                         "_l", sum(qtl_dev$type == "local"), #Local eQTLs
                         "_d", sum(qtl_dev$type == "distant"), #Distant eQTLS
                         ".png")
)

# Explore the different conditions
qtl_hsplot <- qtlCTplot(qtl_hs) + ggtitle("eQTLs found during heat shock")
qtl_hsplot

qtl_recplot <- qtlCTplot(qtl_rec) + ggtitle("eQTLs found during recovery")
qtl_recplot

qtl_all <- bind_rows(
  qtl_dev %>% mutate(condition = "dev"),
  qtl_hs %>% mutate(condition = "hs"),
  qtl_rec %>% mutate(condition = "rec")
)

qtlCompare(qtl_all, compare = "condition") + ggtitle("eQTLs in the three conditions")

qtl_bothplot <- qtlCompBoth(qtl_all, compare = "condition")
qtl_bothplot
ggsave2(plot = qtl_bothplot,
        filename = paste0("bothplot_",
                          "t", nrow(qtl_all), #Total eQTLs
                          "_l", sum(qtl_all$type == "local"), #Local eQTLs
                          "_d", sum(qtl_all$type == "distant"), #Distant eQTLS
                          ".png")
)

## Create combined cis-trans plot of both importance measures
#Recreate the data for with purity instead

qtl_comppur <- qtlGenome(data.pur.dev, cutoff = 4)

qtl_compmse <- mutate(qtl_dev, data = "incmse")
qtl_comppur <- mutate(qtl_ctpur, data = "pur")
#Merge the mse and purity datasets
qtl_msepur <- bind_rows(qtl_comppur, qtl_compmse) 

#Plot them together
qtlCompare(qtl_msepur, compare = "data") + ggtitle("All eQTLs from both measures of importance")


#Find unique QTLs in each dataset
qtl_uniquemse <- anti_join(qtl_compmse, qtl_comppur, by = c('geneid', 'mrkid'))
qtl_uniquepur <- anti_join(qtl_comppur, qtl_compmse, by = c('geneid', 'mrkid'))
qtl_unique <- bind_rows(qtl_uniquemse, qtl_uniquepur)
nrow(qtl_uniquemse)
nrow(qtl_uniquepur)

qtlCompare(qtl_unique, "data") + ggtitle("Unique eQTLs per importance marker")

### Analyse individual eQTLs ------------------------------------------------------------------------------------
## Prepare individual eQTL functions -------------------------
qtlData <- function(marker = qtl_mrk, gene = qtl_gene){
  df <- data.frame(
    Development = qtl_prax,
    Allele = qtl_gm[,marker],
    Expression = as.double(qtl_gexp[gene,])
  )
  return(df)
}

qtlPlot <- function(df, marker = qtl_mrk, gene = qtl_gene){
  ggplot(df, aes(Development, Expression, col=Allele)) +
    geom_smooth(se=F) +
    geom_point()+
    ggtitle(paste(gene, marker))+
    labs(x = "Projection on development axis [arbitrary units]", y = "Expression level") +
    scale_color_manual(values = c("blue","red"))
}

## Prepare and sort RF data ----------------------

      #Use Dev (*.dev.*), HS (*.hs.*) or Rec (*.rec.*) data?
      
      qtl_gexp <-          gexp.dev.nd ; dim(qtl_gexp)
      qtl_rf <-     data.incmse.dev ; dim(qtl_rf)
      qtl_gm <-          use.gm.dev ; dim(qtl_gm)
      qtl_prax <-          prax.dev ; length(qtl_prax)
      
      qtl_sig <-  qtlGenome(qtl_rf) ; dim(qtl_sig)

#Properly label genetic map in preparation for plots      
qtl_gm[qtl_gm == -1 ] <- "CB"
qtl_gm[qtl_gm == 1 ] <- "N2"

#To replace significance: only work w/ high values, choose cutoff manually, based on data distribution
qtl_cutoff <- 0.2 #0.2 for incmse, 4 for purity


## Remake non-cutoff QTL data for exploration if needed (warning: very large)
          # qtl_mrks <- melt(qtl_rf[,-1], varnames = c("geneid", "mrkid"), value.name = "importance") %>%
          #   arrange(desc(importance))
          # qtl_mrks <- left_join(qtl_mrks, rownames_to_column(data.frame(qtl_rf[,1]), "geneid"), by = "geneid") %>%
          #   rename(devImportance = qtl_rf[,1])
          # qtl_mrks$totalImportance <- qtl_mrks$importance/mean(qtl_mrks$importance) + qtl_mrks$devImportance/mean(qtl_mrks$devImportance)
          # 
          # ggplot(data=qtl_mrks, aes(1:length(importance), importance)) +
          #   geom_point(alpha=0.5) +
          #   geom_hline(yintercept=0.2, color="red", linetype="dashed") +
          #   xlab("Index") +
          #   ggtitle("Importance values of best marker of all genes")
          # 
          # qtl_sig <- qtl_mrks[qtl_mrks$importance >= qtl_cutoff,]


## Sort the QTLs by marker importance (default), development importance or combined.
#Sort by importance of development for eQTL
qtl_sig <- qtl_sig[order(-qtl_sig$devImportance),]
#Sort by combined importance
qtl_sig <- qtl_sig[order(-qtl_sig$totalImportance),]
#Resort by importance of marker for eQTL
qtl_sig <- qtl_sig[order(-qtl_sig$importance),]


### Start plotting individual eQTLs ----------------
## Look at one specific gene and its top marker
    #Choose rank of desired gene by QTL importance
    qtl_rank <- 2
qtl_mrk <- as.character(qtl_sig[qtl_rank,"mrkid"])
qtl_gene <- as.character(qtl_sig[qtl_rank,"geneid"])
qtlData(marker = qtl_mrk, gene = qtl_gene) %>%
  qtlPlot()

## Look at n genes with strongest (most important) eQTLs
    
    #Amount of best eQTLs to plot
    n_QTLs <- 9
#Create an empty dataframe to append to
qtl_multidata <- data.frame(Development = double(), 
                               Allele = character(),
                               Expression = double(),
                               Gene = character())
#Get the data for each QTL, label it by gene and append
for (i in 1:n_QTLs){
  qtl_mrk <- as.character(qtl_sig[i,"mrkid"])
  qtl_gene <- as.character(qtl_sig[i,"geneid"])
  qtl_loopdata <- cbind(qtlData(qtl_mrk, qtl_gene), qtl_gene)
  qtl_multidata <- rbind(qtl_multidata, qtl_loopdata)
}

qtl_facetplot <- ggplot(qtl_multidata, aes(Development, Expression, col=Allele))+
  geom_smooth(se=F) +
  geom_point()+
  ggtitle(paste("Top", n_QTLs, "QTLs by importance"))+
  labs(x = "Projection on development axis [arbitrary units]", y = "Expression level") +
  scale_color_manual(values = c("blue","red")) +
  facet_wrap(vars(qtl_gene))
qtl_facetplot


## Check all QTLs of one gene
#Writes all plots to folder named after the gene, doesn't show in Rstudio itself
    #Which gene to study?
    qtl_gene <- "WBGene00008352"
    dir.create("All_QTLs")
    
geneQTLs <- sort(qtl_rf[qtl_gene,-1], decreasing = T) #Sort all markers of one gene by importance
for (i in 1:length(geneQTLs)){     #For each marker, starting with the strongest affinity
  qtl_mrk <- names(geneQTLs)[i]
  
  qtlData(marker = qtl_mrk, gene = qtl_gene) %>%
    qtlPlot() %>%
    ggsave(plot = ., filename = file.path("All_QTLs", paste(qtl_gene, i, ".png", sep="_")))
  print(paste(i, "/", length(geneQTLs)))
}

## Just straight up check all genes with their significant markers
#Writes all plots to folder, doesn't show in Rstudio itself, can take a long time!
    dir.create("All_genes")

for (generank in 1:length(qtl_sig$Importance)){
  qtl_mrk <- as.character(qtl_sig[generank,1])
  qtl_gene <- rownames(qtl_sig[generank,])
  
  qtlData(qtl_mrk, qtl_gene) %>%
    qtlPlot() %>%
    ggsave(plot = ., file.path("All_genes", paste(generank, qtl_gene, qtl_mrk, ".png", sep="_")))
  print(paste(generank, "/", length(qtl_sig$Importance)))
}


## Plot all genes in one big plot
# qtl_bigplot <- ggplot() +
#   labs(x = "Projection on development axis [arbitrary units]", y = "Expression level") +
#   scale_color_manual(values = c("blue","red"))
# for (generank in 1:200){
#   qtl_mrk <- as.character(qtl_sig[generank,1])
#   qtl_gene <- rownames(qtl_sig[generank,])
#   qtl_bigdata <- qtlData(qtl_mrk, qtl_gene)
#   
#   qtl_bigplot <- qtl_bigplot +
#     stat_smooth(qtl_bigdata, geom="line", mapping = aes(Development, Expression, col=Allele), se=F, alpha=0.2)
# 
#   print(paste(generank, "/", length(qtl_sig$Importance)))
# }
# qtl_bigplot



randomForestExplainer::explain_forest(rf.res)

##################### END ##############################################################################################
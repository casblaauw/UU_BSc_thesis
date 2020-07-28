### Import, reshape and explore starting data --------------------------------------------------------------------------------------------
## Libraries ----
library(openxlsx)
library(ggplot2)
library(gplots)
library(cowplot)
library(tidyverse)
library(forcats)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)
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

devColor <- "#00BA38"
hsColor <- "#F8766D"
recColor <- "#619CFF"

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
qtlCTplot <- function(df, color = "#000000"){
  ggplot(df, aes(mrkstart, genestart)) + 
    geom_point(alpha = 0.5, color = color) + 
    facet_grid(genechr~mrkchr) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none") +
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
    scale_colour_manual(values = c(devColor, hsColor, recColor)) +
    labs(x = "Marker position", y = "Gene position")
}

## Hotspot plots
#Function to plot a histogram of distant eQTLs to find hotspots
qtlHotspot <- function(df, color = "grey35"){
  df <- filter(df, type == "distant")
  ggplot(data = df, aes(mrkstart)) + 
    geom_histogram(fill = color, binwidth = 1e6) + #Not sure about best binwidth
    facet_wrap(~mrkchr,nrow = 1) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"), #Mimic c-tplot style
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none") +
    labs(x = "Marker position", y = "#eQTLs") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) #Actually start from 0
}

#Function to plot a distant eQTLs histogram with multiple datasets
qtlCompHotspot <- function(df, compare){
  df <- filter(df, type == "distant")
  ggplot(data = df, aes_string("mrkstart", fill = compare)) + 
    geom_histogram(binwidth = 1e6) + #Not sure about best binwidth
    facet_wrap(~mrkchr,nrow = 1) + 
    theme_bw() +
    theme(panel.spacing = unit(0.1, "lines"), #Mimic c-tplot style
          axis.text = element_blank(), 
          axis.ticks = element_blank()) +
    labs(x = "Marker position", y = "#eQTLs") +
    scale_fill_manual(values = c(devColor, hsColor, recColor)) +
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
                                            legend.title = element_text(size = 11)) +
                                      guides(color = guide_legend(override.aes = list(size = 7, 
                                                                                      alpha = 1, 
                                                                                      shape = "square"))
  ), 
  qtlCompHotspot(df, compare) + theme(panel.grid.minor.x = element_blank(),
                                      plot.margin = unit(c(0,0,0.2,0), "cm"),
                                      strip.background = element_blank(),
                                      strip.text.x = element_blank(),
                                      legend.position = "none"
  ), 
  ncol = 1, 
  align = 'v', 
  axis = 'lr',
  rel_heights = c(4,1))
}


## Create genome-wide data ----------------  
# Set the eQTL cutoff
qtl_cutoff <- 0.2 #0.2 for incmse, 4 for purity

#Get the genome-wide data for all three conditions
qtl_dev <- qtlGenome(data.incmse.dev, cutoff = 0.2)
qtl_hs <- qtlGenome(data.incmse.hs, cutoff = 0.2)
qtl_rec <- qtlGenome(data.incmse.rec, cutoff = 0.2)
head(qtl_dev)
qtl_all <- bind_rows(
  qtl_dev %>% mutate(condition = "dev"),
  qtl_hs %>% mutate(condition = "hs"),
  qtl_rec %>% mutate(condition = "rec")
)

#Count total eQTLs
matrix(c(
  sum(qtl_dev$type == "local"), sum(qtl_hs$type == "local"), sum(qtl_rec$type == "local"),
  sum(qtl_dev$type == "distant"), sum(qtl_hs$type == "distant"), sum(qtl_rec$type == "distant"),
  nrow(qtl_dev), nrow(qtl_hs), nrow(qtl_rec)),
  nrow = 3, dimnames = list(c("Dev", "HS", "Rec"), c("Local", "Distant", "Total"))
)

print(paste("Total:", "dev", nrow(qtl_dev), "hs", nrow(qtl_hs), "rec", nrow(qtl_rec)))
qtl_counts <- qtl_all %>% group_by(geneid, mrkid) %>% summarise(n = n()) %>% group_by(n) %>% summarise(total = n())
print(qtl_counts) #The amount of qtls shared between 1, 2 or 3 conditions
print(paste("dev exclusive:", nrow(anti_join(qtl_dev, rbind(qtl_hs,qtl_rec), by = c("geneid", "mrkid"))),
            "hs exclusive:", nrow(anti_join(qtl_hs, rbind(qtl_dev,qtl_rec), by = c("geneid", "mrkid"))),
            "rec exclusive:", nrow(anti_join(qtl_rec, rbind(qtl_dev,qtl_hs), by = c("geneid", "mrkid")))))
print(paste("dev-hs:", nrow(inner_join(qtl_dev, qtl_hs, by = c("geneid", "mrkid")))-qtl_counts[3,2],
            "dev-rec:", nrow(inner_join(qtl_dev, qtl_rec, by = c("geneid", "mrkid")))-qtl_counts[3,2],
            "hs-rec:",  nrow(inner_join(qtl_hs, qtl_rec, by = c("geneid", "mrkid"))) -qtl_counts[3,2]))

venn.diagram(list(transmute(qtl_dev, paste0(geneid, mrkid))[,1],
                  transmute(qtl_hs, paste0(geneid, mrkid))[,1],
                  transmute(qtl_rec, paste0(geneid, mrkid))[,1]),
             filename = "venn.png",
             category.names = c("Development", "Heat stress", "Recovery"),
             #Adapted from https://www.r-graph-gallery.com/14-venn-diagramm.html
             # Circles
             lwd = 1,
             # lty = 'blank',
             fill = brewer.pal(3, "Pastel2"),
             col = brewer.pal(3, "Set2"),
             # Numbers
             cex = 1.8,
             fontface = "bold",
             fontfamily = "sans",
             # Set names
             cat.cex = 1.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)

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



  #Plot all QTLs for three different conditions side by side
qtl_sidebyside <- plot_grid(qtlCTplot(qtl_dev, devColor) + 
                              theme(plot.margin = unit(rep(0.2,4), "cm",),
                                    axis.title.x = element_blank()), 
                            qtlCTplot(qtl_hs, hsColor) + 
                              theme(plot.margin = unit(rep(0.2,4), "cm",),
                                    axis.title.x = element_blank()), 
                            qtlCTplot(qtl_rec, recColor) + 
                              theme(plot.margin = unit(rep(0.2,4), "cm",),
                                    axis.title.x = element_blank()), 
                            qtlHotspot(qtl_dev, devColor) + 
                              theme(plot.margin = unit(rep(0.2,4), "cm",),
                                    strip.background = element_blank(),
                                    strip.text.x = element_blank()),
                            qtlHotspot(qtl_hs, hsColor) + 
                              theme(plot.margin = unit(rep(0.2,4), "cm",),
                                    strip.background = element_blank(),
                                    strip.text.x = element_blank()),
                            qtlHotspot(qtl_rec, recColor) + 
                              theme(plot.margin = unit(rep(0.2,4), "cm",),
                                    strip.background = element_blank(),
                                    strip.text.x = element_blank()),
                            align = "hv",
                            axis = 'lr',
                            nrow = 2,
                            rel_heights = c(4,1),
                            labels = c("A", "B", "C"))

ggsave2(plot = qtl_sidebyside,
        filename = paste0("sbsplot2_",
                          "d", nrow(qtl_dev),
                          "_h", nrow(qtl_hs),
                          "_r", nrow(qtl_rec),
                          ".png"),
                          width = 18,
                          height = 6)


  # Create a plot of all eQTLs per condition overlaid
qtl_bothplot <- qtlCompBoth(qtl_all, compare = "condition")
qtl_bothplot

ggsave2(plot = qtl_bothplot,
        filename = paste0("bothplot_",
                          "t", nrow(qtl_all), #Total eQTLs
                          "_l", sum(qtl_all$type == "local"), #Local eQTLs
                          "_d", sum(qtl_all$type == "distant"), #Distant eQTLS
                          ".png")
)


    #Find the unique eQTLs per condition
qtl_uniquedev <- anti_join(qtl_dev, rbind(qtl_hs, qtl_rec), by = c('geneid', 'mrkid'))
qtl_uniquehs <- anti_join(qtl_hs, rbind(qtl_dev, qtl_rec), by = c('geneid', 'mrkid'))
qtl_uniquerec <- anti_join(qtl_rec, rbind(qtl_dev, qtl_hs), by = c('geneid', 'mrkid'))
qtl_uniqueall <- bind_rows(
  qtl_uniquedev %>% mutate(condition = "Development"), 
  qtl_uniquehs %>% mutate(condition = "Heat stress"),
  qtl_uniquerec %>% mutate(condition = "Recovery"))

    #Distribution of local and distant amongst the unique eQTLs
qtl_uniqueall %>% group_by(condition) %>% count(type)

qtl_uniqueplot <- qtlCompBoth(qtl_uniqueall, compare = "condition")
ggsave2(plot = qtl_uniqueplot,
        filename = paste0("uniqueplot_",
                          "d", nrow(qtl_uniquedev),
                          "_h", nrow(qtl_uniquehs),
                          "_r", nrow(qtl_uniquerec),
                          ".png"),
        width = 8,
        height = 6)

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
    geom_smooth(se=F, span = 1) +
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
    n_QTLs <- 16
#Create an empty dataframe to append to
qtl_multidata <- data.frame(Development = double(), 
                               Allele = character(),
                               Expression = double(),
                               QTL = character())
#Get the data for each QTL, label it by gene and append
for (i in 1:n_QTLs){
  qtl_mrk <- as.character(qtl_sig[i,"mrkid"])
  qtl_gene <- as.character(qtl_sig[i,"geneid"])
  qtl_loopdata <- cbind(qtlData(qtl_mrk, qtl_gene), 
                        QTL = paste0(
                          substr(qtl_gene, 10, 14),
                          "\n",
                          substr(qtl_mrk, 5, nchar(qtl_mrk))))
  qtl_multidata <- rbind(qtl_multidata, qtl_loopdata)
}

qtl_facetplot <- ggplot(qtl_multidata, aes(Development, Expression, col=Allele))+
  geom_smooth(se=F, span = 1) +
  geom_point()+
  ggtitle(paste("Top", n_QTLs, "QTLs by importance"))+
  labs(x = "Projection on development axis [arbitrary units]", y = "Expression level") +
  scale_color_manual(values = c("blue","red")) +
  facet_wrap(vars(QTL))
qtl_facetplot

# Make a facet plot of the QTLS of your own choice
qtl_choose <- c(72, 93, 355,
                324, 195,317,
                263, 252, 108,
                310, 304, 286)
#Create an empty dataframe to append to
qtl_choosedata <- data.frame(Development = double(), 
                            Allele = character(),
                            Expression = double(),
                            Gene = character(),
                            Marker = character(),
                            QTL = factor())
#Get the data for each QTL, label it by gene and append
for (i in qtl_choose){
  qtl_mrk <- as.character(qtl_sig[i,"mrkid"])
  qtl_gene <- as.character(qtl_sig[i,"geneid"])
  qtl_chooseloop <- cbind(qtlData(qtl_mrk, qtl_gene),
                          Gene = qtl_gene,
                          Marker = qtl_mrk,
                          QTL = factor(paste0(
                            WormGenes[qtl_gene, "name"],
                            "\n",
                            qtl_mrk)))
  qtl_choosedata <- rbind(qtl_chooseloop, qtl_choosedata)
}

qtl_chooseplot <- ggplot(qtl_choosedata, aes(Development, Expression, col=Allele))+
  geom_smooth(se=F, span = 1) +
  geom_point(size = 1) +
  labs(x = "Projection on development axis [arbitrary units]", y = "Expression level") +
  scale_color_manual(values = c("blue","red")) +
  theme_bw() + 
  theme(strip.text = element_text(size = 10, margin = ggplot2::margin(0.1, 0, .1, 0, "cm"))) +
  facet_wrap(vars(QTL), scales = "free_y", nrow = 4, as.table = FALSE)
qtl_chooseplot

ggsave(plot = qtl_chooseplot,
       filename = paste0(c("chooseplot2",
                           qtl_choose,
                           ".png"),
                         collapse = "_"),
       height = 12,
       width = 8)



## Just straight up check all genes with their significant markers
#Writes all plots to folder, doesn't show in Rstudio itself, can take a long time!
    dir.create("All_genes")

for (generank in 1:(length(qtl_sig$importance))){
  qtl_mrk <- as.character(qtl_sig[generank,"mrkid"])
  qtl_gene <- as.character(qtl_sig[generank,"geneid"])
  
  qtlData(qtl_mrk, qtl_gene) %>%
    qtlPlot() %>%
    ggsave(plot = ., file.path("All_genes", paste(generank, qtl_gene, qtl_mrk, ".png", sep="_")))
  print(paste(generank, "/", length(qtl_sig$importance)))
}

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

## Supplemental ---------------------------------------------------
# Figure 1
#Plot the distribution of importance values for incmse and purity
qtl_incmse <- melt(data.incmse.dev[,-1], varnames = c("geneid", "mrkid"), value.name = "importance") %>%
  arrange(desc(importance))
    # Add the extra importance dimensions if needed
    #qtl_incmse <- left_join(qtl_incmse, rownames_to_column(data.frame(qtl_rf[,1]), "geneid"), by = "geneid") %>%
    #     rename(devImportance = qtl_rf[,1])
    #qtl_incmse$totalImportance <- qtl_mrks$importance/mean(qtl_mrks$importance) + qtl_mrks$devImportance/mean(qtl_mrks$devImportance)
    
qtl_pur <- melt(data.pur.dev[,-1], varnames = c("geneid", "mrkid"), value.name = "importance") %>%
  arrange(desc(importance))

qtl_incmseplot <- ggplot(data=qtl_incmse, aes(1:length(importance), importance)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=0.2, color="red", linetype="dashed") +
  xlab("Index") +
  # xlim(0,1000) +
  ggtitle("Distribution of incMSE importance values")

qtl_purplot <- ggplot(data=qtl_pur, aes(1:length(importance), importance)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=4, color="red", linetype="dashed") +
  xlab("Index") +
  # xlim(0,1000) +
  ggtitle("Distribution of purity importance values") 

qtl_incpurplot <- plot_grid(qtl_incmseplot, qtl_purplot, align = "h", labels = "AUTO")
ggsave2(plot = qtl_incpurplot, filename = "add1.png", width = 30, height = 15, units = "cm")

## Create combined cis-trans plot of both importance measures
#Recreate the data for with purity instead

qtl_comppur <- qtlGenome(data.pur.dev, cutoff = 4)

qtl_compmse <- mutate(qtl_dev, data = "incmse")
qtl_comppur <- mutate(qtl_comppur, data = "pur")
#Merge the mse and purity datasets
qtl_msepur <- bind_rows(qtl_comppur, qtl_compmse) 

#Plot them together
qtl_comp1 <- qtlCompare(qtl_msepur, compare = "data") + ggtitle("All eQTLs from both measures of importance")
qtl_comp1

#Find unique QTLs in each dataset
qtl_uniquemse <- anti_join(qtl_compmse, qtl_comppur, by = c('geneid', 'mrkid'))
qtl_uniquepur <- anti_join(qtl_comppur, qtl_compmse, by = c('geneid', 'mrkid'))
qtl_uniqueimp <- bind_rows(qtl_uniquemse, qtl_uniquepur)
nrow(qtl_uniquemse)
nrow(qtl_uniquepur)
qtl_uniquemse %>% count(type)
qtl_uniquepur %>% count(type)

qtl_comp2 <- qtlCompare(qtl_uniqueimp, "data") + ggtitle("Unique eQTLs per importance marker")
qtl_comp2

qtl_compboth <- plot_grid(qtl_comp1, qtl_comp2, align = "h", labels = "AUTO")
ggsave2(qtl_compboth, filename = "add2.png", width = 30, height = 15, units = "cm")

#Find the major hotspot bins

# qtl_hotspotmrks <- qtl_all %>% 
#   group_by(mrkid, mrkstart, mrkchr) %>% 
#   summarise(n = n(), .groups = "keep") %>% 
#   arrange(desc(n)) %>%
#   filter(n > 1)
# 
# qtl_hotspotbins <-  qtl_all %>% 
#   filter(type == "distant") %>%
#   mutate(bin = cut_width(mrkstart, 1e6)) %>%
#   group_by(bin, mrkchr, condition) %>% 
#   summarise(n = n(), .groups = "keep") %>% 
#   arrange(mrkchr) %>%
#   pivot_wider(names_from = condition, values_from = n) 
# qtl_hotspotbins[is.na(qtl_hotspotbins)] <- 0
# qtl_hotspotbins <- qtl_hotspotbins %>%
#   mutate(total = dev + hs + rec) %>%
#   select(bin, mrkchr, dev, hs, rec, total) %>%
#   arrange(desc(total))


qtl_hotspotuni <- qtl_all %>%
  filter(type == "distant") %>%
  mutate(bin = cut_width(mrkstart, 1e6)) %>%
  distinct(bin, genename, .keep_all = TRUE) %>%
  group_by(bin, mrkchr) %>% 
  summarise(total = n(), .groups = "keep")
qtl_hotspotdev <- qtl_dev %>%
  filter(type == "distant") %>%
  mutate(bin = cut_width(mrkstart, 1e6)) %>%
  distinct(bin, genename, .keep_all = TRUE) %>%
  group_by(bin, mrkchr) %>% 
  summarise(dev = n(), .groups = "keep")
qtl_hotspoths <- qtl_hs %>%
  filter(type == "distant") %>%
  mutate(bin = cut_width(mrkstart, 1e6)) %>%
  distinct(bin, genename, .keep_all = TRUE) %>%
  group_by(bin, mrkchr) %>% 
  summarise(hs = n(), .groups = "keep") 
qtl_hotspotrec <- qtl_rec %>%
  filter(type == "distant") %>%
  mutate(bin = cut_width(mrkstart, 1e6)) %>%
  distinct(bin, genename, .keep_all = TRUE) %>%
  group_by(bin, mrkchr) %>% 
  summarise(rec = n(), .groups = "keep")


qtl_hotspottotal <- qtl_hotspotuni %>%
  left_join(qtl_hotspotdev, by = c("bin", "mrkchr")) %>%
  left_join(qtl_hotspoths, by = c("bin", "mrkchr")) %>%
  left_join(qtl_hotspotrec, by = c("bin", "mrkchr")) %>%
  mutate(across(c(total, dev, hs, rec), replace_na, 0)) %>%
  select(bin, chr = mrkchr, dev, hs, rec, total)


write.table(qtl_hotspottotal, file = "hotspotbins2.txt", sep = "\t", quote = FALSE, row.names = F)
write.table(qtl_hotspottotal %>% arrange(desc(total)) %>% head(10), file = "hotspotbinstop2.txt", sep = "\t", quote = FALSE, row.names = F)
randomForestExplainer::explain_forest(rf.res)

##################### END ###########e###################################################################################
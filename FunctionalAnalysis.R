# Script:       FuncationalAnalysis.R
# Description:  In this script, we will explore a differential gene 
#               expression dataset comparing disease vs. control
#               samples, perform pathway enrichment analysis and 
#               build a PPI network
# Version: 2.0
# Last updated: 2026-01-30
# Author: mkutmon


# #############################################
# R INSTRUCTIONS
# #############################################

# * Lines that start with a # are comments
# * You can run a code line by placing the cursor in the line and clicking 
#   CTRL/Command + Enter

# #############################################
# R SETUP
# #############################################

# Here we install and load all required packages. 
if (!("BiocManager" %in% installed.packages())) { install.packages("BiocManager", update=FALSE) }
if (!("rstudioapi" %in% installed.packages())) { BiocManager::install("rstudioapi", update=FALSE) }
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db", update=FALSE) }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr", update=FALSE) }
if (!("EnhancedVolcano" %in% installed.packages())) { BiocManager::install("EnhancedVolcano", update=FALSE) }
if (!("readxl" %in% installed.packages())) { BiocManager::install("readxl", update=FALSE) }
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler", update=FALSE) }
if (!("enrichplot" %in% installed.packages())) { BiocManager::install("enrichplot", update=FALSE) }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3", update=FALSE) }
if (!("msigdbr" %in% installed.packages())) { BiocManager::install("msigdbr",update=FALSE) }
if (!("RColorBrewer" %in% installed.packages())) { BiocManager::install("RColorBrewer",update=FALSE) }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr",update=FALSE) }

library(rstudioapi)
library(org.Hs.eg.db)
library(dplyr)
library(EnhancedVolcano)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(msigdbr)
library(RColorBrewer)
library(readr)

loaded_pkgs <- c("BiocManager", "rstudioapi", "org.Hs.eg.db", "dplyr",
                 "EnhancedVolcano", "readxl", "clusterProfiler",
                 "enrichplot", "RCy3", "msigdbr", "RColorBrewer", "readr")

for (pkg in loaded_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("%s: %s\n", pkg, packageVersion(pkg)))
  } else {
    cat(sprintf("%s: not loaded\n", pkg))
  }
}

# We will set the working directory to the location where the current 
# script is located. This way, we can use relative file path locations. 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# We will create an output folder where all figures and files will be stored
out.folder <- "output/"
dir.create(out.folder)

# #############################################
# DATA EXPLORATION
# #############################################

# Go through the code, add comments and document the code to make
# sure you understand everything

data <- read_excel("GSE239914-differential-analysis.xlsx")
data <- data[,c(8,1,6,2,3,10,11)]

data.pc <- data[data$GeneType=="protein-coding",]
data.pc <- data.pc[!is.na(data.pc$GeneID),]
data.pc <- data.pc[,-6]

# 0.26 = 1.2 | 0.58 = 1.5 | 1 = 2 
log2fc.cutoff <- 1
pvalue.cutoff <- 0.05
degs <- data.pc[abs(data.pc$log2FoldChange) > log2fc.cutoff & data.pc$padj < pvalue.cutoff,]

write.table(degs, file=paste0(out.folder,"degs.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

# How many up and down-regulated genes are there? Try to write the code yourself.
genes.up <- degs[degs$log2FoldChange > log2fc.cutoff,]
genes.down <- degs[degs$log2FoldChange < -log2fc.cutoff,]

pmax <- max(-log10(data.pc$padj))
fcmax <- max(abs(data.pc$log2FoldChange))
EnhancedVolcano(data.pc, title = paste0("Volcanoplot (",nrow(degs), " DEGs)"), lab = data.pc$Symbol, x = "log2FoldChange", y = "padj", pCutoff = pvalue.cutoff, FCcutoff = log2fc.cutoff, labSize = 3, xlim = c(-fcmax,fcmax), ylim=c(0,pmax))

filename <- paste0(out.folder,"volcano-plot.png")
png(filename , width = 2000, height = 1500, res = 150)
EnhancedVolcano(data.pc, title = paste0("Volcanoplot (",nrow(degs), " DEGs)"), lab = data.pc$Symbol, x = "log2FoldChange", y = "padj", pCutoff = pvalue.cutoff, FCcutoff = log2fc.cutoff, labSize = 3, xlim = c(-fcmax,fcmax), ylim=c(0,pmax))
dev.off()

# #############################################
# PATHWAY ENRICHMENT 
# #############################################

# You can start with the WikiPathways collection but you can also look up other pathway 
# collections on MSigDb
genesets.wp <- msigdbr(species = "Homo sapiens", category= "C2", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, entrez_gene)

#genesets.wp <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, entrez_gene)

res.wp <- clusterProfiler::enricher(degs$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 400)
res.wp.df <- as.data.frame(res.wp)
write.table(res.wp.df, file=paste0(out.folder,"WP-Enrichment.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

res.wp.sim <- enrichplot::pairwise_termsim(res.wp)
treeplot(res.wp.sim, label_format = 0.1, showCategory = 80, cluster.params = list(label_words_n = 0))

filename <- paste0(out.folder,"WP_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(treeplot(res.wp.sim, label_format = 0.3, showCategory = nrow(res.wp.df), cluster.params = list(label_words_n = 0)))
dev.off()

# =================

res.wp.up <- clusterProfiler::enricher(genes.down$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 400)
res.wp.up.df <- as.data.frame(res.wp.up)
write.table(res.wp.up.df, file=paste0(out.folder,"WP-Enrichment-Up.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

res.wp.up.sim <- enrichplot::pairwise_termsim(res.wp.up)
treeplot(res.wp.up.sim, label_format = 0.3, showCategory = nrow(res.wp.df), cluster.params = list(label_words_n = 0))

filename <- paste0(out.folder,"WP_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(treeplot(res.wp.sim, label_format = 0.3, showCategory = nrow(res.wp.df), cluster.params = list(label_words_n = 0)))
dev.off()


# #############################################
# PATHWAY VISUALIZATION 
# #############################################

# Check if Cytoscape is running - keep it open and check what is happening while you are 
# running the code
cytoscapePing()

# Check if WikiPathways app is installed
if(!"name: WikiPathways, version: 3.3.10, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("WikiPathways")
}

# Open Pathway of interest - based on the res.wp.df, you can
# select pathways of interest
# Find the associated pathway identifier
# https://www.wikipathways.org/browse/table.html
# Make sure you select the ID of the human pathway
# Example: MAPK Cascade

pw.id <- "WP422"
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',pw.id)) 

toggleGraphicsDetails()

# load the data into Cytoscape (columns get added at the bottom)
loadTableData(data.pc, data.key.column = "EnsemblGeneID", table.key.column = "Ensembl")

# visualize the log2FC as a node fill color gradient
RCy3::setNodeColorMapping(table.column = 'log2FoldChange', mapping.type = 'c', table.column.values = c(-2,0,2), colors = paletteColorBrewerRdBu, default.color = '#FFFFFF', style.name = 'WikiPathways')

# Select significant genes and change border color
x <- RCy3::createColumnFilter('padj', 'padj', 0.05, "LESS_THAN")
RCy3::setNodeBorderColorBypass(x$nodes, new.colors = "#009900")
RCy3::setNodeBorderWidthBypass(x$nodes, new.sizes = 7)
RCy3::clearSelection()


# ==================================================================
# PPI network creation with the stringApp for Cytoscape
# ==================================================================

# make sure Cytoscape is running
RCy3::cytoscapePing()

if(!"name: stringApp, version: 2.0.3, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("stringApp")
}

degs.strict <- degs[abs(degs$log2FoldChange) > 1,]

query <- format_csv(as.data.frame(degs.strict$Symbol), col_names=F, escape = "double", eol =",")
commandsPOST(paste0('string protein query cutoff=0.7 newNetName="PPI network" query="',query,'" limit=0 species="Homo sapiens"'))


# =======================
# if you run into this error: Error:  reason: URI Too Long
# uncomment (remove #) and run the following line - create the network "manually" in Cytoscape, then run the rest of the code normally (check video!!)

# writeClipboard(query)
# =======================

RCy3::analyzeNetwork()
hist(RCy3::getTableColumns(columns = "Degree")$Degree, breaks=100)


# network topology
RCy3::createVisualStyle("centrality")
RCy3::setNodeLabelMapping("display name", style.name = "centrality")
colors <-  c ('#FFFFFF', '#DD8855')
setNodeColorMapping("Degree", c(0,60), colors, style.name = "centrality", default.color = "#C0C0C0")
setNodeSizeMapping("Degree", table.column.values = c(0,60), sizes = c(30,100), mapping.type = "c", style.name = "centrality", default.size = 10)
RCy3::setVisualStyle("centrality")
RCy3::toggleGraphicsDetails()

# data visualization# data visualization
RCy3::mapTableColumn("query term", species="Human", map.from = "HGNC", map.to = "Ensembl")
RCy3::loadTableData(data=data.pc, data.key.column = "EnsemblGeneID", table = "node", table.key.column = "Ensembl")
RCy3::createVisualStyle("log2FC vis")
RCy3::setNodeLabelMapping("display name", style.name = "log2FC vis")
control.points <- c (-5.0, 0.0, 5.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FoldChange", control.points, colors, style.name = "log2FC vis", default.color = "#C0C0C0")
RCy3::setVisualStyle("log2FC vis")
RCy3::lockNodeDimensions("TRUE", "log2FC vis")

write.csv(degs, file=paste0(out.folder,"degs.csv"), row.names = FALSE)

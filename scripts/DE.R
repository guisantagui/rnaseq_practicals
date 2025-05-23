################################################################################
# Bulk RNAseq preprocessing pipeline: Perform differential expression analysis.#
################################################################################

if (!require("BiocManager", quietly = T))
        install.packages("BiocManager")
if (!require("DESeq2", quietly = T))
        BiocManager::install("DESeq2", update = F)
if (!require(ggplot2, quietly = T)){
        install.packages("ggplot2", repos='http://cran.us.r-project.org')
}
if (!require(ggtext, quietly = T)){
        install.packages("ggtext", repos='http://cran.us.r-project.org')
}
if (!require(ggrepel, quietly = T)){
        install.packages("ggrepel", repos='http://cran.us.r-project.org')
}
if(!require(factoextra, quietly = T)){
        install.packages("factoextra", repos='http://cran.us.r-project.org')
}
if(!require(ggpubr, quietly = T)){
        install.packages("ggpubr", repos='http://cran.us.r-project.org')
}
if (!require("argparser", quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
if (!require("devtools", quietly = T)){
        install.packages("devtools", repos='http://cran.us.r-project.org')
}


library(DESeq2)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(factoextra)
library(ggpubr)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Given a directory where aligmnent BAM files are and a gtf annotation file, obtain a counts matrix.")

parser <- add_argument(parser = parser,
                       arg = c("--gene_counts",
                               "--sample_info",
                               "--outDir"),
                       help = c("CSV file with the gene counts.",
                                "CSV file with sample names in one column, and treatment info in another one.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if (!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Obtaining the top N contributor variables
# to the specified components (in our case PC1 and PC2)
getTopContrib <- function(PC, topN = 12, x = "PC1", y = "PC2"){
        contrib <- facto_summarize(PC, 
                                   "var", 
                                   axes = c(as.numeric(gsub("PC", "", x)),
                                            as.numeric(gsub("PC", "", y))))
        contrib <- contrib[order(contrib$contrib, decreasing = T), 
                           c("name", "contrib")]
        topContrib <- as.character(contrib$name[1:topN])
        return(topContrib)
}

# Do and plot a PCA.
plotPCACust <- function(df, scale = F, x = "PC1", y = "PC2", ntop = NA,
                        labs = F, topNFeats = NULL, biplot = F, color,
                        coord_fixed = T){
        if (scale){
                df <- apply(df, 1, function(x) (x - mean(x))/sd(x))
        }else{
                df <- t(df)
        }
        rv <- colVars(df)
        sel <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                        length(rv)))]
        
        if(scale){
                PC <- prcomp(df[, sel], scale. = F, center = F)
        }else{
                PC <- prcomp(df[, sel], scale. = F, center = T)
        }
        
        
        dat <- data.frame(obsnames=row.names(PC$x), PC$x)
        dat <- dat[, c("obsnames", x, y)]
        colDat <- data.frame(dds@colData)
        dat <- cbind.data.frame(dat,
                                colDat[match(dat$obsnames, colDat$sample), !colnames(colDat) %in% c("accession", "sample")])
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        datapc <- data.frame(varnames=rownames(PC$rotation), 
                             PC$rotation)
        mult <- min(
                (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
                (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
        datapc <- transform(datapc,
                            v1 = .7 * mult * (get(x)),
                            v2 = .7 * mult * (get(y))
        )
        datapc$x0 <- rep(0, nrow(datapc))
        datapc$y0 <- rep(0, nrow(datapc))
        if(!is.null(topNFeats)){
                varPlotFilt <- getTopContrib(PC, topN = topNFeats, x = x, y = y)
                datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
        }
        
        pcaPlt <- ggplot(dat, aes(PC1, PC2, color=!!sym(color), label = obsnames)) +
                geom_point(size=3) +
                xlab(sprintf("PC1 (%s %%)", propX)) +
                ylab(sprintf("PC1 (%s %%)", propY)) + 
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if (coord_fixed){
                pcaPlt <- pcaPlt + coord_fixed()
        }
        if (labs){
                pcaPlt <- pcaPlt +
                        geom_text_repel()
        }
        if (biplot){
                pcaPlt <- pcaPlt +
                        geom_text_repel(data=datapc, 
                                        aes(x=v1, y=v2, label=varnames), 
                                        color = "black", 
                                        size = 3,
                                        max.overlaps = 100) + 
                        geom_segment(data = datapc, aes(x=x0, 
                                                        y=y0, 
                                                        xend=v1, 
                                                        yend=v2, 
                                                        label = varnames),
                                     arrow = arrow(length=unit(0.2,"cm"),
                                                   type = "closed",
                                                   angle = 20), 
                                     alpha=0.75, 
                                     color="black", 
                                     size = 0.5)
        }
        return(pcaPlt)
}

# Plot p-value histogram
plot_pval_hist <- function(res_df, what_pval = "pvalue",
                           bins = NULL,
                           ylims = NULL){
        valVec <- res_df[, what_pval]
        valVec <- valVec[!is.na(valVec)]
        df <- data.frame(value = valVec)
        hst <- ggplot(data = df, mapping = aes(x = value)) +
                geom_histogram(bins = bins) +
                xlab(what_pval) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if (!is.null(ylims)){
                hst <- hst +
                        ylim(ylims)
        }
        return(hst)
}

# Obtain a volcano plot
getVolcPlot <- function(fitObjct, alpha = 0.05,
                        useLog2FC4Sign = T,
                        log2FC_thshld = 0.5, outputSignGenes = F,
                        labsThrshld = NULL,
                        topNGenes = NULL,
                        topNGenes_by = NULL){
        #labsThrshld <- c(4, 25)
        #fitObjct <- results
        alphaMinLog10 <- -log10(alpha)
        #compVec <- strsplit(colnames(fitObjct$coefficients), split = "-")[[1]]
        modelResults <- fitObjct
        
        volcPlotDF <- data.frame(proteins = rownames(modelResults),
                                 adj.P.Val = modelResults$padj,
                                 adj.pVal_minLog10 = -log10(modelResults$padj),
                                 log2FC = modelResults$log2FoldChange)
        volcPlotDF <- volcPlotDF[!is.na(volcPlotDF$adj.pVal_minLog10), ]
        
        volcPlotDF$sign <- volcPlotDF$adj.pVal_minLog10 >= alphaMinLog10
        if(useLog2FC4Sign){
                log2FC_bool <- volcPlotDF$log2FC >= log2FC_thshld | 
                        volcPlotDF$log2FC <= -1 *log2FC_thshld
                volcPlotDF$sign <- volcPlotDF$sign & log2FC_bool
        }
        volcPlotDF$sign <- c("not_significant",
                             "significant")[as.numeric(volcPlotDF$sign) + 1]
        
        volcPlotDF$sign <- factor(volcPlotDF$sign)
        volcPlotDF$labs <- volcPlotDF$sign == "significant"
        print(sprintf("There are %s significant genes with significance level of %s.",
                      sum(volcPlotDF$labs),
                      alpha))
        if(!is.null(labsThrshld)){
                volcPlotDF$labs <- volcPlotDF$labs & ((volcPlotDF$log2FC <= -labsThrshld[1] | volcPlotDF$log2FC >= labsThrshld[1]) | volcPlotDF$adj.pVal_minLog10 >= labsThrshld[2])
        }
        if(!is.null(topNGenes) & !is.null(topNGenes_by)){
                signDF <- volcPlotDF[volcPlotDF$sign == "significant", ]
                topGenes <- signDF$proteins[order(abs(signDF[, topNGenes_by]),
                                                  decreasing = T)][1:topNGenes]
                volcPlotDF$labs <- volcPlotDF$labs & (volcPlotDF$proteins %in% topGenes)
        }
        volcPlotDF_labs <- volcPlotDF[volcPlotDF$labs, ]
        volcPlotDF_subSet <- volcPlotDF[volcPlotDF$sign == "significant", ]
        
        #x_label <- sprintf("log2 FC (%s vs %s)", compVec[1], compVec[2])
        #x_label <- gsub("sick", "patients", x_label)
        #x_label <- gsub("ctrl", "controls", x_label)
        
        volcPlot <- ggplot(volcPlotDF, mapping = aes(x = log2FC,
                                                     y = adj.pVal_minLog10,
                                                     color = sign)) +
                geom_point() +
                scale_color_manual(values = c("significant" = "red",
                                              "not_sigificant" = "black")) +
                geom_hline(yintercept = alphaMinLog10, linetype = "dashed") +
                labs(x = "log2 FC (TFGB vs controls)",
                     y ="-log10(adj. p-val)") +
                geom_text_repel(data = volcPlotDF_labs,
                                aes(x = log2FC,
                                    y = adj.pVal_minLog10,
                                    label = proteins)) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_line(colour = "gray"), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, linewidth = 1),
                      legend.position = "none")
        if(useLog2FC4Sign){
                volcPlot <- volcPlot +
                        geom_vline(xintercept = log2FC_thshld,
                                   linetype = "dashed") +
                        geom_vline(xintercept = -1 * log2FC_thshld,
                                   linetype = "dashed")
        }
        if(outputSignGenes){
                signGenes <- volcPlotDF_subSet
                out <- list(significative = signGenes,
                            plot = volcPlot)
        }else{
                out <- volcPlot
        }
        return(out)
}

# Directory stuff
################################################################################
gene_counts_file <- "/Users/guillem.santamaria/Documents/postdoc/teaching_practical/RNAseq_pipe/results/preprocessing/counts/counts.csv"
samp_info_file <- "/Users/guillem.santamaria/Documents/postdoc/teaching_practical/RNAseq_pipe/data/sample_info.csv"
pathway_file <- "/Users/guillem.santamaria/Documents/postdoc/teaching_practical/RNAseq_pipe/data/kegg_pathways.csv"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/teaching_practical/RNAseq_pipe/results/DE/"

gene_counts_file <- parsed$gene_counts
samp_info_file <- parsed$sample_info
outDir <- addSlashIfNot(parsed$outDir)

createIfNot(outDir)

# Load data
################################################################################
gene_counts <- read.csv(gene_counts_file, row.names = 1)
samp_info <- read.csv(samp_info_file, row.names = 1)
pathways <- read.csv(pathway_file, row.names = 1)

# Filter low expressed genes
gene_counts <- gene_counts[rowSums(gene_counts) != 0, ]

print(sprintf("%s genomic features don't have zeros accross all the samples",
              nrow(gene_counts)))

#gene_counts <- gene_counts[, colnames(gene_counts) != "SRR5223570"]
#samp_info <- samp_info[samp_info$sample != "SRR5223570", ]

samp_info$rhl_pheno_2_cats <- factor(samp_info$rhl_pheno_2_cats, levels = c("not_producer", "producer"))

# Create DESeq object
################################################################################

# Change accession names to actual sample names
colnames(gene_counts) <- samp_info$sample[match(colnames(gene_counts),
                                                samp_info$accession)]
samp_info <- samp_info[match(colnames(gene_counts), samp_info$sample), ]

dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = samp_info,
                              design = ~rhl_pheno_2_cats)

# Run DESeq
################################################################################
dds <- DESeq(dds)

# Obtain results dataframe
results <- results(dds)
# Shrink logFC, given that we have small sample size
resLFC <- lfcShrink(dds, coef="rhl_pheno_2_cats_producer_vs_not_producer", type="apeglm")

# Save results dataframes
write.csv(as.data.frame(results),
          file = sprintf("%sDESeqRes.csv",
                         outDir))

write.csv(as.data.frame(resLFC),
          file = sprintf("%sDESeqRes_LFCShrunk.csv",
                         outDir))

pdf(file = sprintf("%sMAPlot.pdf", outDir), height = 5, width = 5)
DESeq2::plotMA(results)
dev.off()

pdf(file = sprintf("%sMAPlot_lfcShrunk.pdf", outDir), height = 5, width = 5)
DESeq2::plotMA(resLFC)
dev.off()

# Do PCA and plot it 
pcaPlot <- plotPCACust(assay(rlog(dds)),
                       scale = F,
                       labs = T,
                       biplot = T,
                       topNFeats = 5,
                       color = "rhl_pheno_2_cats",
                       coord_fixed = F)

ggsave(filename = sprintf("%spca_rlog.pdf", outDir),
       width = 5,
       height = 4,
       plot = pcaPlot)

# Do PCA and plot it 
pcaPlot <- plotPCACust(assay(rlog(dds)),
                       scale = F,
                       labs = T,
                       biplot = T,
                       topNFeats = 5,
                       color = "rhl_pheno_3_cats",
                       coord_fixed = F)

ggsave(filename = sprintf("%spca_rlog_3cats.pdf", outDir),
       width = 5,
       height = 4,
       plot = pcaPlot)

# Plot p-value histograms
pval_hist <- plot_pval_hist(results,
                            what_pval = "pvalue",
                            bins = 50,
                            ylims = c(0, 2500))
padj_hist <- plot_pval_hist(results,
                            what_pval = "padj",
                            bins = 50,
                            ylims = c(0, 2500))

pval_hists_comb <- ggarrange(pval_hist, padj_hist)

ggsave(filename = sprintf("%spval_hists.pdf", outDir),
       height = 4, width = 8,
       plot = pval_hists_comb)

# Do Volcano plots
volcPlot_top10Genes_logFC <- getVolcPlot(fitObjct = resLFC,
                                         topNGenes = 10,
                                         topNGenes_by = "log2FC",
                                         useLog2FC4Sign = F)
volcPlot_top10Genes_pAdj <- getVolcPlot(fitObjct = resLFC,
                                        topNGenes = 10,
                                        topNGenes_by = "adj.pVal_minLog10",
                                        useLog2FC4Sign = F)

volcPlotComb <- ggarrange(volcPlot_top10Genes_logFC,
                          volcPlot_top10Genes_pAdj)

ggsave(filename = sprintf("%svolcPlots.pdf", outDir),
       height = 4, width = 8,
       plot = volcPlotComb)

# Run functional enrichment
################################################################################
DEGs <- rownames(results[results$padj <= 0.05, ])
do_ORA <- function(DEGs, pathway_df, alternative = "greater", padj_method = "BH"){
        DEGs <- DEGs[DEGs %in% pathway_df$locus_tag]
        bg_genes <- unique(pathway_df$locus_tag)
        alternative = "greater"
        pws <- unique(pathway_df$pathway_id)
        pvals <- c()
        for (i in seq_along(pws)){
                pw <- pws[i]
                pw_genes <- unique(pathway_df$locus_tag[pathway_df$pathway_id == pw])
                DEGs_in_pw <- sum(DEGs %in% pw_genes)
                DEGs_not_pw <- sum(!DEGs %in% pw_genes)
                genes_in_pw <- length(unique(pw_genes)) - DEGs_in_pw
                genes_not_pw <- sum((!bg_genes %in% DEGs) & (!bg_genes %in% pw_genes))
                cont_tab <- matrix(c(DEGs_in_pw, DEGs_not_pw, genes_in_pw, genes_not_pw),
                                   nrow = 2)
                pval <- fisher.test(cont_tab, alternative = alternative)$p.value
                pvals <- c(pvals, pval)
        }
        out <- data.frame(pathway_id = pws,
                          pathway_name = pathway_df$pathway_name[match(pws, pathway_df$pathway_id)],
                          p_value = pvals,
                          p_adj = p.adjust(pvals, method = padj_method))
        return(out)
}

ora_res <- do_ORA(DEGs, pathway_df = pathways)
ora_res[ora_res$p_adj <= 0.05, ]

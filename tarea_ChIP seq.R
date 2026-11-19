##################################################
## 1. Direcciones de entrada y salida          ##
##################################################
#Guardar dirección donde se almacenan los datos de entrada
data_folder <- file.path(".")


#Crear carpeta para almacenar las salidas
results_folder <- file.path(".", "resultsABF3")
dir.create(results_folder, showWarnings = FALSE, recursive = TRUE)

##################################################
## 2. A. Anotación genómica de picos (RAMA CENTRAL) ##
##################################################

## 1. Cargar picos y anotación genómica
library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

ABF3.peaks <- readPeakFile(
  peakfile = file.path(data_folder, "GSE80564_AT4G34000_ABF3_ABA_optimal_narrowPeak_p16_reanotado.bed"),
  header = FALSE
)

# Número total de picos
length(ABF3.peaks)
head(ABF3.peaks)

## 2. Asignación a elementos genómicos          

ABF3.peak.annotation <- annotatePeak(
  peak = ABF3.peaks,
  TxDb = txdb,
  tssRegion = c(-1000, 1000))


## 3. QC y caracterización global             

## Distribución por categorías genómicas
png(file.path(results_folder, "ABF3_peak_annotation_pie.png"))
plotAnnoPie(ABF3.peak.annotation)
dev.off()

png(file.path(results_folder, "ABF3_peak_annotation_bar.png"))
plotAnnoBar(ABF3.peak.annotation)
dev.off()

## Distancia de picos al TSS
png(file.path(results_folder, "ABF3_peak_dist_to_TSS.png"))
plotDistToTSS(ABF3.peak.annotation)
dev.off()


## 4. Definición de genes diana                 

## Extraer genes asociados a picos
ABF3.peak.annotation.df <- as.data.frame(ABF3.peak.annotation)

## Filtrar aquellos genes cuyos picos se encuentran en la zona de su promotor
ABF3.target.genes <- ABF3.peak.annotation.df$geneId[ABF3.peak.annotation.df$annotation == "Promoter"]

length(ABF3.target.genes)

## Guardar lista
write.table(
  ABF3.target.genes,
  file = file.path(results_folder, "target_genes_from_ABF3.peaks.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

######################################################
## 2. B. Anclaje manual de picos y genes de ref.    ##
######################################################
ABF3.table <- read.table(file=file.path(data_folder, "GSE80564_AT4G34000_ABF3_ABA_optimal_narrowPeak_p16_reanotado.bed"),
                         header = F)
ABF3.table <- ABF3.table[,1:3]
head(ABF3.table)


ath.genes <- as.data.frame(genes(txdb))
head(ath.genes)


target.genes <- vector(mode = "character",length = nrow(ABF3.table))
peak.annotation <- vector(mode = "character",length = nrow(ABF3.table))


for(i in 1:nrow(ABF3.table))
{
  current.chr <- ABF3.table[i,1]
  current.start <- ABF3.table[i,2]
  current.end <- ABF3.table[i,3]
  
  complete.gene.body <- subset(ath.genes,
                               seqnames == current.chr &
                                 start >= current.start &
                                 end <= current.end)
  
  overlap.tss.pos.strand <- subset(ath.genes,
                                   seqnames == current.chr &
                                     start >= current.start &
                                     start <= current.end &
                                     strand == "+")
  
  overlap.tss.neg.strand <- subset(ath.genes,
                                   seqnames == current.chr &
                                     end >= current.start & end <= current.end &
                                     strand == "-")
  
  overlap.tes.pos.strand <- subset(ath.genes,
                                   seqnames == current.chr &
                                     end >= current.start &
                                     end <= current.end &
                                     strand == "+")
  
  overlap.tes.neg.strand <- subset(ath.genes,
                                   seqnames == current.chr &
                                     start >= current.start &
                                     start <= current.end &
                                     strand == "-")
  
  inside.gene.body <- subset(ath.genes, seqnames == current.chr &
                               start <= current.start &
                               end >= current.end)
  
  
  if(nrow(complete.gene.body) > 0)
  {
    target.genes[i] <- paste(unlist(complete.gene.body$gene_id),collapse = ",")
    peak.annotation[i] <- "gene.body"
  } else if(nrow(overlap.tss.pos.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tss.pos.strand$gene_id),collapse = ",")
    peak.annotation[i] <- "TSS"
  } else if(nrow(overlap.tss.neg.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tss.neg.strand$gene_id),collapse = ",")
    peak.annotation[i] <- "TSS"
  } else if(nrow(overlap.tes.pos.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tes.pos.strand$gene_id),collapse = ",")
    peak.annotation[i] <- "TES"
  } else if(nrow(overlap.tes.neg.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tes.neg.strand$gene_id),collapse = ",")
    peak.annotation[i] <- "TES"
  } else if(nrow(inside.gene.body) > 0)
  {
    target.genes[i] <- paste(unlist(inside.gene.body$gene_id),collapse = ",")
    peak.annotation[i] <- "internal.gene.body"
  } else
  {
    peak.annotation[i] <- "intergenic"
  }
}


freq.peak.annotation <- table(peak.annotation)
png(filename = file.path(results_folder, "freq.peak.annotation.pie.png"))
pie(freq.peak.annotation[c("gene.body",
                           "internal.gene.body",
                           "TSS",
                           "TES",
                           "intergenic")],
    labels = c("Full Gene Body",
               "Gene Body",
               "TSS",
               "TES",
               "Intergenic"),
    clockwise = T,
    border = "black",
    col=c("blue","cyan","green","red","grey"))

dev.off()


ABF3.table.peak.annotation <- cbind(ABF3.table[,1:3],
                                    data.frame(peak.annotation,target.genes))

head(ABF3.table.peak.annotation)

rownames(ABF3.table.peak.annotation) <- apply(X = ABF3.table.peak.annotation[,1:3],
                                              MARGIN = 1,
                                              FUN = function(x) {paste(x,collapse = "_")})


write.table(x = ABF3.table.peak.annotation,
            file = file.path(results_folder, "ABF3.table.peak.annotation.tsv"),
            quote = F,
            sep = "\t",
            row.names = F)


all.ABF3.genes <- unique(unlist(sapply(X =ABF3.table.peak.annotation$target.genes,
                                       FUN = function(x) {strsplit(x,",")[[1]]})))
length(all.ABF3.genes)
length(ABF3.target.genes)

##################################################
## 5. Enriquecimiento funcional (GO)             ##
##################################################

library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)

ABF3.enrich.go <- enrichGO(
  gene = ABF3.target.genes,
  OrgDb = org.At.tair.db,
  keyType       = "TAIR",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

write.table(
  as.data.frame(ABF3.enrich.go),
  file = file.path(results_folder, "GO_ABF3_enrichment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





png(filename = file.path(results_folder, "GO_ABF3_barplot.png"), width = 1200, height = 1000, res = 150)
barplot(ABF3.enrich.go,showCategory = 20)
dev.off()

png(filename = file.path(results_folder, "GO_ABF3_dotplot.png"), width = 1200, height = 1000, res = 150)
dotplot(ABF3.enrich.go,showCategory = 20)
dev.off()

png(filename = file.path(results_folder, "GO_ABF3_emapplot.png"), width = 1200, height = 1000, res = 150)
emapplot(pairwise_termsim(ABF3.enrich.go),showCategory = 20)
dev.off()

png(filename = file.path(results_folder, "GO_ABF3_cnetplot.png"), width = 1200, height = 1000, res = 150)
cnetplot(ABF3.enrich.go,showCategory = 20)
dev.off()

##################################################
## 6. KEGG                                    ##
##################################################

ABF3.enrich.kegg <- enrichKEGG(gene  = ABF3.target.genes,
                               organism = "ath",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05)
df.ABF3.enrich.kegg <- as.data.frame(ABF3.enrich.kegg)
head(df.ABF3.enrich.kegg)

write.csv(df.ABF3.enrich.kegg, 
          file = file.path(results_folder, "KEGG_ABF3_enrichment.tsv"),
          row.names = TRUE)

png(filename = file.path(results_folder, "KEGG_ABF3_barplot.png"), width = 1400, height = 1400, res = 150)
barplot(ABF3.enrich.kegg, showCategory = 20)
dev.off()

png(filename = file.path(results_folder, "KEGG_ABF3_dotplot.png"), width = 1400, height = 1400, res = 150)
dotplot(ABF3.enrich.kegg, showCategory = 20)
dev.off()


library(pathview)
## Vector nombrado: genes del reguloma de ABF3
gene.data <- rep(1, length(ABF3.target.genes))

## Conversión a ENTREZ
gene.conv <- bitr(
  ABF3.target.genes,
  fromType = "TAIR",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)

names(gene.data) <- gene.conv$ENTREZID


ABF3.enrich.kegg@result$ID

pathview(
  gene.data   = gene.data,
  pathway.id  = ABF3.enrich.kegg@result$ID[8],   # Carbon metabolism
  species     = "ath",
  gene.idtype = "entrez",
  kegg.native = TRUE
)



# Seqplotsscript
library(methods)
library()
library(IRanges)
library(BSgenome)
library(digest)
library(rtracklayer)
library(GenomicRanges)
library(seqplots)
#import file into r
#atac_wild_l3_ce11<- import("atac_wild_l3_ce11.bw")
#to change the format of the chr names
#seqlevels(atac_wild_l3_ce11)<-gsub("chr","",seqlevels(atac_wild_l3_ce11))
#loop to change format of the chr names
#filess<-read.table("./convseqname.txt", header=TRUE, stringsAsFactors = FALSE)
#for (i in 1:dim(filess)[1]) {
 # seqlevelchange<-gsub("chr","",seqlevels(i))
#}
#loop for feature import and change of chr names
seqPlotPath="/Users/imac/SeqPlots_data/files/"
featureFiles<-list.files(path =seqPlotPath, pattern = ".gtf")

feature<-read.table("./files/feature.txt", header=TRUE, stringsAsFactors = FALSE)
for (i in 1:length(featureFiles)) {
  data<-import(paste0(seqPlotPath,featureFiles[i]))
  seqlevels(data)<-paste0("chr",seqlevels(data))
  dataName<-gsub(".gtf","",featureFiles[i])
  export(data, paste0(seqPlotPath,dataName, "_UCSC.gtf"), "gtf")
}


#to change format of the chr namess
#genesTSSX<-import("./files/genesTSSx.gtf")
#seqlevels(genesTSSX)<-gsub("x","chrX", seqlevels(genesTSSX))
#genesTSSA<-import("./files/genesTSSA.gtf")
#seqlevels(genesTSSA)<-paste0("chr",seqlevels(genesTSSA))
#export(genesTSSX, "genesTSSchrX.gtf", "gtf")
#export(genesTSSA, "genesTSSchrA.gtf", "gtf")
#
#to load files into seqplots
atac_wild_l3_ce11<- import("atac_wild_l3_ce11.bw")
add_local_file(atac_wild_l3_ce11, name =deparse(substitute(atac_wild_l3_ce11)), file_genome = "ce11", file_user = "Bolaji", file_comment = "atac", root = file.path(path.expand("~"), "SeqPlots_data"))
#to process genomic signals from tracks or motif data
plot2<-getPlotSetArray(tracks = "/Users/imac/SeqPlots_data/files/enrichment.bw",features = c("./files/genesTSSX.gtf", "./files/genesTSSA.gtf"), refgenome = "ce11", bin = 10L, rm0 = FALSE, ignore_strand = TRUE, xmin = 1000L, xmax = 1000L, xanchored = 2000L, type = "af", add_heatmap = TRUE, verbose = TRUE, stat = "median", lvl1m = message, lvl2m = message)
plot(plot2, main='supercoiling around genes', xlab='Relative position[bp]',ylab='log2 bTMP',error.estimates = FALSE,cex.axis = 14, cex.lab = 16, cex.main = 20, cex.legend = 10, colvec = c("darkgreen","blue3"), ylim=c(-0.2,0.4))

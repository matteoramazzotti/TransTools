#!/usr/bin/R
#(c) Matteo Ramazzotti - matteo.ramazzotti_at_unifi.it
args = commandArgs(trailingOnly=TRUE)

fmRNA <- args[1]
fmiRNA <- args[2]
fmRNADE <- args[3]
fmiRNADE <- args[4]
fANN <- args[5]
title<-args[6]

#comment the following lines to run your data
#the following linea are for test data only
fmRNA <- "test/mRNA.data"
fmiRNA <- "test/miRNA.data"
fANN <- "test/Mannot.txt"
fmRNADE <- "test/mrna.de.txt"
fmiRNADE <- "test/mirna.de.txt"
title<- "example"

cat("Running mM_Net.R with arguments:",fmRNA,fmiRNA,fmRNADE,fmiRNADE,fANN,title,"\n",sep="\n");

cat("Loading MirTarBase...")
#mirtarbase
mirdb1<-read.table(pipe("zcat test/mirtarbase-hsa.ok.txt"),sep="\t",quote="")
m<-gsub("-5p","",mirdb1$V2)
m<-gsub("-3p","",m)
mirdb1[,2]<-m
mtb_mix<-unique(paste(mirdb1$V3,mirdb1$V2,sep="@"))
cat(length(mtb_mix),"records\n")

cat("Loading PITA...")
#pita
mirdb3<-read.table(pipe("zcat test/PITA_targets_hg18_3_15_TOP.tab.gz | cut -f2,3 | sort | uniq"),sep="\t",quote="")
m<-gsub("-5p","",mirdb3$V2)
m<-gsub("-3p","",m)
pita_mix<-unique(paste(mirdb3$V1,m,sep="@"))
cat(length(pita_mix),"records\n")

cat("Merging databases...")
library(gplots)
v<-venn(list("pita"=pita_mix,"mirtatbase"=mtb_mix),show.plot=FALSE)
sel<-as.vector(attr(v, "intersections")[[3]])
cat(length(sel),"records remained\n")

cat("Loading annotations...")
#load annotation (mainly entrez2symbol)
ann<-read.table(file=fANN,sep="\t",row.names=1,header=T)
cat(dim(ann)[1],"entries\n")

cat("Loading miRNA data...")
#load mirna data matrix (row names are mirna IDs)
mirna<-read.table(file=fmiRNA,sep="\t",row.names=1,header=T)
#pita does not always report 3p 5p, so they must be striped
m<-gsub("-5p","",rownames(mirna))
m<-gsub("-3p","",m)
#aggregation of possible 3p 5p pairs
a<-aggregate(mirna,by=list(m),FUN=median)
mirna<-a[,-1]
rownames(mirna)<-a[,1]
cat(dim(mirna)[1],"entries\n")

cat("Loading mRNA data...")
#load mrna data matrix (row names are ENTREZ ID)
mrna<-read.table(file=fmRNA,sep="\t",row.names=1,header=T)
#convert rownames of mRNA from entrez to symbol
rownames(mrna)<-ann[as.character(rownames(mrna)),1]
cat(dim(mrna)[1],"entries\n")

cat("Loading mRNA DE...")
#load mRNA DE of exp2 (LTED vs MCFT-)
mrnaDE<-read.table(file=fmRNADE,sep="\t",header=F)
#only significant probes are retained
#mrnaDE<-mrnaDE[mrnaDE[,1]<=0.05,c(2,3,4)]
#genewise aggregation of probes
mrnaDE<-aggregate(mrnaDE[,2],by=list(mrnaDE[,1]),FUN=median)
#then filter for logFC of gene > 1
mrnaDE<-mrnaDE[abs(mrnaDE[,2]) > 1,]
mrnaDE<-data.frame(mrnaDE[,-1],row.names=mrnaDE[,1])
names(mrnaDE)<-"logFC"
#rownames(mrnaDE)<-ann[as.character(rownames(mrnaDE)),1]
rows<-ann[as.character(rownames(mrnaDE)),1]
mrnaDE<-data.frame("logFC"=mrnaDE[!is.na(rows),],row.names=rows[!is.na(rows)])
#names(mrnaDE)<-rows[!is.na(rows)]
cat(dim(mrnaDE)[1],"entries\n")

cat("Loading miRNA DE...")
#load miRNA DE of exp2  (LTED vs MCFT-)
#procedure: read file, filter out probes with adj.p < 0.05 and abs(logFC) < 1, rebuild data.frame with mirnames as rownames
mirnaDE<-read.table(file=fmiRNADE,sep="\t",header=T)
#mirnaDE<-mirnaDE[mirnaDE[,5]<=0.05 & abs(mirnaDE[,6])>= 1,c(1,6)]
#pita does not always report 3p 5p, so they must be striped in DE
m<-gsub("-5p","",mirnaDE[,1])
m<-gsub("-3p","",m)
#aggregation of possible 3p 5p pairs
a<-aggregate(mirnaDE[,2],by=list(m),FUN=median)
mirnaDE<-data.frame("logFC"=a[,2],row.names=a[,1])
cat(dim(mirnaDE)[1],"entries\n")

cat("Building network...\n")
cnt<-0
out<-NULL
for (i in 1:length(sel)) {
	a<-unlist(strsplit(sel[i],"@"))
	#exact match, returns position or false
	Md<-a[1] %in% rownames(mrnaDE)
	md<-a[2] %in% rownames(mirnaDE)
	M<-a[1] %in% rownames(mrna)
	m<-a[2] %in% rownames(mirna)
	if (Md && md && M && m) {
		Md1<-match(a[1],rownames(mrnaDE))
		md1<-match(a[2],rownames(mirnaDE))
		M1<-match(a[1],rownames(mrna))
		m1<-match(a[2],rownames(mirna))
		cnt<-cnt+1
		c<-cor.test(t(mrna[M1,]),t(mirna[m1,]))
		cat(cnt,"/",i,":",a[1],a[2],as.numeric(c$estimate),c$p.value,"\n")
		out<-rbind(out,data.frame(i,a[1],a[2],mrnaDE[Md1,],mirnaDE[md1,],as.numeric(c$estimate),c$p.value))
	}
}
padj<-p.adjust(out[,7],method="BH")
out<-cbind(out,"p.adj"=padj)
names(out)<-c("index","gene","miRNA","Genelog2FC","miRlog2FC","Cor","p.val","adj.p.val")
network<-out[out[,6]<0 & out[,8]<0.05,]
cat(dim(network)[1], "relationships detected\n")
write.table(file=paste(title,"network.txt",sep="."),network[,-1],sep="\t",row.names=F,quote=F)
#write.table(file=paste(title,"nodes2.txt",sep="."), rbind(mrnaDE,mirnaDE),sep="\t")

#### Network with igraph (see http://igraph.org/r/doc/ and http://kateto.net/network-visualization)
library(igraph)

network<-read.table(file=paste(title,"network.txt",sep="."),sep="\t",header=T)

nodes<-c(as.vector(network[,1]),as.vector(network[,2]))
edges<-data.frame(network[,1],gsub("hsa-","",network[,2]))
net <- graph_from_data_frame(d=edges,directed=F)

#node color by FC
Mnodes<-network[!duplicated(network[,1]),c(1,3)]
names(Mnodes)<-"logFC"
mnodes<-network[!duplicated(network[,2]),c(2,4)]
names(mnodes)<-"logFC"
nodes<-rbind(Mnodes,mnodes)
row.names(nodes)<-gsub("hsa-","",nodes[,1])
n<-as_ids(V(net))
fcs<-NULL
for (i in 1:length(n)) {
	fcs<-c(fcs,nodes[match(n[i],rownames(nodes)),2])
}
cr<-colorRamp(c("red","white","green"))
minv=-10
maxv=10
scale<-(fcs-minv)/(maxv-minv)
nodecol<-cr(scale)
nodecol<-rgb(nodecol,maxColorValue = 255)

#edge color by cor
e<-abs(network$Cor)
cr<-colorRamp(c("red","white"))
scale<-(e-min(e))/(max(e)-min(e))
edgecol<-cr(scale)
edgecol<-rgb(edgecol,maxColorValue = 255)

type<-rep(NA,length(as_ids(V(net))))
type[grep("miR",as_ids(V(net)))]<-1
type[grep("miR",as_ids(V(net)),invert=T)]<-0
type[grep("let",as_ids(V(net)))]<-1
V(net)$type<-type

nodesize<-rep(5,length(as_ids(V(net))))

#specifications
vlcol<-rep("azure4",length(as_ids(V(net))))
vlcol[grep("miR-23b",as_ids(V(net)))]<-"black"
vlcol[grep("SLC",as_ids(V(net)))]<-"black"

nodesize<-rep(5,length(as_ids(V(net))))
nodesize[grep("miR-23b",as_ids(V(net)))]<-15
nodesize[grep("SLC6A14",as_ids(V(net)))]<-15

nodelabelsize<-rep(0.7,length(as_ids(V(net))))
nodelabelsize[grep("miR-23b",as_ids(V(net)))]<-1.2
nodelabelsize[grep("SLC6A14",as_ids(V(net)))]<-1.2


pdf(file=paste(title,"network.pdf",sep="."),width=10,height=10)
set.seed(123)
l <- do.call("layout_with_fr", list(net))
plot(net,layout=l,
vertex.color=nodecol,
vertex.frame.color="white",
vertex.size=nodesize,
edge.lty=1,
vertex.label.family = "sans",
edge.color=edgecol
)
dev.off()

#vertex.label.color=vlcol,
#vertex.label.cex=nodelabelsize,





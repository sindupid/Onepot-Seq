#! -S /usr/bin/R

##Seurat version 2.3.4
# library(Seurat)

##Seurat version 3.1.5
library(Seurat)
library(Matrix)
library(gplots)
library(ggrepel)
library(plyr)
library(dplyr)
library(reshape)
library(cowplot)
library(viridis)
theme_set(theme_cowplot(7))



printswitch<-FALSE
RD<-"/Data2/SDJ/200615_NextSeq/rdata/"

######################## Generate Seurat object ########################

totsmpnames<-c()
WD <- "/Path/to/dge/files/"

## sampleprefix in dge files
Samples <- c('pbmc_1','pbmc_2','pbmc_3','pbmc_4')

## New labels for samples
newsmpnames<- c("Hash1","Hash2","Hash3","CD3")
pbmclist<-readSeurat(WD, Samples, newsmpnames, list())
totsmpnames<-c(totsmpnames,newsmpnames)


### data merge
pbmc <- merge(x= pbmclist[[1]],y= pbmclist[2:length(totsmpnames)], add.cell.ids=totsmpnames)


######################## Single cell data QC ########################
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc$orig.ident <-factor(pbmc$orig.ident, levels=totsmpnames)

pbmc_ori<-pbmc



######################## Cell filtering ########################
pbmc<-pbmc_ori
pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 15)
pbmc <- subset(pbmc, subset = percent.mt > 1)

clipr::write_clip(as.data.frame(table(pbmc[["orig.ident"]])))


######################## Normalize & Findvariablegenes ########################
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)

######################## Scaling data (SCTransform) ########################
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


######################## PCA analysis and determine dimensions ########################
pbmc <- RunPCA(object = pbmc)

pbmc<-pbmc_bfc
######################## Clustering & Visualization ########################
res=0.5
dim=1:15

pbmc <- FindNeighbors(pbmc, dims = dim)
pbmc <- FindClusters(pbmc, resolution = res)
pbmc <- RunUMAP(pbmc, dims = dim)

p1<-DimPlot(pbmc, reduction = "umap")
p2<-DimPlot(pbmc, reduction = "umap", group.by="orig.ident")
CombinePlots(list(p1,p2),ncol=2)


pbmc<-subset(pbmc, subset = seurat_clusters!=6)




DimPlot(pbmc, label = TRUE) + NoLegend()


genelist=c("IL7R", "CCR7","S100A4", "CD14","MS4A1","CD8A","FCGR3A","GNLY","FCER1A","CD3E","CD3D", "CD19", "AGR2")


CombinePlots(lapply(FeaturePlot(pbmc, features=genelist, combine=FALSE), function(x){x+NoAxes()+NoLegend()}))

CombinePlots(lapply(FeaturePlot(pbmc, features=genelist, slot="scale.data", combine=FALSE), function(x){x+NoAxes()+NoLegend()}))


CombinePlots(lapply(VizDimLoadings(pbmc, dims = 1:15, reduction = "pca", combine=FALSE), function(x){x+theme_cowplot(10)}))

VlnPlot(pbmc, features=genelist, slot="scale.data")




# CombinePlots(lapply(FeaturePlot(pbmc, features=genelist, combine=FALSE), function(x){x+NoAxes()+NoLegend()}))

genelist=c("IL7R","S100A4", "CD14","MS4A1","CD8A","GNLY","CD3D", "AGR2")
p<-CombinePlots(ncol=4,lapply(FeaturePlot(pbmc, pt.size=0.1,features=genelist, combine=FALSE,
cols=c("lightgrey", "red")),
function(x){x+NoAxes()+
# NoLegend()+
theme(plot.title=element_text(face="plain"), panel.border=element_rect(color="black"))}))

save_plot("/Data2/SDJ/script/Snowmanseq/image/supfig_pbmcexpsamples.png",plot=p,
base_height=4, base_width=10)






# Cluster ID	Markers	Cell Type
# 0	IL7R, CCR7	Naive CD4+ T
# 1	IL7R, S100A4	Memory CD4+
# 2	CD14, LYZ	CD14+ Mono
# 3	MS4A1	B
# 4	CD8A	CD8+ T
# 5	FCGR3A, MS4A7	FCGR3A+ Mono
# 6	GNLY, NKG7	NK
# 7	FCER1A, CST3	DC
# 8	PPBP	Platelet



######################## Save Seurat object ########################

# saveRDS(pbmc, file = paste(RD, "Ovarian_combined:15_0.6.rds", sep=''))





######################## Find all clusters Markers ########################
Idents(pbmc)<-"celltype"


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()+theme(plot.margin = unit(c(0,3,0,0), "cm"))
DoHeatmap(pbmc, features = top10$gene) + NoLegend()+theme(plot.margin = unit(c(0,3,0,0), "cm"), size=3)

######################## cell type assign ########################

# mapfrom<-c('0', '1', '2', '3', '4', '5')
# mapto1<-c('CD4+ T-cell', 'Monocyte', 'CD8+ T-cell', 'NK-cell', 'B-cell', 'MKI67+')
# celltypeorder<-c('CD4+ T-cell', 'CD8+ T-cell', 'NK-cell', 'Monocyte','B-cell', 'MKI67+')


mapfrom<-c('0', '1', '2', '3', '4', '5')
mapto1<-c('CD4+ T-cell', 'CD8+ T-cell', 'H3122', 'Monocyte','NK-cell','B-cell')
celltypeorder<-c('CD4+ T-cell', 'CD8+ T-cell', 'NK-cell', 'Monocyte','B-cell', 'H3122')



pbmc$celltype<-factor(mapvalues(pbmc$seurat_clusters, mapfrom, mapto1), levels=celltypeorder)
# pbmc$detailtype<-factor(mapvalues(pbmc$seurat_clusters, mapfrom, mapto2), levels=mapto2[order(mapto2)])

p<-DimPlot(pbmc, reduction = "umap", group.by="celltype", label=TRUE) + NoLegend()+
scale_color_brewer(palette="Set2")



labeldat<-DimPlot(pbmc, group.by="celltype")$data%>%
group_by(celltype)%>%
summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2))%>%
as.data.frame

p<-ggplot(DimPlot(pbmc, group.by="celltype")$data, aes(UMAP_1,UMAP_2, color=celltype))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
geom_text_repel(data=labeldat, aes(label=celltype),box.padding = 0.5, max.overlaps = Inf, color="black",
size=2.4)+
NoLegend()+NoAxes()+
theme(panel.border = element_rect(colour = "black", fill = NA), plot.title=element_blank())

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3a_umapplot.png",plot=p,
base_height=1.8, base_width=1.8)


p<-ggplot(DimPlot(pbmc, group.by="celltype")$data, aes(UMAP_1,UMAP_2, color=celltype))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
geom_text_repel(data=labeldat, aes(label=celltype),box.padding = 0.8, max.overlaps = Inf, color="black",
size=6,segment.size = 0)+
NoLegend()+NoAxes()

theme(panel.border = element_rect(colour = "black", fill = NA), plot.title=element_blank())

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3a_umapplot.png",plot=p,
base_height=3.6, base_width=3.6)


p<-ggplot(DimPlot(pbmc, group.by="celltype")$data, aes(UMAP_1,UMAP_2, color=celltype))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
geom_text_repel(data=labeldat, aes(label=celltype),box.padding = 1, max.overlaps = Inf, color="black",
size=6,segment.size = 0)+
NoLegend()+NoAxes()

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3a_umapplot.png",plot=p,
base_height=4.3, base_width=4.3)





p<-ggplot(DimPlot(pbmc, group.by="celltype")$data, aes(UMAP_1,UMAP_2, color=celltype))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
geom_text_repel(data=labeldat, aes(label=celltype),box.padding = 1, max.overlaps = Inf, color="black",
size=6,segment.size = 0)+theme_cowplot(15)+theme(panel.border=element_rect(color="black"))+
NoLegend()


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3a_umapplotnew.png",plot=p,
base_height=4.3, base_width=4.3)







p<-ggplot(DimPlot(pbmc, group.by="celltype", split.by="orig.ident")$data %>% filter(is.element(orig.ident, c("1_Hash1","2_Hash","3_Hash"))),
aes(UMAP_1,UMAP_2, color=celltype))+
scale_x_continuous(position = "top", limits=c(-1,11))+
scale_y_continuous(limits=c(-8,12))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
facet_wrap(~orig.ident)+theme(strip.text=element_blank(), strip.background=element_blank(),
axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(),
axis.line.y=element_blank(),panel.spacing = unit(1, "lines"))+NoLegend()


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3b_umapdonorsplit.png",plot=p,
base_height=1, base_width=3)






p<-ggplot(DimPlot(pbmc, group.by="celltype", split.by="orig.ident")$data %>% 
filter(is.element(orig.ident, c("1_Hash1","1_Hash2","1_Hash3"))),
aes(UMAP_1,UMAP_2, color=celltype))+
scale_x_continuous(position = "top", limits=c(-1,11))+
scale_y_continuous(limits=c(-8,12))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
facet_wrap(~orig.ident)+theme(strip.text=element_blank(), strip.background=element_blank(),
axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(),
axis.line.y=element_blank(),panel.spacing = unit(1, "lines"))+NoLegend()


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3b_umapreplicatesplit.png",plot=p,
base_height=1, base_width=3)





pbmc$smpnames<- pbmc$orig.ident

levels(pbmc$smpnames)<-gsub("3_","Donor C ", gsub("2_","Donor B ", gsub("1_","Donor A ",levels(pbmc$smpnames))))



p<-ggplot(DimPlot(pbmc, group.by="celltype", split.by="smpnames")$data,
aes(UMAP_1,UMAP_2, color=celltype))+
geom_point(size=0.1)+scale_color_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
facet_wrap(~smpnames)+theme(legend.title=element_blank(),strip.background=element_rect(fill = NA), panel.border=element_rect(color="black"))


save_plot("/Data2/SDJ/script/Snowmanseq/image/supfig_pbmcfacetsamples.png",plot=p,
base_height=4, base_width=5)















Idents(pbmc)<-"celltype"
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


cd8mark<-FindMarkers(pbmc, ident.1 = "CD8+ T-cell", ident.2="CD4+ T-cell",  min.pct=0.1, logfc.threshold=0.1)
cd8mark_top10 <- rbind(cd8mark %>% top_n(n = 20, wt = -avg_logFC),cd8mark %>% top_n(n = 20, wt = avg_logFC))


library(pheatmap)
library(grid)

expdat<-GetAssayData(object = pbmc, slot = "scale.data")
expdat<-expdat[unique(top10$gene)[is.element(unique(top10$gene), rownames(expdat))],]
expdat<-expdat[,order(pbmc$celltype)]

expdat[expdat>10]<-10
expdat[expdat<(-1)]<--1

p<-pheatmap(expdat,
show_colnames=FALSE,
cluster_rows=FALSE,
cluster_cols=FALSE,
color=magma(n=100),
annotation_col=pbmc@meta.data[c("celltype")],
annotation_colors=list(celltype=c('CD4+ T-cell'="#99B898",'CD8+ T-cell'="#FECEA8",'NK-cell'="#FF847C",'Monocyte'="#E84A5F",'B-cell'="#2A363B",'H3122'="#B2AC88")),
annotation_names_col=FALSE, return=TRUE, annotation_legend=FALSE)

png(filename="/Data2/SDJ/script/Snowmanseq/image/figure3b_celltypeheatmap.png", units="in",res=2000,
height=2, width=4.4)
add.flag(p,kept.labels=c("IL7R", "CCL5","GNLY","GZMB","GZMA","CD14","LYZ", "MS4A1","CCND1","KRT18"),repel.degree=1)
dev.off()







# pbmcsubfilt<-pbmc
# mapfrom<- c("1_Hash1","1_Hash2","1_Hash3","1_CD3","2_CD8","2_Hash","3_CD3","3_CD4","3_CD14","3_CD16","3_CD19","3_Hash","3_Spike1%","3_Spike0.1%")
# mapto<-c("Hash1\tDonor A\t","Hash2\tDonor A\t","Hash3\tDonor A\t","CD3\tDonor A\t","CD8\tDonor B\t","Hash\tDonor B\t","CD3\tDonor C\t","CD4\tDonor C\t","CD14\tDonor C\t","CD16\tDonor C\t","CD19\tDonor C\t","Hash\tDonor C\t","filtered","filtered")
# maporder<-c("Hash1\tDonor A\t","Hash2\tDonor A\t","Hash3\tDonor A\t","Hash\tDonor B\t","Hash\tDonor C\t","CD3\tDonor A\t","CD3\tDonor C\t","CD4\tDonor C\t","CD8\tDonor B\t","CD14\tDonor C\t","CD16\tDonor C\t","CD19\tDonor C\t","filtered")
# pbmcsubfilt$newident<-factor(mapvalues(pbmcsubfilt$orig.ident, mapfrom, mapto), levels=maporder)
# pbmcsubfilt<-subset(x = pbmcsubfilt, subset=orig.ident!='3_Spike1%'&orig.ident!='3_Spike0.1%')
# maporder<-c("Hash1\tDonor A\t","Hash2\tDonor A\t","Hash3\tDonor A\t","Hash\tDonor B\t","Hash\tDonor C\t","CD3\tDonor A\t","CD3\tDonor C\t","CD4\tDonor C\t","CD8\tDonor B\t","CD14\tDonor C\t","CD16\tDonor C\t","CD19\tDonor C\t")
# pbmcsubfilt$newident<-factor(pbmcsubfilt$newident, levels=maporder)



pbmcsubfilt<-pbmc
mapfrom<- c("1_Hash1","1_Hash2","1_Hash3","1_CD3","2_CD8","2_Hash","3_CD3","3_CD4","3_CD14","3_CD16","3_CD19","3_Hash","3_Spike1%","3_Spike0.1%")
mapto<-c("Hash1\tA\t","Hash2\tA\t","Hash3\tA\t","CD3\tA\t","CD8\tB\t","Hash\tB\t","CD3\tC\t","CD4\tC\t","CD14\tC\t","CD16\tC\t","CD19\tC\t","Hash\tC\t","filtered","filtered")
maporder<-c("Hash1\tA\t","Hash2\tA\t","Hash3\tA\t","Hash\tB\t","Hash\tC\t","CD3\tA\t","CD3\tC\t","CD4\tC\t","CD8\tB\t","CD14\tC\t","CD16\tC\t","CD19\tC\t","filtered")
pbmcsubfilt$newident<-factor(mapvalues(pbmcsubfilt$orig.ident, mapfrom, mapto), levels=maporder)
pbmcsubfilt<-subset(x = pbmcsubfilt, subset=orig.ident!='3_Spike1%'&orig.ident!='3_Spike0.1%')
maporder<-c("Hash1\tA\t","Hash2\tA\t","Hash3\tA\t","Hash\tB\t","Hash\tC\t","CD3\tA\t","CD3\tC\t","CD4\tC\t","CD8\tB\t","CD14\tC\t","CD16\tC\t","CD19\tC\t")
pbmcsubfilt$newident<-factor(pbmcsubfilt$newident, levels=maporder)


p<-ggplot(pbmcsubfilt@meta.data, aes(x=newident, fill=celltype))+
geom_bar(stat="count", position="fill", color="black",size=0.2)+
labs(x="",y="Relative frequency", fill="")+
scale_y_continuous(breaks=c(0,1),labels=scales::percent)+
scale_fill_manual(values=c("#99B898","#FECEA8","#FF847C","#E84A5F","#2A363B","#B2AC88"))+
coord_flip()+scale_x_discrete(limits = rev(levels(pbmcsubfilt$newident)))+
theme(legend.position="bottom")+guides(fill=guide_legend(nrow=2,byrow=FALSE))


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3e_celltype composition.png",plot=p,
base_height=1.8, base_width=3)

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3e_celltype composition.png",plot=p,
base_height=2, base_width=2.4)




a<-table(pbmcsubfilt$orig.ident, pbmcsubfilt$celltype)[1:12,1:5]
a<-a/rowSums(a)
b<-rbind("1_Hash"=colMeans(a[1:3,]), a[c(6,12),])
fcdata<-log2(a/b[c(1,1,1,1,2,2,3,3,3,3,3,3),]+1)

# p<-pheatmap(fcdata[c(4,7,8,5,9,10,11),],
# color=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="Purples"))(100),
# cluster_rows=FALSE,
# cluster_cols=FALSE, display_numbers=TRUE,
# legend_breaks=c(0,1,2,3),legend_labels=c("0","1","2","3"),
# number_color="black",border_color="black",
# return=TRUE, angle_col="45", labels_row=c("\tCD3 \tDonor A","\tCD3 \tDonor C","\tCD4 \tDonor C","\tCD8 \tDonor B","\tCD14\tDonor C","\tCD16\tDonor C","\tCD19\tDonor C")
# )

p<-pheatmap(fcdata[c(4,7,8,5,9,10,11),],
color=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="Purples"))(100),
cluster_rows=FALSE,
cluster_cols=FALSE, display_numbers=TRUE,
legend_breaks=c(0,1,2,3),legend_labels=c("0","1","2","3"),
number_color="black",border_color="black",
return=TRUE, angle_col="45", labels_row=c("\tCD3 \tA","\tCD3 \tC","\tCD4 \tC","\tCD8 \tB","\tCD14\tC","\tCD16\tC","\tCD19\tC")
)


png(filename="/Data2/SDJ/script/Snowmanseq/image/figure3b_fcheatmap.png", units="in",res=300,
height=2.5, width=4)
p
dev.off()




a<-table(pbmc$orig.ident, pbmc$celltype)[13:14,]
a<-a/rowSums(a)
a<-log2(a[,"H3122"]/c(0.01,0.001)+1)
a<-data.frame(sample=names(a),value=as.vector(a),newlabel=c("1.0%","0.1%"))

p<-ggplot(a, aes(x=newlabel,y=value))+geom_bar(stat="identity", color="black",fill="#353a41")+coord_flip()+
labs(x="",y=expression("Enrichment (log"[2]*"(1+Fold change))"))

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3e_h3122enrich.png",plot=p,
base_height=0.6, base_width=1.8)








p1<-DimPlot(pbmc, pt.size=0.1,reduction = "umap", label=FALSE, sizes.highlight=0.5,
cells.highlight=colnames(subset(x = pbmc, subset=orig.ident=="3_Spike1%")), cols.highlight="#CF0029") + NoLegend()+NoAxes()

p2<-DimPlot(pbmc, pt.size=0.1,reduction = "umap", label=FALSE, sizes.highlight=0.5,
cells.highlight=colnames(subset(x = pbmc, subset=orig.ident=="3_Spike0.1%")), cols.highlight="#CF0029") + NoLegend()+NoAxes()

p<-CombinePlots(list(p1,p2))

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure3e_h3122highlight.png",plot=p,
base_height=3.6, base_width=7.2)










p<-ggplot(substatdat, aes(x=sample,y=value,fill=Var.2))+
geom_bar(stat="identity", position="fill",color="black",size=0.2)+
scale_fill_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
labs(x="",y="Relative frequency",fill="")+
scale_y_continuous(breaks=c(0,1),labels=scales::percent)+
theme(strip.background = element_blank(), axis.ticks.x=element_blank(),legend.position="top")+
coord_flip()











CombinePlots(lapply(FeaturePlot(pbmc, features=genelist, combine=FALSE), function(x){x+NoAxes()+NoLegend()}))







DimPlot(pbmc, reduction = "umap", group.by="celltype", split.by = "orig.ident")+
scale_color_brewer(palette="Set2")


VlnPlot(pbmc, features=genelist, group.by="celltype", assay="SCT")








p<-DimPlot(pbmc, rduction = "umap", group.by="celltype",label=FALSE, sizes.highlight=4,
cells.highlight=colnames(subset(x = pbmc, subset = celltype == "Tcell")), cols.highlight="black") + NoLegend()


p<-DimPlot(pbmc, reduction = "umap", group.by="celltype",label=FALSE, sizes.highlight=4,
cells.highlight=colnames(subset(x = pbmc, subset = celltype == "Tcell")), cols.highlight="black") + NoLegend()

ggplot(p$data)+geom_point(aes(UMAP_1, UMAP_2), color="black",size=3)+
geom_point(aes(UMAP_1,UMAP_2,color=celltype))+scale_color_brewer(palette="Set2")+labs(color="")+
theme(legend.position = c(1, 1),legend.justification = c("right", "top"))+
guides(color = guide_legend(override.aes = list(size = 5)))





DoHeatmap(pbmc, features = top10$gene, group.by="celltype") + NoLegend()


Idents(object = pbmc_test) <- "newcluster"



tempmetdat<-pbmc@meta.data[is.element(pbmc@meta.data$orig.ident, c("1_Hash1","2_Hash","3_Hash","1_CD3","3_CD3","3_CD4","3_CD8","3_CD14","3_CD19","3_CD16","3_Spike1%","3_Spike0.1%")),]



clipr::write_clip(table(pbmc@meta.data$orig.ident, pbmc@meta.data$celltype))




p<-ggplot(tempmetdat, aes(x=orig.ident, fill=celltype))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
scale_fill_brewer(palette="Set2")+
RotatedAxis()







p<-ggplot(pbmc@meta.data, aes(x=orig.ident, fill=celltype))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
scale_fill_brewer(palette="Set2")+
RotatedAxis()



p<-ggplot(pbmc@meta.data, aes(x=patientid, fill=detailtype))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
RotatedAxis()

p<-ggplot(pbmc@meta.data[pbmc$celltype=="Epithelial",], aes(x=patientid, fill=detailtype))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
scale_fill_brewer(palette="Set1")+
RotatedAxis()



p<-ggplot(pbmc@meta.data, aes(x=drugfn, fill=celltype))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
scale_fill_brewer(palette="Set2")+
RotatedAxis()+
facet_wrap(~patientid)







######################## cell type subclustering ########################

# pbmc_epi<-subset(pbmc, subset=celltype=="Epithelial")
# pbmc_epi<-subset(pbmc, subset=celltype=="Fibroblast")
pbmc_epi<-subset(pbmc, subset=celltype=="Tcell")

pbmc_epi <- FindVariableFeatures(object = pbmc_epi)

VariableFeaturePlot(pbmc_epi)

pbmc_epi <- SCTransform(pbmc_epi, vars.to.regress = "percent.mt", verbose = FALSE)

pbmc_epi <- RunPCA(object = pbmc_epi)


VizDimLoadings(pbmc_epi, dims = 1:8, reduction = "pca")
DimHeatmap(pbmc_epi, dims = 1:15, cells = 500, balanced = TRUE)


######################## Clustering & Visualization ########################
res=0.5
dim=c(1,3,4,5,6)

pbmc_epi <- FindNeighbors(pbmc_epi, dims = dim)
pbmc_epi <- FindClusters(pbmc_epi, resolution = res)


pbmc_epi <- RunUMAP(pbmc_epi, dims = dim)

p1<-DimPlot(pbmc_epi, reduction = "umap")
p2<-DimPlot(pbmc_epi, reduction = "umap", group.by="patientid")
p3<-DimPlot(pbmc_epi, reduction = "umap", group.by="drugfn")
CombinePlots(list(p1,p2,p3),ncol=3)



Idents(pbmc_epi)<-"detailtcell"
pbmc.markers <- FindAllMarkers(pbmc_epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc_epi, features = top10$gene) + NoLegend()+theme(plot.margin = unit(c(0,3,0,0), "cm"))




mapfrom<-c('0', '1', '2', '3', '4')
mapto<-c('CD8 Tcell', 'IL7R+ Tcell', 'Activated Tcell', 'NKcell', 'MEG3+ Tcell')
celltypeorder<-mapto

pbmc_epi$detailtcell<-factor(mapvalues(pbmc_epi$seurat_clusters, mapfrom, mapto), levels=celltypeorder)

p<-ggplot(pbmc_epi@meta.data, aes(x=patientid, fill=detailtcell))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
scale_fill_brewer(palette="Set1")+
RotatedAxis()


p<-ggplot(pbmc_epi@meta.data, aes(x=drugfn, fill=detailtcell))+
geom_bar(stat="count", position="fill", color="black")+
labs(x="",y="Cell type composition", fill="")+
scale_fill_brewer(palette="Set1")+
RotatedAxis()+
facet_wrap(~patientid)



DimPlot(pbmc_epi, reduction = "umap", group.by="detailtcell")+scale_color_brewer(palette="Set1")



genelist=c("GNLY", "CD3D", "CD3E","CD8A", "CD4", "IL2RA","IL7R", "FOXP3", "IFNG", "TOX")

CombinePlots(lapply(FeaturePlot(pbmc_epi, features=genelist, combine=FALSE), function(x){x+NoAxes()+NoLegend()}))










FeaturePlot(pbmc, features=c("COL3A1"))

+facet_wrap(~patientid)


CombinePlots(lapply(VlnPlot(pbmc, features=c("HP", "CCR2","CCR3","CCR4","CCR5","CCR6","CCR7"), group.by="orig.ident", combine=FALSE)
, function(x){x+NoAxes()+facet_wrap(~patientid)}))

VlnPlot(pbmc, features=c("HP", "CCR2","CCR3","CCR4","CCR5","CCR6","CCR7"), group.by="orig.ident", combine=FALSE)

VlnPlot(pbmc, features="HP", group.by="patientid")

genelist<-c("GNLY")


assaydat<-GetAssayData(object = pbmc, slot = "data")
plotdat<-t(assaydat[genelist,,drop=FALSE])
plotdat<-cbind(pbmc[[]],plotdat)
plotdat<-melt(plotdat, id.vars=colnames(pbmc[[]]))

ggplot(plotdat, aes(x=celltype,y=value, fill=orig.ident,color=orig.ident))+geom_violin()+
geom_jitter(size=1,position=position_jitterdodge())+RotatedAxis()+
scale_fill_viridis(discrete = TRUE, option = "D")+
scale_color_viridis(discrete = TRUE, option = "D")+
labs(x="",y="",fill="",color="")+
ggtitle(genelist)


ggplot(plotdat, aes(x=celltype,y=value, fill=drugfn))+geom_boxplot()+
facet_wrap(~patientid)+RotatedAxis()+
scale_fill_viridis(discrete = TRUE, option = "D")+
labs(x="",y="",fill="",color="")+
ggtitle(genelist)



ggplot(plotdat, aes(x=celltype,y=value, fill=drugfn))+geom_boxplot()+
facet_wrap(~patientid)+RotatedAxis()+
scale_fill_viridis(discrete = TRUE, option = "D")+
labs(x="",y="",fill="",color="")+
ggtitle(genelist)















######################## Cell specific drug DEG ########################
pbmcsub <- subset(pbmc, subset = orig.ident!="CD3"&celltype!="MKI67+")

Idents(pbmcsub)<-"celltype"
bardat<-data.frame()
sum.markers<-data.frame()
sum.avgexp<-data.frame()

for (subquery in setdiff(celltypeorder,"MKI67+")){
pbmc_sub <- subset(pbmcsub, idents = subquery)
Idents(pbmc_sub) <- "orig.ident"
celltype.avg <- log1p(AverageExpression(pbmc_sub, verbose = FALSE)$RNA)

for (v in c("Hash2","Hash3")){

temp.markers<-data.frame()
temp.markers <- FindMarkers(pbmc_sub, ident.1 = v, ident.2="Hash1",  min.pct=0.1, logfc.threshold=0.1)
temp.markers$gene<-rownames(temp.markers)
temp.markers$sample<-v
temp.markers$celltype<-subquery
sum.markers<-rbind(sum.markers, temp.markers)

temp.markers <- temp.markers %>% filter(p_val_adj<0.05)

bardat<-rbind(bardat, data.frame(
Positive=nrow(temp.markers[temp.markers$avg_logFC>0,]),
Negative=nrow(temp.markers[temp.markers$avg_logFC<0,]),
drug=v,
celltype=subquery))


topn <- temp.markers %>% top_n(n = 3, wt = avg_logFC)
temp.avg<-celltype.avg
temp.avg$gene <- rownames(temp.avg)
temp.avg$Label<-"non-significant"
temp.avg$topn<-FALSE
temp.avg$celltype<-subquery
if(nrow(temp.avg[temp.markers$gene,])!=0){temp.avg[temp.markers$gene,]$Label<-"Significant"}
if(nrow(temp.avg[topn$gene,])!=0){temp.avg[topn$gene,]$topn<-TRUE}
print(v)

tempmelt<-melt(temp.avg, id.vars=c("Hash1", "gene","Label","topn","celltype"))
tempmelt<-tempmelt %>% filter(variable==v)
sum.avgexp<-rbind(sum.avgexp, tempmelt)
}
}


# write.table(sum.markers, file=paste(RD, nameplot, "_unfiltered_deg.txt", sep=''), quote=FALSE, row.names=FALSE)
# sum.markers.filtered<-sum.markers %>% filter(p_val_adj<0.05)
# write.table(sum.markers.filtered, file=paste(RD, nameplot, "_filtered_deg.txt", sep=''), quote=FALSE, row.names=FALSE)

## xy-average expression plot
p<-ggplot(sum.avgexp, aes(Hash1, value, color=Label))+
geom_point(alpha=0.2)+scale_color_manual(values=c("black","red"))+
geom_text_repel(data=sum.avgexp[sum.avgexp$topn,], aes(label=gene))+
facet_grid(variable~celltype)+ NoLegend()+labs(x="",y="")+
theme(strip.background = element_blank())




save_plot(paste(RD, nameplot, "_Averageexpression.png",sep=''),plot=p,
base_width=10, base_height=5)


## Number of DEG barplot
p<-ggplot(melt(bardat), aes(x=celltype, y=value, fill=variable))+
geom_bar(stat="identity")+facet_wrap(~drug)+ scale_fill_brewer(palette="Set1")+
labs(x="", y="Differentially expressed gene (FDR < 0.05)", fill="")+
RotatedAxis()

save_plot(paste(RD, nameplot, "_DEGnum.png",sep=''),plot=p,
base_height=6)









pbmc[["newident"]]<-paste(pbmc@meta.data$celltype,pbmc@meta.data$orig.ident,sep="_")

Idents(pbmc) <- "newident"

pbmcsub <- subset(pbmc, subset = orig.ident!="CD3"&celltype!="MKI67+")
celltype.avg <- log1p(AverageExpression(pbmcsub, verbose = FALSE)$RNA)

cordata<-cor(celltype.avg,method="pearson")


pheatmap(cordata[order(rownames(cordata)),order(rownames(cordata))], show_colnames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE,
treeheight_row=0)



genelist=c("IL7R", "CCR7","S100A4", "CD14","MS4A1","CD8A","FCGR3A","GNLY","FCER1A","CD3E","CD3D", "CD19")
pheatmap(celltype.avg[genelist,order(colnames(celltype.avg))], cluster_rows=FALSE, cluster_cols=FALSE,
scale="row",treeheight_row=0)




cusorder<-heatmapdata$tree_col$labels[heatmapdata$tree_col$order]
cordata_order<-cordata[cusorder,]
cordata_order<-cordata_order[,cusorder]

usrcol=colorRampPalette(c("blue", "white", "red"))(1024)
usrcol=colorRampPalette(c("midnightblue","navy", "white", "firebrick3", "firebrick4"))(1024)


pheatmap(cordata_order, color=usrcol, show_colnames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE,
treeheight_row=0, annotation_colors=ann_colors, annotation_legend=FALSE)










######################## GO analysis ########################
library(topGO)
library(ALL)

sample="Carboplatin"
celltype="Epithelial"
posneg="Negative"

if (posneg=="Positive"){
geneList.inter<-sum.markers$gene[sum.markers$sample==sample&sum.markers$celltype==celltype&sum.markers$avg_logFC>0]
}else if(posneg=="Negative"){
geneList.inter<-sum.markers$gene[sum.markers$sample==sample&sum.markers$celltype==celltype&sum.markers$avg_logFC<0]
}


# geneList.inter<-pbmc.markers[pbmc.markers$cluster=="Myeloid",]$gene

geneList<-rownames(pbmc)
genetf<-is.element(geneList, geneList.inter)
genetf<-factor(as.integer(genetf))
names(genetf)<-geneList


GOdata <- new("topGOdata", ontology = "BP", allGenes = genetf,
annot = annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher") 

# allRes$Term <- factor(allRes$Term, levels=rev(allRes$Term))
ggplot(allRes, aes(reorder(Term,-log(as.numeric(classicFisher))) , -log(as.numeric(classicFisher), base=10)))+geom_bar(stat='identity')+coord_flip() +
 labs(y="- log p-value")+theme(axis.title.y=element_blank())+ggtitle(paste(sample,celltype,posneg,length(geneList.inter),sep="_"))


# goID="GO:0010941"
# goID="GO:0002252"
# geneList.inter[is.element(geneList.inter, genesInTerm(GOdata, goID)[[1]])]







######################## Accumulative curve ########################
for (Samplename in Samples)
{
tmpdata <-read.table(paste(WD, Samplename, "/summary/", Samplename, "_used_read_ratio.txt", sep=''),
 header=FALSE, stringsAsFactors=FALSE, comment.char = "U")
# tmpdata$sample<-Samplename
tmpdata$sample<-factor(mapvalues(Samplename, Samples, newsmpnames),levels=newsmpnames)

if (Samplename==Samples[1]){
plotdata<-tmpdata
}
else{
plotdata<-rbind(plotdata, tmpdata)
}
}

p <- ggplot(plotdata, aes(V1, V2, color=sample))+facet_wrap(~sample)
p + geom_line(size=1) + theme_bw() + xlab("Cell Number") + ylab("Accumulative fraction (%)") +
ylim(0,100) +
scale_x_continuous(minor_breaks=seq(0, 5000, by = 500), limits=c(0,5000))+
theme(axis.text.x = element_text(colour="grey20", size=12), axis.text.y = element_text(colour="grey20", size=12), text=element_text(size=16, family="Arial"))

######################## Samfile summary ########################
for (Samplename in Samples)
{
# tmpdata <-read.table(paste(WD, Samplename, "/summary/", Samplename, "_samfile_summary.txt", sep=''),
 # header=FALSE, stringsAsFactors=FALSE,  sep="\t")[1:5,]

tmpdata <-read.table(paste(WD, Samplename, "/summary/", Samplename, "_samfile_summary_total.txt", sep=''),
 header=FALSE, stringsAsFactors=FALSE,  sep="\t")[1:7,]
 
tmpdata$sample<-mapvalues(Samplename, Samples, newsmpnames)

if (Samplename==Samples[1]){
plotdata<-tmpdata
}
else{
plotdata<-rbind(plotdata, tmpdata)
}
}
plotdata<-melt(plotdata, id.vars=c("V1","sample"))

plotdata$V1<-factor(plotdata$V1, levels=rev(c("Coding reads", "UTR reads", "Intergenic reads", "Intronic reads", "Multimapped reads","Unmapped reads","Dimer reads")))

ggplot(plotdata, aes(x=sample, y=value, fill=V1))+geom_bar(stat="identity", position="stack")

ggplot(plotdata, aes(x=sample, y=value, fill=V1))+geom_bar(stat="identity", position="fill")


RColorBrewer::brewer.pal(name="Set1", n=3)
















######################### Functions ############################
readSubsmp<-function(WD, Samples, newsmpnames, smporder){
for (Samplename in Samples){
tmpdata <- read.table(paste(WD,Samplename,"/subsample/",Samplename,"_subsample_summary.txt", sep="")
, header=TRUE, stringsAsFactors=FALSE)
if (nrow(tmpdata)!=0){
tmpdata$sample<-factor(mapvalues(Samplename, Samples, newsmpnames, warn_missing=FALSE),levels=newsmpnames)
if (Samplename==Samples[1]){
mergedat<-tmpdata
}
else{
mergedat<-rbind(mergedat, tmpdata)
}
}
}
mergedat$sample<-factor(mergedat$sample, levels=smporder)
return(mergedat)
}


Downsmpplot_fitted<-function(subdownsmpdat){
p<-ggplot(subdownsmpdat, aes(Rawreads, UMI,color=sample))+geom_point()+
labs(x="Median raw reads per cell", y="Median UMIs",color="")+
xlim(0,2*10^5)+ylim(0,4000)+
geom_smooth(method = "nls",
method.args = list(formula = y~a+(b/(x+c)),start = c(a=2000, b=-2*10^7, c=10000)),
data=subdownsmpdat,
se=FALSE,
fullrange=TRUE,
linetype=2,
size=0.5
)+
theme(legend.position = c(.95, .95),
legend.justification = c("right", "top"),
plot.margin = unit(c(0,1,0,0), "cm"))+
scale_color_viridis(discrete = TRUE, option = "D")

return(p)
}



readSeurat<-function(WD, Samples, newsmpnames, pbmclist){

for (Samplename in Samples)
{
rawdata <- as(data.matrix(read.table(paste(WD, Samplename, "/dge/", Samplename, "_hg_dge.txt", sep='')
, header=TRUE
, sep="\t"
, stringsAsFactors=TRUE
, row.names=1 )), "sparseMatrix")

tempbmc <- CreateSeuratObject(counts = rawdata[rownames(rawdata)!="RNA28S5",], min.cells = 3, min.features = 0, project = mapvalues(Samplename, Samples, newsmpnames))
pbmclist<-append(pbmclist, tempbmc)

}

return(pbmclist)
}

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- grid::segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                   grobs = new.flag,
                                   t = 4, 
                                   l = 4
  )

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid.newpage()
  grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


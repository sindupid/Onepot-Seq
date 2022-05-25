library(Seurat)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(plyr)
library(dplyr)
library(reshape)
library(viridis)
theme_set(theme_cowplot(7))

library(wesanderson)
# wes_palette(n=4, name="Chevalier1"))










########################## Figure 1e #######################################
plotdata<-
metadata%>%
filter(SequencingDate=="210303")%>%
filter(is.element(sample, c('Snowman_1', 'Snowman_2', 'Snowman_8','Snowman_4','Snowman_5')))


plotdata$y=ifelse(plotdata$SequencingDate=="210217",
plotdata$mx500_n/(plotdata$mm500_n+plotdata$hs200_n+plotdata$mx500_n),
plotdata$mx500_n/(plotdata$mm500_n+plotdata$hs500_n+plotdata$mx500_n))
plotdata$cy=ifelse(plotdata$SequencingDate=="210217",
(plotdata$mm500_n+plotdata$hs200_n+plotdata$mx500_n)/(plotdata$HumanInput+plotdata$MouseInput),
(plotdata$mm500_n+plotdata$hs500_n+plotdata$mx500_n)/(plotdata$HumanInput+plotdata$MouseInput))
# plotdata$hcy=ifelse(plotdata$SequencingDate=="210217",
# plotdata$hs200_n/plotdata$HumanInput,
# plotdata$hs500_n/plotdata$HumanInput)
# plotdata$mcy=ifelse(plotdata$SequencingDate=="210217",
# plotdata$mm500_n/plotdata$MouseInput,
# plotdata$mm500_n/plotdata$MouseInput)


plotdata$LysisTime<-factor(plotdata$LysisTime, levels=c(100,300,900,1800,3600))

p<-ggplot(plotdata, aes(x=LysisTime, y=cy*100, fill=LysisTime))+
geom_bar(stat="identity")+
scale_fill_manual(values=wes_palette("Zissou1", n = 5))+
labs(x="Lysis time (sec)",y="Cell Yield (%)", fill="", color="")+
RotatedAxis()+NoLegend()

save_plot("figure1e_lysistime.png",plot=p,
base_height=1, base_width=1.2)


WD <- "/Data2/SDJ/210303_NextSeq/"
Samples <-c('Snowman_1', 'Snowman_2', 'Snowman_3','Snowman_4','Snowman_5','Snowman_6','Snowman_7','Snowman_8','Snowman_9','Snowman_10','Snowman_11')
newsmpnames<-Samples
smporder<-Samples
dir.create(paste(WD, "rdata", sep=''))

hUMIcut<-500
mUMIcut<-500
mergedat<-readHMplot(WD, Samples, newsmpnames, smporder, hUMIcut, mUMIcut)
mergedat<-mergedat%>%filter(Total<50000)

subsamples <-c('Snowman_1', 'Snowman_2', 'Snowman_8','Snowman_4','Snowman_5')
newsubsmpnames<-c('100','300','900','1800','3600')
ncol=length(subsamples)

subdat<-mergedat[is.element(mergedat$sample, subsamples),]
subdat$sample<-factor(mapvalues(subdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)

p<-ggplot(subdat[is.element(subdat$ident, c("Human","Mouse")),], aes(x=sample, y=spepury, color=ident))+
geom_boxplot(outlier.shape=NA,lwd=0.3)+
# geom_violin(position=position_dodge(width=0.8))+
# geom_jitter(size=0.1,alpha=0.1,position=position_jitterdodge())+
labs(x="Lysis time (sec)",y="Species purity",color="")+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(legend.position=c(0.05,0.4))+
ylim(0.89,1)+
RotatedAxis()

save_plot("figure1e_speciespuri.png",plot=p,
base_height=1, base_width=1.27)




########################## Figure 1f #######################################

plotdata<-
metadata%>%
filter(Well!="Drop-seq")%>%
filter(Well=="12well")%>%
filter(SequencingDate=="210217"|SequencingDate=="210303")%>%
filter(sample!="Bsnowman_1_cDNA"&sample!="Bsnowman_2_cDNA")%>%
filter(Staining=="None")%>%
filter(CentriRTLysis=="Normal")%>%
filter(Wash=="0")

plotdata$y=ifelse(plotdata$SequencingDate=="210217",
plotdata$mx500_n/(plotdata$mm500_n+plotdata$hs200_n+plotdata$mx500_n),
plotdata$mx500_n/(plotdata$mm500_n+plotdata$hs500_n+plotdata$mx500_n))

p<-ggplot(plotdata, aes(x=(HumanInput+MouseInput)/WellArea, y=y*100,
fill=as.character(SequencingDate), shape=Well))+
stat_function(fun=function(x) mult_wobind(x,90), aes(color="r=90"), linetype=2)+
stat_function(fun=function(x) mult_wobind(x,60), aes(color="r=60"), linetype=2)+
# stat_function(fun=function(x) mult_wobind(x,30), aes(color="r=30"), linetype=2)+
scale_color_manual(labels=c(expression("r"[e]*"=60µm"),expression("r"[e]*"=90µm")),
values = wes_palette("Moonrise1", n = 2))+
geom_point(size=2,shape=21)+
scale_fill_manual(labels=c("Direct lysis", "Diffusive lysis"), values=c("#2C5F2D","#97BC62FF"))+
guides(fill = guide_legend(override.aes=list(shape=21)))+
labs(x=expression("Cell density (cells/"*cm^2*")"),y="Multiplets (%)", fill="", color="")+
theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin(),
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5))


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure1_multipletcelldensity.png",plot=p,
base_height=1.9, base_width=1.6)

# save_plot("/Data2/SDJ/script/Snowmanseq/image/figure1_multipletcelldensity.png",plot=p,
# base_height=2.1, base_width=1.8)



plotdata$cy=ifelse(plotdata$SequencingDate=="210217",
(plotdata$mm500_n+plotdata$hs200_n+plotdata$mx500_n)/(plotdata$HumanInput+plotdata$MouseInput),
(plotdata$mm500_n+plotdata$hs500_n+plotdata$mx500_n)/(plotdata$HumanInput+plotdata$MouseInput))

p<-ggplot(plotdata, aes(x=(HumanInput+MouseInput)/WellArea, y=cy*100,
fill=as.character(SequencingDate), shape=Well))+
geom_hline(aes(yintercept=yield_wobind(5714.29,90),color="r=90"),linetype=2)+
geom_hline(aes(yintercept=yield_wobind(5714.29,60),color="r=60"),linetype=2)+
# geom_hline(aes(yintercept=yield_wobind(5714.29,30),color="r=30"),linetype=2)+
scale_color_manual(labels=c(expression("r"[e]*"=60µm"),expression("r"[e]*"=90µm")),
values = wes_palette("Moonrise1", n = 2))+
geom_point(size=2, shape=21)+
ylim(0,100)+
scale_fill_manual(labels=c("Direct lysis", "Diffusive lysis"), values=c("#2C5F2D","#97BC62FF"))+
guides(fill = guide_legend(override.aes=list(shape=21)))+
labs(x=expression("Cell density (cells/"*cm^2*")"),y="Cell yield (%)", fill="", color="")+
theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin(),
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5))

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure1_yieldcelldensity.png",plot=p,
base_height=1.9, base_width=1.6)

# save_plot("/Data2/SDJ/script/Snowmanseq/image/figure1_yieldcelldensity.png",plot=p,
# base_height=2.1, base_width=1.8)



######################## Figure 2c ##################################


facsdat<-read.table("/Data2/SDJ/script/Snowmanseq/facsdata.txt",stringsAsFactors=FALSE,header=TRUE)
facsdat$value=facsdat$Beadcell/(facsdat$Beadcell+facsdat$Onlybead)

facsdatstat<- facsdat %>% 
filter(ratio=="1_1")%>%
filter(time!="60min-96")%>%
group_by(time)%>%
summarise(mean= mean(value*100),sd=sd(value*100))

facsdatstat$time<-factor(facsdatstat$time, levels=c("10min","30min","60min","120min"))

p<-ggplot(facsdatstat,aes(x=time, y=mean))+
geom_errorbar(aes(x=time, ymin=mean-sd, ymax=mean+sd), width=0.2)+
geom_bar(stat="identity",fill="#353a41")+
geom_hline(yintercept=5.6, linetype=2)+
labs(x="Incubation time",y="Bead-cell doublets \n/ Total beads (%)")+RotatedAxis()


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure2_facsbinding.png",plot=p,
base_height=1.4, base_width=1.3)








######################## Figure 2d ##################################

plotdata<-
metadata%>%
filter(Well!="Drop-seq")%>%
filter(Well=="12well")%>%
filter(SequencingDate=="210217")%>%
filter(sample!="Bsnowman_1_cDNA"&sample!="Bsnowman_2_cDNA")%>%
filter(CentriRTLysis=="Normal")%>%
filter(HumanInput<=4000)%>%
filter(Wash=="0")

plotdata$y=ifelse(plotdata$SequencingDate=="210217",
plotdata$mx500_n/(plotdata$mm500_n+plotdata$hs200_n+plotdata$mx500_n),
plotdata$mx500_n/(plotdata$mm500_n+plotdata$hs500_n+plotdata$mx500_n))

plotdata$cy=ifelse(plotdata$SequencingDate=="210217",
(plotdata$mm500_n+plotdata$hs200_n+plotdata$mx500_n)/(plotdata$HumanInput+plotdata$MouseInput),
(plotdata$mm500_n+plotdata$hs500_n+plotdata$mx500_n)/(plotdata$HumanInput+plotdata$MouseInput))

plotdata$cd<-(plotdata$HumanInput+plotdata$MouseInput)/plotdata$WellArea

plotdatasub1<- plotdata %>% filter(Staining=="None")
plotdatasub2<- plotdata %>% filter(Staining=="Both")

plotdatasub<-merge(plotdatasub1, plotdatasub2, by="cd")



p<-ggplot(plotdata, aes(x=(HumanInput+MouseInput)/WellArea, y=cy*100))+
geom_point(size=2, shape=21,aes(fill=as.character(Staining)))+
ylim(0,100)+xlim(200,2300)+
scale_fill_manual(labels=c("Stained", "Not stained"), values=c("#00203FFF","#ADEFD1FF"))+
guides(fill = guide_legend(override.aes=list(shape=21)))+
geom_segment(data=plotdatasub, aes(x = cd,xend=cd,y=cy.x*100, yend = cy.y*100),
arrow=arrow(length=unit(0.12,"cm"),angle=30),color='orange', size=0.5)+
labs(x=expression("Cell density (cells/"*cm^2*")"),y="Cell yield (%)", fill="", color="")+
theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin(),
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5))

save_plot("/Data2/SDJ/script/Snowmanseq/image/figure2_abeffcellyield.png",plot=p,
base_height=1.5, base_width=1.5)


p<-ggplot(plotdata, aes(x=(HumanInput+MouseInput)/WellArea, y=y*100))+
geom_point(size=2, shape=21,aes(fill=as.character(Staining)))+
xlim(200,2300)+ylim(0,40)+
scale_fill_manual(labels=c("Stained", "Not stained"), values=c("#00203FFF","#ADEFD1FF"))+
guides(fill = guide_legend(override.aes=list(shape=21)))+
geom_segment(data=plotdatasub, aes(x = cd,xend=cd,y=y.x*100, yend = y.y*100),
arrow=arrow(length=unit(0.12,"cm"),angle=30),color='orange', size=0.5)+
labs(x=expression("Cell density (cells/"*cm^2*")"),y="Multiplets (%)", fill="", color="")+
theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin(),
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5))


save_plot("/Data2/SDJ/script/Snowmanseq/image/figure2_abeffmultiplet.png",plot=p,
base_height=1.5, base_width=1.45)









################### Function ################################



Downsmpplot<-function(subdownsmpdat){
p<-ggplot(subdownsmpdat, aes(Rawreads, UMI, color=sample))+geom_point()+geom_line()+
labs(x="Median raw reads per cell", y="Median UMIs",color="")+
xlim(0,10^5)+ylim(0,4000)
return(p)
}


Downsmpplot_fitted<-function(subdownsmpdat){
p<-ggplot(subdownsmpdat, aes(Rawreads, UMI,color=sample))+geom_point()+
labs(x="Median raw reads per cell", y="Median UMIs",color="")+
xlim(0,10^5)+ylim(0,4000)+
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
scale_color_brewer(palette=3, type="qual")

return(p)
}

hmDownsmpplot_fitted<-function(subdownsmpdat){
p1<-ggplot(subdownsmpdat, aes(hRawreads, hUMI,color=sample))+geom_point()+
labs(x="Median raw reads per cell", y="Median UMIs",color="")+
xlim(0,10^5)+ylim(0,6000)+
geom_smooth(method = "nls",
method.args = list(formula = y~a+(b/(x+c)),start = c(a=2000, b=-2*10^7, c=10000)),
data=subdownsmpdat,
se=FALSE,
fullrange=TRUE,
linetype=2,
size=0.5
)+
theme(legend.position = "none",
plot.margin = unit(c(0,1,0,0), "cm"))+
scale_color_viridis(discrete = TRUE, option = "D")+
ggtitle("Human cells")

p2<-ggplot(subdownsmpdat, aes(mRawreads, mUMI,color=sample))+geom_point()+
labs(x="Median raw reads per cell", y="Median UMIs",color="")+
xlim(0,10^5)+ylim(0,6000)+
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
scale_color_viridis(discrete = TRUE, option = "D")+
ggtitle("Mouse cells")


p<-CombinePlots(list(p1,p2))
return(p)
}





HMplot<-function(subdat){
p<-ggplot(subdat, aes(Human_transcripts, Mouse_transcripts, color=ident))+geom_point(alpha = 0.3)+
labs(x="Human transcripts", y="Mouse transcripts",color="")+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(strip.background = element_blank(), axis.text.x=element_text(angle=30, hjust=1), legend.position="none")
return(p)
}

HMpurityplot<-function(subdat){
p<-ggplot(subdat[is.element(subdat$ident, c("Human","Mouse")),], aes(x=ident, y=spepury, color=ident))+
geom_boxplot(outlier.shape=NA)+geom_jitter(size=0.3)+
labs(x="",y="Species purity",color="")+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(strip.background = element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(), legend.position="none")
return(p)
}
HMtranscriptplot<-function(subdat){
p<-ggplot(subdat[is.element(subdat$ident, c("Human","Mouse")),], aes(x=ident, y=Total, color=ident))+
geom_violin()+geom_boxplot(width=0.5,outlier.shape=NA)+
# geom_jitter(size=0.7)+
labs(x="",y="Number of transcripts",color="")+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(strip.background = element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(),legend.position="none")
return(p)
}

HMnumberplot<-function(substatdat){
p<-ggplot(substatdat, aes(x="",y=value,fill=Var.2))+
geom_bar(stat="identity", position=position_stack(reverse=TRUE))+
scale_fill_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
labs(x="",y="Number of cells",fill="")+
theme(strip.background = element_blank(), axis.ticks.x=element_blank(),legend.position="none")
return(p)
}


HMaccplot<-function(subplotdat){
p <- ggplot(subplotdat, aes(V1, V2, color=sample))+
geom_line(size=1) + xlab("Cell Number") + ylab("Accumulative fraction (%)") +
ylim(0,100) + xlim(0,5000)+
# theme_bw()+
# scale_x_continuous(minor_breaks=seq(0, 10000, by = 500))+
theme(strip.background = element_blank(),axis.text.x = element_text(angle=30, hjust=1), legend.position="none",
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5),
axis.title.y=element_text(size=12),
axis.line=element_blank()
)
return(p)
}




HMaccplot_ident<-function(subdat){
p<-ggplot(subdat, aes(bcdorder, cumsumrat, color=ident))+
geom_point(size=0.5, alpha=0.5) + labs(x="Ordered cell barcodes",y="Accumulative fraction (%)", color="") +
ylim(0,1) + xlim(0,10000)+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(strip.background = element_blank(),axis.text.x = element_text(angle=30, hjust=1),
legend.position="none",
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5),
axis.title.y=element_text(size=12),
axis.line=element_blank()
)
return(p)
}


HMaccplot_spepury<-function(subdat){
p<-ggplot(subdat, aes(bcdorder, cumsumrat, color=spepury))+
geom_point(size=0.5,alpha=0.5) + labs(x="Ordered cell barcodes", y="Accumulative fraction (%)", color="") +
ylim(0,1) + xlim(0,10000)+
theme(strip.background = element_blank(),axis.text.x = element_text(angle=30, hjust=1),
legend.position="none",
panel.grid.major = element_line(colour = "grey92"),
panel.grid.minor = element_line(colour = "grey92",size = rel(0.5)),
panel.border = element_rect(fill = NA,colour = "black", linetype=1,size=0.5),
axis.title.y=element_text(size=12),
axis.line=element_blank()
)
return(p)
}



mini_HMpurityplot<-function(subdat){
p<-ggplot(subdat[is.element(subdat$ident, c("Human","Mouse")),], aes(x=sample, y=spepury, color=ident))+
geom_boxplot(outlier.shape=NA)+
labs(x="",y="Species purity",color="")+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(strip.background = element_blank(), legend.position="none")+
coord_flip()

return(p)
}

mini_HMnumberplot<-function(substatdat){
p<-ggplot(substatdat, aes(x=sample,y=value,fill=Var.2))+
geom_bar(stat="identity", position=position_stack(reverse=TRUE))+
scale_fill_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
labs(x="",y="Number of cells",fill="")+
theme(strip.background = element_blank(), legend.position="none")+
coord_flip()
return(p)
}

mini_HMtranscriptplot<-function(subdat){
p<-ggplot(subdat[is.element(subdat$ident, c("Human","Mouse")),], aes(x=sample, y=Total, color=ident))+
geom_violin(position=position_dodge(width=0.8))+geom_boxplot(width=0.5,outlier.shape=NA,position=position_dodge(width=0.8))+
# geom_jitter(size=0.7)+
labs(x="",y="Number of transcripts",color="")+
scale_color_manual(values=c("#446455","#FDD262","#C7B19C","#D3DDDC"))+
theme(strip.background = element_blank(),legend.position="none")+
coord_flip()
return(p)
}


HMreadUMIplot<-function(subdat){
a_mean <- subdat[subdat$ident!="Emptybeads",] %>% group_by(sample) %>% dplyr::summarize(slope_avg = sum(Total)/sum(rawreads), intercept=0)
p<-ggplot(subdat[subdat$ident!="Emptybeads",], aes(rawreads, Total, color=sample))+geom_point(alpha=0.3)+
theme(strip.background = element_blank(), axis.text.x=element_text(angle=30, hjust=1), legend.position="none")+
labs(x="Raw reads per cell", y="UMIs per cell",color="")+
geom_abline(data= a_mean, linetype="dashed",aes(slope=slope_avg,color=sample, intercept=intercept))+
geom_text(data=a_mean, mapping = aes(x = -Inf, y = Inf, label = round(slope_avg, digits=4)),hjust=-0.5, vjust=3)+
facet_wrap(~sample,ncol=ncol)
return(p)
}

readSubsmp<-function(WD, Samples, newsmpnames, smporder){
for (Samplename in Samples){
print(Samplename)
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

readhmSubsmp<-function(WD, Samples, newsmpnames, smporder){
for (Samplename in Samples){
tmpdata <- read.table(paste(WD,Samplename,"/subsample/",Samplename,"_subsample_summary_hm.txt", sep="")
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




readHMplot<-function(WD, Samples, newsmpnames, smporder, hUMIcut, mUMIcut){
	for (Samplename in Samples){

	tmpdata <- read.table(paste(WD,Samplename,"/dge/",Samplename,"_mix2_dge.summary.txt", sep="")
	, header=TRUE, stringsAsFactors=FALSE)
	tmpdata[,"HtoM"] <-tmpdata$Human_transcripts/(tmpdata$Mouse_transcripts+tmpdata$Human_transcripts)
	tmpdata[,"Total"] <-tmpdata$Mouse_transcripts+tmpdata$Human_transcripts
	tmpdata$spepury<-ifelse(tmpdata$HtoM>0.5, tmpdata$HtoM, 1-tmpdata$HtoM)
	# tmpdata<-head(tmpdata[order(tmpdata$Total, decreasing=TRUE),],n=200)

	tmpdata$sample<-factor(mapvalues(Samplename, Samples, newsmpnames, warn_missing=FALSE),levels=newsmpnames)
	tmpdata2 <- read.table(paste(WD,Samplename,"/dge/",Samplename,"_cell_readcounts.txt", sep="")
	, header=FALSE, stringsAsFactors=FALSE)
	colnames(tmpdata2)<-c("rawreads","Barcode")
	tmpdata<-merge(tmpdata, tmpdata2, by="Barcode",all.x=TRUE)
    
    tmpdata<-tmpdata[order(tmpdata$rawreads, decreasing=TRUE),]
    tmpdata$cumsumrat<-cumsum(tmpdata$rawreads)/sum(tmpdata2$rawreads)
    tmpdata$bcdorder<-1:nrow(tmpdata)
       
	if (Samplename==Samples[1]){
	mergedat<-tmpdata
	}
	else{
	mergedat<-rbind(mergedat, tmpdata)
	}
	}

	mergedat$sample<-factor(mergedat$sample, levels=smporder)

	mergedat$ident<-"Emptybeads"
	mergedat[mergedat$Total>max(hUMIcut,mUMIcut),]$ident<-"Mixed"
	mergedat[mergedat$Total>hUMIcut&mergedat$HtoM>=0.9,]$ident<-"Human"
	mergedat[mergedat$Total>mUMIcut&mergedat$HtoM<0.1,]$ident<-"Mouse"
	mergedat$ident<-factor(mergedat$ident, levels=c("Human","Mouse", "Mixed","Emptybeads"))
	return(mergedat)
}

readACCplot<-function(WD, Samples, newsmpnames, smporder){
for (Samplename in Samples)
{
tmpdata <-read.table(paste(WD, Samplename, "/summary/", Samplename, "_used_read_ratio.txt", sep=''),
 header=FALSE, stringsAsFactors=FALSE, comment.char = "U")
# tmpdata$sample<-Samplename
tmpdata$sample<-mapvalues(Samplename, Samples, newsmpnames, warn_missing=FALSE)
if (Samplename==Samples[1]){
plotdata<-tmpdata
}
else{
plotdata<-rbind(plotdata, tmpdata)
}
}
plotdata$sample<-factor(plotdata$sample, levels=smporder)
return(plotdata)
}





save_SMplots<-function(WD,subsamples,newsubsmpnames,ncol,nameplot,mergedat,statdat,plotdata,downsmpdat){

subdat<-mergedat[is.element(mergedat$sample, subsamples),]
subdat$sample<-factor(mapvalues(subdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)
substatdat<-statdat[is.element(statdat$Var.1, subsamples),]
substatdat$sample<-factor(mapvalues(substatdat$Var.1, subsamples, newsubsmpnames),levels=newsubsmpnames)
substatdat$Var.2<-factor(substatdat$Var.2, levels=c("Human","Mouse","Mixed"))
subplotdat<-plotdata[is.element(plotdata$sample, subsamples),]
subplotdat$sample<-factor(mapvalues(subplotdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)
subdownsmpdat<-downsmpdat[is.element(downsmpdat$sample, subsamples),]
subdownsmpdat$sample<-factor(mapvalues(subdownsmpdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)


#### Human-Mouse plot ###
p<-HMplot(subdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_HMplot.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)

### Species purity ###
# p<-HMpurityplot(subdat)+facet_wrap(~sample,ncol=ncol)
# save_plot(paste(WD, "rdata/", nameplot, "_SPplot.png",sep=''),plot=p,ncol=ncol,
# base_height=2.5, base_aspect_ratio=0.6)
p<-mini_HMpurityplot(subdat)
save_plot(paste(WD, "rdata/", nameplot, "_SPplot.png",sep=''),plot=p,
base_width=4, base_height=0.67+0.33*ncol)


### Number of transcripts ###
# p<-HMtranscriptplot(subdat)+facet_wrap(~sample,ncol=ncol)
# save_plot(paste(WD, "rdata/", nameplot, "_NTplot.png",sep=''),plot=p,ncol=ncol,
# base_height=2.5, base_aspect_ratio=0.6)
p<-mini_HMtranscriptplot(subdat)
save_plot(paste(WD, "rdata/", nameplot, "_NTplot.png",sep=''),plot=p,
base_width=4, base_height=0.67+0.33*ncol)

#### Number of cells ###
# p<-HMnumberplot(substatdat)+facet_wrap(~sample,ncol=ncol)
# save_plot(paste(WD, "rdata/", nameplot, "_NCplot.png",sep=''),plot=p,ncol=ncol,
# base_height=2.5, base_aspect_ratio=0.6)
p<-mini_HMnumberplot(substatdat)
save_plot(paste(WD, "rdata/", nameplot, "_NCplot.png",sep=''),plot=p,
base_width=4, base_height=0.67+0.33*ncol)


### Accumulative curve ###
p<-HMaccplot(subplotdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_ACplot.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)


p<-HMaccplot_ident(subdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_ACplot_ident.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)

p<-HMaccplot_spepury(subdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_ACplot_spepury.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)



### rawreads vs UMI plot ####
p<-HMreadUMIplot(subdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_RUplot.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)


p<-Downsmpplot_fitted(subdownsmpdat)
save_plot(paste(WD, "rdata/", nameplot, "_downsmpplot.png",sep=''),plot=p)
}



save_SMplots_wods<-function(WD,subsamples,newsubsmpnames,ncol,nameplot,mergedat,statdat,plotdata){

subdat<-mergedat[is.element(mergedat$sample, subsamples),]
subdat$sample<-factor(mapvalues(subdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)
substatdat<-statdat[is.element(statdat$Var.1, subsamples),]
substatdat$sample<-factor(mapvalues(substatdat$Var.1, subsamples, newsubsmpnames),levels=newsubsmpnames)
substatdat$Var.2<-factor(substatdat$Var.2, levels=c("Human","Mouse","Mixed"))
subplotdat<-plotdata[is.element(plotdata$sample, subsamples),]
subplotdat$sample<-factor(mapvalues(subplotdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)


#### Human-Mouse plot ###
p<-HMplot(subdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_HMplot.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)

### Species purity ###
# p<-HMpurityplot(subdat)+facet_wrap(~sample,ncol=ncol)
# save_plot(paste(WD, "rdata/", nameplot, "_SPplot.png",sep=''),plot=p,ncol=ncol,
# base_height=2.5, base_aspect_ratio=0.6)
p<-mini_HMpurityplot(subdat)
save_plot(paste(WD, "rdata/", nameplot, "_SPplot.png",sep=''),plot=p,
base_width=4, base_height=0.67+0.33*ncol)


### Number of transcripts ###
# p<-HMtranscriptplot(subdat)+facet_wrap(~sample,ncol=ncol)
# save_plot(paste(WD, "rdata/", nameplot, "_NTplot.png",sep=''),plot=p,ncol=ncol,
# base_height=2.5, base_aspect_ratio=0.6)
p<-mini_HMtranscriptplot(subdat)
save_plot(paste(WD, "rdata/", nameplot, "_NTplot.png",sep=''),plot=p,
base_width=4, base_height=0.67+0.33*ncol)

#### Number of cells ###
# p<-HMnumberplot(substatdat)+facet_wrap(~sample,ncol=ncol)
# save_plot(paste(WD, "rdata/", nameplot, "_NCplot.png",sep=''),plot=p,ncol=ncol,
# base_height=2.5, base_aspect_ratio=0.6)
p<-mini_HMnumberplot(substatdat)
save_plot(paste(WD, "rdata/", nameplot, "_NCplot.png",sep=''),plot=p,
base_width=4, base_height=0.67+0.33*ncol)


### Accumulative curve ###
p<-HMaccplot(subplotdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_ACplot.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)


### rawreads vs UMI plot ####
p<-HMreadUMIplot(subdat)+facet_wrap(~sample,ncol=ncol)
save_plot(paste(WD, "rdata/", nameplot, "_RUplot.png",sep=''),plot=p,
base_height=2.5, base_width=0.95+1.51*ncol)

}







save_customplots<-function(WD,subsamples,newsubsmpnames,ncol,nameplot,hmdownsmpdat){
subhmdownsmpdat<-hmdownsmpdat[is.element(hmdownsmpdat$sample, subsamples),]
subhmdownsmpdat$sample<-factor(mapvalues(subhmdownsmpdat$sample, subsamples, newsubsmpnames),levels=newsubsmpnames)
p<-hmDownsmpplot_fitted(subhmdownsmpdat)
save_plot(paste(WD, "rdata/", nameplot, "_hmdownsmpplot.png",sep=''),plot=p, ncol=2, base_height=3, base_aspect_ratio=1.2)
}






writeHMbarcodes<-function(WD, Samples){
for (Samplename in Samples){
for (species in c("Human","Mouse")){
write.table(mergedat[mergedat$sample==Samplename&mergedat$ident==species,]$Barcode,
file=paste(WD,Samplename,"/dge/",species,"_barcode.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
}
}
}


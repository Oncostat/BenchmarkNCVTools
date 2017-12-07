

ARGS<-commandArgs(trailingOnly=TRUE)
minRec<-max(2,as.integer(ARGS[1]))
maxRec<-as.integer(ARGS[2])
ARGS<-ARGS[-c(1:2)]

x<-c("Clinvar","COSMIC","HGMD")[sapply(c("Clinvar","COSMIC","HGMD"),function(x) length(grep(x,ARGS))==0)]
if(length(x)!=0){stop(paste0(paste(x,collapse=" and ")," file",ifelse(length(x)==1,"was","s were")," not found."))}


stroke_width<-1.5
colors<-c(1,1,2,"darkgreen",4,4,4,5,rep("purple",4),rep("yellow",4))
typeline<-c(1,2,1,1,1,2,3,1,1,2,3,4,1,2,3,4)
scores<-c("CADD","DANN","FATHMM_MKL","Funseq2","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","snp","somcll","somliver","somlung","sommel")




clin<-read.delim(ARGS[grep("Clinvar",ARGS)],h=T,stringsAsFactors=F)
cosmic<-read.delim(ARGS[grep("COSMIC",ARGS)],h=T,stringsAsFactors=F)
hgmd<-read.delim(ARGS[grep("HGMD",ARGS)],h=T,stringsAsFactors=F)
colnames(clin)[ncol(clin)]<-colnames(cosmic)[ncol(cosmic)]<-colnames(hgmd)[ncol(hgmd)]<-"No_GENCODE_feature"
clin<-subset(clin,select=-c(gene))
cosmic<-subset(cosmic,select=-c(gene))
hgmd<-subset(hgmd,select=-c(gene))

clin$intron<-clin$transcript-clin$exon
hgmd$intron<-hgmd$transcript-hgmd$exon
cosmic$intron<-cosmic$transcript-cosmic$exon
clin<-clin[,c(colnames(clin)[c(1:2)],sort(setdiff(colnames(clin)[3:ncol(clin)],"No_GENCODE_feature")),"No_GENCODE_feature")]
hgmd<-hgmd[,c(colnames(hgmd)[c(1:2)],sort(setdiff(colnames(hgmd)[3:ncol(hgmd)],"No_GENCODE_feature")),"No_GENCODE_feature")]
cosmic<-cosmic[,c(colnames(cosmic)[c(1:5)],sort(setdiff(colnames(cosmic)[6:ncol(cosmic)],"No_GENCODE_feature")),"No_GENCODE_feature")]


#############
#Make the table and the plot of frequencies for the paper
freq<-data.frame(Features=colnames(clin)[-c(1:2)])
freq$Clinvar<-apply(clin[,-c(1:2)],2,function(x){100*sum(x>0)/length(x)})
freq$HGMD<-apply(hgmd[,-c(1:2)],2,function(x){100*sum(x>0)/length(x)})
freqConf<-freqUnkn<-freq
countRecConf<-countRecUnkn<-countRecTotal<-NULL
for(i in minRec:maxRec){
	nbRec<-sum(cosmic$recTotal>=i)
	countRecTotal<-c(countRecTotal,nbRec)
	if(nbRec>0){freq<-within(freq,{assign(paste0("Rec",i),apply(cosmic[cosmic$recTotal>=i,-c(1:5)],2,function(x){100*sum(x>0)/length(x)}))})}
}
write.table(freq,"Results/featDistTotal.txt",row.names=F,quote = F,sep='\t')

freq$Features<-gsub("[[:punct:]]"," ",freq$Features)
freq$Features[freq$Features=="No_GENCODE_feature"]<-"no_GENCODE_feature"
freq$Features[freq$Features=="Selenocysteine"]<-"selenocysteine"

legInset<-c(-0.55,-0.7,-0.3,-0.55)
postscript("FigS1_Feature_distribution.eps",width=1120,height=560)
par(mfrow=c(2,2),mar=c(5.1, 4.1, 4.1, 12.1))
stockRows<-NULL
rec<-minRec:(length(grep("Rec",colnames(freq)))+1)
for(i in c("rna","pseudo","gene","other")){
	titl<-ifelse(i=="rna","RNAs",ifelse(i=="pseudo","Pseudogenes",ifelse(i=="gene","Genes","Other features")))
	if(i=="other"){
		rows<-c(1:nrow(freq))[!c(1:nrow(freq))%in%stockRows]
	}else{
		if(i=="gene"){
			rows<-grep(i,tolower(freq$Features))
			rows<-rows[!rows%in%grep("pseudo",tolower(freq$Features))]
			stockRows<-c(stockRows,rows)
		}else{
			rows<-grep(i,tolower(freq$Features))
			stockRows<-c(stockRows,rows)
		}
	}
	data<-freq[rows,]
	rows2rm<-NULL
	for(j in 1:nrow(data)){if(sum(data[j,-1])==0){rows2rm<-c(rows2rm,j)}}
	if(length(rows2rm)!=0){data<-data[-rows2rm,]}
	plot(0,0,type='l',xlim=c(minRec,max(rec)),ylim=c(0,max(data[,-1])),xlab="Recurrence",ylab="Frequency (%)",lwd=3,col="white",main=titl)  #/max(data[,1])     
	abline(h=data[,2],col=rep(1:8,2),lty=sort(rep(1:2,8)))
	for(j in 1:nrow(data)){lines(rec,data[j,4:ncol(data)],col=j*I(j<8)+j*I(j>=8),lty=1*(1+I(j>8)),lwd=3)}
	par(xpd=TRUE)
	legend("topright",inset=c(legInset[grep(i,c("rna","pseudo","gene","other"))],0),legend=data$Features,col=rep(1:8,3),bty="n",lty=1*(1+I(1:nrow(data)>8)),cex=0.75)
	par(xpd=FALSE)
}
dev.off()



#############
#out the most rec features distribution  (>10)
mostRecTotal<-(maxRec:minRec)[which(rev(countRecTotal)>=10)[1]]

data<-cosmic[cosmic[,"recTotal"]>=get("mostRecTotal"),]
out<-data[,c("chr","pos","recTotal")]
out$f<-unlist(lapply(apply(data,1,function(x) colnames(data)[-c(1:5)][x[-c(1:5)]>0]),function(x) paste(x,collapse=", ")))
colnames(out)<-c("Chromosome","Position","Recurrence","Features")
out<-out[order(out$Chromosome,out$Position),]
out<-out[order(out$Recurrence,decreasing = T),]
write.table(out,paste0("Results/mostRecTotal.txt"),row.names=F,quote = F,sep='\t')


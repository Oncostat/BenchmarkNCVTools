

ARGS<-commandArgs(trailingOnly=TRUE)
minRec<-max(2,as.integer(ARGS[1]))
maxRec<-as.integer(ARGS[2])
ARGS<-ARGS[-c(1:2)]

x<-c("Clinvar","COSMIC","HGMD")[sapply(c("Clinvar","COSMIC","HGMD"),function(x) length(grep(x,ARGS))==0)]
if(length(x)!=0){stop(paste0(paste(x,collapse=" and ")," files ",ifelse(length(x)==1,"was","were")," not found."))}

#ARGS<-paste0("Databases/annot",c("Clinvar","COSMIC","HGMD"),".txt")

clin<-read.delim(ARGS[grep("Clinvar",ARGS)],h=T,stringsAsFactors=F)
cosmic<-read.delim(ARGS[grep("COSMIC",ARGS)],h=T,stringsAsFactors=F)
hgmd<-read.delim(ARGS[grep("HGMD",ARGS)],h=T,stringsAsFactors=F)
colnames(clin)[ncol(clin)]<-colnames(cosmic)[ncol(cosmic)]<-colnames(hgmd)[ncol(hgmd)]<-"No_GENCODE_feature"
clin<-subset(clin,select=-c(gene))
cosmic<-subset(cosmic,select=-c(gene))
hgmd<-subset(hgmd,select=-c(gene))

#############
#make the table and the plot of frequencies for the paper
freq<-data.frame(Features=colnames(clin)[-c(1:2)])
freq$Clinvar<-apply(clin[,-c(1:2)],2,function(x){100*sum(x>0)/length(x)})
freq$HGMD<-apply(hgmd[,-c(1:2)],2,function(x){100*sum(x>0)/length(x)})
freqConf<-freqUnkn<-freq
countRecConf<-countRecUnkn<-countRecTotal<-NULL
for(i in minRec:maxRec){
	nbRec<-sum(cosmic$recConf>=i)
	countRecConf<-c(countRecConf,nbRec)
	if(nbRec>0){freqConf<-within(freqConf,{assign(paste0("Rec",i),apply(cosmic[cosmic$recConf>=i,-c(1:5)],2,function(x){100*sum(x>0)/length(x)}))})}
	nbRec<-sum(cosmic$recUnkn>=i)
	countRecUnkn<-c(countRecUnkn,nbRec)
	if(nbRec>0){freqUnkn<-within(freqUnkn,{assign(paste0("Rec",i),apply(cosmic[cosmic$recUnkn>=i,-c(1:5)],2,function(x){100*sum(x>0)/length(x)}))})}
	nbRec<-sum(cosmic$recTotal>=i)
	countRecTotal<-c(countRecTotal,nbRec)
	if(nbRec>0){freq<-within(freq,{assign(paste0("Rec",i),apply(cosmic[cosmic$recTotal>=i,-c(1:5)],2,function(x){100*sum(x>0)/length(x)}))})}
}
write.table(freqConf,"Results/featDistConf.txt",row.names=F,quote = F,sep='\t')
write.table(freqUnkn,"Results/featDistUnkn.txt",row.names=F,quote = F,sep='\t')
write.table(freq,"Results/featDistTotal.txt",row.names=F,quote = F,sep='\t')


jpeg(paste0("Results/FeaturesDistConf.jpeg"),width=7,height=6,units="in",res=350)
par(mfrow=c(2,2))
stockRows<-NULL
rec<-minRec:(length(grep("Rec",colnames(freqConf)))+1)
for(i in c("rna","pseudo","gene","other")){
	titl<-ifelse(i=="rna","RNAs",ifelse(i=="pseudo","Pseudogenes",ifelse(i=="gene","Genes","Other features")))
	if(i=="other"){
		rows<-c(1:nrow(freqConf))[!c(1:nrow(freqConf))%in%stockRows]
	}else{
		if(i=="gene"){
			rows<-grep(i,tolower(freqConf$Features))
			rows<-rows[!rows%in%grep("pseudo",tolower(freqConf$Features))]
			stockRows<-c(stockRows,rows)
		}else{
			rows<-grep(i,tolower(freqConf$Features))
			stockRows<-c(stockRows,rows)
		}
	}
	data<-freqConf[rows,]
	rows2rm<-NULL
	for(j in 1:nrow(data)){if(sum(data[j,-1])==0){rows2rm<-c(rows2rm,j)}}
	if(length(rows2rm)!=0){data<-data[-rows2rm,]}
	plot(0,0,type='l',xlim=c(minRec,max(rec)),ylim=c(0,max(data[,-1])),xlab="Recurrence",ylab="Frequency (%)",lwd=3,col="white",main=titl)  #/max(data[,1])     
	abline(h=data[,2],col=rep(1:8,2))
	for(j in 1:nrow(data)){lines(rec,data[j,4:ncol(data)],col=j*I(j<8)+j*I(j>=8),lty=1*(1+I(j>8)),lwd=3)}
	legend("topright",legend=data$Features,col=rep(1:8,3),bty="n",lty=1*(1+I(1:nrow(data)>8)))
}
dev.off()


jpeg(paste0("Results/FeaturesDistUnkn.jpeg"),width=7,height=6,units="in",res=350)
par(mfrow=c(2,2))
stockRows<-NULL
rec<-minRec:(length(grep("Rec",colnames(freqUnkn)))+1)
for(i in c("rna","pseudo","gene","other")){
	titl<-ifelse(i=="rna","RNAs",ifelse(i=="pseudo","Pseudogenes",ifelse(i=="gene","Genes","Other features")))
	if(i=="other"){
		rows<-c(1:nrow(freqUnkn))[!c(1:nrow(freqUnkn))%in%stockRows]
	}else{
		if(i=="gene"){
			rows<-grep(i,tolower(freqUnkn$Features))
			rows<-rows[!rows%in%grep("pseudo",tolower(freqUnkn$Features))]
			stockRows<-c(stockRows,rows)
		}else{
			rows<-grep(i,tolower(freqUnkn$Features))
			stockRows<-c(stockRows,rows)
		}
	}
	data<-freqUnkn[rows,]
	rows2rm<-NULL
	for(j in 1:nrow(data)){if(sum(data[j,-1])==0){rows2rm<-c(rows2rm,j)}}
	if(length(rows2rm)!=0){data<-data[-rows2rm,]}
	plot(0,0,type='l',xlim=c(minRec,max(rec)),ylim=c(0,max(data[,-1])),xlab="Recurrence",ylab="Frequency (%)",lwd=3,col="white",main=titl)  #/max(data[,1])     
	abline(h=data[,2],col=rep(1:8,2))
	for(j in 1:nrow(data)){lines(rec,data[j,4:ncol(data)],col=j*I(j<8)+j*I(j>=8),lty=1*(1+I(j>8)),lwd=3)}
	legend("topright",legend=data$Features,col=rep(1:8,3),bty="n",lty=1*(1+I(1:nrow(data)>8)))
}
dev.off()

#jpeg(paste0("Results/FeaturesDist.jpeg"),width=7,height=6,units="in",res=350)
freq$Features<-gsub("[[:punct:]]"," ",freq$Features)
freq$Features[freq$Features=="No_GENCODE_feature"]<-"no_GENCODE_feature"
freq$Features[freq$Features=="Selenocysteine"]<-"selenocysteine"
postscript("Results/FeaturesDist.eps",width=7,height=6)
postscript("../../Article/FigS1_Feature_distribution.eps",width=1120,height=960)
par(mfrow=c(2,2))
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
	abline(h=data[,2],col=rep(1:8,2))
	for(j in 1:nrow(data)){lines(rec,data[j,4:ncol(data)],col=j*I(j<8)+j*I(j>=8),lty=1*(1+I(j>8)),lwd=3)}
	if(i=="other"){
		coords1<-1:(floor(length(data$Features)/2))
		coords2<-(floor(length(data$Features)/2+1)):length(data$Features)
		legend("topleft", legend=data$Features[coords1],col=coords1,lwd=stroke_width,bty="n",lty=1,cex=0.75)
		legend("topright", legend=data$Features[coords2],col=coords2,lwd=stroke_width,bty="n",lty=2,cex=0.75)
	}else{
		legend("topleft",legend=data$Features,col=rep(1:8,3),bty="n",lty=1*(1+I(1:nrow(data)>8)),cex=0.75)
	}
}
dev.off()



#############
#out the most rec features distribution  (>10?)
mostRecConf<-(maxRec:minRec)[which(rev(countRecConf)>=10)[1]]
mostRecUnkn<-(maxRec:minRec)[which(rev(countRecUnkn)>=10)[1]]
mostRecTotal<-(maxRec:minRec)[which(rev(countRecTotal)>=10)[1]]

for(rec in c("Conf","Unkn","Total")){
	data<-cosmic[cosmic[,paste0("rec",rec)]>=get(paste0("mostRec",rec)),]
	out<-data[,c("chr","pos",paste0("rec",rec))]
	out$f<-unlist(lapply(apply(data,1,function(x) colnames(data)[-c(1:5)][x[-c(1:5)]>0]),function(x) paste(x,collapse=", ")))
	colnames(out)<-c("Chromosome","Position","Recurrence","Features")
	out<-out[order(out$Chromosome,out$Position),]
	out<-out[order(out$Recurrence,decreasing = T),]
	write.table(out,paste0("Results/mostRec",rec,".txt"),row.names=F,quote = F,sep='\t')
}

#feats<-c(read.delim("Databases/feat_gencode.txt",h=F,stringsAsFactors=F)) -> use colnames...









#abc<-read.delim("../Databases/most_rec_features.txt",h=T,stringsAsFactors=F)
#abc[,4]<-gsub("_"," ",abc[,4],fixed=T)
#abc[abc[,4]=="notingencode",4]<-"Unknown"
#print(xtable(abc),include.rownames=F,sanitize.text.function=force)

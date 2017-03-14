# TODO: Add comment
# 
# Author: drubay
###############################################################################
ARGS<-commandArgs(trailingOnly=TRUE)
if(length(ARGS)<4){
	stop("At least the number of cores or one input file (Clinvar, 1000 genomes or COSMIC) is missing. Please check the input.")
}
nbcores<-ARGS[1]
ARGS<-ARGS[-1]
for(i in ARGS){
	if(!file.exists(i)){
		stop(paste0("The file ",i," does not exists."))
	}
}

ARGS<-c("Databases/Clinvar.txt","Databases/COSMIC.txt","Databases/data1KG.unmatched","Databases/data1KG.regionclin",
	apply(expand.grid("Databases/data1KG.regioncosrec",c("Conf","Unkn","Total"),as.character(1:20),".txt"),1,paste,collapse=""))

logFile<-NULL

args<-tolower(ARGS)
tmp<-grep("conf",args)
if(length(tmp)==0){
	conf<-FALSE
}else{
	conf<-TRUE
}
rm(tmp)

thresholds<-unique(na.omit(as.numeric(unlist(strsplit(unlist(lapply(strsplit(unlist(args[grep("cos",args)]),".",fixed=T),function(vec) vec[2])), "[^0-9]+")))))
minRec<-max(2,min(thresholds))
maxRec<-max(thresholds)
if(min(thresholds)==1){
	logFile<-"The recurrence threshold 1 was not considered for the benchmark to save memory and time. If you need to do it, you can use the benchmarkOne.R script in the Codes folder."
	write.table(logFile,"benchmarkResults.log",row.names=F,quote=F)
}

#initialize
library(PRROC)
scores<-c("CADD","DANN","FATHMM_MKL","Funseq2","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","snp","somcll","somliver","somlung","sommel")
unmatched<-read.delim(ARGS[grep("unmatched",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(unmatched)<-c("chr","pos","ref","alt","AF","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
unmatched<-subset(unmatched,select=-c(AF))

results<-function(dataPathos,controls,folder=NULL){ #allows to parallelize ?
	if(is.null(folder) | (folder!="Clinvar" & folder!="COSMIC")){stop("Please indicate an output folder (Clinvar or COSMIC only)")}
	eff<-data.frame(controlType=controls,
		pathos=rep(nrow(get(dataPathos)),length(controls)),
		controls=as.integer(unlist(sapply(controls,function(vec){nrow(get(vec,envir=.GlobalEnv))}))))
	savetab<-data.frame(scores=scores,AUC_unmatched=rep(NA,length(scores)),AUC_region=rep(NA,length(scores)))
	savetabpr<-data.frame(scores=scores,PRAUC_unmatched=rep(NA,length(scores)),PRAUC_region=rep(NA,length(scores)))
	direction<-data.frame(unmatched=scores,region=scores)
	for(k in controls){
		sc<-scores
		for(j in 1:length(scores)){
			colPatho<-match(scores[j],colnames(get(dataPathos)))
			colControl<-match(scores[j],colnames(get(k)))
			assign(paste0("roc_",k,"_",scores[j]),roc.curve(scores.class0 = get(dataPathos)[,colPatho], scores.class1 = get(k)[,colControl],curve=T))
			if(get(paste0("roc_",k,"_",scores[j]))$auc<0.5){
				assign(paste0("roc_",k,"_",scores[j]),roc.curve(scores.class0 = -get(dataPathos)[,colPatho], scores.class1 = -get(k)[,colControl],curve=T))
				assign(paste0("pr_",k,"_",scores[j]),pr.curve(scores.class0 = -get(dataPathos)[,colPatho], scores.class1 = -get(k)[,colControl],curve=T))
				sc[j]<-paste0("Negative_",scores[j])
			}else{
				assign(paste0("pr_",k,"_",scores[j]),pr.curve(scores.class0 = get(dataPathos)[,colPatho], scores.class1 = get(k)[,colControl],curve=T))
			}

			roc<-get(paste0("roc_",k,"_",scores[j]))
			savetab[j,grep(k,colnames(savetab))]<-roc$auc
			curvetable<-roc$curve
			colnames(curvetable)<-c("FPR","Sensitivity","Threshold")
			write.table(curvetable,paste0("Results/",folder,"/Figures/ROC_table_",dataPathos,"_",k,"_",scores[j],".txt"),row.names=F,quote=F)
			pr<-get(paste0("pr_",k,"_",scores[j]))
			savetabpr[j,grep(k,colnames(savetabpr))]<-pr$auc.davis.goadrich
			curvetable<-pr$curve
			colnames(curvetable)<-c("Recall","Precision","Threshold")
			write.table(curvetable,paste0("Results/",folder,"/Figures/PR_table_",dataPathos,"_",k,"_",scores[j],".txt"),row.names=F,quote=F)
		}
		direction[,grep(k,colnames(direction))]<-sc
	}
	out<-list(savetab,savetabpr,direction,eff)
	names(out)<-c(paste0("Results/",folder,"/Tables/ROCAUC_",dataPathos,".txt"),
		paste0("Results/",folder,"/Tables/PRAUC_",dataPathos,".txt"),
		paste0("Results/",folder,"/Tables/direction_",dataPathos,".txt"),
		paste0("Results/",folder,"/Tables/eff_",dataPathos,".txt"))
	return(out)
}


Clinvar<-read.delim(ARGS[grep("clinvar",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(Clinvar)<-c("chr","pos","ref","alt","id","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
Clinvar<-subset(Clinvar,select=-c(id))
region<-read.delim(ARGS[grep("regionclin",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(region)<-c("chr","pos","ref","alt","AF","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
region<-subset(region,select=-c(AF))
out<-results("Clinvar",c("unmatched","region"),folder="Clinvar")
for(i in 1:length(out)){write.table(out[[i]],names(out)[i],row.names=F,quote=F)}
rm(Clinvar,region,out)
gc()

logFile<-c(logFile,"Succeed in Clinvar benchmark.")
write.table(logFile,"benchmarkResults.log",row.names=F,quote=F)


COSMIC<-read.delim(ARGS[grep("cosmic",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(COSMIC)<-c("chr","pos","ref","alt","project","recConf","recUnkn","recTotal","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
COSMIC<-subset(COSMIC,select=-c(project))
out<-list()
for(thres in minRec:maxRec){
	#print(thres)
	for(rec in c("recConf","recUnkn","recTotal")){
		#print(rec)
		dataPathos<-paste0("COSMIC",rec,thres)
		assign(dataPathos,COSMIC[COSMIC[,rec]>=thres,])
		if(nrow(get(dataPathos))>0){
			region<-read.delim(ARGS[grep(paste0("regioncos",tolower(rec),thres,".txt"),args)],h=F,stringsAsFactor=F,row.names=NULL)
			if(nrow(region)>0){
				colnames(region)<-c("chr","pos","ref","alt","AF","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
				region<-subset(region,select=-c(AF))
				out<-c(out,results(dataPathos,c("unmatched","region"),folder="COSMIC"))
				rm(region)
				gc()
			}
		}
		rm(list=dataPathos)
		gc()
	}
	logFile<-c(logFile,paste0("Succeed in COSMIC rec ",rec," benchmark."))
	write.table(logFile,"benchmarkResults.log",row.names=F,quote=F)
}
for(i in 1:length(out)){
	write.table(out[[i]],names(out)[i],row.names=F,quote=F)
}

if(minRec==maxRec){
}else{
	for(i in 1:length(out)){
		for(j in 1:length(out[[i]])){
			write.table(out[[i]][[j]],names(out[[i]])[j],row.names=F,quote=F)
		}
	}	
}
rm(COSMIC,out)
gc()


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
library(parallel)
ciAuc<-function(patho,control,type="ROC",nBoot=2000,alpha=0.05,cl=4){
	n0<-length(control)
	n1<-length(patho)
	n<-n0+n1
	aucs<-unlist(mclapply(1:nBoot,function(i){
		samp0<-sample(1:n0,n0,replace=TRUE)
		samp1<-sample(1:n1,n1,replace=TRUE)
		if(type=="ROC"){
			return(as.numeric(roc.curve(scores.class0 = patho[samp1], scores.class1 = control[samp0])$auc))
		}else{
			return(as.numeric(pr.curve(scores.class0 = patho[samp1], scores.class1 = control[samp0])$auc.davis.goadrich))
		}
	},mc.cores = cl))
	aucs<-sort(aucs)
	if(type=="ROC"){
		auc<-as.numeric(roc.curve(scores.class0 = patho, scores.class1 = control)$auc)
	}else{
		auc<-as.numeric(pr.curve(scores.class0 = patho, scores.class1 = control)$auc.davis.goadrich)
	}
	return(paste0(sprintf("%0.3f",auc)," [",sprintf("%0.3f",aucs[nBoot*alpha/2]),"; ",sprintf("%0.3f",aucs[nBoot*(1-alpha/2)]),"]"))
}



results<-function(dataPathos,controls,folder=NULL){ #allows to parallelize ?
	if(is.null(folder) | (folder!="Clinvar" & folder!="COSMIC")){stop("Please indicate an output folder (Clinvar or COSMIC only)")}
	if(length(grep(".",controls,fixed=TRUE))>0){controlName<-strsplit(controls[1],".",fixed=TRUE)[[1]][1]}else{controlName<-""}
	un<-reg<-FALSE
	if(length(grep("unmatched",tolower(controls)))!=0){un<-TRUE}
	if(length(grep("region",tolower(controls)))!=0){reg<-TRUE}
	eff<-data.frame(controlType=controls,
		pathos=rep(nrow(get(dataPathos)),length(controls)),
		controls=as.integer(unlist(sapply(controls,function(vec){nrow(get(vec,envir=.GlobalEnv))}))))
	if(un & reg){
		direction<-data.frame(unmatched=scores,region=scores)
		savetab<-data.frame(scores=scores,AUC_unmatched=rep(NA,length(scores)),AUC_region=rep(NA,length(scores)))
		savetabpr<-data.frame(scores=scores,PRAUC_unmatched=rep(NA,length(scores)),PRAUC_region=rep(NA,length(scores)))
	}else{
		if(un & !reg){
			direction<-data.frame(unmatched=scores)
			savetab<-data.frame(scores=scores,AUC_unmatched=rep(NA,length(scores)))
			savetabpr<-data.frame(scores=scores,PRAUC_unmatched=rep(NA,length(scores)))
		}else{
			if(!un & reg){
				direction<-data.frame(region=scores)
				savetab<-data.frame(scores=scores,AUC_region=rep(NA,length(scores)))
				savetabpr<-data.frame(scores=scores,PRAUC_region=rep(NA,length(scores)))
			}
		}
	}
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
			if(sc[j]==paste0("Negative_",scores[j])){
				savetab[j,grep(strsplit(k,".",fixed=TRUE)[[1]][2],colnames(savetab))]<-ciAuc(-get(dataPathos)[,colPatho],-get(k)[,colControl],type="ROC")
			}else{
				savetab[j,grep(strsplit(k,".",fixed=TRUE)[[1]][2],colnames(savetab))]<-ciAuc(get(dataPathos)[,colPatho],get(k)[,colControl],type="ROC")
			}
			curvetable<-roc$curve
			colnames(curvetable)<-c("FPR","Sensitivity","Threshold")
			write.table(curvetable,paste0("Results/",folder,"/Figures/ROC_table_",dataPathos,"_",k,"_",scores[j],".txt"),row.names=F,quote=F)
			pr<-get(paste0("pr_",k,"_",scores[j]))
			if(sc[j]==paste0("Negative_",scores[j])){
				savetabpr[j,grep(strsplit(k,".",fixed=TRUE)[[1]][2],colnames(savetabpr))]<-ciAuc(-get(dataPathos)[,colPatho],-get(k)[,colControl],type="PR")
			}else{
				savetabpr[j,grep(strsplit(k,".",fixed=TRUE)[[1]][2],colnames(savetabpr))]<-ciAuc(get(dataPathos)[,colPatho],get(k)[,colControl],type="PR")
			}
			curvetable<-pr$curve
			colnames(curvetable)<-c("Recall","Precision","Threshold")
			write.table(curvetable,paste0("Results/",folder,"/Figures/PR_table_",dataPathos,"_",k,"_",scores[j],".txt"),row.names=F,quote=F)
		}
		direction[,grep(strsplit(k,".",fixed=TRUE)[[1]][2],colnames(direction))]<-sc
	}
	out<-list(savetab,savetabpr,direction,eff)
	names(out)<-c(paste0("Results/",folder,"/Tables/ROCAUC_",dataPathos,"_",controlName,".txt"),
		paste0("Results/",folder,"/Tables/PRAUC_",dataPathos,"_",controlName,".txt"),
		paste0("Results/",folder,"/Tables/direction_",dataPathos,"_",controlName,".txt"),
		paste0("Results/",folder,"/Tables/eff_",dataPathos,"_",controlName,".txt"))
	return(out)
}


scores<-c("CADD","DANN","FATHMM_MKL","Funseq2","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","snp","somcll","somliver","somlung","sommel")
data1KG.unmatched<-read.delim(ARGS[grep("unmatched",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(data1KG.unmatched)<-c("chr","pos","ref","alt","AF","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
data1KG.unmatched<-subset(data1KG.unmatched,select=-c(AF))


#Clinvar
Clinvar<-read.delim(ARGS[grep("clinvar",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(Clinvar)<-c("chr","pos","ref","alt","id","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
Clinvar<-subset(Clinvar,select=-c(id))
benign.unmatched<-read.delim("Databases/benign.txt",h=F,stringsAsFactor=F,row.names=NULL)
colnames(benign.unmatched)<-c("chr","pos","ref","alt","id","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
benign.unmatched<-subset(benign.unmatched,select=-c(id))

benign.region<-read.delim(ARGS[grep("benign.regionclin",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(benign.region)<-c("chr","pos","ref","alt","id","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
benign.region<-subset(benign.region,select=-c(id))
out<-results("Clinvar",c("benign.unmatched","benign.region"),folder="Clinvar")
for(i in 1:length(out)){write.table(out[[i]],names(out)[i],row.names=F,quote=F,sep='\t')}
rm(benign.region,out)

data1KG.region<-read.delim(ARGS[grep("data1kg.regionclin",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(data1KG.region)<-c("chr","pos","ref","alt","AF","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
data1KG.region<-subset(data1KG.region,select=-c(AF))
out<-results("Clinvar",c("data1KG.unmatched","data1KG.region"),folder="Clinvar")
for(i in 1:length(out)){write.table(out[[i]],names(out)[i],row.names=F,quote=F,sep='\t')}
rm(Clinvar,data1KG.region,out)
gc()

logFile<-c(logFile,"Succeed in Clinvar benchmark.")
write.table(logFile,"benchmarkResults.log",row.names=F,quote=F)


#COSMIC
COSMIC<-read.delim(ARGS[grep("cosmic",args)],h=F,stringsAsFactor=F,row.names=NULL)
colnames(COSMIC)<-c("chr","pos","ref","alt","project","recConf","recUnkn","recTotal","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
COSMIC<-subset(COSMIC,select=-c(project))
write.table(COSMIC$recTotal[COSMIC$recTotal>=2],"Databases/cosRec.txt",sep="\t",row.names=F,quote=F)

#COSMIC vs 1KG controls
out<-list()
for(thres in maxRec:minRec){
	print(thres)
		dataPathos<-paste0("COSMICrecTotal",thres)
		assign(dataPathos,COSMIC[COSMIC[,"recTotal"]>=thres,])
		if(nrow(get(dataPathos))>0){
			data1KG.region<-read.delim(ARGS[grep(tolower(paste0("data1KG.regioncosrecTotal",thres,".txt")),args)],h=F,stringsAsFactor=F,row.names=NULL)
			if(nrow(data1KG.region)>0){
				colnames(data1KG.region)<-c("chr","pos","ref","alt","AF","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
				data1KG.region<-subset(data1KG.region,select=-c(AF))
				out<-c(out,results(dataPathos,c("data1KG.unmatched","data1KG.region"),folder="COSMIC"))
				rm(data1KG.region)
				gc()
			}
		}
		rm(list=dataPathos)
		gc()

	logFile<-c(logFile,paste0("Succeed in COSMIC rec ",thres," benchmark."))
	write.table(logFile,"benchmarkResults.log",row.names=F,quote=F)
}

for(i in 1:length(out)){
	write.table(out[[i]],names(out)[i],row.names=F,quote=F,sep='\t')
}

rm(COSMIC,out)
gc()



#COSMIC vs ClinVar benign controls
out<-list()
for(thres in maxRec:minRec){
	dataPathos<-paste0("COSMICrecTotal",thres)
	assign(dataPathos,COSMIC[COSMIC[,"recTotal"]>=thres,])
	if(nrow(get(dataPathos))>0){
		if(file.info(ARGS[grep(paste0("databases/benign.regioncosrecTotal",thres,".txt"),args)])$size>0){
			benign.region<-read.delim(ARGS[grep(paste0("databases/benign.regioncosrecTotal",thres,".txt"),args)],h=F,stringsAsFactor=F,row.names=NULL)
			if(nrow(benign.region)>0){
				colnames(benign.region)<-c("chr","pos","ref","alt","id","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","CADD","DANN","FATHMM_MKL","Funseq2","snp","somcll","somliver","somlung","sommel")
				benign.region<-subset(benign.region,select=-c(id))
				out<-c(out,results(dataPathos,c("benign.unmatched","benign.region"),folder="COSMIC"))
				rm(benign.region)
				gc()
			}
		}else{
			out<-c(out,results(dataPathos,"benign.unmatched",folder="COSMIC"))
		}
	}
	rm(list=dataPathos)
	gc()
}

for(i in 1:length(out)){
	write.table(out[[i]],names(out)[i],row.names=F,quote=F,sep='\t')
}

rm(COSMIC,out)
gc()








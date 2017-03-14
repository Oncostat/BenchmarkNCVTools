



ARGS<-commandArgs(trailingOnly=TRUE)
#nbCores<-ARGS
minRec<-max(2,ARGS[1])
maxRec<-ARGS[2]

#ARGS<-ARGS[-c(1,2)]

legCex<-1
labAxisCex<-1.5
stroke_width<-1.5
colors<-c(1,1,2,"darkgreen",4,4,4,5,rep("purple",4))
typeline<-c(1,2,1,1,1,2,3,1,1,2,3,4)
scores<-c("CADD","DANN","FATHMM_MKL","Funseq2","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","snp","somcll","somliver","somlung","sommel")


for(folder in c("Clinvar","COSMIC")){ #constitute the list of the files and make checks
	listFiles<-dir(path=paste0("Results/",folder,"/Figures/"))
	if(folder=="COSMIC"){
		conf<-I(length(grep("Conf",listFiles))!=0)
		unkn<-I(length(grep("Unkn",listFiles))!=0)
		total<-I(length(grep("Total",listFiles))!=0)
	}
	if(length(grep("unmatched",listFiles))!=0){
		assign(paste0(folder,"_controls"),"unmatched")
		if(length(grep("region",listFiles))!=0){
			assign(paste0(folder,"_controls"),c(get(paste0(folder,"_controls")),"region"))
		}else{
			warning("Only the results for the unmatched control set were found for the analysis of the ",folder," set (not those for the region control set)")
		}
	}else{
		if(length(grep("region",listFiles))!=0){
			assign(paste0(folder,"_controls"),"region")
			warning("Only the results for the region control set were found for the analysis of the ",folder," set (not those for the unmatched control set)")
		}else{
			stop("No results were found for the analysis of the ",folder," set. Check that the results tables were generated with unmatched and/or region word in their file names")
		}
	}
}


legCex<-1
labAxisCex<-1.5
stroke_width<-1.5
colors<-c(1,1,2,"darkgreen",4,4,4,5,rep("purple",4))
typeline<-c(1,2,1,1,1,2,3,1,1,2,3,4)
scores<-c("CADD","DANN","FATHMM_MKL","Funseq2","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","snp","somcll","somliver","somlung","sommel")

#Clinvar
if(file.exists("Results/Clinvar/Tables/ROCAUC_Clinvar.txt") & file.exists("Results/Clinvar/Tables/PRAUC_Clinvar.txt") & 
	file.exists("Results/Clinvar/Tables/direction_Clinvar.txt") & file.exists("Results/Clinvar/Tables/eff_Clinvar.txt")){
	direction<-read.delim(paste0("Results/Clinvar/Tables/direction_Clinvar.txt"),sep=" ",stringsAsFactors=F)
	eff<-read.delim(paste0("Results/Clinvar/Tables/eff_Clinvar.txt"),sep=" ",stringsAsFactors=F)
	trainingSet<-colnames(direction)
	#paste0(toupper(substr(tolower(trainingSet),1,1)),substr(tolower(trainingSet),2,nchar(trainingSet)))
	for(j in trainingSet){
		for(r in c("ROC","PR")){
			auc<-read.delim(paste0("Results/Clinvar/Tables/",r,"AUC_Clinvar.txt"),sep=" ",stringsAsFactors=F)
			npathos<-eff[grep(j,eff[,1]),grep("pathos",colnames(eff))]
			ncontrols<-eff[grep(j,eff[,1]),grep("controls",colnames(eff))]
			colauc<-grep(j,colnames(auc))
			jpeg(paste0("Results/Clinvar/Figures/",r,"_",j,"_Clinvar.jpeg"),width=7,height=6,units="in",res=350)
			par(mar=c(4.1,4.1,0.1,2.1))
			for(n in 1:length(scores)){
				rowauc<-grep(scores[n],auc$scores)
				tmp<-paste0("Results/Clinvar/Figures/",r,"_table_Clinvar_",j,"_",scores[n],".txt")
				if(!file.exists(tmp)){
					print(paste0(tmp," does not exist"))
					next
				}
				curvesc<-read.delim(tmp,sep=" ")
				if(n==1){ #change labels according to PR or ROC
					plot(curvesc[,1],curvesc[,2],type="l",
						xlab=ifelse(r=="ROC","False Positive Rate (1 - Specificity)","Sensitivity (Recall)"),
						ylab=ifelse(r=="ROC","Sensitivity (Recall)","Precision"),
						font.lab=2,cex.axis=0.85,cex.lab=labAxisCex,
						#main=paste0("Controls: ",m),
						lwd=stroke_width,
						xlim=c(0,1),ylim=c(0,1))
					leg<-paste0(direction[n,grep(j,colnames(direction))]," (AUC = ",sprintf("%.2f",round(auc[n,grep(j,colnames(auc))],2)),")")
				}else{
					lines(curvesc[,1],curvesc[,2],col=colors[n],lty=typeline[n],lwd=stroke_width)
					leg<-c(leg,paste0(direction[n,grep(j,colnames(direction))]," (AUC = ",sprintf("%.2f",round(auc[n,grep(j,colnames(auc))],2)),")"))
				}
				
			}
			posleg<-ifelse(r=="ROC","bottomright",ifelse(npathos/(ncontrols+npathos)>0.4,"bottomright","topright"))
			if(r=="ROC"){
				curve(1*x,add=T,col="lightgray",lwd=2)
				if(max(auc[,grep(j,colnames(auc))])>0.85){
					legend("bottomleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}else{
					legend("topleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}				
			}else{
				leg<-c(leg,paste0("prevalence = ",sprintf("%.2f",npathos/(ncontrols+npathos))))
				colors<-c(colors,"lightgray")
				abline(h=npathos/(ncontrols+npathos),col="lightgray",lwd=2)
				if(npathos/(ncontrols+npathos)>0.4){
					legend("bottomleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}else{
					legend("topleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}				
			}
			legend(posleg, legend=leg,col=colors,lwd=stroke_width,bty="n",cex=legCex,lty=typeline)
			dev.off()
		}
	}
}


#COSMIC
for(i in c("recConf","recUnkn","recTotal")){
	for(rec in minRec:maxRec){
		if(file.exists(paste0("Results/COSMIC/Tables/ROCAUC_COSMIC",i,rec,".txt")) & file.exists(paste0("Results/COSMIC/Tables/PRAUC_COSMIC",i,rec,".txt")) & 
			file.exists(paste0("Results/COSMIC/Tables/direction_COSMIC",i,rec,".txt")) & file.exists(paste0("Results/COSMIC/Tables/eff_COSMIC",i,rec,".txt"))){
			direction<-read.delim(paste0("Results/COSMIC/Tables/direction_COSMIC",i,rec,".txt"),sep=" ",stringsAsFactors=F)
			eff<-read.delim(paste0("Results/COSMIC/Tables/eff_COSMIC",i,rec,".txt"),sep=" ",stringsAsFactors=F)
			trainingSet<-colnames(direction)
			#paste0(toupper(substr(tolower(trainingSet),1,1)),substr(tolower(trainingSet),2,nchar(trainingSet)))
			for(j in trainingSet){
				for(r in c("ROC","PR")){
					auc<-read.delim(paste0("Results/COSMIC/Tables/",r,"AUC_COSMIC",i,rec,".txt"),sep=" ",stringsAsFactors=F)
					npathos<-eff[grep(j,eff[,1]),grep("pathos",colnames(eff))]
					ncontrols<-eff[grep(j,eff[,1]),grep("controls",colnames(eff))]
					colauc<-grep(j,colnames(auc))
					jpeg(paste0("Results/COSMIC/Figures/",r,"_",j,"_COSMIC",i,rec,".jpeg"),width=7,height=6,units="in",res=350)
					par(mar=c(4.1,4.1,0.1,2.1))
					for(n in 1:length(scores)){
						rowauc<-grep(scores[n],auc$scores)
						tmp<-paste0("Results/COSMIC/Figures/",r,"_table_COSMIC",i,rec,"_",j,"_",scores[n],".txt")
						if(!file.exists(tmp)){
							print(paste0(tmp," does not exist"))
							next
						}
						curvesc<-read.delim(tmp,sep=" ")
						if(n==1){ #change labels according to PR or ROC
							plot(curvesc[,1],curvesc[,2],type="l",
								xlab=ifelse(r=="ROC","False Positive Rate (1 - Specificity)","Sensitivity (Recall)"),
								ylab=ifelse(r=="ROC","Sensitivity (Recall)","Precision"),
								font.lab=2,cex.axis=0.85,cex.lab=labAxisCex,
								#main=paste0("Controls: ",m),
								lwd=stroke_width,
								xlim=c(0,1),ylim=c(0,1))
							leg<-paste0(direction[n,grep(j,colnames(direction))]," (AUC = ",sprintf("%.2f",round(auc[n,grep(j,colnames(auc))],2)),")")
						}else{
							lines(curvesc[,1],curvesc[,2],col=colors[n],lty=typeline[n],lwd=stroke_width)
							leg<-c(leg,paste0(direction[n,grep(j,colnames(direction))]," (AUC = ",sprintf("%.2f",round(auc[n,grep(j,colnames(auc))],2)),")"))
						}
						
					}
					posleg<-ifelse(r=="ROC","bottomright",ifelse(npathos/(ncontrols+npathos)>0.4,"bottomright","topright"))
					if(r=="ROC"){
						curve(1*x,add=T,col="lightgray",lwd=2)
						if(max(auc[,grep(j,colnames(auc))])>0.85){
							legend("bottomleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
						}else{
							legend("topleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
						}
					}else{
						leg<-c(leg,paste0("prevalence = ",sprintf("%.2f",npathos/(ncontrols+npathos))))
						colors<-c(colors,"lightgray")
						abline(h=npathos/(ncontrols+npathos),col="lightgray",lwd=2)
						if(npathos/(ncontrols+npathos)>0.4){
							legend("bottomleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
						}else{
							legend("topleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
						}
					}
					legend(posleg, legend=leg,col=colors,lwd=stroke_width,bty="n",cex=legCex,lty=typeline)
					dev.off()
				}
			}
		}
	}	
}



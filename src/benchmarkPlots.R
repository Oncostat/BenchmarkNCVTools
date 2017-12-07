
ARGS<-commandArgs(trailingOnly=TRUE)
minRec<-max(2,ARGS[1])
maxRec<-ARGS[2]

legCex<-0.75
labAxisCex<-1.5
stroke_width<-1.5
colors<-c(1,1,2,"darkgreen",4,4,4,5,rep("purple",4))
typeline<-c(1,2,1,1,1,2,3,1,1,2,3,4)
scores<-c("CADD","DANN","FATHMM_MKL","Funseq2","GWAVA_region","GWAVA_TSS","GWAVA_unmatched","snp","somcll","somliver","somlung","sommel")


#########################################
#1KG controls
#########################################
for(folder in c("Clinvar","COSMIC")){ #constitute the list of the files and make checks
	listFiles<-dir(path=paste0("Results/",folder,"/Figures/"))
	if(folder=="COSMIC"){
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


#Clinvar
if(file.exists("Results/Clinvar/Tables/ROCAUC_Clinvar_data1KG.txt") & file.exists("Results/Clinvar/Tables/PRAUC_Clinvar_data1KG.txt") & 
	file.exists("Results/Clinvar/Tables/direction_Clinvar_data1KG.txt") & file.exists("Results/Clinvar/Tables/eff_Clinvar_data1KG.txt")){
	direction<-read.delim(paste0("Results/Clinvar/Tables/direction_Clinvar_data1KG.txt"),stringsAsFactors=F)
	for(d in 1:ncol(direction)){
		direction[,d]<-gsub("Negative_","Decreasing ",direction[,d])
		direction[,d]<-gsub("FATHMM_MKL","FATHMM-MKL",direction[,d])
		direction[,d]<-gsub("snp","SNP",direction[,d])
		direction[,d]<-gsub("somcll","SOM[CLL]",direction[,d])
		direction[,d]<-gsub("somliver","SOM[Liver]",direction[,d])
		direction[,d]<-gsub("somlung","SOM[Lung]",direction[,d])
		direction[,d]<-gsub("sommel","SOM[Melanoma]",direction[,d])
		direction[,d]<-gsub("GWAVA_","GWAVA[",direction[,d])
		direction[grep("GWAVA",direction[,d]),d]<-paste0(direction[grep("GWAVA",direction[,d]),d],"]")
		direction[,d]<-gsub(" ","~",direction[,d])
	}
	eff<-read.delim(paste0("Results/Clinvar/Tables/eff_Clinvar_data1KG.txt"),stringsAsFactors=F)
	trainingSet<-colnames(direction)
	for(j in trainingSet){
		for(r in c("ROC","PR")){
			auc<-read.delim(paste0("Results/Clinvar/Tables/",r,"AUC_Clinvar_data1KG.txt"),stringsAsFactors=F)
			npathos<-eff[grep(j,eff[,1]),grep("pathos",colnames(eff))]
			ncontrols<-eff[grep(j,eff[,1]),grep("controls",colnames(eff))]
			colauc<-grep(j,colnames(auc))
			postscript(paste0("Results/Clinvar/Figures/",r,"_",j,"_Clinvar.eps"),width=7,height=6)
			par(mar=c(4.1,4.1,0.1,2.1))
			for(n in 1:length(scores)){
				rowauc<-grep(scores[n],auc$scores)
				tmp<-paste0("Results/Clinvar/Figures/",r,"_table_Clinvar_data1KG.",j,"_",scores[n],".txt")
				if(!file.exists(tmp)){
					print(paste0(tmp," does not exist"))
					next
				}
				curvesc<-read.delim(tmp,sep=" ")
				if(n==1){
					plot(curvesc[,1],curvesc[,2],type="l",
						xlab=ifelse(r=="ROC","False Positive Rate (1 - Specificity)","Sensitivity (Recall)"),
						ylab=ifelse(r=="ROC","Sensitivity (Recall)","Precision"),
						font.lab=2,cex.axis=0.85,cex.lab=labAxisCex,
						lwd=stroke_width,
						xlim=c(0,1),ylim=c(0,1))
					leg<-paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'")
				}else{
					lines(curvesc[,1],curvesc[,2],col=colors[n],lty=typeline[n],lwd=stroke_width)
					leg<-c(leg,paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'"))
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
				leg<-c(leg,paste0("'prevalence = ",sprintf("%.3f",npathos/(ncontrols+npathos)),"'"))
				colors<-c(colors,"lightgray")
				abline(h=npathos/(ncontrols+npathos),col="lightgray",lwd=2)
				if(npathos/(ncontrols+npathos)>0.4){
					legend("bottomleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}else{
					legend("topleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}				
			}
			legend(posleg, legend=parse(text=leg),col=colors,lwd=stroke_width,bty="n",cex=legCex,lty=typeline)
			dev.off()
		}
	}
}


#COSMIC
for(rec in minRec:maxRec){
	if(file.exists(paste0("Results/COSMIC/Tables/ROCAUC_COSMICrecTotal",rec,"_data1KG.txt")) & file.exists(paste0("Results/COSMIC/Tables/PRAUC_COSMICrecTotal",rec,"_data1KG.txt")) & 
		file.exists(paste0("Results/COSMIC/Tables/direction_COSMICrecTotal",rec,"_data1KG.txt")) & file.exists(paste0("Results/COSMIC/Tables/eff_COSMICrecTotal",rec,"_data1KG.txt"))){
		direction<-read.delim(paste0("Results/COSMIC/Tables/direction_COSMICrecTotal",rec,"_data1KG.txt"),stringsAsFactors=F)
		for(d in 1:ncol(direction)){
			direction[,d]<-gsub("Negative_","Decreasing ",direction[,d])
			direction[,d]<-gsub("FATHMM_MKL","FATHMM-MKL",direction[,d])
			direction[,d]<-gsub("snp","SNP",direction[,d])
			direction[,d]<-gsub("somcll","SOM[CLL]",direction[,d])
			direction[,d]<-gsub("somliver","SOM[Liver]",direction[,d])
			direction[,d]<-gsub("somlung","SOM[Lung]",direction[,d])
			direction[,d]<-gsub("sommel","SOM[Melanoma]",direction[,d])
			direction[,d]<-gsub("GWAVA_","GWAVA[",direction[,d])
			direction[grep("GWAVA",direction[,d]),d]<-paste0(direction[grep("GWAVA",direction[,d]),d],"]")
			direction[,d]<-gsub(" ","~",direction[,d])
		}
		eff<-read.delim(paste0("Results/COSMIC/Tables/eff_COSMICrecTotal",rec,"_data1KG.txt"),stringsAsFactors=F)
		trainingSet<-colnames(direction)
		for(j in trainingSet){
			for(r in c("ROC","PR")){
				auc<-read.delim(paste0("Results/COSMIC/Tables/",r,"AUC_COSMICrecTotal",rec,"_data1KG.txt"),stringsAsFactors=F)
				npathos<-eff[grep(j,eff[,1]),grep("pathos",colnames(eff))]
				ncontrols<-eff[grep(j,eff[,1]),grep("controls",colnames(eff))]
				colauc<-grep(j,colnames(auc))
				postscript(paste0("Results/COSMIC/Figures/",r,"_",j,"_COSMICrecTotal",rec,".eps"),width=7,height=6)
				par(mar=c(4.1,4.1,0.1,2.1))
				for(n in 1:length(scores)){
					rowauc<-grep(scores[n],auc$scores)
					tmp<-paste0("Results/COSMIC/Figures/",r,"_table_COSMICrecTotal",rec,"_",j,"_",scores[n],".txt")
					if(!file.exists(tmp)){
						print(paste0(tmp," does not exist"))
						next
					}
					curvesc<-read.delim(tmp,sep=" ")
					if(n==1){
						plot(curvesc[,1],curvesc[,2],type="l",
							xlab=ifelse(r=="ROC","False Positive Rate (1 - Specificity)","Sensitivity (Recall)"),
							ylab=ifelse(r=="ROC","Sensitivity (Recall)","Precision"),
							font.lab=2,cex.axis=0.85,cex.lab=labAxisCex,
							lwd=stroke_width,
							xlim=c(0,1),ylim=c(0,1))
						leg<-paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'")
					}else{
						lines(curvesc[,1],curvesc[,2],col=colors[n],lty=typeline[n],lwd=stroke_width)
						leg<-c(leg,paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'"))
					}
					
				}
				posleg<-ifelse(r=="ROC","bottomright",ifelse(npathos/(ncontrols+npathos)>0.4,"bottomright","topright"))
				if(r=="ROC"){
					curve(1*x,add=T,col="lightgray",lwd=2)
					if(max(auc[,grep(j,colnames(auc))])>0.85){
						legend("bottomleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
					}else{
						legend("topleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
					}
				}else{
					leg<-c(leg,paste0("'prevalence = ",sprintf("%.3f",npathos/(ncontrols+npathos)),"'"))
					colors<-c(colors,"lightgray")
					abline(h=npathos/(ncontrols+npathos),col="lightgray",lwd=2)
					if(npathos/(ncontrols+npathos)>0.4){
						legend("bottomleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
					}else{
						legend("topleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
					}
				}
				legend(posleg, legend=parse(text=leg),col=colors,lwd=stroke_width,bty="n",cex=legCex,lty=typeline)
				dev.off()
			}
		}
	}
}	







#########################################
#ClinVar benign controls
#########################################
for(folder in c("Clinvar","COSMIC")){ #constitute the list of the files and make checks
	listFiles<-dir(path=paste0("Results/",folder,"/Figures/"))
	if(folder=="COSMIC"){
		cont<-unkn<-FALSE
		total<-I(length(grep("Total",listFiles))!=0)
	}
	if(length(grep("benign.unmatched",listFiles))!=0){
		assign(paste0(folder,"_controls"),"unmatched")
		if(length(grep("region",listFiles))!=0){
			assign(paste0(folder,"_controls"),c(get(paste0(folder,"_controls")),"region"))
		}else{
			warning("Only the results for the unmatched control set were found for the analysis of the ",folder," set (not those for the region control set)")
		}
	}else{
		if(length(grep("benign.region",listFiles))!=0){
			assign(paste0(folder,"_controls"),"region")
			warning("Only the results for the region control set were found for the analysis of the ",folder," set (not those for the unmatched control set)")
		}else{
			stop("No results were found for the analysis of the ",folder," set. Check that the results tables were generated with unmatched and/or region word in their file names")
		}
	}
}


#Clinvar
if(file.exists("Results/Clinvar/Tables/ROCAUC_Clinvar_benign.txt") & file.exists("Results/Clinvar/Tables/PRAUC_Clinvar_benign.txt") & 
	file.exists("Results/Clinvar/Tables/direction_Clinvar_benign.txt") & file.exists("Results/Clinvar/Tables/eff_Clinvar_benign.txt")){
	direction<-read.delim(paste0("Results/Clinvar/Tables/direction_Clinvar_benign.txt"),stringsAsFactors=F)
	for(d in 1:ncol(direction)){
		direction[,d]<-gsub("Negative_","Decreasing ",direction[,d])
		direction[,d]<-gsub("FATHMM_MKL","FATHMM-MKL",direction[,d])
		direction[,d]<-gsub("snp","SNP",direction[,d])
		direction[,d]<-gsub("somcll","SOM[CLL]",direction[,d])
		direction[,d]<-gsub("somliver","SOM[Liver]",direction[,d])
		direction[,d]<-gsub("somlung","SOM[Lung]",direction[,d])
		direction[,d]<-gsub("sommel","SOM[Melanoma]",direction[,d])
		direction[,d]<-gsub("GWAVA_","GWAVA[",direction[,d])
		direction[grep("GWAVA",direction[,d]),d]<-paste0(direction[grep("GWAVA",direction[,d]),d],"]")
		direction[,d]<-gsub(" ","~",direction[,d])
	}
	eff<-read.delim(paste0("Results/Clinvar/Tables/eff_Clinvar_benign.txt"),stringsAsFactors=F)
	trainingSet<-colnames(direction)
	for(j in trainingSet){
		for(r in c("ROC","PR")){
			auc<-read.delim(paste0("Results/Clinvar/Tables/",r,"AUC_Clinvar_benign.txt"),stringsAsFactors=F)
			npathos<-eff[grep(j,eff[,1]),grep("pathos",colnames(eff))]
			ncontrols<-eff[grep(j,eff[,1]),grep("controls",colnames(eff))]
			colauc<-grep(j,colnames(auc))
			postscript(paste0("Results/Clinvar/Figures/",r,"_",j,"_Clinvar_benign.eps"),width=7,height=6)
			par(mar=c(4.1,4.1,0.1,2.1))
			for(n in 1:length(scores)){
				rowauc<-grep(scores[n],auc$scores)
				tmp<-paste0("Results/Clinvar/Figures/",r,"_table_Clinvar_benign.",j,"_",scores[n],".txt")
				if(!file.exists(tmp)){
					print(paste0(tmp," does not exist"))
					next
				}
				curvesc<-read.delim(tmp,sep=" ")
				if(n==1){
					plot(curvesc[,1],curvesc[,2],type="l",
						xlab=ifelse(r=="ROC","False Positive Rate (1 - Specificity)","Sensitivity (Recall)"),
						ylab=ifelse(r=="ROC","Sensitivity (Recall)","Precision"),
						font.lab=2,cex.axis=0.85,cex.lab=labAxisCex,
						lwd=stroke_width,
						xlim=c(0,1),ylim=c(0,1))
					leg<-paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'")
				}else{
					lines(curvesc[,1],curvesc[,2],col=colors[n],lty=typeline[n],lwd=stroke_width)
					leg<-c(leg,paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'"))
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
				leg<-c(leg,paste0("'prevalence = ",sprintf("%.3f",npathos/(ncontrols+npathos)),"'"))
				colors<-c(colors,"lightgray")
				abline(h=npathos/(ncontrols+npathos),col="lightgray",lwd=2)
				if(npathos/(ncontrols+npathos)>0.4){
					legend("bottomleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}else{
					legend("topleft",legend=paste0(c("n ClinVar = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
				}				
			}
			legend(posleg, legend=parse(text=leg),col=colors,lwd=stroke_width,bty="n",cex=legCex,lty=typeline)
			dev.off()
		}
	}
}


#COSMIC
for(rec in minRec:maxRec){
	if(file.exists(paste0("Results/COSMIC/Tables/ROCAUC_COSMICrecTotal",rec,".txt")) & file.exists(paste0("Results/COSMIC/Tables/PRAUC_COSMICrecTotal",rec,".txt")) & 
		file.exists(paste0("Results/COSMIC/Tables/direction_COSMICrecTotal",rec,".txt")) & file.exists(paste0("Results/COSMIC/Tables/eff_COSMICrecTotal",rec,".txt"))){
		direction<-read.delim(paste0("Results/COSMIC/Tables/direction_COSMICrecTotal",rec,"_benign.txt"),stringsAsFactors=F)
		for(d in 1:ncol(direction)){
			direction[,d]<-gsub("Negative_","Decreasing ",direction[,d])
			direction[,d]<-gsub("FATHMM_MKL","FATHMM-MKL",direction[,d])
			direction[,d]<-gsub("snp","SNP",direction[,d])
			direction[,d]<-gsub("somcll","SOM[CLL]",direction[,d])
			direction[,d]<-gsub("somliver","SOM[Liver]",direction[,d])
			direction[,d]<-gsub("somlung","SOM[Lung]",direction[,d])
			direction[,d]<-gsub("sommel","SOM[Melanoma]",direction[,d])
			direction[,d]<-gsub("GWAVA_","GWAVA[",direction[,d])
			direction[grep("GWAVA",direction[,d]),d]<-paste0(direction[grep("GWAVA",direction[,d]),d],"]")
			direction[,d]<-gsub(" ","~",direction[,d])
		}
		if(file.exists(paste0("Results/COSMIC/Tables/eff_COSMICrecTotal",rec,"_benign.txt"))){
			eff<-read.delim(paste0("Results/COSMIC/Tables/eff_COSMICrecTotal",rec,"_benign.txt"),stringsAsFactors=F)
			trainingSet<-colnames(direction)
			for(j in trainingSet){
				for(r in c("ROC","PR")){
					nextCurve<-FALSE
					if(file.exists(paste0("Results/COSMIC/Tables/",r,"AUC_COSMICrecTotal",rec,"_benign.txt"))){
						auc<-read.delim(paste0("Results/COSMIC/Tables/",r,"AUC_COSMICrecTotal",rec,"_benign.txt"),stringsAsFactors=F)
						npathos<-eff[grep(j,eff[,1]),grep("pathos",colnames(eff))]
						ncontrols<-eff[grep(j,eff[,1]),grep("controls",colnames(eff))]
						colauc<-grep(j,colnames(auc))
						postscript(paste0("Results/COSMIC/Figures/",r,"_",j,"_COSMICrecTotal",rec,"_benign.eps"),width=7,height=6)
						par(mar=c(4.1,4.1,0.1,2.1))
						for(n in 1:length(scores)){
							rowauc<-grep(scores[n],auc$scores)
							tmp<-paste0("Results/COSMIC/Figures/",r,"_table_COSMICrecTotal",rec,"_benign.",j,"_",scores[n],".txt")
							curvesc<-read.delim(tmp,sep=" ")
							if(n==1){ 
								plot(curvesc[,1],curvesc[,2],type="l",
									xlab=ifelse(r=="ROC","False Positive Rate (1 - Specificity)","Sensitivity (Recall)"),
									ylab=ifelse(r=="ROC","Sensitivity (Recall)","Precision"),
									font.lab=2,cex.axis=0.85,cex.lab=labAxisCex,
									lwd=stroke_width,
									xlim=c(0,1),ylim=c(0,1))
								leg<-paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'")
							}else{
								lines(curvesc[,1],curvesc[,2],col=colors[n],lty=typeline[n],lwd=stroke_width)
								leg<-c(leg,paste0(direction[n,grep(j,colnames(direction))],"~'(AUC = ",auc[n,grep(j,colnames(auc))],")'"))
							}
						}
						posleg<-ifelse(r=="ROC","bottomright",ifelse(npathos/(ncontrols+npathos)>0.4,"bottomright","topright"))
						if(r=="ROC"){
							curve(1*x,add=T,col="lightgray",lwd=2)
							if(max(auc[,grep(j,colnames(auc))])>0.85){
								legend("bottomleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
							}else{
								legend("topleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
							}
						}else{
							leg<-c(leg,paste0("'prevalence = ",sprintf("%.3f",npathos/(ncontrols+npathos)),"'"))
							colors<-c(colors,"lightgray")
							abline(h=npathos/(ncontrols+npathos),col="lightgray",lwd=2)
							if(npathos/(ncontrols+npathos)>0.4){
								legend("bottomleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
							}else{
								legend("topleft",legend=paste0(c("n COSMIC = ","n controls = "),c(npathos,ncontrols)),bty='n',cex=legCex)
							}
						}
						legend(posleg, legend=parse(text=leg),col=colors,lwd=stroke_width,bty="n",cex=legCex,lty=typeline)
						dev.off()
					}else{
						print(paste0("There is no file for ",rec," recurrences for the ",j," benchmark"))
						nextCurve<-TRUE
						break
					}
				}
			}
		}
	}
}	


###################################
#Variables initialization
###################################
#Variant input data labels
COSMICsubsets = ["Conf","Unkn","Total"]
#Recurrence threshold range
minrec="2"
maxrec="10"
#minimal allele frequency to consider a 1KG variant as non pathogenic
AF="0.01"
#number of random variant for the random control set
rdmN="100000"
#seed for random sampling
seed="123"

###################################
#Import
###################################
rule CADD_DL:
	output:
		"Scores/CADD_scores.tsv"
	shell:
		"wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz -O Scores/CADD_scores.tsv.gz;"
		"zcat -c Scores/CADD_scores.tsv.gz > {output};"
		"rm Scores/CADD_scores.tsv.gz"

rule DANN_DL:
	output:
		"Scores/DANN_scores.tsv"
	shell:
		"wget -c https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz -O Scores/DANN_scores.tsv.gz;"
		"zcat -c  Scores/DANN_scores.tsv.gz > {output};"
		"rm Scores/DANN_scores.tsv.gz"

#source = https://github.com/HAShihab/fathmm-MKL
rule FATHMM_MKL_DL:
	output:
		"Scores/FATH_scores.tsv"
	shell:
		"wget -c http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current_zerobased.tab.gz -O Scores/FATH_scores.tab.gz;"
		"zcat -c Scores/FATH_scores.tab.gz > {output}"
		"rm Scores/FATH_scores.tab.gz"

#source = http://funseq2.gersteinlab.org/downloads (May 2016)
rule Funseq2_DL:
	output:
		"Scores/Funseq2_scores.tsv"
	shell:
		"wget -c http://archive.gersteinlab.org/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz -O Scores/Funseq2_scores.tsv.gz;"
		"zcat -c  Scores/Funseq2_scores.tsv.gz > {output}"
		"rm Scores/Funseq2_scores.tsv.gz"


rule SNPSOM_DL:
	output:
		"Scores/SNPSOM.bed"
	shell:
		"wget http://biodev.cea.fr/98drivers/Genomewise_SNP_and_SOM_scores.zip -O Scores/SNPSOM.zip;"
		"unzip Scores/SNPSOM.zip -d Scores;"
		"mv Scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_CLL_RF_1Kb.bed Scores/CLL.bed;"
		"mv Scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_liver_cancer_RF_1Kb.bed Scores/liver.bed;"
		"mv Scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_lung_cancer_RF_1Kb.bed Scores/lung.bed;"
		"mv Scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_melanoma_RF_1Kb.bed Scores/melanoma.bed;"
		"tail -n+2 Scores/Genomewise_SNP_and_SOM_scores/Human_genome_Uniform_SNP_score_nc.bed > Scores/tmpsnp.bed;"
		"rm -r Scores/Genomewise_SNP_and_SOM_scores;"
		"rm -r Scores/__MACOSX;"
		"cut -f 4 Scores/liver.bed | paste Scores/CLL.bed - > Scores/tmpsom.bed;"
		"cut -f 4 Scores/lung.bed | paste Scores/tmpsom.bed - > Scores/tmpsnpsom.bed;"
		"cut -f 4 Scores/melanoma.bed | paste Scores/tmpsnpsom.bed -  | tail -n+2 > Scores/tmpsom.bed;"
		"rm Scores/tmpsnpsom.bed Scores/CLL.bed Scores/liver.bed Scores/lung.bed Scores/melanoma.bed;"
		"gawk -i inplace -F'\t' '{{$2-=1 FS $2;}}1' OFS='\t' Scores/tmpsnp.bed;"
		"gawk -i inplace -F'\t' '{{$2-=1 FS $2;}}1' OFS='\t' Scores/tmpsom.bed;"
		"./Codes/coordfusion Scores/tmpsnp.bed > Scores/snp.bed;"
		"./Codes/coordfusion Scores/tmpsom.bed > Scores/som.bed;"
		"rm Scores/tmpsnp.bed Scores/tmpsom.bed;"
		"sed -i 's/chr//g' Scores/som.bed;"
		"echo $'Y\t0\t59373567\tNA\tNA\tNA\tNA' >> Scores/som.bed;"
		"sed -i 's/chr//g' Scores/snp.bed;"
		"./Codes/snpsomDM Scores/snp.bed Scores/som.bed Scores/tmpSNPSOM.bed;"
		"rm Scores/snp.bed Scores/som.bed;"
		"bedtools complement -i Scores/tmpSNPSOM.bed -g Codes/grch37p13 > Scores/SNPSOMcomplement.bed;"
		"sed -i 's/$/\tNA\tNA\tNA\tNA\tNA/g' Scores/SNPSOMcomplement.bed;"
		"sed -i s/X/23/g Scores/tmpSNPSOM.bed;"
		"sed -i s/Y/24/g Scores/tmpSNPSOM.bed;"
		"sed -i s/X/23/g Scores/SNPSOMcomplement.bed;"
		"sed -i s/Y/24/g Scores/SNPSOMcomplement.bed;"
		"cat Scores/tmpSNPSOM.bed Scores/SNPSOMcomplement.bed | sort -nk1,1 -nk2,2 > {output};"
		"rm Scores/tmpSNPSOM.bed Scores/SNPSOMcomplement.bed;"

rule GWAVA_DL:
	output:
		"Scores/GWAVA"
	shell:
		"wget -r -P {output} -np -nH --cut-dirs 5  ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/ ;"

rule KG_DL:
	output:
		"Databases/KG.vcf"
	shell:
		"wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -O Databases/KG.vcf.gz;"
		"zcat Databases/KG.vcf.gz > {output}"

rule Clinvar_DL:	
	output:
		"Databases/Clinvar.vcf"
	shell:
		"wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20160531.vcf.gz -O Databases/Clinvar.vcf.gz;"
		"zcat Databases/Clinvar.vcf.gz > {output};"


rule GENCODE_annotations_DL:
	output:
		"Databases/Gencode_annots.gtf","Databases/annot_gencode.txt","Databases/feat_gencode.txt"
	shell:
		"wget -c ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz -O {output[0]}.gz;"
		"zcat {output[0]}.gz > {output[0]};"
		"tail -n+6 {output[0]} > Databases/tmp.txt;"
		"cut -f 1,4,5 Databases/tmp.txt > Databases/tmp2.txt;"
		"cut -f 3 Databases/tmp.txt | paste Databases/tmp2.txt - > Databases/gencodepos.txt;"
		"rm Databases/tmp2.txt;"
		"sed -i 's/chr//g' Databases/gencodepos.txt;"
		"sed -i 's/X/23/g' Databases/gencodepos.txt;"
		"sed -i 's/Y/24/g' Databases/gencodepos.txt;"
		"sed -i 's/M/25/g' Databases/gencodepos.txt;"
		"cut -f 9 Databases/tmp.txt > Databases/tmp2.txt;"
		"cut -d ';' -f 3,6 Databases/tmp2.txt > Databases/gencodeannots.txt;"
		"rm Databases/tmp2.txt;"
		"sed -i 's/; /\t/g' Databases/gencodeannots.txt;"
		"sed -i 's/gene_type //g' Databases/gencodeannots.txt;"
		"sed -i 's/transcript_type //g' Databases/gencodeannots.txt;"
		"""sed -i 's/"//g' Databases/gencodeannots.txt;"""
		"paste Databases/gencodepos.txt Databases/gencodeannots.txt > Databases/tmp.txt;"
		"sort -n -k1,1 -k2,2 -k3,3  Databases/tmp.txt > {output[1]};"
		"cut -f 4 {output[1]} | sort -u > Databases/tmp.txt;"
		"cut -f 5 {output[1]} | sort -u >> Databases/tmp.txt;"
		"cut -f 6 {output[1]} | sort -u >> Databases/tmp.txt;"
		"sed -i 's/ //g' Databases/tmp.txt;"
		"sort -u Databases/tmp.txt > {output[2]};"
		"rm Databases/tmp.txt;"


rule Non_coding_GENCODE:
	input:
		"Databases/Gencode_annots.gtf" 
	output:
		"Databases/NCS.bed" 
	shell:	
		"sed -i 's/^chr//g' {input};"	
		"egrep $'\tCDS|\tstart|\tstop' {input} | bedtools complement -i - -g Codes/grch37p13 > {output};"  #grch37p13 is a data with chr and their length of the last release (13) of grch37 + 1 to ensure to take all positions with 1-based format
		"awk -i inplace '{$1;$2+=1;$3}1' OFS='\t' {output}"

rule COSMIC_DM:
	output:
		"Databases/COSMIC.raw"
	shell:
		"zcat Databases/CosmicNCV.tsv.gz > Databases/CosmicNCV.tsv;"
		"cut -f 1,6,7,8,9,15 Databases/CosmicNCV.tsv > Databases/CosmicNCV.txt;"
		"""egrep $'\ty' Databases/CosmicNCV.txt | awk -F $'\t'  '{{ print $2,$4,$5,$1,$3 }}' OFS=$'\t'  | sort -u > Databases/cosTmp.txt;""" 
		"./Codes/COSMICDM Databases/cosTmp.txt Databases/cosmicWGS.txt;"
		"rm Databases/cosTmp.txt;"
		"sort -nk1,1 -nk2,2 -k3,4 Databases/cosmicWGS.txt > {output};"
		"rm Databases/cosmicWGS.txt;"
 

rule KG_DM:
	input:
		"Databases/KG.vcf"
	output:
		"Databases/data1KG.raw"
	shell:
		"./Codes/rmUseless {input} > Databases/KG.tmp;"
		"./Codes/rmMultiDupl Databases/KG.tmp > {output};"
		"rm  Databases/KG.tmp {input};"
		
rule Clinvar_DM:
	input:
		"Databases/Clinvar.vcf"
	output:
		"Databases/Clinvar.raw"
	shell:
		"a=$(egrep $'#' {input} | wc -l);"
		"b=$(bc <<< $a+1);"
		"tail -n+$b {input} > Databases/clitmp.txt;"
		"egrep CLNSIG=5 Databases/clitmp.txt | cut -f 1,2,4,5,8 > Databases/certain.txt;"
		"sed -i 's/X/23/' Databases/certain.txt;"
		"sed -i 's/Y/24/' Databases/certain.txt;"
		"egrep -v $'^MT\t' Databases/certain.txt | sort -nk1,1 -nk2,2 -k3,3 -k4,4 - > Databases/clitmp.txt;"
		"./Codes/clinvarDM Databases/clitmp.txt {output};"
		"sed -i 's/^23\t/X\t/g' {output};"
		"sed -i 's/^24\t/Y\t/g' {output};"
		"rm Databases/certain.txt Databases/clitmp.txt;"
		"unset a b;"

rule benign:
	input:
		"Databases/Clinvar.vcf"
	output:
		"Databases/benign.raw"
	shell:
		"a=$(egrep $'#' {input} | wc -l);"
		"b=$(bc <<< $a+1);"
		"tail -n+$b {input} > Databases/clitmp.txt;"
		"egrep CLNSIG=2 Databases/clitmp.txt | cut -f 1,2,4,5,8 > Databases/certain.txt;"
		"sed -i 's/X/23/' Databases/certain.txt;"
		"sed -i 's/Y/24/' Databases/certain.txt;"
		"egrep -v $'^MT\t' Databases/certain.txt | sort -nk1,1 -nk2,2 -k3,3 -k4,4 - > Databases/clitmp.txt;"
		"./Codes/clinvarDM Databases/clitmp.txt {output};"
		"sed -i 's/^23\t/X\t/g' {output};"
		"sed -i 's/^24\t/Y\t/g' {output};"
		"rm Databases/certain.txt Databases/clitmp.txt;"
		"unset a b;"

rule HGMD_DM:
	input:
		"Databases/NCS.bed"
	output:
		"Databases/hgmd.bed","Databases/nCDShgmd.bed"
	shell:
		"tail -n+2 Databases/HGMD_Advanced_Substitutions.bed | cut -f 1,2 > Databases/hgmd.tmp;"
		"sed -i 's/chr//g' Databases/hgmd.tmp;"
		"sort -nk1,1 -nk2,2 Databases/hgmd.tmp > {output[0]};"
		"rm Databases/hgmd.tmp;"
		"gawk -i inplace -F'\t' '{{$2=$2 FS $2;}}1' OFS='\t' {output[0]};"
		"bedtools intersect -a {output[0]} -b {input} > {output[1]}"



###############################################################
#Common DM
rule Remove_CDS_positions:
	input:
		"Databases/{file}.raw","Databases/NCS.bed"
	output:
		"Databases/{file}.bed"
	shell:
		"gawk -i inplace -F'\t' '{{$2=$2 FS $2;}}1' OFS='\t' {input[0]};"
		"bedtools intersect -a {input[0]} -b {input[1]} > {output}"

rule Remove_HGMD_variants:
	input:
		"Databases/{file}.bed","Databases/hgmd.bed"
	output:
		"Databases/{file}.clean"
	shell: 
		"bedtools subtract -a {input[0]} -b {input[1]} > {input[0]}.tmp;"
		"""if [ "$(echo {input[0]})" != "Databases/Clinvar.bed" ]; then bedtools subtract -a {input[0]}.tmp -b Databases/Clinvar.bed > {output}; else cp {input[0]}.tmp {output}; fi;"""
		"rm {input[0]}.tmp;"



###############################################################
#Scoring
rule GWAVA_scoring:
	input:
		"Databases/{file}.clean"
	output:
		"Databases/{file}.gw"
	shell:
		"a={input};"
		"b=(${{a//\// }});"
		"b=${{b[1]}};"
		"a=(${{b//./ }});"
		"a=${{a[0]}};"
		"./Codes/exportGWAVA {input} > Scores/GWAVA/src/$a.tmp;"
		"cd Scores/GWAVA/src;"
		"split -l 1000000 $a.tmp $a.x;" 	#split file because too many lines for RAM saving; ~ 2hours / subdata
		"""myFiles=$(ls | grep "$a.x");""" 
		"for f in $myFiles;" 
		"do "
		"mv $f $f.bed;"
		"python gwava_annotate.py $f.bed annots$a.csv;"
		"python gwava.py region annots$a.csv region$a.bed;"
		"python gwava.py tss annots$a.csv tss$a.bed;"
		"python gwava.py unmatched annots$a.csv unmatched$a.bed;"
		"cat region$a.bed | sort -k4n,4n > regionS$a.bed;"
		"cat tss$a.bed | sort -k4n,4n > tssS$a.bed;"
		"cat unmatched$a.bed | sort -k4n,4n > unmatchedS$a.bed;"
		"rm annots$a.csv region$a.bed tss$a.bed unmatched$a.bed;"
		"cut -f 5 tssS$a.bed | paste regionS$a.bed - > tmp$a.bed;"
		"cut -f 5 unmatchedS$a.bed | paste tmp$a.bed - > tmp2$a.bed;"
		"cat tmp2$a.bed >> ../../../{input}.tmp;"  
		"rm $f.bed tmp$a.bed tmp2$a.bed regionS$a.bed tssS$a.bed unmatchedS$a.bed;"
		"done;"
		"rm $a.tmp;"
		"cd ../../../;"
		"egrep -v NA {input}.tmp > {input}.gw;"
		"rm {input}.tmp ;"
		"./Codes/getGWAVAsc {input} {input}.gw> {output};"
		"rm {input}.gw;"
		"sed -i -r 's/(\s+)?\S+//2' {output};"
		"unset a b myFiles f;"

rule MakeIndex:
	input:
		"Scores/CADD_scores.tsv",
		"Scores/DANN_scores.tsv",
		"Scores/FATH_scores.tsv",
		"Scores/Funseq2_scores.tsv"
	output:
		"Scores/initScores","Scores/index"
	shell:
		"./Codes/makeIndex"



rule Scoring:
	input:
		"Databases/{file}.gw",
		"Scores/SNPSOM.bed",
		"Scores/initScores",
		"Scores/index"
	output:
		"Databases/{file}.txt"
	shell:
		"./Codes/variantScoring {input[0]} {input[0]}.tmp;"
		"egrep -v NA {input[0]}.tmp | sort -nk1,1 -nk2,2 -k4,4 > {input[0]}.temp;" #Remove_unscorable
		"rm {input[0]}.tmp;"
		"./Codes/mergeSnpsom {input[0]}.temp {input[0]}.out;"
		"rm {input[0]}.temp;"
		"egrep -v NA {input[0]}.out | sort -nk1,1 -nk2,2 -k4,4 > {output};" #Remove_unscorable
		"rm {input[0]}.out;"

###############################################################
#Finalization of the DM		
rule COSMIC_recurrence:
	input:
		"Databases/COSMIC.txt"
	output:
		['Databases/COSMICrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)]
	shell:
		"./Codes/filterRec {input} {minrec} {maxrec};"


#Modify KG_subset for region of all COSMIC subset
rule KG_subsets:
	input:
		"Databases/data1KG.txt","Databases/Clinvar.txt",['Databases/COSMICrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)]
	output:
		"Databases/data1KG.unmatched","Databases/data1KG.regionclin",['Databases/data1KG.regioncosrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)]
	shell:
		"./Codes/filterAF {input[0]} {AF} > Databases/KG1pc.txt;"	
		"./Codes/sampling Databases/KG1pc.txt {rdmN} {seed} > {output[0]};"
		"./Codes/region Databases/KG1pc.txt {input[1]} > {output[1]};"
		"thresholds=$(seq {minrec} {maxrec});"
		"for f in $thresholds;" 
		"do "
		"./Codes/region Databases/KG1pc.txt Databases/COSMICrecUnkn$f.txt Databases/COSMICrecConf$f.txt > Databases/data1KG.regioncosrecConf$f.txt;"
		"./Codes/region Databases/KG1pc.txt Databases/COSMICrecUnkn$f.txt > Databases/data1KG.regioncosrecUnkn$f.txt;"
		"./Codes/region Databases/KG1pc.txt Databases/COSMICrecTotal$f.txt > Databases/data1KG.regioncosrecTotal$f.txt;"
		"done;"

rule benign_subsets:
	input:
		"Databases/benign.txt","Databases/Clinvar.txt",['Databases/COSMICrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)]#map(lambda x: 'Databases/COSMICrec'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
	output:
		"Databases/benign.regionclin",['Databases/benign.regioncosrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)] #map(lambda x: 'Databases/data1KG.regioncos'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))#,,map(lambda x: 'Databases/data1KG.regioncoscont'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
	shell:
		"./Codes/region Databases/benign.txt {input[1]} > {output[0]};"
		"thresholds=$(seq {minrec} {maxrec});"
		"for f in $thresholds;" 
		"do "
		"./Codes/region Databases/benign.txt Databases/COSMICrecUnkn$f.txt Databases/COSMICrecConf$f.txt > Databases/benign.regioncosrecConf$f.txt;"
		"./Codes/region Databases/benign.txt Databases/COSMICrecUnkn$f.txt > Databases/benign.regioncosrecUnkn$f.txt;"
		"./Codes/region Databases/benign.txt Databases/COSMICrecTotal$f.txt > Databases/benign.regioncosrecTotal$f.txt;"
		"done;"


###############################################################
#Description
rule Format2annot:
	input:
		"Databases/Clinvar.txt","Databases/COSMIC.txt","Databases/nCDShgmd.bed"
	output:
		"Databases/Clinvar.annot","Databases/COSMIC.annot","Databases/HGMD.annot"
	shell:
		"cut -f 1,2 {input[0]} > {output[0]};"
		"cut -f 1,2 {input[1]} > {output[1]};"
		"cut -f 1,2 {input[2]} > {output[2]};"

rule Annotation:
	input:
		"Databases/{descFile}.annot"
	output:
		"Databases/annot{descFile}.txt"
	shell:
		"gawk -F'\t' -i inplace '{{$NF=$2 FS $NF;}}1' OFS='\t' {input};"
		"bedtools intersect -loj -a {input} -b Databases/annot_gencode.txt > {input}.tmp;"
		"cut -f 3-6 --complement {input}.tmp > {input};"
		"rm {input}.tmp;"
		"./Codes/uniquePos {input} {output};"
		"""if [ "$(echo {input})" == "Databases/COSMIC.annot" ]; then
		echo $'chr\tpos\trecConf\trecUnkn\trecTotal' > tmp
		cut -f 1,2,6,7,8 Databases/COSMIC.txt >> tmp
		cat {output} > tmp2
		awk 'FNR==NR{{a[$1,$2]=$3;for(i=4;i<=NF;++i) a[$1,$2] = a[$1,$2] "\t" $i;next}}{{ print $0 "\t" a[$1,$2]}}' OFS='\n' tmp2 tmp > {output}
		rm tmp tmp2;
		fi;"""

rule Descriptive:
	input:
		expand("Databases/annot{descFile}.txt",descFile=['Clinvar','COSMIC','HGMD'])
	output:
		"Results/FeaturesDist.jpeg"
	shell:
		"Rscript --vanilla Codes/descriptiveResults.R {minrec} {maxrec} {input};" 


###############################################################
#Analyze
rule Analyze:
	input:
		"Databases/Clinvar.txt","Databases/COSMIC.txt","Databases/data1KG.unmatched","Databases/benign.txt",
		"Databases/data1KG.regionclin",['Databases/data1KG.regioncosrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)], #map(lambda x: 'Databases/data1KG.regioncos'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
		"Databases/benign.regionclin",['Databases/benign.regioncosrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)]
	output:
		"benchmarkResults.log" #written by benchmark.R
	threads: 
		1000
	shell:
		"Rscript --vanilla Codes/benchmark.R {threads} {input};"  #{threads} will be the min between 1000 and the "--cores" option
		"Rscript --vanilla Codes/benchmarkPlots.R {minrec} {maxrec};"

###############################################################

rule pipEnd:
	input:
		"benchmarkResults.log","Results/FeaturesDist.jpeg"
	output:
		"pipeline.log"
	shell:
		"cp {input[0]} > {output};"

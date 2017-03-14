rule all:
	input:
		"benchmarkresults.log","results/FeaturesDist.jpeg"




#snakemake --dag | dot | display



###################################
#Variables initialization
###################################
#COSMIC subset according to their somatic status confidence
COSMICsubsets = ["Conf","Unkn","Total"]
#Recurrence threshold range
minrec="1"
maxrec="20"
#minimal allele frequency to consider a 1KG variant as non pathogenic
AF="0.01"
#number of random variant for the random control set
rdmN="100000"
#seed for random sampling
seed="123456"



###################################
#Import
###################################
rule CADD_DL:
	output:
		"scores/CADD_scores.tsv"
	shell:
		"wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz -O scores/CADD_scores.tsv.gz;"
		"zcat -c scores/CADD_scores.tsv.gz > {output};"
		"rm scores/CADD_scores.tsv.gz"

rule DANN_DL:
	output:
		"scores/DANN_scores.tsv"
	shell:
		"wget -c https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz -O scores/DANN_scores.tsv.gz;"
		"zcat -c  scores/DANN_scores.tsv.gz > {output};"
		"rm scores/DANN_scores.tsv.gz"

#source = https://github.com/HAShihab/fathmm-MKL
rule FATHMM_MKL_DL:
	output:
		"scores/FATH_scores.tsv"
	shell:
		"wget -c http://fathmm.biocompute.org.uk/databases/fathmm-MKL_Current_zerobased.tab.gz -O scores/FATH_scores.tab.gz;"
		"zcat -c scores/FATH_scores.tab.gz > {output}"
		"rm scores/FATH_scores.tab.gz"

#source = http://funseq2.gersteinlab.org/downloads (May 2016)
rule Funseq2_DL:
	output:
		"scores/Funseq2_scores.tsv"
	shell:
		"wget -c http://archive.gersteinlab.org/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz -O scores/Funseq2_scores.tsv.gz;"
		"zcat -c  scores/Funseq2_scores.tsv.gz > {output}"
		"rm scores/Funseq2_scores.tsv.gz"


rule SNPSOM_DL:
	output:
		"scores/SNPSOM.bed"
	shell:
		"wget http://biodev.cea.fr/98drivers/Genomewise_SNP_and_SOM_scores.zip -O scores/SNPSOM.zip;"
		"unzip scores/SNPSOM.zip -d scores;"
		"mv scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_CLL_RF_1Kb.bed scores/CLL.bed;"
		"mv scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_liver_cancer_RF_1Kb.bed scores/liver.bed;"
		"mv scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_lung_cancer_RF_1Kb.bed scores/lung.bed;"
		"mv scores/Genomewise_SNP_and_SOM_scores/Human_genome_som_melanoma_RF_1Kb.bed scores/melanoma.bed;"
		"tail -n+2 scores/Genomewise_SNP_and_SOM_scores/Human_genome_Uniform_SNP_score_nc.bed > scores/tmpsnp.bed;"
		"rm -r scores/Genomewise_SNP_and_SOM_scores;"
		"rm -r scores/__MACOSX;"
		"cut -f 4 scores/liver.bed | paste scores/CLL.bed - > scores/tmpsom.bed;"
		"cut -f 4 scores/lung.bed | paste scores/tmpsom.bed - > scores/tmpsnpsom.bed;"
		"cut -f 4 scores/melanoma.bed | paste scores/tmpsnpsom.bed -  | tail -n+2 > scores/tmpsom.bed;"
		"rm scores/tmpsnpsom.bed scores/CLL.bed scores/liver.bed scores/lung.bed scores/melanoma.bed;"
		"gawk -i inplace -F'\t' '{{$2-=1 FS $2;}}1' OFS='\t' scores/tmpsnp.bed;"
		"gawk -i inplace -F'\t' '{{$2-=1 FS $2;}}1' OFS='\t' scores/tmpsom.bed;"
		"./src/coordfusion scores/tmpsnp.bed > scores/snp.bed;"
		"./src/coordfusion scores/tmpsom.bed > scores/som.bed;"
		"rm scores/tmpsnp.bed scores/tmpsom.bed;"
		"sed -i 's/chr//g' scores/som.bed;"
		"echo $'Y\t0\t59373567\tNA\tNA\tNA\tNA' >> scores/som.bed;"
		"sed -i 's/chr//g' scores/snp.bed;"
		"./src/snpsomDM scores/snp.bed scores/som.bed scores/tmpSNPSOM.bed;"
		"rm scores/snp.bed scores/som.bed;"
		"bedtools complement -i scores/tmpSNPSOM.bed -g src/grch37p13 > scores/SNPSOMcomplement.bed;"
		"sed -i 's/$/\tNA\tNA\tNA\tNA\tNA/g' scores/SNPSOMcomplement.bed;"
		"sed -i s/X/23/g scores/tmpSNPSOM.bed;"
		"sed -i s/Y/24/g scores/tmpSNPSOM.bed;"
		"sed -i s/X/23/g scores/SNPSOMcomplement.bed;"
		"sed -i s/Y/24/g scores/SNPSOMcomplement.bed;"
		"cat scores/tmpSNPSOM.bed scores/SNPSOMcomplement.bed | sort -nk1,1 -nk2,2 > {output};"
#		"sed -i s/^23/X/g {output};"
#		"sed -i s/^24/Y/g {output};"
		"rm scores/tmpSNPSOM.bed scores/SNPSOMcomplement.bed;"

rule GWAVA_DL:
	output:
		"scores/GWAVA"
	shell:
		"wget -r -P {output} -np -nH --cut-dirs 5  ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/ ;"
#		"a=$(pip freeze | grep scikit-learn);"
#		"b=(${{a//=/ }});"
#		"b=${{b[1]}};"
#		"./updatePython $b;"

rule KG_DL:
	output:
		"databases/KG.vcf"
	shell:
		"wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -O databases/KG.vcf.gz;"
		"zcat databases/KG.vcf.gz > {output}"

rule ClinVar_DL:	
	output:
		"databases/Clinvar.vcf"
	shell:
		"wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20160531.vcf.gz -O databases/Clinvar.vcf.gz;"
		"zcat databases/Clinvar.vcf.gz > {output};"


rule GENCODE_annotations_DL:
	output:
		"databases/Gencode_annots.gtf","databases/annot_gencode.txt","databases/feat_gencode.txt"
	shell:
		"wget -c ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz -O {output[0]}.gz;"
		"zcat {output[0]}.gz > {output[0]};"
		"tail -n+6 {output[0]} > databases/tmp.txt;"
		"cut -f 1,4,5 databases/tmp.txt > databases/tmp2.txt;"
		"cut -f 3 databases/tmp.txt | paste databases/tmp2.txt - > databases/gencodepos.txt;"
		"rm databases/tmp2.txt;"
		"sed -i 's/chr//g' databases/gencodepos.txt;"
		"sed -i 's/X/23/g' databases/gencodepos.txt;"
		"sed -i 's/Y/24/g' databases/gencodepos.txt;"
		"sed -i 's/M/25/g' databases/gencodepos.txt;"
		"cut -f 9 databases/tmp.txt > databases/tmp2.txt;"
		"cut -d ';' -f 3,6 databases/tmp2.txt > databases/gencodeannots.txt;"
		"rm databases/tmp2.txt;"
		"sed -i 's/; /\t/g' databases/gencodeannots.txt;"
		"sed -i 's/gene_type //g' databases/gencodeannots.txt;"
		"sed -i 's/transcript_type //g' databases/gencodeannots.txt;"
		"""sed -i 's/"//g' databases/gencodeannots.txt;"""
		"paste databases/gencodepos.txt databases/gencodeannots.txt > databases/tmp.txt;"
		"sort -n -k1,1 -k2,2 -k3,3  databases/tmp.txt > {output[1]};"
		"cut -f 4 {output[1]} | sort -u > databases/tmp.txt;"
		"cut -f 5 {output[1]} | sort -u >> databases/tmp.txt;"
		"cut -f 6 {output[1]} | sort -u >> databases/tmp.txt;"
		"sed -i 's/ //g' databases/tmp.txt;"
		"sort -u databases/tmp.txt > {output[2]};"
		"rm databases/tmp.txt;"


rule Non_coding_GENCODE_positions:
	input:
		"databases/Gencode_annots.gtf" 
	output:
		"databases/NCS.bed" 
	shell:	
		"sed -i 's/^chr//g' {input};"	
		"egrep $'\tCDS|\tstart|\tstop' {input} | bedtools complement -i - -g src/grch37p13 > {output}"  #grch37p13 is a data with chr and their length of the last release (13) of grch37 + 1 to ensure to take all positions with 1-based format

rule cosmicDM:
	output:
		"databases/COSMIC.raw"#,"databases/COSMICconf.raw"
	shell:
		"zcat databases/CosmicNCV.tsv.gz > databases/CosmicNCV.tsv;"
		"cut -f 1,6,7,8,9,15 databases/CosmicNCV.tsv > databases/CosmicNCV.txt;"
		"""egrep $'\ty' databases/CosmicNCV.txt | awk -F $'\t'  '{{ print $2,$4,$5,$1,$3 }}' OFS=$'\t'  | sort -u > databases/cosTmp.txt;""" #sort -u to remove duplicated lines (some inconsistencies for the ID_NCV, pubmed_PMID, SNP status, zygosity, FATHMM-MKL score or group... -> if same sample/individual ID, position and mutated allele, keep only one)
		"./src/cosmicDM databases/cosTmp.txt databases/cosmicWGS.txt;"
		"rm databases/cosTmp.txt;"
		"sort -nk1,1 -nk2,2 -k3,4 databases/cosmicWGS.txt > {output};" #{output[0]}
		"rm databases/cosmicWGS.txt;"
 

rule KG_DM:
	input:
		"databases/KG.vcf"
	output:
		"databases/data1KG.raw"
	shell:
		"./src/rmUseless {input} > databases/KG.tmp;" #remove indels (but keep chrM & no filter for AF, see below)
		"./src/rmMultiDupl databases/KG.tmp > {output};"		#remove duplicata for multiallelic sites (keep only the AF indicated by the multiallelic line)
		"rm  databases/KG.tmp {input};"
		
rule ClinVar_DM:
	input:
		"databases/Clinvar.vcf"
	output:
		"databases/Clinvar.raw"
	shell:
		"a=$(egrep $'#' {input} | wc -l);"
		"b=$(bc <<< $a+1);"
		"tail -n+$b {input} > databases/clitmp.txt;"
		"egrep CLNSIG=5 databases/clitmp.txt | cut -f 1,2,4,5,8 > databases/certain.txt;"
		"sed -i 's/X/23/' databases/certain.txt;"
		"sed -i 's/Y/24/' databases/certain.txt;"
		"egrep -v $'^MT\t' databases/certain.txt | sort -nk1,1 -nk2,2 -k3,3 -k4,4 - > databases/clitmp.txt;"
		"./src/clinvarDM databases/clitmp.txt {output};"
		"sed -i 's/^23\t/X\t/g' {output};"
		"sed -i 's/^24\t/Y\t/g' {output};"
		"rm databases/certain.txt databases/clitmp.txt;"
		"unset a b;"

rule HGMD_DM:
	output:
		"databases/hgmd.bed"
	shell:
		"tail -n+2 databases/HGMD_Advanced_Substitutions.bed | cut -f 1,2 > databases/hgmd.tmp;"
		"sed -i 's/chr//g' databases/hgmd.tmp;"
		"sort -nk1,1 -nk2,2 databases/hgmd.tmp > {output};"
		"rm databases/hgmd.tmp;"
		"gawk -i inplace -F'\t' '{{$2=$2 FS $2;}}1' OFS='\t' {output};"



###############################################################
#Common DM
rule Remove_CDS_positions:
	input:
		"databases/{file}.raw","databases/NCS.bed"
	output:
		"databases/{file}.bed"
	shell:
		"gawk -i inplace -F'\t' '{{$2=$2 FS $2;}}1' OFS='\t' {input[0]};"
		"bedtools intersect -a {input[0]} -b {input[1]} > {output}"

rule Remove_HGMD_variants:
	input:
		"databases/{file}.bed","databases/hgmd.bed"
	output:
		"databases/{file}.clean"
	shell: 
		"bedtools subtract -a {input[0]} -b {input[1]} > {input[0]}.tmp;"
		"""if [ "$(echo {input[0]})" != "databases/Clinvar.bed" ]; then bedtools subtract -a {input[0]}.tmp -b databases/Clinvar.bed > {output}; else cp {input[0]}.tmp {output}; fi;"""
		"rm {input[0]}.tmp;"



###############################################################
#Scoring
rule GWAVA_scoring:
	input:
		"databases/{file}.clean"
	output:
		"databases/{file}.gw"
	shell:
		"a={input};"
		"b=(${{a//\// }});"
		"b=${{b[1]}};"
		"a=(${{b//./ }});"
		"a=${{a[0]}};"
		"./src/export_GWAVA {input} > scores/GWAVA/src/$a.tmp;"
		"cd scores/GWAVA/src;"
		"split -l 1000000 $a.tmp $a.x;" 	#split file because too many lines for RAM saving; ~ 2hours / subdata
		"""myFiles=$(ls | grep "$a.x");"""  # not a comment, neither string tips error, do not remove !
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
		"./src/getGWAVAsc {input} {input}.gw> {output};"
		"rm {input}.gw;"
		"sed -i -r 's/(\s+)?\S+//2' {output};"
		"unset a b myFiles f;"
#		"sed -i 's/X\t/23\t/g' {output};"
#		"sed -i 's/Y\t/24\t/g' {output};"

#larger is better, for the 2,281,292 lines of COSMIC, 4 days using files of 1000, 4-5 hours using files of 1000000 lines
#but be careful to the RAM !

rule MakeIndex:
	input:
		"scores/CADD_scores.tsv",
		"scores/DANN_scores.tsv",
		"scores/FATH_scores.tsv",
		"scores/Funseq2_scores.tsv"
	output:
		"scores/initscores","scores/index"
	shell:
		"./src/makeIndex"



#rule split1KG:
#	input:
#		"databases/data1KG.gw"
#	output: 
##		map(lambda x: 'databases/data1KG'+str(x)+'gw',range(0,len(filesToScore)+1))
#		dynamic("databases/data1KG{num}.gw")
#	shell:
#		"src/split {input} 1000000 .gw"
##	run:
##		call("./src/split","databases/data1KG.gw","1000000","gw")
##		filesToScore = [s for s in [f for f in os.listdir('.') if os.path.isfile(f)] if ".gw" in s]
##		filesToScore = [s for s in filesToScore if "data1KG" in s]
##		if "data1KG.gw" in filesToScore: filesToScore.remove("data1KG.gw")




rule Scoring:
	input:
		#"databases/data1KG{num}.gw",
		#["databases/Clinvar.gw","databases/COSMIC.gw","databases/data1KG{num}.gw"],
		"databases/{file}.gw",
		#"scores/CADD_scores.tsv",
		#"scores/DANN_scores.tsv",
		#"scores/FATH_scores.tsv",
		#"scores/Funseq2_scores.tsv",
		"scores/SNPSOM.bed",
		"scores/initscores",
		"scores/index"
	output:
		"databases/{file}.txt"
		#"databases/data1KG{num}.txt"
		#["databases/Clinvar.txt","databases/COSMIC.txt","databases/data1KG{num}.txt"]
		#["databases/Clinvar.txt","databases/COSMIC.txt",map(lambda x: 'databases/data1KG'+str(x)+'gw',range(0,len(filesToScore)+1))],
	shell:
		"src/variantScoring {input[0]} {input[0]}.tmp;"
		"egrep -v NA {input[0]}.tmp | sort -nk1,1 -nk2,2 -k4,4 > {input[0]}.temp;" #Remove_unscorable
		"rm {input[0]}.tmp;"
		"src/mergeSnpsom {input[0]}.temp {input[0]}.out;"
		"rm {input[0]}.temp;"
		"egrep -v NA {input[0]}.out | sort -nk1,1 -nk2,2 -k4,4 > {output};" #Remove_unscorable
		"rm {input[0]}.out;"

#rule cat1KG:
#	input:Conf
#		dynamic("databases/data1KG{num}.txt")		
#	output:
#		"databases/data1KG.txt"
#	shell:
#		"cat {input} >> {output}"

###############################################################
#finalization of the DM		

rule COSMIC_recurrence:
	input:
		"databases/COSMIC.txt"#,"databases/COSMICconf.txt"
	output:
		['databases/COSMICrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)] #map(lambda x: 'databases/COSMICrec'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))#,map(lambda x: 'databases/COSMICconfrec'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
	shell:
		"./src/filterRec {input} {minrec} {maxrec};"


# modify KG_subset for region of all COSMIC subset...
#for R code, keep only the large dataset (unnecessary to load each dataset...) -> COSMIC1.final, COSMIC2.final,..., can be removed after the KG_subset
rule KG_subsets:
	input:
		"databases/data1KG.txt","databases/Clinvar.txt",['databases/COSMICrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)]#map(lambda x: 'databases/COSMICrec'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
	output:
		"databases/data1KG.unmatched","databases/data1KG.regionclin",['databases/data1KG.regioncosrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)] #map(lambda x: 'databases/data1KG.regioncos'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))#,,map(lambda x: 'databases/data1KG.regioncoscont'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
	shell:
		"./src/filterAF {input[0]} {AF} > databases/KG1pc.txt;"	#AF cutoff (1pc = 1 percent) & specific disease (default = all)
		"./src/sampling databases/KG1pc.txt {rdmN} {seed} > {output[0]};"
		"./src/region databases/KG1pc.txt {input[1]} > {output[1]};"
		"thresholds=$(seq {minrec} {maxrec});"
		"for f in $thresholds;" 
		"do "
		"./src/region databases/KG1pc.txt databases/COSMICrecUnkn$f.txt databases/COSMICrecConf$f.txt > databases/data1KG.regioncosrecConf$f.txt;"
		"./src/region databases/KG1pc.txt databases/COSMICrecUnkn$f.txt > databases/data1KG.regioncosrecUnkn$f.txt;"
		"./src/region databases/KG1pc.txt databases/COSMICrecTotal$f.txt > databases/data1KG.regioncosrecTotal$f.txt;"
		"rm databases/COSMICrecConf$f.txt databases/COSMICrecUnkn$f.txt databases/COSMICrecTotal$f.txt;"
		"done;"
		#"rm databases/KG1pc.txt"


###############################################################
#Description
rule Format2annot:
	input:
		"databases/Clinvar.txt","databases/COSMIC.txt","databases/hgmd.bed"
	output:
		"databases/Clinvar.annot","databases/COSMIC.annot","databases/HGMD.annot"
	shell:
		"cut -f 1,2 {input[0]} > {output[0]};"
		"cut -f 1,2 {input[1]} > {output[1]};"
		"cut -f 1,2 {input[2]} > {output[2]};"

rule Annotation:
	input:
		"databases/{descFile}.annot"
	output:
		"databases/annot{descFile}.txt"
	shell:
		"gawk -F'\t' -i inplace '{{$NF=$2 FS $NF;}}1' OFS='\t' {input};"
		"bedtools intersect -loj -a {input} -b databases/annot_gencode.txt > {input}.tmp;"
		"cut -f 3-6 --complement {input}.tmp > {input};"
		"rm {input}.tmp;"
		"./src/uniquePos {input} {output};"
		#"rm {input};"
		"""if [ "$(echo {input})" == "databases/COSMIC.annot" ]; then
		echo $'recConf\trecUnkn\trecTotal' > tmp;
		cut -f 6,7,8 databases/COSMIC.txt >> tmp;
		gawk -F'\t' -i inplace 'FNR==NR{{a[NR]=$1;b[NR]=$2;c[NR]=$3;next}}{{$3=a[FNR] FS $3;$3=b[FNR] FS $3;$3=c[FNR] FS $3;}}1' OFS='\t' tmp {output};
		rm tmp;
		fi;"""

rule Descriptive:
	input:
		expand("databases/annot{descFile}.txt",descFile=['Clinvar','COSMIC','HGMD'])
		#"databases/annotClinvar.txt","databases/annotCOSMIC.txt","databases/annotHGMD.txt"
	output:
		"results/FeaturesDist.jpeg"
	shell:
		"Rscript --vanilla src/descriptiveresults.R {minrec} {maxrec} {input};" 


###############################################################
#Analyze
rule Analyze:
	input:
		"databases/Clinvar.txt","databases/COSMIC.txt","databases/data1KG.unmatched",
		"databases/data1KG.regionclin",['databases/data1KG.regioncosrec'+x+str(y)+'.txt' for x in COSMICsubsets for y in range(int(minrec),int(maxrec)+1)] #map(lambda x: 'databases/data1KG.regioncos'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
		#,map(lambda x: 'databases/data1KG.regionconfcos'+str(x)+'.txt',range(int(minrec),int(maxrec)+1))
		#dynamic("databases/data1KG{num}.txt")
		#["databases/Clinvar.txt","databases/COSMIC.txt",dynamic("databases/data1KG{num}.txt")]
	output:
		"benchmarkresults.log" #write by benchmark.R
	threads: 
		1000
	shell:
		"Rscript --vanilla src/benchmark.R {threads} {input};"  #{threads} will be the min between 1000 and the "--cores" option
		"Rscript --vanilla src/benchmarkPlots.R {minrec} {maxrec};"

###############################################################

#rule all:
#	input:
#		#"descriptiveresults.txt","benchmarkresults.txt"
#		"benchmarkresults.log","results/FeaturesDist.jpeg"
#	output:
#		"report.txt"
#	shell:
#		"bc <<< 2+2 > {output};"
		#"Rscript report.R"

#!/share/apps/R-3.6.1/bin/Rscript 

# script to produce ID/Phenotype file for PrevAD
# and APOE allele doses
# to use third wave of ADSP case-control data

library(stringr)

ToDropFile="/home/rejudcu/ADSP2024/downloaded/R4WGS_recommendation_list.csv"
PhenoFile="/home/rejudcu/ADSP2024/downloaded/ADSPCaseControlPhenotypes_DS_2022.08.18.v4_ALL.csv"
SampleManifestFile="/home/rejudcu/ADSP2024/downloaded/SampleManifest_DS_2022.08.18.v3_ALL.csv"
Phase1IDFile="/home/rejudcu/ADSP2024/PhenoFiles/origIDPrevAD.txt"
Phase2IDFile="/home/rejudcu/ADSP2024/PhenoFiles/newIDPrevAD.txt"

# Get list of Sample IDs from VCFFile
#  zcat /cluster/project9/bipolargenomes/ADSP/ADSP2024/vcf/gcad.qc.r2.wes.chr22.20504.GATK.2020.06.26.v2.biallelic.genotypes.ALL.vcf.gz | head -n 5000 | grep '#CHROM' | while IFS=`\t` read line; do for w in $line; do echo $w >> sampleIDs.txt ; done; done
SampleIDs=data.frame(read.table("sampleIDs.txt",header=FALSE,stringsAsFactors=FALSE))

# Read manifest and retain samples in VCFFile
SampleManifest=data.frame(read.csv(SampleManifestFile,header=TRUE,stringsAsFactors=FALSE))
print(nrow(SampleManifest))
colnames(SampleIDs)[1]="SampleID"
SampleManifest=merge(SampleManifest,SampleIDs,by="SampleID")
print(nrow(SampleManifest))

# Read phenotype file and add Sample IDs of retained samples
Phenos=data.frame(read.csv(PhenoFile,header=TRUE,stringsAsFactors=FALSE))
print(nrow(Phenos))
Phenos=merge(Phenos,SampleManifest,by="SUBJID")
print(nrow(Phenos))

# Read list of samples to be dropped and exlude them
ToDrop=data.frame(read.csv(ToDropFile,header=TRUE,stringsAsFactors=FALSE))
ToDrop=ToDrop[ToDrop$ADSP_Recommendation!="Keep",]
Phenos=Phenos[!(Phenos$SampleID %in% ToDrop$SampleID),]
print(nrow(Phenos))

# Use reported APOE genotype if not available from WGS and exclude samples with neither
Phenos$APOE=Phenos$APOE_WGS
Phenos$APOE[is.na(Phenos$APOE)]=Phenos$APOE_reported[is.na(Phenos$APOE)]
Phenos=Phenos[!is.na(Phenos$APOE),]
print(nrow(Phenos))

# Calculate doses for APOE 3 and 4 alleles
# https://stackoverflow.com/questions/12427385/how-to-calculate-the-number-of-occurrence-of-a-given-character-in-each-row-of-a
Phenos$DoseAPOE3=str_count(Phenos$APOE,"3")
Phenos$DoseAPOE4=str_count(Phenos$APOE,"4")
write.table(Phenos,"ADSP2024.phenos.txt",row.names=FALSE,quote=FALSE,sep="\t")

# Get samples for first cohort and second cohorts
FirstCohort=data.frame(read.table(Phase1IDFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))
# Subject IDs, not sample IDs
FirstCohort=merge(FirstCohort,SampleManifest,by="SUBJID")
SecondCohort=data.frame(read.table(Phase2IDFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))
OldSamples=c(FirstCohort$SampleID,SecondCohort$SampleID)

# Save phenotype files for PrevAD, APOE3Dose and APOE3Dose
PrevAD=Phenos[,c("SampleID","PrevAD")]
write.table(PrevAD[!(PrevAD$SampleID %in% OldSamples),],"PrevAD.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(PrevAD[PrevAD$SampleID %in% OldSamples,],"OldPrevAD.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(PrevAD,"AllPrevAD.txt",row.names=FALSE,quote=FALSE,sep="\t")

DoseAPOE3=Phenos[,c("SampleID","DoseAPOE3")]
write.table(DoseAPOE3,"DoseAPOE3.txt",row.names=FALSE,quote=FALSE,sep="\t")
DoseAPOE4=Phenos[,c("SampleID","DoseAPOE4")]
write.table(DoseAPOE4,"DoseAPOE4.txt",row.names=FALSE,quote=FALSE,sep="\t")





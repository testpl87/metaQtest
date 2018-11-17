
###### QTest Example ######
source("QTest_metaQTest_code.R")
anno0<-read.table("wave3_anno_score12.txt",header=T)[,c(2,3,5)]
colnames(anno0)<-c("CHR","gene","snp")
anno0$gene<-as.character(anno0$gene);anno0$snp<-as.character(anno0$snp);anno0$CHR<-as.numeric(as.matrix(anno0$CHR))

#Use plink files #
library(MultiPhen)
allgeno0<-read.plink("chr21_9960") ## read binary plink file
allgeno0[allgeno0==0]<-4;allgeno0[allgeno0==2]<-0;allgeno0[allgeno0==4]<-2;  ## appropriate 012 coding
bim<-read.table("chr21_9960.bim")
fam<-read.table("chr21_9960.fam")
allgeno0<-data.frame(ID=fam[,2],allgeno0)
colnames(allgeno0)<-c("ID",as.character(bim[,2]))
colnames(pheno)[2]<-"ID"
rslt<-QTest.all(phe.cova=pheno,allgeno=allgeno0,anno=anno0,yname="HEIGHT_ALL_RAW_RES_INV",covaname=NULL,cut.r2=0.05,STT=0.2,min.mac=5,maf.cut=0.05,outname="HEIGHT_ALL_RAW_RES_INV_weightF")

## Use genotype file per each gene
fam<-read.table("chr21_9960.fam")
pheno<-pheno[match(fam[,2],pheno[,2]),]
genename<-read.table("a.txt")[,1];genename1<-apply(as.matrix(genename),1,function(v)strsplit(v,split=".geno")[[1]])
setwd("/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat")
rslt<-QTest.bibs(phe.cova=pheno,yname="HDLexcmed_ALL_RAW_RES_INV",genename=genename,STT=0.2,min.mac=5,maf.cut=0.05,weight=TRUE,preprocess=FALSE,outname="HDLexcmed_ALL_RAW_RES_INV_weightT")
rslt.skat<-SKAT.bibs(phe.cova=pheno,yname="HDLexcmed_ALL_RAW_RES_INV",genename=genename,min.mac=5,maf.cut=0.05,weight="yes",outname="HDLexcmed_ALL_RAW_RES_INV")
GS.QTest(pheno,yname="LDLcalexcmed_ALL_RAW_RES_INV",set.anno=set.anno0,genename=genename1,genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",minSetsize=200,STT=0.2,maf.cut=0.05,min.mac=5,weight=FALSE,outname="_LDLcalexcmed_ALL_RAW_RES_INV_c2_weightF")->rslt

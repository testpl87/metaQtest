##~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#  R code for QTest & meta QTest            O O
#                                
#  by Jaehoon Lee   & Jieun Ka           #
#        (2013.04.03)          #
#                              #
     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(SKAT)
###################################
##    Anno file                  ##  
##colnames:  "CHR" "gene","snp"  ##
###################################
#anno<-read.table("...",header=T)

############################################################
##    Pheno & covariate file                              ##  
##colnames: "ID" "pheno1" "pheno2" ..."cova1" "cova2",..  ##
############################################################

#phe.cova<-read.table("...",header=T)
#colnames(phe.cova)[1]<-"ID"

########################################
##    Genotype file for all SNPs      ##  
##colnames:   "snp1" "snp2" ...       ##
########################################
#allgeno<-read.table("...",header=T)


#############################################################
## get.clumped.geno(geno,r2)                               ##
##                                                         ##
## Function for preprocessing genotype based on r^2 cutoff ##
#############################################################
get.clumped.geno<-function(origeno,cut.r2){
  if(ncol(origeno)>1){
       temp<-origeno
       cor.temp<-cor(temp,use="pairwise.complete.obs")
       cor.temp[lower.tri(cor.temp,diag=TRUE)==TRUE]<-0
       cor.temp[is.na(cor.temp)==TRUE]<-0
       idx.row<-which((cor.temp)^2>cut.r2)%%ncol(origeno)
       idx.col<-ceiling(which((cor.temp)^2>cut.r2)/ncol(origeno))
       if(length(idx.row)>1){
           corr.mat<-rbind(idx.row,idx.col)
           edge<-list();edge[[1]]<-as.numeric(corr.mat[,1])
           for(i in 2:length(idx.row)){
                 edg.idx<-which(lapply(edge,function(vec0)length(which(vec0%in%corr.mat[,i])==TRUE))>0)
                 if(length(edg.idx)==1){edge[[edg.idx]]<-unique(c(as.numeric(corr.mat[,i]),as.numeric(edge[[edg.idx]])))}
                 if(length(edg.idx)>1){temp1<-c(unlist(edge[edg.idx]),as.numeric(corr.mat[,i]))
		                       edge[[edg.idx[1]]]<-temp1;edge[edg.idx[2:length(edg.idx)]]<-NULL}
                 if(length(edg.idx)==0){edge[[length(edge)+1]]<-as.numeric(corr.mat[,i])}
           }
			edge<-lapply(edge,unique)
			collapse.z<-lapply(edge,function(vec1){
                            c.z<-apply(temp[,vec1],1,function(v)sum(na.omit(v)))/length(vec1)
							return(c.z)
			})
	   	   
			names(collapse.z)<-unlist(lapply(edge,function(vec2)paste("collapse",paste("_",vec2,collapse="",sep=""),sep="")))
			newgeno<-cbind(temp[,-unlist(edge)],data.frame(collapse.z));
			if(length(c(1:ncol(temp))[-unlist(edge)])==1){colnames(newgeno)[1]<-colnames(temp)[-unlist(edge)]}
			return(newgeno)
       }
       if(length(idx.row)==1){
           collapse.z<-list(apply(temp[,c(idx.row,idx.col)],1,sum))
           names(collapse.z)<-paste("collapse_",idx.row,"_",idx.col,sep="")
           newgeno<-cbind(temp[,-c(idx.row,idx.col)],data.frame(collapse.z))
	   if(length(c(1:ncol(temp))[-c(idx.row,idx.col)])==1){colnames(newgeno)[1]<-colnames(temp)[-c(idx.row,idx.col)]};return(newgeno)}
       
       if(length(idx.row)==0){return(temp)}
   }
   if(ncol(origeno)==1){return(origeno)}
 }

## Empirical null distribution for Q3 test
null.dist.Q3<-read.table("n.log10.minp.1e09.txt",header=T)
totlength<-sum(null.dist.Q3[,2])
get.finalp<-function(pval){(sum(null.dist.Q3[null.dist.Q3[,1]>-log10(pval),2])+1)/(totlength+1)}


## Get a given STT
get.stt<-function(L,a,STT){
	1-pgamma(L*qgamma(1-STT,a,1),L*a,1)
}


get.a<-function(L,STT){
	aa<-diff<-seq(0,1,length=200)
	for(i in 1:length(aa)){
	diff[i]<-abs(get.stt(L,aa[i],STT)-STT)
	}
	return(aa[which.min(diff)])
}

get.imputed.S<-function(v){
	na.v<-which(is.na(v)==TRUE)
	if(length(na.v)>0){v[na.v]<-mean(na.omit(v))}
	return(v)
}


###############################################################################################
## QTest.one(phe.cova,geno,yname, r2)  #QTest for single gene                                ##
##                                                                                           ##
## y: trait, newgeno: preprocessed genotype matrix                                           ##
## cut.r2: r^2 cutoff, n.perm: # of pumutations                                              ##
## n.perm: number of permutation for GM method                                               ##
## a: shape parameter for GM method                                                          ##
## use.GM: whether Gamma method is used                                                      ##
## EV.perm=c(E.GM.perm,V.GM.perm,E.IG.perm,V.IG.perm): estimated parameters from permutation ##
###############################################################################################

library(CompQuadForm)
QTest.one<-function(y,covadat=NULL,newgeno,STT,weight=FALSE,weight.type="beta"){
	if(length(covadat)!=0){resid<-try(resid(glm(y~.,data=data.frame(covadat),na.action=na.exclude)),TRUE)}
	if(length(covadat)==0){resid<-try(resid(glm(y~1,na.action=na.exclude)),TRUE)}
	S<-try(as.matrix(newgeno),TRUE)
#	S<-apply(S,2,get.imputed.S)
	fit<-try(glm(resid~.,data=data.frame(S)),TRUE)
	na.S<-try(which(is.na(coef(fit)[-1])==TRUE),TRUE)
	if(length(na.S)>0){
		S<-try(as.matrix(S[,-na.S]),TRUE)
		fit<-try(glm(resid~.,data=data.frame(S)),TRUE)
	}
	coef<-try(coef(summary(fit))[-1,1:2],TRUE)
	
	
	if(mode(fit)=="character"){length(coef)<-0}
	if(length(coef)!=0){
		if(length(coef)!=2){beta1<-coef[,1];se1<-coef[,2]}
		if(length(coef)==2){beta1<-coef[1];se1<-coef[2]}
		SS<-cbind(1,S); n<-length(na.omit(y));maf.S<-apply(S,2,function(v)mean(na.omit(v))/2)
		if(weight.type=="inv.sd"){w.S0<-1/sqrt(maf.S*(1-maf.S))}
		if(weight.type=="beta"){w.S0<-qbeta(maf.S,1,25,lower.tail=F)}
		WS<-diag(w.S0);if(length(beta1)==1){WS<-w.S0}
		vv<-vcov(fit)[-1,-1];alpha<-(1/(se1^2))
		#vv<-solve(t(SS)%*%SS)[-1,-1]*var(na.omit(y))*(n-1)/n;alpha<-(1/(se1^2))
		
		##Q1##
		if(weight==FALSE){
			var.pool0<-t(alpha)%*%vv%*%alpha
			beta.pool0<-t(alpha)%*%beta1
			z.score0<-beta.pool0/sqrt(var.pool0)
			Q1<-z.score0^2
			p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
		}
		if(weight==TRUE){		
			var.pool<-t(alpha)%*%WS%*%vv%*%WS%*%alpha
			beta.pool<-t(alpha)%*%WS%*%beta1
			z.score<-beta.pool/sqrt(var.pool)
			Q1<-z.score^2
			p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
		}
		##Q2##
		#y.new<-y-mean(y)
		#Q2<-t(beta1)%*%solve(vv)%*%beta1 ## =t(y.new)%*%S%*%solve(t(SS)%*%SS)[-1,-1]%*%t(S)%*%y.new/var(y)*n/(n-1)
		#p.Q2<-pchisq(Q2,df=ncol(S),lower.tail=FALSE)
		
		Q2.eigen<-eigen(vv);U2<-Q2.eigen$vectors; l2<-Q2.eigen$values;na.l2<-which(l2/mean(l2)<0.01)
		if(length(na.l2)>0){l2<-l2[-na.l2];U2<-U2[,-na.l2]; maf.S<-maf.S[-na.l2]}
		p2<-pchisq((t(U2)%*%beta1)^2/l2,df=1,lower.tail=FALSE)
		a<-get.a(length(p2),STT)
		q2<-2*(qgamma(p2,a,1,lower.tail=FALSE))
		Q2<-sum(q2)
		p.Q2<-pchisq(Q2,df=2*a*length(q2),lower.tail=FALSE)
		if(length(beta1)==1){p.Q2<-p.Q1}



		##Q3 ##

		if(weight==FALSE){
			c(t(alpha)%*%vv)->cov.beta
			b.star<-beta1-beta.pool0*cov.beta/var.pool0[1]
			vv.star<-vv-(cov.beta%*%t(cov.beta))/var.pool0[1]
		}
		if(weight==TRUE){
			w.vv<-WS%*%vv%*%WS
			c(t(alpha)%*%w.vv)->cov.beta
			b.star<-WS%*%beta1-beta.pool*cov.beta/var.pool[1]
			vv.star<-w.vv-(cov.beta%*%t(cov.beta))/var.pool[1]
		}

		if(length(beta1)!=1){
			Q3.eigen<-eigen(vv.star); U3<-Q3.eigen$vectors; l3<-Q3.eigen$values;na.l3<-which(l3/mean(l3)<0.001)
			if(length(na.l3)>0){l3<-l3[-na.l3];U3<-U3[,-na.l3]}
			L3<-try(diag(l3),TRUE);if(length(l3)==1){L3<-l3}
			q2.proj<-(t(U3)%*%b.star)^2/l3; 
			p2.1<-pchisq(q2.proj,df=1,lower.tail=FALSE)
			a<-get.a(length(p2.1),STT)
			Q2.proj<-sum(2*qgamma(p2.1,a,1,lower.tail=FALSE))
			p.Q2.proj<-pchisq(Q2.proj,df=2*a*length(l3),lower.tail=FALSE)
			Q2.1<-qchisq(p.Q2.proj,df=1,lower.tail=FALSE)
			
			pi0<-seq(0,1,by=0.1);p.Q3.can<-rep(1,11)
			p.Q3.can[1]<-p.Q2.proj
			p.Q3.can[11]<-p.Q1
			for(h in 2:10){
				Q3<-pi0[h]*Q1+(1-pi0[h])*Q2.1
				p.Q3.can[h]<-davies(Q3,c(pi0[h],(1-pi0[h])),c(1,1))$Qq
				if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-imhof(Q3,c(pi0[h],(1-pi0[h])),c(1,1))$Qq}
				if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-liu(Q3,c(pi0[h],(1-pi0[h])),c(1,1))[1]}
			}
			
			
			p.Q3<-get.finalp(min(p.Q3.can))
		}

		if(length(beta1)==1){p.Q3<-p.Q1}
		#size<-ncol(as.matrix(newgeno))
		#burden.maf<-length(which(apply(as.matrix(S),1,function(g)sum(na.omit(g)))>0))/(2*nrow(S))
		#rslt<-data.frame(size,burden.maf,p.Q1,p.Q2,p.Q3)
		rslt<-data.frame(p.Q1,p.Q2,p.Q3)
		

		
	}

	if(length(coef)==0){rslt<-NA}
	options (warn=-1)
	return(rslt)	
}

QTest.one.bin<-function(y,covadat=NULL,newgeno,STT,weight=FALSE,weight.type="beta"){
	S<-try(as.matrix(newgeno),TRUE)
	fit<-try(glm(y~.,family=binomial,data=data.frame(S)),TRUE)
	na.S<-try(which(is.na(coef(fit)[-1])==TRUE),TRUE)
	if(length(na.S)>0){
		S<-try(as.matrix(S[,-na.S]),TRUE)
		fit<-try(glm(y~.,family=binomial,data=data.frame(S)),TRUE)
	}
	coef<-try(coef(summary(fit))[-1,1:2],TRUE)
	
	
	if(mode(fit)=="character"){length(coef)<-0}
	if(length(coef)!=0){
		if(length(coef)!=2){beta1<-coef[,1];se1<-coef[,2]}
		if(length(coef)==2){beta1<-coef[1];se1<-coef[2]}
		SS<-cbind(1,S); n<-length(na.omit(y));maf.S<-apply(S,2,function(v)mean(na.omit(v))/2)
		if(weight.type=="inv.sd"){w.S0<-1/sqrt(maf.S*(1-maf.S))}
		if(weight.type=="beta"){w.S0<-qbeta(maf.S,1,25,lower.tail=F)}
		WS<-diag(w.S0);if(length(beta1)==1){WS<-w.S0}
		vv<-vcov(fit)[-1,-1];alpha<-(1/(se1^2))
		#vv<-solve(t(SS)%*%SS)[-1,-1]*var(na.omit(y))*(n-1)/n;alpha<-(1/(se1^2))
		
		##Q1##
		if(weight==FALSE){
			var.pool0<-t(alpha)%*%vv%*%alpha
			beta.pool0<-t(alpha)%*%beta1
			z.score0<-beta.pool0/sqrt(var.pool0)
			Q1<-z.score0^2
			p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
		}
		if(weight==TRUE){		
			var.pool<-t(alpha)%*%WS%*%vv%*%WS%*%alpha
			beta.pool<-t(alpha)%*%WS%*%beta1
			z.score<-beta.pool/sqrt(var.pool)
			Q1<-z.score^2
			p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
		}
		##Q2##
		#y.new<-y-mean(y)
		#Q2<-t(beta1)%*%solve(vv)%*%beta1 ## =t(y.new)%*%S%*%solve(t(SS)%*%SS)[-1,-1]%*%t(S)%*%y.new/var(y)*n/(n-1)
		#p.Q2<-pchisq(Q2,df=ncol(S),lower.tail=FALSE)
		
		Q2.eigen<-eigen(vv);U2<-Q2.eigen$vectors; l2<-Q2.eigen$values;na.l2<-which(l2/mean(l2)<0.01)
		if(length(na.l2)>0){l2<-l2[-na.l2];U2<-U2[,-na.l2]; maf.S<-maf.S[-na.l2]}
		p2<-pchisq((t(U2)%*%beta1)^2/l2,df=1,lower.tail=FALSE)
		a<-get.a(length(p2),STT)
		q2<-2*(qgamma(p2,a,1,lower.tail=FALSE))
		Q2<-sum(q2)
		p.Q2<-pchisq(Q2,df=2*a*length(q2),lower.tail=FALSE)
		if(length(beta1)==1){p.Q2<-p.Q1}



		##Q3 ##

		if(weight==FALSE){
			c(t(alpha)%*%vv)->cov.beta
			b.star<-beta1-beta.pool0*cov.beta/var.pool0[1]
			vv.star<-vv-(cov.beta%*%t(cov.beta))/var.pool0[1]
		}
		if(weight==TRUE){
			w.vv<-WS%*%vv%*%WS
			c(t(alpha)%*%w.vv)->cov.beta
			b.star<-WS%*%beta1-beta.pool*cov.beta/var.pool[1]
			vv.star<-w.vv-(cov.beta%*%t(cov.beta))/var.pool[1]
		}

		if(length(beta1)!=1){
			Q3.eigen<-eigen(vv.star); U3<-Q3.eigen$vectors; l3<-Q3.eigen$values;na.l3<-which(l3/mean(l3)<0.001)
			if(length(na.l3)>0){l3<-l3[-na.l3];U3<-U3[,-na.l3]}
			L3<-try(diag(l3),TRUE);if(length(l3)==1){L3<-l3}
			q2.proj<-(t(U3)%*%b.star)^2/l3; 
			p2.1<-pchisq(q2.proj,df=1,lower.tail=FALSE)
			a<-get.a(length(p2.1),STT)
			Q2.proj<-sum(2*qgamma(p2.1,a,1,lower.tail=FALSE))
			p.Q2.proj<-pchisq(Q2.proj,df=2*a*length(l3),lower.tail=FALSE)
			Q2.1<-qchisq(p.Q2.proj,df=1,lower.tail=FALSE)
			
			pi0<-seq(0,1,by=0.1);p.Q3.can<-rep(1,11)
			p.Q3.can[1]<-p.Q2.proj
			p.Q3.can[11]<-p.Q1
			for(h in 2:10){
				Q3<-pi0[h]*Q1+(1-pi0[h])*Q2.1
				p.Q3.can[h]<-davies(Q3,c(pi0[h],(1-pi0[h])),c(1,1))$Qq
				if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-imhof(Q3,c(pi0[h],(1-pi0[h])),c(1,1))$Qq}
				if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-liu(Q3,c(pi0[h],(1-pi0[h])),c(1,1))[1]}
			}
			
			
			p.Q3<-get.finalp(min(p.Q3.can))
		}

		if(length(beta1)==1){p.Q3<-p.Q1}
		#size<-ncol(as.matrix(newgeno))
		#burden.maf<-length(which(apply(as.matrix(S),1,function(g)sum(na.omit(g)))>0))/(2*nrow(S))
		#rslt<-data.frame(size,burden.maf,p.Q1,p.Q2,p.Q3)
		rslt<-data.frame(p.Q1,p.Q2,p.Q3)
		

		
	}

	if(length(coef)==0){rslt<-NA}
	options (warn=-1)
	return(rslt)	
}


###############################################################################################
## QTest.preprocess(allgeno,anno)  ## preprocess step                                        ##
###############################################################################################

QTest.preprocess<-function(allgeno,anno){
	geno.set<-tapply(anno$snp,anno$gene,function(v)allgeno[,colnames(allgeno)%in%v==TRUE])	
	geno.set1<-geno.set
	for(k in 1:length(geno.set1)){geno.set1[[k]]<-get.clumped.geno(as.matrix(geno.set[[k]]),0.05)}
	#;cat("preprocessing step for gene ",k,"/",length(geno.set1),"\n");flush.console();}
	options (warn=-1)
	return(geno.set1)
}




###################################################################################################
## QTest.all(phe.cova,geno,yname,covaname, r2)                                                   ##
##                                                                                               ##
###################################################################################################
QTest.all<-function(phe.cova,allgeno,anno,yname,covaname=NULL,STT=0.2,cut.r2=0.05,min.mac=5,weight=FALSE,maf.cut=0.05,preprocess=FALSE,outname){
	#preprocessing
	cat("spliting genotype data for each gene..","\n")
#	allgeno<-allgeno[match(phe.cova$ID,allgeno$ID),]
#	allgeno<-allgeno[,-1]
	where.y<-which(colnames(phe.cova)==yname)
	where.cova<-which(colnames(phe.cova)%in%covaname==TRUE)
	y<-phe.cova[,where.y]
	covadat<-phe.cova[,where.cova]
	na.y<-which(is.na(y)==TRUE)	
	if(length(na.y)>0){
		y<-y[-na.y]
		if(length(covaname)!=0){covadat<-covadat[-na.y,]}
		allgeno<-allgeno[-na.y,]
	}
	allgeno<-allgeno[,apply(allgeno,2,function(v)length(which(v>0)))>=min.mac]
	if(length(covaname)==0){covadat<-NULL}
	frq<-apply(allgeno,2,get.maf)
	allgeno<-allgeno[,frq<maf.cut]
	anno<-anno[anno$snp%in%colnames(allgeno)==TRUE,]
	anno$gene<-as.character(anno$gene);anno$snp<-as.character(anno$snp);anno$CHR<-as.character(anno$CHR)
	na.gene<-names(table(anno$gene)[table(anno$gene)==1])
	anno<-anno[anno$gene%in%na.gene==FALSE,]
	geno.set1<-tapply(anno$snp,anno$gene,function(v)allgeno[,colnames(allgeno)%in%v==TRUE])	
	if(preprocess==TRUE){
		cat("preprocessing step for gene..","\n")
		for(k in 1:length(geno.set1)){geno.set1[[k]]<-get.clumped.geno(as.matrix(geno.set1[[k]]),0.05)
		}
	}

	rslt<-list()
	for(i in 1:length(geno.set1)){
		rslt[[i]]<-QTest.one(y,newgeno=geno.set1[[i]],covadat,STT,weight)
		cat("QTest for gene ",i,"/",length(geno.set1),"\n");flush.console();
	}
	rslt1<-do.call(rbind,rslt)

	maf.set<-lapply(geno.set1,function(v){v1<-na.omit(apply(v,1,sum));length(which(v1>0))/(2*length(v1))})
	gene.size<-lapply(geno.set1,ncol)
	genename<-names(geno.set1)
	rslt.dat<-data.frame(gene=genename,CHR=anno$CHR[match(genename,anno$gene)],size=unlist(gene.size),burden.MAF=unlist(maf.set),rslt1)

	options (warn=-1)
	write.table(rslt.dat,paste(outname,"_QTest.rslt",sep=""),quote=FALSE,row.names=F)
	return(rslt.dat)
}

get.maf <-function(vec){
	(length(which(vec==1))+2*length(which(vec==2)))/(2*length(vec))
}


SKAT.one<-function(y,covadat,geno,weight){
	if(length(covadat)==0){obj.jh<-SKAT_Null_Model(y~1,out_type="C")}
	if(length(covadat)!=0){obj.jh<-SKAT_Null_Model(y~.,data=data.frame(covadat),out_type="C")}

#	size<-ncol(as.matrix(geno))
#	burden.maf<-length(which(apply(as.matrix(geno),1,function(g)sum(na.omit(g)))>0))/(2*nrow(geno))

	if(weight=="yes"){
		p.skat<-try(SKAT(as.matrix(geno), obj.jh,weights.beta=c(1,25))$p.value,TRUE)		
		p.skat.o<-try(SKAT(as.matrix(geno), obj.jh,method="optimal",weights.beta=c(1,25))$p.value,TRUE)
		return(data.frame(p.skat,p.skat.o))
	}
	if(weight=="no"){
		p.skat.u<-try(SKAT(as.matrix(geno), obj.jh,weights.beta=c(1,1))$p.value,TRUE)
		p.skat.o.u<-try(SKAT(as.matrix(geno), obj.jh,method="optimal",weights.beta=c(1,1))$p.value,TRUE)
		return(data.frame(p.skat.u,p.skat.o.u))
	}
	if(weight=="both"){
		p.skat<-try(SKAT(as.matrix(geno), obj.jh,weights.beta=c(1,25))$p.value,TRUE)		
		p.skat.o<-try(SKAT(as.matrix(geno), obj.jh,method="optimal",weights.beta=c(1,25))$p.value,TRUE)
		p.skat.u<-try(SKAT(as.matrix(geno), obj.jh,weights.beta=c(1,1))$p.value,TRUE)
		p.skat.o.u<-try(SKAT(as.matrix(geno), obj.jh,method="optimal",weights.beta=c(1,1))$p.value,TRUE)
		return(data.frame(p.skat,p.skat.u,p.skat.o,p.skat.o.u))
	}


	
#	return(c(p.skat))
}


SKAT.all<-function(phe.cova,allgeno,anno,yname,covaname=NULL,min.mac=5,maf.cut=0.05,weight="yes",outname){
#	allgeno<-allgeno[match(phe.cova$ID,allgeno$ID),]
#	allgeno<-allgeno[,-1]	


	where.y<-which(colnames(phe.cova)==yname)
	where.cova<-which(colnames(phe.cova)%in%covaname==TRUE)
	y<-phe.cova[,where.y]
	covadat<-phe.cova[,where.cova]
	if(length(covaname)==0){covadat<-NULL}
	allgeno<-allgeno[,apply(allgeno,2,function(v)length(which(v>0)))>=min.mac]
	frq<-apply(allgeno,2,get.maf)
	allgeno<-allgeno[,frq<maf.cut]


	anno<-anno[anno$snp%in%colnames(allgeno)==TRUE,]
	anno$gene<-as.character(anno$gene);anno$snp<-as.character(anno$snp)
	na.gene<-names(table(anno$gene)[table(anno$gene)==1])
	anno<-anno[anno$gene%in%na.gene==FALSE,]
	
	
	
	geno.set<-tapply(anno$snp,anno$gene,function(v)allgeno[,colnames(allgeno)%in%v==TRUE])	
	rslt<-list()
	for(i in 1:length(geno.set)){
		rslt[[i]]<-SKAT.one(y,covadat,geno.set[[i]],weight=weight)
		cat("SKAT for gene ",i,"/",length(geno.set),"\n");flush.console();
	}
	rslt1<-do.call(rbind,rslt)
	
	maf.set<-lapply(geno.set,function(v){v1<-na.omit(apply(v,1,sum));length(which(v1>0))/(2*length(v1))})
	gene.size<-lapply(geno.set,ncol)
	genename<-names(geno.set)

	rslt.dat<-data.frame(gene=genename,CHR=anno$CHR[match(genename,anno$gene)],size=unlist(gene.size),burden.MAF=unlist(maf.set),rslt1)
	options(warn=-1)
	write.table(rslt.dat,paste(outname,"_SKAT",".rslt",sep=""),quote=F,row.names=F)
	return(rslt.dat)
}

QTest.bibs<-function(phe.cova=pheno,yname,genename,STT=0.2,min.mac=5,weight=TRUE,weight.type="beta",maf.cut=0.05,preprocess=TRUE,outname){
	where.y<-which(colnames(phe.cova)==yname)
	y<-phe.cova[,where.y]
	na.y<-which(is.na(y)==TRUE)	
	if(length(na.y)>0){
		y<-y[-na.y]
	}
	rslt<-list()
	maf.set<-size<-rep(10,length(genename))
	for(i in 1:length(genename)){
		geno<-read.table(as.character(genename[i]),header=T)
		geno<-as.matrix(geno[,apply(geno,2,function(v)length(which(v>0)))>=min.mac])
		frq<-apply(geno,2,get.maf)
		if(ncol(geno)>1){geno<-as.matrix(geno[,frq<maf.cut])}
		if(ncol(geno)>1){
			if(length(na.y)>0){geno<-geno[-na.y,]}
			if(preprocess==TRUE){geno<-get.clumped.geno(geno,0.05)}
			v1<-na.omit(apply(geno,1,sum))
			maf.set[i]<-length(which(v1>0))/(2*length(v1))
			size[i]<-ncol(geno)
			rslt[[i]]<-QTest.one(y,newgeno=geno,covadat=NULL,STT=STT,weight=weight,weight.type=weight.type)
		}
		if(ncol(geno)<2){rslt[[i]]<-maf.set[i]<-size[i]<-NA}
		
		cat("QTest for gene ",i,"/",length(genename),"\n");flush.console();
	}
	genename1<-tapply(as.character(genename),1:length(genename),function(v)strsplit(v,split=".gen")[[1]][1])
	rslt1<-do.call(rbind,rslt)
	rslt.dat<-data.frame(gene=genename1,size=size,burden.MAF=maf.set,rslt1)
	options (warn=-1)
	write.table(rslt.dat,paste(outname,"_QTest",".rslt",sep=""),quote=FALSE,row.names=F)
	return(rslt.dat)
}








SKAT.bibs<-function(phe.cova=pheno,yname,genename,min.mac=5,maf.cut=0.05,weight="yes",outname){
	where.y<-which(colnames(phe.cova)==yname)
	y<-phe.cova[,where.y]
	na.y<-which(is.na(y)==TRUE)	
	if(length(na.y)>0){
		y<-y[-na.y]
	}
	rslt<-list()
	maf.set<-size<-rep(10,length(genename))
	for(i in 1:length(genename)){
		geno<-as.matrix(read.table(as.character(genename[i]),header=T))
		geno<-as.matrix(geno[,apply(geno,2,function(v)length(which(v>0)))>=min.mac])
		frq<-apply(geno,2,get.maf)
		if(ncol(geno)>1){geno<-as.matrix(geno[,frq<maf.cut])}
		if(ncol(as.matrix(geno))>1){
			if(length(na.y)>0){geno<-geno[-na.y,]}
			v1<-na.omit(apply(geno,1,sum))
			maf.set[i]<-length(which(v1>0))/(2*length(v1))
			size[i]<-ncol(geno)
			rslt[[i]]<-SKAT.one(y,covadat=NULL,geno=geno,weight)
		}
		if(ncol(as.matrix(geno))<2){rslt[[i]]<-maf.set[i]<-size[i]<-NA}
		
		cat("SKAT for gene ",i,"/",length(genename),"\n");flush.console();
	}
	genename1<-tapply(as.character(genename),1:length(genename),function(v)strsplit(v,split=".gen")[[1]][1])
	rslt1<-do.call(rbind,rslt)
	rslt.dat<-data.frame(gene=genename1,size=size,burden.MAF=maf.set,rslt1)
	options (warn=-1)
	write.table(rslt.dat,paste(outname,"_SKAT",".rslt",sep=""),quote=FALSE,row.names=F)
	return(rslt.dat)

}

GS.QTest<-function(pheno,yname,set.anno,anno0,genename,genefolder,minSetsize=200,STT=0.2,maf.cut=0.05,min.mac=5,weight=FALSE,outname){
	rslt.set<-list()
	cat("Preprocessing...","\n");flush.console();
	apply(set.anno[,-1],1,function(v)as.character(unlist(as.matrix(v))))->set.list
	set.list[set.list==""]<-NA
	set.list<-apply(set.list,2,function(v)as.character(na.omit(v)))
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%genename==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<minSetsize&geneleng1>0)
	y=pheno[,colnames(pheno)%in%yname==TRUE]
	setwd(paste("/",genefolder,"/",sep=""))
	for(i in idx.accept){
		conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)
		setgeno<-apply(as.matrix(set.list1[[conv.idx]]),1,function(v)read.table(paste(v,".geno",sep=""),header=T))
		setgeno1<-do.call(cbind,setgeno)
		snp.maf<-apply(setgeno1,2,get.maf)
		snp.mac<-apply(setgeno1,2,function(v)length(which(v>0)))
		filter.idx<-which(snp.maf<maf.cut&snp.mac>=min.mac)
		unique.snpname<-unique(names(filter.idx));unique.snpname<-unique.snpname[unique.snpname!="x"]
		if(length(unique.snpname)>1){
			setgeno2<-as.matrix(setgeno1[,match(unique.snpname,colnames(setgeno1))])		
			setgeno3<-get.clumped.geno(setgeno2,cut.r2=0.05)
			new.anno<-data.frame(snp=colnames(setgeno3))
			new.anno$gene<-as.character(anno0$gene[match(new.anno$snp,anno0$snp)]);new.anno$snp<-as.character(new.anno$snp)
			naidx<-which(is.na(new.anno$gene)==TRUE)
			new.anno0<-new.anno[-naidx,]
			setgeno.s<-tapply(new.anno0$snp,new.anno0$gene,function(v)setgeno3[,colnames(setgeno3)%in%v==TRUE])
			if(length(naidx)>0){
				for(j in 1:length(naidx)){
					collap.idx<-as.numeric(strsplit(as.character(new.anno[naidx[j],]$snp),split="_")[[1]][-1])
					tag.snp<-colnames(setgeno2)[collap.idx[which.min(apply(setgeno2[,collap.idx],2,function(v)as.numeric(try(coef(summary(lm(y~v)))[2,4],TRUE))))]]
					if(length(tag.snp)>0){which.tag<-which(names(setgeno.s)==anno0$gene[which(anno0$snp==tag.snp)[1]])	
						if(length(which.tag)>0){setgeno.s[[which.tag]]<-cbind(setgeno.s[[which.tag]],setgeno3[,naidx[j]])}
						if(length(which.tag)==0){setgeno.s[[length(setgeno.s)+1]]<-setgeno3[,naidx[j]];names(setgeno.s)[length(setgeno.s)]<-anno0$gene[anno0$snp==tag.snp]}
					}
				}
			}
			if(length(naidx)==0){setgeno.s<-tapply(new.anno$snp,new.anno$gene,function(v)setgeno3[,colnames(setgeno3)%in%v==TRUE])}
			rslt<-do.call(rbind,lapply(setgeno.s,function(v)QTest.one(y=y,newgeno=as.matrix(v),STT=STT,weight=weight)))
			rslt.skat<-do.call(rbind,lapply(setgeno.s,function(v)SKAT.one(y=y,geno=as.matrix(v),covadat=NULL,weight="no")))
			
			a<-get.a(length(setgeno.s),STT)
			q1<-2*(qgamma(rslt[,1],a,1,lower.tail=FALSE))
			q2<-2*(qgamma(rslt[,2],a,1,lower.tail=FALSE))
			q3<-2*(qgamma(rslt[,3],a,1,lower.tail=FALSE))
			s1<-2*(qgamma(rslt.skat[,1],a,1,lower.tail=FALSE))
			s2<-2*(qgamma(rslt.skat[,2],a,1,lower.tail=FALSE))
		
			Q1<-sum(q1);Q2<-sum(q2);Q3<-sum(q3);S1<-sum(s1);S2<-sum(s2)
			p.Q1<-pchisq(Q1,df=2*a*length(q1),lower.tail=FALSE)
			p.Q2<-pchisq(Q2,df=2*a*length(q1),lower.tail=FALSE)
			p.Q3<-pchisq(Q3,df=2*a*length(q1),lower.tail=FALSE)
			p.S1<-pchisq(S1,df=2*a*length(q1),lower.tail=FALSE)
			p.S2<-pchisq(S2,df=2*a*length(q1),lower.tail=FALSE)
			rslt.set[[i]]<-data.frame(set=names(set.list1[conv.idx]),n.gene=geneleng1[i],n.SNP=ncol(setgeno2),p.Q1,p.Q2,p.Q3,p.S1,p.S2)
			}
		cat("Test for set",i,"\n");flush.console()
	}
	rslt.set1<-do.call(rbind,rslt.set)
	write.table(rslt.set1,paste(outname,"_GS.QTest",".rslt",sep=""),quote=F,row.names=F)
	return(rslt.set1)
}

SKAT.geneset.OneStepTest<-function(pheno,yname,set.anno,genename,minSetsize=100,maf.cut=0.05,min.mac=5,outname){
	rslt.skat<-list()
	cat("Preprocessing...","\n");flush.console();
	apply(set.anno[,-1],1,function(v)na.omit(as.character(unlist(v))))->set.list
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%genename==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<100&geneleng1>0)
	y=pheno[,colnames(pheno)%in%yname==TRUE]
	for(i in idx.accept){
		conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)
		setgeno<-apply(as.matrix(set.list1[[conv.idx]]),1,function(v)read.table(paste(v,".geno",sep=""),header=T))
		setgeno1<-do.call(cbind,setgeno)
		snp.maf<-apply(setgeno1,2,get.maf)
		snp.mac<-apply(setgeno1,2,function(v)length(which(v>0)))
		filter.idx<-which(snp.maf<0.05&snp.mac>=min.mac)
		setgeno2<-setgeno1[,filter.idx]
		rslt.skat[[i]]<-SKAT.one(y=y,geno=setgeno2,covadat=NULL,weight="yes")
		cat("Test for set",i,"\n");flush.console()
	}
	rslt.OneSkat<-do.call(rbind,rslt.skat)
	return(rslt.OneSkat)
	#write.table(rslt.set1,paste("/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/geneset/GS.QTest",outname,".rslt",sep=""),quote=F,row.names=F)
	write.table(rslt.set1,paste(outname,"SKAT_GStest.rslt",sep=""),quote=F,row.names=F)
}



GS.BTest<-function(pheno,yname,set.anno,genename,genefolder,minSetsize=200,STT=0.2,maf.cut=0.05,min.mac=5,weight1=TRUE,weight2=TRUE,weight.type="beta",outname){
	rslt.set<-list()
	cat("Preprocessing...","\n");flush.console();
	apply(set.anno[,-1],1,function(v)as.character(unlist(as.matrix(v))))->set.list
	set.list[set.list==""]<-NA
	set.list<-apply(set.list,2,function(v)as.character(na.omit(v)))
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%genename==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<minSetsize&geneleng1>0)
	y0=pheno[,colnames(pheno)%in%yname==TRUE]
	na.y<-which(is.na(y0)==TRUE)
	setwd(paste("/",genefolder,"/",sep=""))
	for(i in idx.accept){
		conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)
		setgeno<-apply(data.frame(set.list1[[conv.idx]]),1,function(v)read.table(paste(v,".geno",sep=""),header=T))
		names(setgeno)<-data.frame(set.list1[[conv.idx]])[,1]
		setleng0<-unlist(lapply(setgeno,ncol))
		idx.l1<-which(setleng0==1)
		if(length(idx.l1)>0){for(k in 1:length(idx.l1)){names(setgeno[[idx.l1[k]]])<-anno0$snp[anno0$gene==names(setgeno)[idx.l1[k]]]}}
		mafS<-lapply(setgeno,function(v0)apply(as.matrix(v0),2,function(v)mean(na.omit(v))/2))
		macS<-lapply(setgeno,function(v0)apply(as.matrix(v0),2,function(v)length(which(v>0))))
		naS<-mapply(function(v1,v2)which(v1>=maf.cut|v2<min.mac),mafS,macS)
		set.snpname<-lapply(setgeno,names)
		setgeno2<-setgeno1<-mapply(function(v1,v2){if(length(v2)>0)return(v1[,-v2]);if(length(v2)==0)return(v1)},setgeno,naS)
		set.snpname1<-mapply(function(v1,v2){if(length(v2)>0)return(v1[-v2]);if(length(v2)==0)return(v1)},set.snpname,naS)
		if(mode(setgeno2)!="list"){	setgeno2<-setgeno1<-mapply(function(v1,v2){if(length(v2)>0)return(list(v1[,-v2]));if(length(v2)==0)return(v1)},setgeno,naS)}
		if(mode(set.snpname1)!="list"){	set.snpname1<-mapply(function(v1,v2){if(length(v2)>0)return(list(v1[-v2]));if(length(v2)==0)return(v1)},set.snpname,naS)}
		setleng1<-lapply(setgeno2,function(v)ncol(as.matrix(v)))
		if(length(which(setleng1==0))>0){
			setgeno2<-setgeno1[-which(setleng1==0)]
			set.snpname1<-set.snpname1[-which(setleng1==0)]
			setleng1<-lapply(setgeno2,function(v)ncol(as.matrix(v)))
		}
		if(length(setgeno2)>0){	
			if(weight1==TRUE){
				mafS1<-lapply(setgeno2,function(v0)apply(as.matrix(v0),2,function(v)mean(na.omit(v))/2))
				if(weight.type=="inv.sd"){wS<-lapply(mafS1,function(v)1/sqrt(v*(1-v)))}
				if(weight.type=="beta"){wS<-lapply(mafS1,function(v)qbeta(v,1,25,lower.tail=F))}
				wS<-lapply(wS,function(v)v/sum(v))
				colS<-mapply(function(v1,v2)apply(as.matrix(v1),1,function(v0)sum(v2*v0)),setgeno2,wS)
			}
			if(weight1==FALSE){colS<-do.call(cbind,lapply(setgeno2,function(v1)apply(as.matrix(v1),1,function(v0)sum(v0))))}
			if(length(na.y)>0){y<-y0[-na.y];colS<-colS[-na.y,]}
			if(length(na.y)==0){y<-y0}
			rslt<-try(QTest.one(y,covadat=NULL,newgeno=colS,STT,weight=weight2,weight.type=weight.type),TRUE);if(mode(rslt)=="character"){rslt<-QTest.one(y,covadat=NULL,newgeno=colS,STT,weight=FALSE)}
			if(mode(rslt)=="character"){rslt<-NA}
			if(is.na(rslt)==TRUE){rslt<-matrix(c(NA,NA,NA),1,3);colnames(rslt)<-c("p.Q1","p.Q2","p.Q3")}
			rslt1<-try(SKAT.one(y,covadat=NULL,geno=colS,weight="yes"),TRUE)
			if(mode(rslt1)=="character"){rslt1<-NA}
			if(is.na(rslt1)==TRUE){rslt1<-matrix(c(NA,NA),1,2);colnames(rslt1)<-c("p.skat","p.skat1")}
			rslt=c(rslt,rslt1)
			rslt.set[[i]]<-data.frame(set=names(set.list1[conv.idx]),n.gene=length(setgeno2),n.SNP=sum(unlist(setleng1)),rslt)

		}	
		cat("Test for set",i,"\n");flush.console()
	}
	rslt.set1<-do.call(rbind,rslt.set)
	write.table(rslt.set1,paste(outname,"_GS.BTest",".rslt",sep=""),quote=F,row.names=F)
	return(rslt.set1)
}




GS.BTest.bin<-function(pheno,yname,set.anno,genename,genefolder,minSetsize=200,STT=0.2,maf.cut=0.05,min.mac=5,weight1=TRUE,weight2=TRUE,weight.type="beta",outname){
	rslt.set<-list()
	cat("Preprocessing...","\n");flush.console();
	apply(set.anno[,-1],1,function(v)as.character(unlist(as.matrix(v))))->set.list
	set.list[set.list==""]<-NA
	set.list<-apply(set.list,2,function(v)as.character(na.omit(v)))
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%genename==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<minSetsize&geneleng1>0)
	y0=pheno[,colnames(pheno)%in%yname==TRUE]
	na.y<-which(is.na(y0)==TRUE)
	setwd(paste("/",genefolder,"/",sep=""))
	for(i in idx.accept){
		conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)
		setgeno<-apply(data.frame(set.list1[[conv.idx]]),1,function(v)read.table(paste(v,".geno",sep=""),header=T))
		names(setgeno)<-data.frame(set.list1[[conv.idx]])[,1]
		setleng0<-unlist(lapply(setgeno,ncol))
		idx.l1<-which(setleng0==1)
		if(length(idx.l1)>0){for(k in 1:length(idx.l1)){names(setgeno[[idx.l1[k]]])<-anno0$snp[anno0$gene==names(setgeno)[idx.l1[k]]]}}
		mafS<-lapply(setgeno,function(v0)apply(as.matrix(v0),2,function(v)mean(na.omit(v))/2))
		macS<-lapply(setgeno,function(v0)apply(as.matrix(v0),2,function(v)length(which(v>0))))
		naS<-mapply(function(v1,v2)which(v1>=maf.cut|v2<min.mac),mafS,macS)
		set.snpname<-lapply(setgeno,names)
		setgeno2<-setgeno1<-mapply(function(v1,v2){if(length(v2)>0)return(v1[,-v2]);if(length(v2)==0)return(v1)},setgeno,naS)
		set.snpname1<-mapply(function(v1,v2){if(length(v2)>0)return(v1[-v2]);if(length(v2)==0)return(v1)},set.snpname,naS)
		if(mode(setgeno2)!="list"){	setgeno2<-setgeno1<-mapply(function(v1,v2){if(length(v2)>0)return(list(v1[,-v2]));if(length(v2)==0)return(v1)},setgeno,naS)}
		if(mode(set.snpname1)!="list"){	set.snpname1<-mapply(function(v1,v2){if(length(v2)>0)return(list(v1[-v2]));if(length(v2)==0)return(v1)},set.snpname,naS)}
		setleng1<-lapply(setgeno2,function(v)ncol(as.matrix(v)))
		if(length(which(setleng1==0))>0){
			setgeno2<-setgeno1[-which(setleng1==0)]
			set.snpname1<-set.snpname1[-which(setleng1==0)]
			setleng1<-lapply(setgeno2,function(v)ncol(as.matrix(v)))
		}
		if(length(setgeno2)>0){	
			if(weight1==TRUE){
				mafS1<-lapply(setgeno2,function(v0)apply(as.matrix(v0),2,function(v)mean(na.omit(v))/2))
				if(weight.type=="inv.sd"){wS<-lapply(mafS1,function(v)1/sqrt(v*(1-v)))}
				if(weight.type=="beta"){wS<-lapply(mafS1,function(v)qbeta(v,1,25,lower.tail=F))}
				wS<-lapply(wS,function(v)v/sum(v))
				colS<-mapply(function(v1,v2)apply(as.matrix(v1),1,function(v0)sum(v2*v0)),setgeno2,wS)
			}
			if(weight1==FALSE){colS<-do.call(cbind,lapply(setgeno2,function(v1)apply(as.matrix(v1),1,function(v0)sum(v0))))}
			if(length(na.y)>0){y<-y0[-na.y];colS<-colS[-na.y,]}
			if(length(na.y)==0){y<-y0}
			rslt<-try(QTest.one.bin(y,covadat=NULL,newgeno=colS,STT,weight=weight2,weight.type=weight.type),TRUE);if(mode(rslt)=="character"){rslt<-QTest.one.bin(y,covadat=NULL,newgeno=colS,STT,weight=FALSE)}
			if(mode(rslt)=="character"){rslt<-NA}
			if(is.na(rslt)==TRUE){rslt<-matrix(c(NA,NA,NA),1,3);colnames(rslt)<-c("p.Q1","p.Q2","p.Q3")}
			#rslt.skat<-try(SKAT.one(y,covadat=NULL,geno=colS,weight="yes"),TRUE)
			rslt.set[[i]]<-data.frame(set=names(set.list1[conv.idx]),n.gene=length(setgeno2),n.SNP=sum(unlist(setleng1)),rslt)
		}	
		cat("Test for set",i,"\n");flush.console()
	}
	rslt.set1<-do.call(rbind,rslt.set)
	write.table(rslt.set1,paste(outname,"_GS.BTest",".rslt",sep=""),quote=F,row.names=F)
	return(rslt.set1)
}

get.gamma.p<-function(ori.p,STT){
	a<-get.a(length(ori.p),STT)
	q2<-2*(qgamma(ori.p,a,1,lower.tail=FALSE))
	Q2<-sum(q2)
	p.Q2<-pchisq(Q2,df=2*a*length(ori.p),lower.tail=FALSE)
	return(p.Q2)
}

get.fisher.p<-function(ori.p,STT){
	q2<-2*(qgamma(ori.p,1,1,lower.tail=FALSE))
	Q2<-sum(q2)
	p.Q2<-pchisq(Q2,df=2*length(ori.p),lower.tail=FALSE)
	return(p.Q2)
}

GS.QTest.pval<-function(set.anno,pvaldat,minSetsize=200,STT=0.2,outname){
	gamma.p<-list()
	apply(set.anno[,-1],1,function(v)as.character(unlist(as.matrix(v))))->set.list
	set.list[set.list==""]<-NA
	set.list<-apply(set.list,2,function(v)as.character(na.omit(v)))
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%pvaldat[,1]==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<minSetsize&geneleng1>0)
	for(i in idx.accept){
		gamma.p[[i]]<-apply(pvaldat[pvaldat[,1]%in%set.list1[[i]]==TRUE,2:13],2,function(v)get.gamma.p(v,STT))
		if(length(na.omit(gamma.p[[i]]))==0){cat("a;sldkfjasdf","\n");flush.console()}
	}
	gamma.pdat<-do.call(rbind,gamma.p)
	rslt<-data.frame(set=set.anno[idx.accept,1],gamma.pdat)
}		
		
GS.QTest.pval1<-function(set.anno,pvaldat,minSetsize=200,STT=0.2,outname){
	gamma.p<-list()
	apply(set.anno[,-1],1,function(v)as.character(unlist(as.matrix(v))))->set.list
	set.list[set.list==""]<-NA
	set.list<-apply(set.list,2,function(v)as.character(na.omit(v)))
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%pvaldat[,1]==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<minSetsize&geneleng1>0)
	for(i in idx.accept){
		gamma.p[[i]]<-apply(pvaldat[pvaldat[,1]%in%set.list1[[i]]==TRUE,2:13],2,function(v)get.fisher.p(na.omit(v)))
		if(length(na.omit(gamma.p[[i]]))==0){cat("a;sldkfjasdf","\n");flush.console()}
	}
	gamma.pdat<-do.call(rbind,gamma.p)
	rslt<-data.frame(set=set.anno[idx.accept,1],gamma.pdat)
}		
		
		
		
#######  Analysis Example 1###############
#source("QTest_code_Final.R")
#null.dist.Q3<-read.table("n.log10.minp.1e09.txt",header=T)
#library(MultiPhen)
#anno0<-read.table("/home/shared/T2Dgenes/project1/2013_March_release_wave3/annotations/wave3_anno_score12.txt",header=T)[,c(2,3,5)]
#colnames(anno0)<-c("CHR","gene","snp")
#pheno<-read.table("Wave3_28feb2013_transformations_allexomes_cleanVCF_exclethnicrelated_n9960_alltraits_w_mega10PCs.ped",header=T)
#allgeno0<-read.plink("chr21_9960_minMAC5") ## read binary plink file
#allgeno0[allgeno0==0]<-4;allgeno0[allgeno0==2]<-0;allgeno0[allgeno0==4]<-2;  ## appropriate 012 coding
#rslt<-QTest.all(phe.cova=pheno,allgeno=allgeno0,anno=anno0,yname="RES.ALT",covaname=NULL,cut.r2=0.05,a=0.1,min.mac=5,maf.cut=0.01,preprocess=TRUE,weight=TRUE,outname="RES.ALT.chr1")
#rslt1<-QTest.all(phe.cova=pheno,allgeno=allgeno0,anno=anno0,yname="LDLcaladj_ALL_RAW_RES_INV",covaname=c("age","sex"),cut.r2=0.05,a=0.2,min.mac=5,maf.cut=0.05,outname="LDLcaladj_ALL_RAW_RES_INV_chr1")
#rslt.SKAT<-SKAT.all(phe.cova=pheno,allgeno=allgeno0,anno=anno0,yname="BMI_ALL_RAW_RES_INV",covaname=paste("PC",1:10,sep=""),min.mac=5,maf.cut=0.01,weight="yes",outname="BMI_ALL_RAW_RES_INV_chr1")


#######  Analysis Example 2 ###############
#source("QTest_code_APR2013.R")
#null.dist.Q3<-read.table("n.log10.minp.1e09.txt",header=T)
#pheno<-read.table("Wave3_28feb2013_transformations_allexomes_cleanVCF_exclethnicrelated_n9960_alltraits_w_mega10PCs.ped",header=T)
#genename<-read.table("a.txt")[,1]
#setwd("/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat")
#rslt<-QTest.bibs(phe.cova=pheno,yname="BMI_ALL_RAW_INV",genename=genename,STT=0.05,min.mac=5,maf.cut=0.05,preprocess=TRUE,weight=TRUE,outname="BMI_ALL_RAW_INV_STT0.2_weightT_preprocessT")
#rslt.skat<-SKAT.bibs(phe.cova=pheno,yname="BMI_ALL_RAW_INV",genename=genename,min.mac=5,maf.cut=0.05,weight="yes",outname="BMI_ALL_RAW_INV")

###### Gene set Analysis NIDDK Example ######

### Under /home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat,  
### you need "QTest_code_APR2013.R", "n.log10.minp.1e09.txt","Wave3_pheno_BIBS.txt", "c2.cp.v3.1.symbols_1.gmt", "wave3_anno_score12_v2.txt"
### Under same folder, you also need "GENENAME.geno" for each gene ("a.txt"  has these file lists)

#	setwd("/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat")
#	source("QTest_code_APR2013.R")
#	null.dist.Q3<-read.table("n.log10.minp.1e09.txt",header=T)
#	pheno<-read.table("Wave3_pheno_BIBS.txt",header=T)
#	set.anno0<-read.table("c2.cp.v3.1.symbols_1.gmt",header=T);set.anno0[,1]<-as.character(set.anno0[,1])
#	anno0<-read.table("wave3_anno_score12_v2.txt",header=T)
#	anno0$gene<-as.character(anno0$gene);anno0$snp<-as.character(anno0$snp)
#	genename<-read.table("a.txt");genename1<-apply(as.matrix(genename),1,function(v)strsplit(v,split=".geno")[[1]])

# for traits 
# "HEIGHT_ALL_RAW_RES_INV",  "BMI_ALL_RAW_RES_INV" ,"HDLexcmed_ALL_RAW_RES_INV","LDLcalexcmed_ALL_RAW_RES_INV", "CHOLexcmed_ALL_RAW_RES_INV","TGexcmed_ALL_LOG_RES_INV"       

# Height#
	# MEGA # GS.QTest(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a0=0.05,a=0.2,maf.cut=0.05,min.mac=5,outname="HEIGHT_ALL_RAW_RES_INV")
	# MEGA with weight # GS.QTest.w(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,genefolder="/home/lgjlh/T2Dp1/plink/Each_Gene_dat/",genename=genename1,maxSetsize=200,a0=0.05,a=0.2,maf.cut=0.05,min.mac=5,outname="HEIGHT_ALL_RAW_RES_INV_weight")
	# cohortID=c("AJ",   "AW",   "EK",   "ES",   "HA"  , "HS" ,  "SL",   "SS" ,  "UA"   ,"UM" ) 
		# for(i in 1:10){rslt<-GS.QTest1(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,cohort=cohortID[i],genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a0=0.05,a=0.2,maf.cut=0.05,min.mac=5,outname=paste("HEIGHT_ALL_RAW_RES_INV_",cohortID[i],sep=""))}
		# for(i in 1:10){rslt<-GS.QTest1.w(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,cohort=cohortID[i],genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a0=0.05,a=0.2,maf.cut=0.05,min.mac=5,outname=paste("HEIGHT_ALL_RAW_RES_INV_",cohortID[i],sep=""))}
	# ethnic.grp<-c("African-American", "East-Asian" ,"European","Hispanic","South-Asian") 
		# for(i in 1:5){rslt<-GS.QTest2(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,cohort=ethnic.grp[i],genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a0=0.05,a=0.2,maf.cut=0.05,min.mac=5,outname=paste("HEIGHT_ALL_RAW_RES_INV_",ethnic.grp[i],sep=""))}
		# for(i in 1:5){rslt<-GS.QTest2.w(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,cohort=ethnic.grp[i],genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a0=0.05,a=0.2,maf.cut=0.05,min.mac=5,outname=paste("HEIGHT_ALL_RAW_RES_INV_",ethnic.grp[i],sep=""))}

# BMI#
# HDL#
# LDL#
# CHOL#
# TG#


#rslt<-GS.BTest(pheno,yname="HEIGHT_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a=0.2,maf.cut=0.05,min.mac=5,outname="HEIGHT_ALL_RAW_RES_INV")
#rslt<-GS.BTest(pheno,yname="BMI_ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a=0.2,maf.cut=0.05,min.mac=5,outname="BMI_ALL_RAW_RES_INV")
#rslt1<-GS.BTest.w(pheno,yname="LDLcalexcmed__ALL_RAW_RES_INV",set.anno=set.anno0,gene.anno=anno0,genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,a=0.2,maf.cut=0.05,min.mac=5,outname="LDLcalexcmed__ALL_RAW_RES_INV")


#setwd("/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat")
#source("QTest_code_Final2013.R")
#set.anno0<-read.table("c2.all.v4.0.symbols_1.gmt",header=F);set.anno0[,1]<-as.character(set.anno0[,1])
#set.anno0<-read.table("c5.all.v4.0.symbols_1.gmt",header=F);set.anno0[,1]<-as.character(set.anno0[,1])
#pheno<-read.table("Wave3_pheno_BIBS.txt",header=T)
#anno0<-read.table("wave3_anno_score12_v2.txt",header=T)
#anno0$gene<-as.character(anno0$gene);anno0$snp<-as.character(anno0$snp)
#genename<-read.table("a.txt");genename1<-apply(as.matrix(genename),1,function(v)strsplit(v,split=".geno")[[1]])

#rslt<-GS.BTest(pheno,yname="HDLexcmed_ALL_RAW_RES_INV_w_PC",set.anno=set.anno0,gene.anno=anno0,genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat",genename=genename1,maxSetsize=200,STT=0.2,maf.cut=0.05,min.mac=5,weight1=TRUE,weight2=TRUE,outname="HDLexcmed_ALL_RAW_RES_INV_w_PC")

#yname="LDLcalexcmed_ALL_RAW_RES_INV"
##set.anno=set.anno0
#gene.anno=anno0
#genefolder="/home/shared/T2Dgenes/project1/2013_March_release_wave3/VCFs_dosages/Each_Gene_dat"
#genename=genename1
#maxSetsize=200
#STT=0.2
#maf.cut=0.05
#min.mac=5
##weight1=TRUE
#weight2=TRUE
#weight.type="beta"


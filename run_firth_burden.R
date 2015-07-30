## Load data to use for calculation of method overlap
args <- commandArgs(trailingOnly = TRUE)
library(reshape)
library(stringi)
library(logistf)

d2=read.table(args[1],header=T)
d1=read.table(args[2],header=TRUE)
d3=read.table(args[3], header=T)
d1[d1==-9] <- NA
d1[d1==1] <- 0
d1[d1==2] <- 1
data <- merge(d1,d2,by=c('PatientID'))
print(head(data),2)
data <- data[,-(2:3),drop=FALSE]
print(head(data),2)
data <- merge(data,d3,by=c('PatientID'))
print(head(data),2)
data <- data[,-(1),drop=FALSE]
index1=1
print(head(data),2)
write.table(paste("pheno index time", sep=""), file=paste("run_firth_burden.log",sep=""), col.names=F, append=FALSE, row.names=F, quote=FALSE)

for(nam in colnames(data)[grep("_pheno", colnames(data))]) {
  index2=1
  data$phenosum=0
  data$controlsum<-0
  data$phenosum[data[[nam]]==1]=1
  data$controlsum[data[[nam]]==0]=1
  x=sum(data$phenosum)
  print(paste(""),delim="")
  print(paste("# ON PHENOTYPE: ",nam),delim="")
  if(x >= args[4]) {
    print(paste("...ok, at least ",args[4]," cases...",x," to be exact", delim=""))
    system(paste("mkdir -p results"))
    system(paste("mkdir -p results/",nam,sep=""))
    write.table(paste(nam,index1,Sys.time(), sep=" "), file=paste("run_firth_burden.log",sep=""), col.names=F, append=TRUE, row.names=F, quote=FALSE)
    for(gene in colnames(data)[grep("_gene", colnames(data))]) {
      print(paste("...compare phenotype",nam,"to gene",gene,"...pheno index",index1,"...gene index",index2))
      y=sum(data[[gene]])
      z=sum(data$controlsum)
      case_carriers=length(which(data[[nam]]==1 & data[[gene]]>=1))
      control_carriers=length(which(data[[nam]]==0 & data[[gene]]>=1))
      print(paste("...case carriers:",case_carriers,"Cont carriers:",control_carriers))

      if (case_carriers >= args[5]) {
	  print(paste("...ok, at least ",args[5]," case carriers...",case_carriers," to be exact", delim=""))
          cf_cases=case_carriers/x
          cf_controls=control_carriers/z
          #print(paste("...gene: ",gene," Phenotype: ",nam,"...pheno index:",index1,"...gene index:",index2,sep=""))
          print(paste("run..."),delim="")
          n1=sum(data$geno)
          data$ynvar1=1
          data$ynvar1[data$geno==0]=0
          n2=sum(data$ynvar1)
          #print(table(data$phenosum))
          #print(Sys.time())
   	  data$geno[data[[gene]]==0]=0
          data$geno[data[[gene]]==1]=1
          data$geno[data[[gene]]==2]=2
          data$geno[data[[gene]]==3]=3
          data$geno[data[[gene]]==4]=4
          data$geno[data[[gene]]==5]=5

#####################################################################################
#if(FALSE) {
      # Firth regression
          glm<-logistf(data=data, phenosum~geno+SEX+AGE+AGESQ+PC1)
          print(glm)
          p=glm$prob[2]
          cil=glm$ci.lower[2]
          ciu=glm$ci.upper[2]
          coef=glm$coefficients[2]
          write.table(paste(nam,x,z,case_carriers,control_carriers,cf_cases,cf_controls,coef,cil,ciu,exp(coef),exp(cil),exp(ciu),p,n1,n2,gene,sep=" "), file = paste("results/",nam,"/LoF_test.",nam,
".",gene,".txt",sep=""), col.names=F, append=FALSE, row.names=F, quote=FALSE)
#}

#####################################################################################
  if (FALSE) {
      # Logistic regression as an alternative
          glm<-summary(glm(phenosum~geno+SEX+AGE+AGESQ+PC1, data=data, family=binomial() ))
          #print(glm)
          cc <- coef(glm,na.action = na.pass)
          #print(cc)
	  list <- c(gene)
	  for (h in list) {
            p=try(cc["geno","Pr(>|z|)"])
            c=try(cc["geno","Estimate"])
	    if (class(p) == "try-error") { p="NA" }
            if (class(c) == "try-error") { c="NA" }
            if (grepl("Error", p) == "TRUE") { p="NA" }
            if (grepl("Error", c) == "TRUE") { c="NA" }
            write.table(paste(p,c,n1,n2,case_carriers,control_carriers,cf_cases,cf_controls,nam,gene,x,sep=" "), file = paste("results/",nam,"/LoF_test.",nam,".",gene,".txt",sep=""), col.names=F, a
ppend=FALSE, row.names=F, quote=FALSE)
          }
  }
#####################################################################################
      }
     index2<-index2+1
    }
  }
index1<-index1+1
}

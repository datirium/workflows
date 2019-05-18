#!/usr/bin/env Rscript

library(MASS)
library(affy)
library(R.basic)

common_peak_count_read1 <- read.table("common_peak_count_read1", header=FALSE)
common_peak_count_read2 <- read.table("common_peak_count_read2", header=FALSE)
#peak_count_read1 <- read.table("peak_count_read1",header=FALSE)
#peak_count_read2 <- read.table("peak_count_read2",header=FALSE)

merge_common_peak_count_read1 <- read.table("merge_common_peak_count_read1", header=FALSE)
merge_common_peak_count_read2 <- read.table("merge_common_peak_count_read2", header=FALSE)

#table_MA <-read.table("MAnorm.bed",header=FALSE)

table_merge_MA <- read.table("MAnorm_merge.bed", header=FALSE)


M <- log2((common_peak_count_read1+1)/(common_peak_count_read2+1))
A <- 0.5*log2((common_peak_count_read1+1)*(common_peak_count_read2+1))
M <- as.matrix(M)
A <- as.matrix(A)

#linear<-lm(M~A)$coefficients
#b<-lm(M~A)$coefficients
#b<-robustRegBS(M,A,beta=linear)
b <- rlm(M~A)$coefficients
#png('MAplot_before_rescaling.png')
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
#ma.plot(A,M,cex=1,main="MA plot before rescaling (common peaks)")
#abline(b,col="green")
#dev.off()

#cat("M = b[1] + b[2] * A\n")
#log2_peak_count_read1 <- log2(peak_count_read1 + 1)
#log2_peak_count_read2 <- log2(peak_count_read2 + 1)
#log2_peak_count_read1_rescaled <- (2-b[2])*log2_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);
#M_rescaled <- (log2_peak_count_read1_rescaled - log2_peak_count_read2);
#A_rescaled <- (log2_peak_count_read1_rescaled + log2_peak_count_read2)/2;

#png('MAplot_after_rescaling.png')
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
#ma.plot(as.matrix(A_rescaled),as.matrix(M_rescaled),cex=1,main=" MA plot after rescaling (all peaks)")
#dev.off ()


log2_merge_common_peak_count_read1 <- log2(merge_common_peak_count_read1 + 1)
log2_merge_common_peak_count_read2 <- log2(merge_common_peak_count_read2 + 1)
log2_merge_common_peak_count_read1_rescaled <- (2-b[2])*log2_merge_common_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2])
merge_M_rescaled <- (log2_merge_common_peak_count_read1_rescaled - log2_merge_common_peak_count_read2)
merge_A_rescaled <- (log2_merge_common_peak_count_read1_rescaled + log2_merge_common_peak_count_read2)/2





# function for calculating pvalue
pval <- function(x, y){
	if (x+y<20) { # x + y is small
		p1<- nChooseK(x+y,x) * 2^-(x+y+1)
		p2<- nChooseK(x+y,y) * 2^-(x+y+1)
	} else { # if x+y is large, use approximation
		log_p1 <- (x+y)*log(x+y) - x*log(x) - y*log(y) - (x+y+1)*log(2)
		p1<-exp(log_p1)
		p2<-p1
		#log_p2 <- (x+y)*log(x+y) - x*log(x) - y*log(y) - (x+y+1)*log(2);
		#p2<- exp(log_p2);
	}
	pvalue=max(p1,p2)
	return(pvalue)
}
	
#table_MA[,5] <- peak_count_read1
#table_MA[,6] <- peak_count_read2
#table_MA[,7] <- M_rescaled
#table_MA[,8] <- A_rescaled
#table_MA <-as.data.frame(table_MA)
#table_MA[,9] <- 0
#log2_peak_count_read1_rescaled <- as.matrix(log2_peak_count_read1_rescaled)
#peak_count_read2 <- as.matrix(peak_count_read2)
#for (n in c(1:nrow(table_MA))) {
#        cat(n,'\t',round(2^log2_peak_count_read1_rescaled[n]),'\t',peak_count_read2[n],'\n')
#        table_MA[n,9]<--log10(pval(round(2^(log2_peak_count_read1_rescaled[n])),peak_count_read2[n]))
#}


#colnames(table_MA)[1] = "chr"
#colnames(table_MA)[2] = "start"
#colnames(table_MA)[3] = "end"
#colnames(table_MA)[4] = "description"
#colnames(table_MA)[5] = "#raw_read_1"
#colnames(table_MA)[6] = "#raw_read_2"
#colnames(table_MA)[7] = "M_value_rescaled"
#colnames(table_MA)[8] = "A_value_rescaled"
#colnames(table_MA)[9] = "-log10(p-value)"

#write.table(table_MA,"MAnorm_result.xls",sep="\t",quote=FALSE,row.names=FALSE)

# table_merge
table_merge_MA[,5] <- merge_common_peak_count_read1
table_merge_MA[,6] <- merge_common_peak_count_read2
table_merge_MA[,7] <- merge_M_rescaled
table_merge_MA[,8] <- merge_A_rescaled
table_merge_MA <-as.data.frame(table_merge_MA)
table_merge_MA[,9] <- 0
log2_merge_common_peak_count_read1_rescaled <- as.matrix(log2_merge_common_peak_count_read1_rescaled)
merge_common_peak_count_read2 <- as.matrix(merge_common_peak_count_read2)
for (n in c(1:nrow(table_merge_MA))) {
#        cat(n,'\t',round(2^log2_merge_common_peak_count_read1_rescaled[n]),'\t',merge_common_peak_count_read2[n],'\n')
        table_merge_MA[n,9]<--log10(pval(round(2^(log2_merge_common_peak_count_read1_rescaled[n])),merge_common_peak_count_read2[n]))
}


colnames(table_merge_MA)[1] = "chr"
colnames(table_merge_MA)[2] = "start"
colnames(table_merge_MA)[3] = "end"
colnames(table_merge_MA)[4] = "description"
colnames(table_merge_MA)[5] = "#raw_read_1"
colnames(table_merge_MA)[6] = "#raw_read_2"
colnames(table_merge_MA)[7] = "M_value_rescaled"
colnames(table_merge_MA)[8] = "A_value_rescaled"
colnames(table_merge_MA)[9] = "-log10(p-value)"

write.table(table_merge_MA,"MAnorm_result_commonPeak_merged.xls",sep="\t",quote=FALSE,row.names=FALSE)


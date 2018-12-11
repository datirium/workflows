options(warn=-1)
args <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(DBI))
suppressMessages(library(RMySQL))
suppressMessages(library(BiocParallel))
register(MulticoreParam(4))

DRV <- dbDriver("MySQL")
con <- dbConnect(DRV, user=args[1], password=args[2], dbname=args[3], host=args[4], client.flag = CLIENT_MULTI_STATEMENTS)
EMS <- args[5]
T1 <- args[6]
T2 <- args[7]
RTYPE <- args[8]

T1T <- dbGetQuery(con,paste("select tableName,name from ",EMS,".genelist where leaf=1 and (parent_id like '",T1,"' or id like '",T1,"')",sep=""))
T2T <- dbGetQuery(con,paste("select tableName,name from ",EMS,".genelist where leaf=1 and (parent_id like '",T2,"' or id like '",T2,"')",sep=""))

print("T1T")
print(T1T)
print("T2T")
print(T2T)

T1C <- dim(T1T)[1]
T2C <- dim(T2T)[1]

DESeqA=1
if(T1C>1 && T2C>1)
DESeqA=2


if(DESeqA==2) {
    suppressMessages(library(DESeq2))
} else {
    suppressMessages(library(DESeq))
}


tblEnd="";
if(RTYPE=="1")
tblEnd="_isoforms";
if(RTYPE=="2")
tblEnd="_genes";
if(RTYPE=="3")
tblEnd="_common_tss";


condition <- c(rep("untreated",T1C),rep("treated",T2C))
DF <- data.frame(condition, row.names=c(paste("R_u", 1:T1C, sep=""), paste("R_t", 1:T2C, sep="")))
DF$condition <- factor(DF$condition,levels=c("untreated","treated"))


names<-c()

SQL1 <- "SELECT a1.refseq_id, a1.gene_id, a1.chrom, a1.txStart, a1.txEnd, a1.strand"
SQL2 <- paste(" from `",paste(T1T[1,1],tblEnd,"` a1",sep=""),sep="")
SQL3 <- " where a1.chrom not like 'control' "

for(i in 1:T1C) {
    SQL1 <- paste(SQL1, ",a", i, ".TOT_R_0 as R_u", i, sep="")
    names <- append(names,c(paste("R_u", i, sep="")))
    if(i>1) {
        SQL2 <- paste(SQL2, ", `", paste(T1T[i,1], tblEnd, "` a", i, sep=""), sep="")
        SQL3 <- paste(SQL3, " and a1.refseq_id=a", i, ".refseq_id and a1.gene_id=a", i, ".gene_id and a1.chrom=a", i, ".chrom and a1.txStart=a", i, ".txStart and a1.txEnd=a", i, ".txEnd and a1.strand=a", i, ".strand ", sep="")
    }
}

for(i in 1:T2C) {
    names <- append(names, c(paste("R_t", i, sep="")))
    SQL1 <- paste(SQL1, ",b", i, ".TOT_R_0 as R_t", i, sep="")
    SQL2 <- paste(SQL2, ", `", paste(T2T[i,1], tblEnd, "` b", i, sep=""), sep="")
    SQL3 <- paste(SQL3, " and a1.refseq_id=b", i, ".refseq_id and a1.gene_id=b", i, ".gene_id and a1.chrom=b", i, ".chrom and a1.txStart=b", i, ".txStart and a1.txEnd=b", i, ".txEnd and a1.strand=b", i, ".strand ", sep="")
}

if(T1C==1) {
    SQL1 <- paste(SQL1, ",a1.RPKM_0 ", sep="")
    T1NAME <- T1T[1,2]
} else {
    T1Th <- dbGetQuery(con, paste("select tableName, name from ", EMS, ".genelist where leaf=0 and id like '", T1, "'", sep=""))
    tblName <- paste(T1Th[1,1], tblEnd, sep="")
    T1NAME <- T1Th[1,2]
    SQL1 <- paste(SQL1, ",a", T1C+1, ".RPKM_0 ", sep="")
    SQL2 <- paste(SQL2, ", `", tblName, "` a", T1C+1, sep="")
    i <- T1C+1
    SQL3 <- paste(SQL3, " and a1.refseq_id=a", i, ".refseq_id and a1.gene_id=a", i, ".gene_id and a1.chrom=a", i, ".chrom and a1.txStart=a", i, ".txStart and a1.txEnd=a", i, ".txEnd and a1.strand=a", i, ".strand ", sep="")
}
if(T2C==1) {
    SQL1 <- paste(SQL1, ",b1.RPKM_0 ", sep="")
    T2NAME <- T2T[1,2]
} else {
    T2Th <- dbGetQuery(con, paste("select tableName, name from ", EMS, ".genelist where leaf=0 and id like '", T2, "'", sep=""))
    tblName <- paste(T2Th[1,1], tblEnd, sep="")
    T2NAME <- T2Th[1,2]
    SQL1 <- paste(SQL1, ",b", T2C+1, ".RPKM_0 ", sep="")
    SQL2 <- paste(SQL2, ", `", tblName, "` b", T2C+1, sep="")
    i <- T2C+1
    SQL3 <- paste(SQL3, " and a1.refseq_id=b", i, ".refseq_id and a1.gene_id=b", i, ".gene_id and a1.chrom=b", i, ".chrom and a1.txStart=b", i, ".txStart and a1.txEnd=b", i, ".txEnd and a1.strand=b", i, ".strand ", sep="")
}


SQL <- paste(SQL1, SQL2, SQL3)


fullData <- dbGetQuery(con, SQL)


colnames(fullData) <- c("refseq_id", "gene_id", "chrom", "txStart", "txEnd", "strand", names, T1NAME, T2NAME)

dataDimention <- dim(fullData)[2]
totReadsIndex <- seq(7, dataDimention-2)
totRPKMIndex <- seq(dataDimention-1, dataDimention)

if(DESeqA==2) {
    dse <- DESeqDataSetFromMatrix(countData=fullData[,totReadsIndex], colData=DF, design =~condition)

    if(T1C==1|T2C==1) {
        dsq <- DESeq(dse,fitType="local")
    } else {
        dsq <- DESeq(dse)
    }

    DESeqRes <- as.data.frame(results(dsq)[,c(2,5,6)])
    DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0;
    DESeqRes[is.na(DESeqRes)] <- 1;
} else {
    cds <- newCountDataSet( fullData[,totReadsIndex] , condition)
    cdsF <- estimateSizeFactors( cds )

    if(T1C==1|T2C==1) {
        cdsD <- estimateDispersions(cdsF, method="blind", sharingMode="fit-only", fitType="local")
    } else {
        cdsD <- estimateDispersions(cdsF)
    }

    DESeqRes <- nbinomTest(cdsD, "untreated", "treated" )
    isinfl <- is.infinite(DESeqRes$log2FoldChange)
    DESeqRes$log2FoldChange[isinfl] <- log2((DESeqRes$baseMeanB[isinfl]+0.1)/(DESeqRes$baseMeanA[isinfl]+0.1))
    DESeqRes <- DESeqRes[,c(6,7,8)]
    DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0;
    DESeqRes[is.na(DESeqRes)] <- 1;
}

final <- data.frame(cbind(fullData[,c(1:6,(dataDimention-1):dataDimention)], DESeqRes), check.names=F, check.rows=F)
colnames(final) <- c("refseq_id", "gene_id", "chrom", "txStart", "txEnd", "strand", paste("RPKM", T1NAME), paste("RPKM", T2NAME), "LOGR", "pvalue", "padj")
write.table(format(final[with(final, order(chrom, gene_id, refseq_id, txStart, txEnd, strand)),], digits=8), file=paste(args[9], ".tsv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

dbDisconnect(con)

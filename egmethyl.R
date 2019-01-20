library(data.table)

library(matlib)


library(bigmemory)
library(mvtnorm)
library(MASS)
source("aSPU.R")
source("aSPUO.R")
source("daSPU.R")
source("dist_support.R")

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))

option_list = list(
make_option("--sumstats", action="store", default=NA, type='character',
help="summary statistics (rds file and must have SNP and Z column headers) [required]"),
make_option("--out", action="store", default=NA, type='character',
help="Path to output files [required]"),
make_option("--out_name", action="store", default=NA, type='character',
help="output file name [required]"),
make_option("--gene_list", action="store", default=NA, type='character',
help="Gene list to be analyzed [required]"),
make_option("--tissue", action="store", default=NA, type='character',
help="Tissue to use. Either MCF7 or Hippo. [required]"),
make_option("--test_method", action="store", default="NA", type='character',
help="gene-based methods to be used. [required]")
)
opt = parse_args(OptionParser(option_list=option_list))


if(opt$test_method=="aSPU") {
    library(Rcpp)
    library(RcppArmadillo)
    sourceCpp("GMaSPU_support.cpp")
}


#opt = list(sumstats = "scz1.txt", out = "./out",out_name = "res",tissue = "Hippo",gene_list = "gene_list.txt",test.method = "SSU")

if(opt$tissue == "Hippo"){
    weights = "./WEIGHTS/pos_hip_all.rds"
}

if(opt$tissue == "MCF7"){
    weights = "./WEIGHTS/pos_hui2_all.rds"
}
outd =opt$out
system(paste("mkdir -p ", opt$out, sep=""))
sumstats = opt$sumstats

gene.list.input = fread(opt$gene_list)
gene.list.input = as.data.frame(gene.list.input)

wgtlist = readRDS(weights)

wgtlist0 = wgtlist[wgtlist[,5] %in% gene.list.input[,1],]
colnames(wgtlist0) = c("WGT","CHR","P0","P1","ID")

#wgtlist0 = wgtlist0[1:20,]
sumstat.orgin =  fread(sumstats)
sumstat.orgin = as.data.frame(sumstat.orgin)

# only save related information
snp.inf = list()

for(i in 1:nrow(wgtlist0)) {
    tryCatch({
        wgt.dat = readRDS(wgtlist0$WGT[i])
        
        snp.inf[[i]] = unique( wgt.dat$SNP)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

snp.inf = unlist(snp.inf)
snp.inf = unique(snp.inf)

m = match(snp.inf, sumstat.orgin$SNP)
## For each wgt file:
sumstat.orgin = sumstat.orgin[m,]

rm(m)
rm(snp.inf)

gene.list = read.table("glist-hg19.txt")
gene.candidate = gene.list[order(gene.list[,1]),]

out.fun <-function(w) {
    tryCatch({
        # Load in summary stats
        genename = wgtlist0$ID[w]
        
        if (opt$tissue == "Hippo") {
            pliot = readRDS("hip_encoderoadmap.rds")
        }
        
        if (opt$tissue == "MCF7") {
            pliot = readRDS("hui_allenhancer.rds")
        }
        
            coding.region = gene.candidate[gene.candidate[,4] %in% genename,1:3]
            coding.region[1,2] = coding.region[1,2] - 5 * 100
            coding.region[1,3] = coding.region[1,3] + 5 * 100
            
            colnames(coding.region) = c("chr","P0","P1")
        
        if(sum(pliot[,1] %in% genename) <1) {
            enchancer = coding.region
            if(is.null(enchancer)) {
                next
            }
        } else {
            pliot.tmp = pliot[pliot[,1] %in% genename,]
            enchancer = pliot.tmp[,2:4]
            colnames(enchancer) = c("chr","P0","P1")
            enchancer = rbind(enchancer,coding.region)
        }
        
        sumstat = sumstat.orgin
        # Load in reference data

        wgt.dat = readRDS(wgtlist0$WGT[w])
        
        wgt.dat.2 = NULL
        for (i in 1:dim(enchancer)[1]) {
            wgt.dat.tmp = wgt.dat[wgt.dat[,"BP"] > enchancer[i,2] & wgt.dat[,"BP"] < enchancer[i,3],]
            wgt.dat.2 = rbind(wgt.dat.2,wgt.dat.tmp)
        }
        
        wgt.dat = wgt.dat.2
        wgt.dat = wgt.dat[!duplicated(wgt.dat),]
        
        
        gene.body.data = wgt.dat[wgt.dat[,"BP"] > coding.region[1,2] & wgt.dat[,"BP"] < coding.region[1,3],]
        
        
        tmp.snp.list = wgt.dat[,1]
        tmp.snp.filename = paste0(outd,"/methy_snp",genename,".txt")
        write.table(tmp.snp.list,tmp.snp.filename, row.names = FALSE, col.names = FALSE,quote = FALSE)
        
        system(paste0("plink --bfile ./LDRF/1000G.EUR.METHY2.CHR",wgtlist0$CHR[w],
        " --extract ",tmp.snp.filename," --maf 0.05 --make-bed --out ",outd,"/test_",genename))
        
        # Load in reference data
        genos = read_plink(paste0(outd,"/test_",genename),impute="avg")
        system(paste0("rm ",outd,"/test_",genename,"*"))
        system(paste0("rm ",tmp.snp.filename))
        
        genos.bim.all = genos$bim
        genos.bim.all = genos.bim.all[nchar(genos.bim.all[,5]) ==1 & nchar(genos.bim.all[,6]) ==1,]
        
        genos.bim = genos.bim.all[genos.bim.all[,2] %in% wgt.dat[,1],]
        genos.bed = genos$bed[,genos.bim[,2]]
        
        # Match summary data to input, record NA where summary data is missing
        m = match( genos.bim[,2] , sumstat$SNP )
        sum.missing = is.na(m)
        sumstat = sumstat[m,]
        sumstat$SNP = genos.bim[,2]
        sumstat$A1[ sum.missing ] = genos.bim[sum.missing,5]
        sumstat$A2[ sum.missing ] = genos.bim[sum.missing,6]
        
        # QC / allele-flip the input and output
        qc = allele.qc( sumstat$A1 , sumstat$A2 , genos.bim[,5] , genos.bim[,6] )
        
        # Flip Z-scores for mismatching alleles
        sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
        sumstat$A1[ qc$flip ] = genos.bim[qc$flip,5]
        sumstat$A2[ qc$flip ] = genos.bim[qc$flip,6]
        
        # Remove strand ambiguous SNPs (if any)
        if ( sum(!qc$keep) > 0 ) {
            genos.bim = genos.bim[qc$keep,]
            genos.bed = genos.bed[,qc$keep]
            sumstat = sumstat[qc$keep,]
        }
        
        cpg.name = unique(wgt.dat[,2])
        
        wgt.matrix.list = list()
        snp.list = list()
        
        cpg.name.used = NULL
        if(length(cpg.name)<1) next
        
        for (i in cpg.name) {
            tryCatch({
                wgt.matrix = wgt.dat[wgt.dat[,2] ==i,]
                wgt.matrix = wgt.matrix[nchar(wgt.matrix[,"A1"])==1,]
                wgt.matrix = wgt.matrix[nchar(wgt.matrix[,"A2"])==1,]
                
                m = match( wgt.matrix[,"SNP"] , genos.bim[,2] )
                m.keep = !is.na(m)
                wgt.matrix = wgt.matrix[m.keep,]
                
                rownames(wgt.matrix) = wgt.matrix[,"SNP"]
                cur.genos = genos.bed[,m[m.keep]]
                
                geno.bim.tmp = genos.bim[m[m.keep],]
                qc = allele.qc( wgt.matrix[,"A1"] , wgt.matrix[,"A2"] , geno.bim.tmp[,5] , geno.bim.tmp[,6] )
                
                # Flip Z-scores for mismatching alleles
                wgt.matrix$beta[ qc$flip ] = -1 * wgt.matrix$beta[ qc$flip ]
                wgt.matrix$A1[ qc$flip ] = geno.bim.tmp[qc$flip,5]
                wgt.matrix$A2[ qc$flip ] = geno.bim.tmp[qc$flip,6]
                
                if( class(cur.genos)=="numeric") {
                    next
                } else if (dim(cur.genos)[2] < 2) {
                    next
                }
                tmp.cor <- cor(cur.genos)
                tmp.cor <- tmp.cor^2
                tmp.cor[upper.tri(tmp.cor)] <- 0
                diag(tmp.cor) <- 0
                
                cur.genos <- cur.genos[,!apply(tmp.cor,2,function(x) any(x > 0.95))]
                if(is.null(dim(cur.genos)) ) next

                ref.tmp = cov(cur.genos)
                wgt.matrix = wgt.matrix[rownames(ref.tmp),]
                weight.tmp = abs(qnorm(wgt.matrix[,"p-value"]/2)) #* sign(wgt.matrix[,"beta"])
                
                # weight.tmp = weight.tmp/ (sum(abs(weight.tmp)) + 0.00001)
                weight.tmp = as.matrix(weight.tmp)
                rownames(weight.tmp) = rownames(ref.tmp)
                
                wgt.matrix.list[[as.character(i)]] = weight.tmp
                snp.list[[as.character(i)]] = rownames(weight.tmp)
                cpg.name.used = c(cpg.name.used,i)
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
        
        unique.snps <- unique(unlist(snp.list))
        
        cur.FAIL = FALSE
        if (length(unique.snps) <1) {
            cat( "WARNING : " , unlist(wgtlist0[w,]) , " do not have non-zero weights\n")
            cur.FAIL = TRUE
        }
        
        m = match( unique.snps , genos.bim[,2] )
        m.keep = !is.na(m)
        unique.snps = unique.snps[m.keep]
        cur.genos = scale(genos.bed[,m[m.keep]])
        
        
        tmp.cor <- cor(cur.genos)^2
        tmp.cor[upper.tri(tmp.cor)] <- 0
        diag(tmp.cor) <- 0
        
        cur.genos <- cur.genos[,!apply(tmp.cor,2,function(x) any(x > 0.95))]

        cur.bim = genos.bim[genos.bim[,2] %in% colnames(cur.genos), ]
        cur.genos = cur.genos[,cur.bim[,2]]
        
        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , sumstat$SNP)
        cur.Z = sumstat$Z[m]
        
        # Compute LD matrix
        cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
        
        cur.miss = is.na(cur.Z)
        # Impute missing Z-scores
        if ( sum(cur.miss) != 0 ) {
            if ( sum(!cur.miss) == 0 ) {
                cat( "WARNING : " , unlist(wgtlist0[w,]) , " had no overlapping GWAS Z-scores\n")
                cur.FAIL = TRUE
            } else {
                cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
            }
        }
        
        len.cpg = length(cpg.name.used)
        weight_mat = matrix(0,dim(cur.bim)[1],len.cpg)
        rownames(weight_mat) = cur.bim[,2]
        colnames(weight_mat) = cpg.name.used
        
        for (i in 1:len.cpg) {
            tmp.weight = wgt.matrix.list[[as.character(cpg.name.used[i])]]
            weight_mat[rownames(tmp.weight),i] = tmp.weight
        }
        
        wgt.matrix = weight_mat
        wgt.matrix.final = matrix(NA,dim(wgt.matrix)[1],3)
        
        wgt.matrix.final[,1] = apply(wgt.matrix,1,sum) # sum it up
        wgt.matrix.final[,2] = as.numeric(wgt.matrix.final[,1]!=0) # non-zero 0-1 weight
        wgt.matrix.final[,3] = apply(wgt.matrix,1,max.abs) # take the maximum; we tried several ways to combine the weights, take 0-1 weight is the best and thus has been used in E + G + Methyl
        wgt.matrix = wgt.matrix.final
        pval.res = NULL
        
        mod.best = 2
        
        non_zero_ind <- (wgt.matrix[,mod.best] != 0)
        
        ### standardize weights
        diag_element <- as.vector(wgt.matrix[non_zero_ind,mod.best])
        diag_sd <- diag_element/sum(abs(diag_element))
        weight_diag <- diag(diag_sd,nrow = length(diag_sd))
        
        Zstat.w <- weight_diag %*% cur.Z[non_zero_ind]
        corSNP.w <- weight_diag %*% cur.LD[non_zero_ind, non_zero_ind] %*% t(weight_diag)
        
        pSum = Sum(U=Zstat.w, CovS=corSNP.w)
        pSSU = SumSqU(U=Zstat.w, CovS=corSNP.w)
        #pUminP = UminPd(U=Zstat.w, CovS=corSNP.w)
        
        U = Zstat.w
        V = corSNP.w
        weight = rep(1,length(U))
        
        pval.res = c(pSum,pSSU)
        if(opt$test_method == "aSPU") {
        PaSPU.w <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e4)$pvs
        
        if(min(PaSPU.w)< 5e-4) {
            PaSPU.w <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e5)$pvs
        }
        
        if(min(PaSPU.w) < 5e-5) {
            PaSPU.w <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e6)$pvs
        }
        
        if(min(PaSPU.w) < 5e-6) {
            PaSPU.w <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e7)$pvs
        }
        
        pval.res = c(pval.res,PaSPU.w)
        }
        # sink(paste(outd,"/out_",job,".txt",sep=""),append=T)
        tmp.res = c(wgtlist0$CHR[w],as.character(wgtlist0$ID[w]),wgtlist0$P0[w], wgtlist0$P1[w],dim(weight_mat)[1],dim(weight_mat)[2]) #7
        
        tmp.res =  c(tmp.res,pval.res) #6
        
        #pUminP = UminPd(U=Zstat.w, CovS=corSNP.w)
        #tmp.res = c(tmp.res, pUminP) #14
        
        cat(w,";")
        if(is.null(tmp.res)) {
            return(rep(NA,5))
        } else {
            return(tmp.res)
        }
        
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

out.res = as.data.frame(matrix(NA,nrow(wgtlist0),16))
tmp.chr = unique(wgtlist0$CHR)

colnames(out.res) = c("CHR","ID","P0","P1","nSNP","nCpG","SPU(1)","SPU(2)","SPU(1)-perm","SPU(2)-perm","SPU(3)","SPU(4)","SPU(5)","SPU(6)","SPU(Inf)","aSPU")
out.i = 1
for(i in 1:dim(wgtlist0)[1]) {
    res = out.fun(i)
    if(!is.null(res)) out.res[out.i,1:length(res)] = res
    out.i = out.i + 1
}

out.res = out.res[!is.na(out.res[,1]),]
out.res = out.res[,colSums(is.na(out.res)) <0.9 * dim(out.res)[1]]

out.file = paste(outd,"/",opt$out_name,".txt",sep="")
write.table(out.res,out.file,quote=F,row.names=F)

system(paste("rm ",outd,"/test_*",sep=""))
system(paste("rm ",outd,"/methy_*",sep=""))

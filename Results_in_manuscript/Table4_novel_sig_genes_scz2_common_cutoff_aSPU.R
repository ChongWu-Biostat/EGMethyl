######################################################
# Multiple Imputation for cell type proportion    ####
# Version 1.00                                    ####
# Feb 2, 2015                                     ####
# Author: Chong Wu, Weihua Guan                   ####
######################################################

setwd("/Users/uniquechong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/Codes/Methy4/ResAna/")

res.out = as.data.frame(matrix(NA,3,11))

cutoff = 0.05/ 10000

file.name = "enhancer_MCF_scz2"
#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

mcf.gene.list = data[,2]


rownames(data) = paste0(data[,2],data[,3],data[,4])
cis = readRDS("summary_best_hui2_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])
data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)


data.tmp = data[data[,6] < cutoff | data[,7] < cutoff,]

data.tmp1 = data.tmp[data.tmp[12] > 5e-8,]



file.name = "enhancer_hip_scz2"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

hip.gene.list = data[,2]


rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]
cis2 = readRDS("summary_best_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)


data.tmp = data[data[,6] < cutoff | data[,7] < cutoff,]
data.tmp2 = data.tmp[data.tmp[12] > 5e-8,]

data.EG = rbind(data.tmp1,data.tmp2)



res = as.data.frame(matrix(NA,1,5))

library(data.table)

gene.inf = fread("gwas_catalog_v1.0.tsv")
gene.inf = as.data.frame(gene.inf)

gene.inf = gene.inf[ grepl("Schizophrenia", gene.inf[,8]),]


for(i in 1:dim(gene.inf)[1]) {
    gene.inf.tmp = gene.inf[i,]
    tmp.name = gene.inf.tmp["MAPPED_GENE"]
    tmp.name = as.character(tmp.name)
    tmp.name = unlist(strsplit(tmp.name, " - "))
    tmp.name = unlist(strsplit(tmp.name, ","))
    if(is.null(tmp.name)) next
    tmp.name = unlist(strsplit(tmp.name, ";"))
    tmp.name = trimws(tmp.name)
    tmp.name1 = tmp.name
    
    tmp.name = gene.inf.tmp[14]
    tmp.name = as.character(tmp.name)
    tmp.name = unlist(strsplit(tmp.name, " - "))
    tmp.name = unlist(strsplit(tmp.name, ", "))
    tmp.name = unlist(strsplit(tmp.name, ","))
    if(is.null(tmp.name)) next
    tmp.name = unlist(strsplit(tmp.name, ";"))
    tmp.name = trimws(tmp.name)
    tmp.name = c(tmp.name,tmp.name1)
    tmp.name = unique(tmp.name)
    
    res.tmp = as.data.frame(matrix(NA,length(tmp.name),5))
    res.tmp[,2:5] = gene.inf.tmp[c(28,8,6,7)] #c("P.VALUE","DISEASE.TRAIT","LINK",Study
    res.tmp[,1] = tmp.name
    res = rbind(res,res.tmp)
}

res = res[-1,]
gene.list = unique(res[,1])


output  = data.EGM[,c(2,1,3,4,5,6,10,11,18,19)]
final.output = as.data.frame(matrix(-1,dim(output)[1],11))
final.output[,1:10] = output[,1:10]
gene.name = output[,1]
tmp.index = NULL
for(index in 1:dim(output)[1]) {
    if(gene.name[index] %in% gene.list) {
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,11] = tmp[4]
    }
}

write.table(final.output,"tmp.txt")

file.name = "std_scz2"

#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_hui2_cis_scz1.rds")
#cis = readRDS("summary_best_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
#cis2 = readRDS("summary_best_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)


data.tmp = data[data[,6] < cutoff | data[,7] < cutoff,]

data.STD = data.tmp[data.tmp[,12] > 5e-8,]



###
file.name = "summary_scz2_YFS"
data = readRDS(paste0(file.name,".rds"))


SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`


fina.SPU1 = SPU1
fina.SPU2 = SPU2

fina.aSPU = SPU1
YFS = readRDS("scz2_YFS.rds")

fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU1 = fina.aSPU[fina.aSPU[,18] > 5e-8,]

fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU2 = fina.aSPU[fina.aSPU[,18] > 5e-8,]

data.TWAS = rbind(fina.aSPU1,fina.aSPU2)


res.out[3,6] = dim(YFS)[1]

file.name = "summary_scz2_NTR"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

fina.aSPU = rbind(fina.aSPU,aSPU)
fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

YFS = readRDS("scz2_NTR.rds")


fina.aSPU = SPU1
fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]


fina.aSPU1 = fina.aSPU[fina.aSPU[,18] > 5e-8,]


fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU2 = fina.aSPU[fina.aSPU[,18] > 5e-8,]
data.TWAS2 = rbind(fina.aSPU1,fina.aSPU2)
data.TWAS = rbind(data.TWAS,data.TWAS2)

file.name = "summary_scz2_METSIM"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

fina.aSPU = rbind(fina.aSPU,aSPU)
fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

fina.aSPU = SPU1

YFS = readRDS("scz2_METSIM.rds")

fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU1 = fina.aSPU[fina.aSPU[,18] > 5e-8,]


fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU2 = fina.aSPU[fina.aSPU[,18] > 5e-8,]
data.TWAS2 = rbind(fina.aSPU1,fina.aSPU2)
data.TWAS = rbind(data.TWAS,data.TWAS2)

file.name = "summary_scz2_CMC2"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

fina.aSPU = rbind(fina.aSPU,aSPU)
fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

fina.aSPU = SPU1

YFS = readRDS("scz2_CMC2.rds")

fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU1 = fina.aSPU[fina.aSPU[,18] > 5e-8,]


fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU2 = fina.aSPU[fina.aSPU[,18] > 5e-8,]
data.TWAS2 = rbind(fina.aSPU1,fina.aSPU2)


data.TWAS = rbind(data.TWAS,data.TWAS2)

#data.STD

#data.EGM

#data.EG

gene.name.other = c(data.TWAS[,2], data.STD[,2], data.EG[,2])
gene.name.other = unique(gene.name.other)


#################################
### New aSPU results
#################################


file.name = "enhancer_MCF_scz2"
#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)

mcf.gene.list = data[,2]


file.name = "new_scz2_ALL_hui2"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]

data = data[data[,2] %in% mcf.gene.list,] # if there is no information in enhancer, there is no point to conduct test.
dim(data)
#data[,3] = data[,3] +500
#data[,4] = data[,4] -500

rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_hui2_cis_scz1.rds")
rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],18))

data1[,1:14] = data[,c(1:7,10,11,8,9,16,17)]
data1[,15:16] = cis[,c(6,5)]
data1[,17:18] = cis2[,c(6,5)]

data = data1
data[,8:18] = apply(data[,8:18], 2, FUN=as.numeric)

data$weight = "MCF"
data.mcf = data

colSums(data< cutoff)
gene.name = data[,2]

#######################
### Hip data        ###
#######################
file.name = "enhancer_hip_scz2"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

rownames(data) = paste0(data[,2],data[,3],data[,4])

hip.gene.list = data[,2]


file.name = "new_scz2_ALL_hip"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]

data = data[data[,2] %in% hip.gene.list,] # if there is no information in enhancer, there is no point to conduct test.
dim(data)

data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_cis_scz1.rds")
rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],18))

data1[,1:14] = data[,c(1:7,10,11,8,9,16,17)]
data1[,15:16] = cis[,c(6,5)]
data1[,17:18] = cis2[,c(6,5)]

data = data1
data[,8:18] = apply(data[,8:18], 2, FUN=as.numeric)

data$weight = "Hip"
data.hip = data


data.aSPU = rbind(data.mcf,data.hip)

data.aSPU = data.aSPU[data.aSPU[,18]> 5 * 10^{-8},]

aSPU = data.aSPU[data.aSPU[,13]<cutoff,]

aSPU = aSPU[!aSPU[,2]%in% gene.name.other, ]
aSPU = aSPU[aSPU[,10]>=cutoff &aSPU[,11]>=cutoff ,]





res = as.data.frame(matrix(NA,1,5))

library(data.table)

gene.inf = fread("gwas_catalog_v1.0.tsv")
gene.inf = as.data.frame(gene.inf)

gene.inf = gene.inf[ grepl("Schizophrenia", gene.inf[,8]),]


for(i in 1:dim(gene.inf)[1]) {
    gene.inf.tmp = gene.inf[i,]
    tmp.name = gene.inf.tmp["MAPPED_GENE"]
    tmp.name = as.character(tmp.name)
    tmp.name = unlist(strsplit(tmp.name, " - "))
    tmp.name = unlist(strsplit(tmp.name, ","))
    if(is.null(tmp.name)) next
    tmp.name = unlist(strsplit(tmp.name, ";"))
    tmp.name = trimws(tmp.name)
    tmp.name1 = tmp.name
    
    tmp.name = gene.inf.tmp[14]
    tmp.name = as.character(tmp.name)
    tmp.name = unlist(strsplit(tmp.name, " - "))
    tmp.name = unlist(strsplit(tmp.name, ", "))
    tmp.name = unlist(strsplit(tmp.name, ","))
    if(is.null(tmp.name)) next
    tmp.name = unlist(strsplit(tmp.name, ";"))
    tmp.name = trimws(tmp.name)
    tmp.name = c(tmp.name,tmp.name1)
    tmp.name = unique(tmp.name)
    
    res.tmp = as.data.frame(matrix(NA,length(tmp.name),5))
    res.tmp[,2:5] = gene.inf.tmp[c(28,8,6,7)] #c("P.VALUE","DISEASE.TRAIT","LINK",Study
    res.tmp[,1] = tmp.name
    res = rbind(res,res.tmp)
}

res = res[-1,]

output  = aSPU[,c(1:6,10,11,13,18,19)]
final.output = as.data.frame(matrix(-1,dim(output)[1],12))
final.output[,1:11] = output[,1:11]
gene.name = output[,1]
tmp.index = NULL
for(index in 1:dim(output)[1]) {
    if(gene.name[index] %in% gene.list) {
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,12] = tmp[4]
    }
}

final.output2 = final.output[,c(2,1,5,6,7,8,9,10,11,12)]
library(xtable)

print(xtable(final.output2,digits = c(0,0,0,0,0,-1,-1,-1,-1,0,0)), math.style.exponents = TRUE,include.rownames=FALSE)



#res.out = res.out[,1:6]
xtable(res.out)

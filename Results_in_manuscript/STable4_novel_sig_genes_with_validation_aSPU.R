setwd("/Users/uniquechong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/Codes/Methy4/ResAna/")


cutoff = 0.05/10000

file.name = "enhancer_MCF_scz1"
#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)

mcf.gene.list = data[,2]


file.name = "new_scz1_ALL_hui2"
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

data1 = as.data.frame(matrix(NA,dim(data)[1],13))

data1[,1:9] = data[,c(1:6,8,9,17)]
data1[,10:11] = cis[,c(6,5)]
data1[,12:13] = cis2[,c(6,5)]

data = data1
data[,3:13] = apply(data[,3:13], 2, FUN=as.numeric)

data$weight = "MCF"
data.mcf = data

colSums(data< cutoff)
gene.name = data[,2]

#######################
### Hip data        ###
#######################
file.name = "enhancer_hip_scz1"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

rownames(data) = paste0(data[,2],data[,3],data[,4])

hip.gene.list = data[,2]


file.name = "new_scz1_ALL_hip"
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

data1 = as.data.frame(matrix(NA,dim(data)[1],13))

data1[,1:9] = data[,c(1:6,8,9,17)]
data1[,10:11] = cis[,c(6,5)]
data1[,12:13] = cis2[,c(6,5)]

data = data1
data[,3:13] = apply(data[,3:13], 2, FUN=as.numeric)

data$weight = "Hippo"
data.hip = data

gene.name2 = data[,2]
gene.name = c(gene.name,gene.name2)


output = rbind(data.mcf,data.hip)



output = output[,c(1:9,13,14)]
output[,1] = as.numeric(output[,1])

output[,3:10] = apply(output[,3:10],2,as.numeric)

#output = output[output[,7]< cutoff |output[,8] < cutoff,]

output = output[output[,9] < cutoff,]
output = output[order(output[,1],output[,3]),]


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


sum(unique(output[,2]) %in% gene.list )
sum(unique(gene.name) %in% gene.list)
length(gene.name)

length(unique(output[,2]))


length(output[,2])

phyper(sum(output[,2] %in% gene.list )-1,sum(gene.name %in% gene.list), length(gene.name) - length(output[,2]), length(output[,2]) , lower.tail = FALSE)




final.output = as.data.frame(matrix(-1,dim(output)[1],12))
final.output[,1:11] = output[,1:11]
gene.name = output[,2]
tmp.index = NULL
for(index in 1:dim(output)[1]) {
    if(gene.name[index] %in% gene.list) {
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,12] = tmp[4]
    }
}
sum(gene.name %in% gene.list)

write.table(final.output,"/Users/uniquechong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/Codes/Methy4/ResAna/Res_manuscript_final/EGM_results_final/S4.txt")
#tmp.index
#library(xtable)
#final.output2 = final.output[,c(2,1,5,6,7,8,10,11,12)]

#print(xtable(final.output2,digits = c(0,0,0,0,0,-1,-1,-1,0,0)), math.style.exponents = TRUE,include.rownames=FALSE)

#xtable(final.output,rownames= F)
#res.out = res.out[,1:9]
#
#res.out = res.out[,1:6]
#xtable(res.out)




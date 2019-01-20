setwd("/Users/uniquechong/Desktop/egmtl/WEIGHTS")

files = list.files()


for(i in 1:length(files)) {
    tryCatch({
        
        input.name = files[i]
        tmp = readRDS(input.name)
        tmp = tmp[,c(7:10,1:5)]
        
        tmp.name = colnames(tmp)
        tmp.name[6] = "CpG"
        colnames(tmp) = tmp.name
        
        input.name = gsub(".rds","",input.name)
        
        if(grepl('hip_all', input.name)) {
            input.name = gsub("hip_all.","",input.name)
            save.name = paste("/Users/uniquechong/Desktop/egmtl/WEIGHTS2/Hippo_",input.name,".bed",sep="")

        } else if (grepl('hui2_all', input.name)) {
            input.name = gsub("hui2_all.","",input.name)
            save.name = paste("/Users/uniquechong/Desktop/egmtl/WEIGHTS2/MCF7_",input.name,".bed",sep="")
        } else {
            next
        }
        
        write.table(tmp,save.name,row.names=F,quote=F,sep=" ")

    }, error=function(e) {
        cat("ERROR :",conditionMessage(e), "\n")
        
    }
    
    )

}

tmp =readRDS("pos_hui_all.rds")

# enhancer information

tmp = readRDS("/Users/uniquechong/Desktop/egmtl/hui_allenhancer.rds")
tmp = tmp[,c(2,3,4,1)]
tmp[,1] = as.numeric(tmp[,1])
tmp = tmp[order(tmp[,1]),]

tmp[is.na(tmp[,1]),1] = "X"
colnames(tmp) = c("CHR","enhancer_start","enhancer_end","gene_ID")
write.table(tmp,"/Users/uniquechong/Desktop/egmtl/MCF7_enhancer.bed",row.names=F,quote=F,sep=" ")


tmp = readRDS("/Users/uniquechong/Desktop/egmtl/hip_encoderoadmap.rds")
tmp = tmp[,c(2,3,4,1)]
tmp[,1] = as.numeric(tmp[,1])
tmp = tmp[order(tmp[,1]),]
tmp[is.na(tmp[,1]),1] = "X"

colnames(tmp) = c("CHR","enhancer_start","enhancer_end","gene_ID")
write.table(tmp,"/Users/uniquechong/Desktop/egmtl/Hippo_enhancer.bed",row.names=F,quote=F,sep=" ")

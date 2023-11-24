
setwd("/Users/huachunyin/desktop/Xinqiao/Panx1/2022年8月/")

library(stringr)
library(tidyverse)
library(tidyr)
library(RColorBrewer)
df<-openxlsx::read.xlsx("expressiondata.xlsx")
#count<-openxlsx::read.xlsx("fpkm.xlsx")
#protein-coding


# 2) Preprocessing.
dim(df)
 

# Unique duplicate genes.
df1 <- aggregate(.~df$gene_name,data = count[,2:15], mean)
#dim(df1)
 

rownames(df1) <- df1[,1]
df1 <- df1[,-1]

#filtering the low expresion 
df1 = df1[rowMeans(df1)>10,] 
dim(df1)

# Check the data.
#boxplot(log(df1+1))


# Get sample name.
get_sample_name <- function(count_matrix){
  sample_name <- colnames(df1)
  sample_name
  for(i in 1:length(sample_name)){
    sample_name_tmp <- strsplit(sample_name[i], split = '_', fixed = T)[[1]][2]
    sample_name[i] <- sample_name_tmp
  }
  sample_name
  colnames(df1) <- sample_name
  
  return(sample_name)
}

group <- get_sample_name(count_matrix = df1)


pheatmap::pheatmap(log10(df1+0.1),
                   show_rownames = F,
                   show_colnames = T,
                   cluster_cols = T,
                   cluster_rows=T,
                   fontsize_row=8, #??????????????????
                   border_color = "NA",#FFFFFF",
                   scale = "row",
                   angle_col=90, #??????????????????
                   #annotation_row = anno,
                   #annotation_colors = ann_colors,
                   color =colorRampPalette(c("#0225A2", "white","#FF2400"))(100)
                   #color =colorRampPalette(c("#5BB93B","#030303", "#C22D1D"))(100)
)



################MSPJ

source("Chunk01-R script for RNA-seq data simulation based on improved compcodeR.R")


### Step-01. Preparing the colors used in this study. 

mypal1 <- terrain.colors(14)

mypal2 <- pal_npg("nrc", alpha = 0.7)(10)

### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Selecting the dataset, for DNA microarray and RNA-seq. 

### @@@@@@@@@@@@@@ Selecting the dataset @@@@@@@@@@@@@@ ###

eset <- df1

eset[eset<0]=NA

eset<-na.omit(eset)

print(eset[1:6, 1:6])

# Visualizing the gene expression matrix. 
boxplot(log2(eset), 
        col = mypal1, 
        main = "RNA-sequencing data",
        #xaxt = "n", 
        names=group,outline = FALSE)


### ------------------------------------------------------------------------ ###
### Step-03. Generating multiple sub-groups based resampling for primary study. 

# ord.gene: which gene you focused on. 

eset <- log2(eset+0.1)

colnames(eset) <- c(rep("Control",7),rep("Experimental",7))

eset <- as.matrix(eset)
sample.sets <- generateSubGroup(eset, 
                                set.n = 30, # the number of sub-groups
                                size.min = 8, # the lower limit of sample size,default Value:10
                                size.max = 14#11
                                ) # the maximum sample size???default Value:20

print(sample.sets[[1]][1:5, 1:5])

### End of Step-03.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-04. Preprocessing for DNA microarray or RNA-seq data, alternatively. 

# For all sub-datasets.  

for (d in 1:length(sample.sets)) {
  
  data <- sample.sets[[d]]
  
  tmp <- normalize.quantiles(data)
  
  rownames(tmp) <- rownames(data)
  
  colnames(tmp) <- colnames(data)
  
  data <- tmp
  
  #DGElist <- DGEList(counts = data)
  
  #DGElist <- calcNormFactors(DGElist, method = "upperquartile")
  
  #. boxplot(log2(DGElist$count))
  
  # plotMDS(DGElist)
  
  #data <- DGElist$count  

  sample.sets[[d]] <- data
  
}

print(sample.sets[[1]][1:5, 1:5])


####################MSPJ########################################################

### ------------------------------------------------------------------------ ###
### Step-01. Computing the statistics used for meta-analysis.

set.n <- length(sample.sets)

num.gene <- nrow(eset)

# 1) Identifying the differentially expressed genes by meta-analysis, in bulk. 

deg <- batchMeta(data.list = sample.sets, 
                 cutoff = 1, 
                 g.start = 1, 
                 g.end = num.gene)
#deg1, cutoff=0

deg.meta <- c(rownames(eset)[deg$UpGene], 
              rownames(eset)[deg$DownGene])

#deg.meta1 <- c(rownames(eset)[deg1$UpGene], rownames(eset)[deg1$DownGene])
#
# -- Generate the meta object. 
gene.ind <- c(deg$UpGene,deg$DownGene)

degmeta <- matrix(NA,nrow=length(deg.meta),ncol=2)

colnames(degmeta) <- c("geneID","SMD")
  
for (g1 in 1:length(deg.meta)) {
  
  res <- singleMeta(data.list = sample.sets, 
                    ord.gene = gene.ind[g1])
  
  SMD<- res$TE.fixed
  degmeta[g1,1] <- deg.meta[g1]
  degmeta[g1,2] <- SMD
  
}

View(degmeta)
### -----------------------SVMRFE-------------------------------------------- ###
### Step-01. Data preprocessing for whole dataset, seq.matrix or mcr.matrix. 

#-- seq.matrix.RData: simulated RNA-seq data; 
#-- mcr.matrix.RData: simulated microarray data. 

eset.mat <- as.data.frame(t(eset))

sam.lab <- as.factor(colnames(eset))

input <- cbind(sam.lab, eset.mat)

input <- as.data.frame(input)

#save(input, file = "input.RData")

### End of Step-02.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-04. Gene selection using multiple SVM-RFE algorithm.

# 1) Generating multiple sub-groups based resampling, and setting up cross validation.

if (nrow(input)<13) {
  
  results <- svmRFE(input)
  
  deg.svm <- colnames(input)[results]
  
  deg.svm<-deg.svm[1:floor(nrow(eset)*0.1)]
  
  
} else{
  
  set.seed(1)
  
  nfold <- 5
  
  nrows <- nrow(input)
  
  folds <- rep(1:nfold, len=nrows)[sample(nrows)]
  
  folds <- lapply(1:nfold, function(x) which(folds == x))
  
  # 2) Perform feature ranking on all training sets
  
  
  feedback<-try(results <- lapply(folds, svmRFE.wrap, input, k = 5, halve.above = 1000),silent = TRUE)
  
  if("try-error"%in%class(feedback)){
    
    set.seed(2)
    
    folds <- rep(1:nfold, len=nrows)[sample(nrows)]
    
    folds <- lapply(1:nfold, function(x) which(folds == x))
    
    results <- lapply(folds, svmRFE.wrap, input, k = 5, halve.above = 1000)
  }
  
  length(results)
  
  # 3) Obtain top features across ALL folds
  
  top.features <- WriteFeatures(results, input, save = FALSE)
  
}

head(top.features)

# 4) Selecting the top n genes as the DEGs identified by SVM-RFE method, eg. 500. 
topn <- 1000
deg.svm <- top.features$FeatureName[1:topn]#floor(nrow(eset)*0.1)

deg.svm <- as.character(deg.svm)

deg.svm


### End of Step-04.
### ------------------------------------------------------------------------ ###

#########################Permutation test#######################################

all.genes <- rownames(eset)

gene.count <- length(all.genes)


deg.count <- NULL

i <- 0

pb <- txtProgressBar(min = 0, max = gene.count, style = 3, char = "+")

for (g in all.genes) {
  
  i <- i + 1
  
  setTxtProgressBar(pb, i)
  
  # Display the progress bar!
  
  tmp <- Summarize(get(g) ~ sam.lab, data = input, digits = 3)
  
  #print(tmp)
  #boxplot(get(g) ~ sam.lab, data = input)
  
  deg.per <- try(independence_test(get(g) ~ sam.lab, 
                                   data = input,
                                   distribution = approximate(nresample = 1000)), 
                 silent = FALSE)
  
  # deg.Z <- deg.per@statistic@teststatistic
  deg.p <- deg.per@distribution@pvalue(deg.per@statistic@teststatistic)
  
  flush.console()
  
  # This variable, deg.count, stores all the differentially expressed genes.
  
  if (!is.na(deg.p) & deg.p < 0.05) deg.count <- c(deg.count, g,deg.p) else next
  
}

close(pb)


deg.count1<-matrix(deg.count,ncol=2,byrow = T)

colnames(deg.count1)<-c("FeatureName","Pvalue")

deg.count1<-as.data.frame(deg.count1)

deg.count1$Pvalue<-as.vector(deg.count1$Pvalue)

deg.count1$Pvalue<-as.numeric(deg.count1$Pvalue)

fdrp=p.adjust(deg.count1$Pvalue, "BH")
fdrp <-data.frame(FeatureName=deg.count1$FeatureName,
                  Pvalue=deg.count1$Pvalue,
                  adjP=fdrp)


top.features1 <- merge(top.features[which(top.features$FeatureName%in%deg.svm),],fdrp,by="FeatureName")

top.features1 <-  merge(top.features1,degmeta,
                        by.x="FeatureName",by.y="geneID")
#############################commongene#########################


setlist <- list(Meta.analysis = deg.meta, 
                SVM.RFE = deg.svm, 
                Permutation = fdrp$FeatureName)

#save(setlist, file = "setlist.RData")

venn(setlist, 
     lty = 0, 
     col = "navyblue", 
     zcolor = 1:3, 
     lwd = 2, 
     box = FALSE)

commongene <- intersect(deg.svm,intersect(deg.meta,deg.per))

degsvm <- top.features[1:topn,c(1,3)]

com.gene<- merge(degmeta,degsvm,by.x="geneID",by.y="FeatureName")


com.gene <- merge(com.gene,deg.per.data,by.x="geneID",by.y="name")

##################################################

e.p <- which(group%in%"OE")
c.p <- which(group%in%"NC")

FC <- apply(eset, 1, function(x) {mean(x[e.p])/mean(x[c.p])})

deg.FC<-data.frame(FC=FC,log2FC=log2(FC))

deg.FC$FeatureName <- rownames(deg.FC)

com.gene1 <- merge(top.features1,deg.FC,by="FeatureName")

openxlsx::write.xlsx(com.gene1,"14samplefpkmDEG.xlsx")

hp <- log2(eset[com.gene1$geneID,]+0.1)

pheatmap::pheatmap(eset[com.gene1$geneID,],
                   show_rownames = T,
                   show_colnames = T,
                   cluster_cols = T,
                   cluster_rows=T,
                   fontsize_row=8, #??????????????????
                   border_color = "NA",#FFFFFF",
                   scale = "row",
                   angle_col=90, #??????????????????
                   #annotation_row = anno,
                   #annotation_colors = ann_colors,
                   color =colorRampPalette(c("#0225A2", "white","#FF2400"))(100)
                   #color =colorRampPalette(c("#5BB93B","#030303", "#C22D1D"))(100)
)


library(GSVA)

library("KEGGREST")

keggGetclass<-function(x){
  
  gs<-keggGet(x)
  
  class<-gs[[1]]$CLASS
  
  if(length(class)>0) {
    
    a <- try(unlist(strsplit(class,"; ")))
    
    b<-a[2]
    
  } else {
    
    b<-NA}
  
  return(b)
  
  Sys.sleep(1)
}

keggGetcompound<-function(x){
  
  gs<-keggGet(x)
  
  COM<-names(gs[[1]]$COMPOUND)
  
}

keggGetgene<-function(x){
  
  gs<-keggGet(x)
  
  Sys.sleep(1)
  
  Gene<-gs[[1]]$GENE
  
  if(length(Gene)>0) {
    
    gene <- Gene[str_which(Gene,";")] #get the gene symbol 
    
    a <- try(unlist(strsplit(gene,";")))

    b<-a[-str_which(a,":")]
    
    } else {
      b<-NA}
 
 
  
  return(b)
  
}

keggGetGC<-function(x){
  
  gs<-keggGet(x)
  
  COM<-names(gs[[1]]$COMPOUND)
  
  gene.index<-seq(1,length(gs[[1]]$GENE),by=2)
  
  Gene<-gs[[1]]$GENE[gene.index]
  
  map<-c(COM,Gene)
}


KEGGlist<-keggList("pathway","rno")

keggnames<-names(KEGGlist)

keggnames<-substring(keggnames,"6")

keggnames_list <- list()

keggnames_list <- c(keggnames_list, keggnames)

KEGGclass <- lapply(keggnames_list, keggGetclass)#obtain the class of  pathway

KEGGg <- lapply(keggnames_list, keggGetgene)#obtain the genes in metabolism pathway

KEGGcom <- lapply(keggnames_list, keggGetcompound)#obtain the genes in metabolism pathway


KEGGlist1<-as.data.frame(KEGGlist)

names(KEGGg)<-KEGGlist1$KEGGlist

names(KEGGcom)<-KEGGlist1$KEGGlist

KEGGg1=KEGGg[!sapply(KEGGg,is.null)]#delet the null element  

KEGGg2 <- KEGGg1[-which(is.na(KEGGg1)==TRUE)]

KEGGcom1=KEGGcom[!sapply(KEGGcom,is.null)]#delet the null element  

#get the metabolites
KEGGmetagene <- KEGGg[names(KEGGcom1)]
  
gsva_es <- GSVA::gsva(eset, KEGGg2, method="plage")
#method="gsva","ssgsea", "zscore", "plage".  kcdf=c("Gaussian", "Poisson", "none")
View(gsva_es)

rownames(gsva_es) <- unlist(strsplit(rownames(gsva_es), 
                                     split = " - Rattus norvegicus (rat)", 
                                     fixed = T))


colnames(gsva_es)<- paste0(colnames(gsva_es),1:14)

#########

annorow <- data.frame(Sample=factor(c(rep("Control",7),rep("OE",7))))#

rownames(annorow) <- colnames(gsva_es)

KEGGclass3 <- unlist(KEGGclass2)

annocol <- as.data.frame(KEGGclass3)

annocol <- cbind(annocol,KEGGlist1)

rownames(annocol) <- annocol$KEGGlist

rownames(annocol) <- unlist(strsplit(rownames(annocol), 
                                     split = " - Rattus norvegicus (rat)", 
                                     fixed = T))

annocol <- annocol %>% 
  mutate(KEGGclass=dplyr::case_when(
    str_detect(KEGGclass3,"metabolism") ~ "Metabolism",
    
    str_detect(KEGGclass3,"Metabolism")  ~ "Metabolism",
    
    str_detect(KEGGclass3,"metabolites")  ~ "Metabolism",
    
    str_detect(KEGGclass3,"Immune")  ~ "Immune",
    
    str_detect(KEGGclass3,"metabolic") ~ "Metabolism")
  )

annocol <- annocol[-which(is.na(annocol$KEGGclass)),]

annocol1 <- as.data.frame(annocol[,-2])
rownames(annocol1) <-rownames(annocol)
colnames(annocol1)[1] <- "Class"


get_anno_for_heatmap2<-function(annocol,annorow=NULL,color=NULL,only.color=F){
  require(plyr)
  require(stringr)
  if(is.null(color)){
    require(RColorBrewer)
    color=c(brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),
            brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"),
            brewer.pal(9,"Pastel1"), brewer.pal(8,"Pastel2"))
  }
  
  annocolor=do.call(as.list,list(x=annocol))
  annocolor=lapply(annocolor,function(x){if(is.factor(x)){x=levels(x);a=color[1:length(x)];names(a)=x;return(a)}else{x=unique(x);a=color[1:length(x)];names(a)=x;return(a)}})
  if(!is.null(annorow)){
    annocolor.row<-do.call(as.list,list(x=annorow))
    annocolor.row=lapply(annocolor.row,function(x){if(is.factor(x)){x=levels(x);a=color[1:length(x)];names(a)=x;return(a)}else{x=unique(x);a=color[1:length(x)];names(a)=x;return(a)}})
  }else{annocolor.row=NULL}
  annocolor=c(annocolor,annocolor.row)
  annocolor_col<-as.list(annocol)
  annocolor_row<-as.list(annorow)
  annocolor<-c(annocol,annorow)
  annocolor<-lapply(annocolor,function(x){if(is.factor(x)){x=levels(x);return(x)}else{x=unique(x);return(x)}})
  annocolor<-do.call(c,annocolor)
  annocolor<-data.frame(var_name=as.factor(stringr::str_replace(names(annocolor),"[0-9]{1,}$","")),
                        var=annocolor,
                        color=color[1:length(annocolor)])
  annocolor<-split(annocolor,annocolor$var_name)
  annocolor<-lapply(annocolor,function(x){a=x$var;b=as.character(x$color);names(b)=a;return(b)})
  

  
  
}

col_colors<-get_anno_for_heatmap2(annorow,annocol1)


gsva_M <- gsva_es[rownames(annocol1),]
gsva_M <- gsva_M[order(gsva_M[,11]),]
pheatmap::pheatmap(gsva_M,
         show_rownames = T,
         show_colnames = F,
         cluster_cols = T,
         cluster_rows=T,
         fontsize_row=7, #??????????????????
         border_color = "NA",#FFFFFF",
         scale = "none",
         angle_col=90, #??????????????????
         cutree_rows =3, 
         treeheight_row = 30,
         annotation_row = annocol1,
         annotation_col = annorow,
         annotation_colors = col_colors,
         color =colorRampPalette(c("#0225A2", "white","#FF2400"))(100)
         #color =colorRampPalette(c("#5BB93B","#030303", "#C22D1D"))(100)
)
###########volcano plot########
top.features1$SMD <- as.numeric(as.character(top.features1$SMD))

com.gene1$SMD <- as.numeric(as.character(com.gene1$SMD))

top.features2 <- top.features1 %>% mutate(
  state=case_when(
    SMD > 1 & as.character(top.features1$FeatureName)%in%as.character(com.gene1$geneID) ~ 'up',
    
    SMD < -1 & as.character(top.features1$FeatureName)%in%as.character(com.gene1$geneID)~ 'down',
    
    TRUE ~ 'none'
  )
)
GM_gene <- com.gene1[order(abs(com.gene1$SMD),decreasing = T),][1:20,]

GM_gene <- GM_gene %>% mutate(
  state=case_when(
    SMD > 1  ~ 'up',
    
    SMD < -1 ~'down'
  )
)
GM_gene$FeatureName <- as.character(GM_gene$FeatureName)
library(ggrepel)


ggplot(top.features2,
       aes(x=SMD,y=-log10(Pvalue),
       color=state))+
  geom_point()+
  scale_color_manual(values = c("ins"="grey",
                                "up"="#EDD3D6",
                                "down"="#CFDEF4"))+
  geom_text_repel(data = GM_gene,
                aes(label = FeatureName),
                 size = 3,segment.color = "black",
                  max.overlaps = 100)+
  
  theme_classic()+
  ylab('-log10 (P value)')+
  xlab('SMD')+
  #ylim(0,5)+
  xlim(-7,7)+
  geom_vline(xintercept=c(-1,1),
             lty=3,
             col="black",
             lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),
             lty=3,
             col="black",
             lwd=0.5)+
  theme(text = element_text(family="serif",size = 8))


##################enrichment

KEGGor <- openxlsx::read.xlsx("KEGG.xlsx")

KEGG <- KEGGor[,c(3,5,11)]


KEGG1 <- strsplit(KEGG$Hit.in.Query.List,",")

names(KEGG1) <- KEGG$Name


KEGG1 <- as.data.frame(unlist(KEGG1))


KEGG1$pathway <- unlist(strsplit(rownames(KEGG1), "[0-9]+$"))
  

KEGGlinks <- merge(KEGG1,KEGG[,1:2],by.x="pathway",by.y="Name")

colnames(KEGGlinks)[2] <- "Gene"

KEGGlinks$`p-value` <- abs(log10(KEGGlinks$`p-value`))
#openxlsx::write.xlsx(KEGGlinks,"KEGGlinks.xlsx")
KEGGNodes <- KEGGor[,c(3,9)]
#openxlsx::write.xlsx(KEGGNodes,"KEGGNodes.xlsx")
KEGGNodes<- openxlsx::read.xlsx("KEGGNodes.xlsx")

library(networkD3)

#simpleNetwork(KEGGlinks[,1:2],
#              fontFamily="serif",#字体设置如"华文行楷" 等
#              fontSize = 20, #节点文本标签的数字字体大小（以像素为单位）。
#              linkColour="grey",#连线颜色,black,red,blue,  
#              nodeColour="#FF4500",#节点颜色,red，蓝色blue,cyan,yellow等
#              charge = -100,#数值表示节点排斥强度（负值）或吸引力（正值）  
#              opacity=0.5,#透明度,1及以上为不透明，0为完全透明  
#              zoom=TRUE #可缩放
#              )
#

KEGGlinks2 <- KEGGlinks
KEGGNodes2 <- KEGGNodes

KEGGlinks2$pathway=as.numeric(factor(KEGGlinks2$pathway,levels = KEGGNodes2$Name))-1

a <- KEGGlinks2 %>% 
  dplyr::group_by(Gene) %>% 
  count_()

gene <- as.character(a$Gene)
class_gene <- c()
for(i in gene){
  index <- which(i == KEGGlinks$Gene)
  pathway_check <- KEGGlinks$pathway[index]
  if(length(pathway_check) != 1){
    pathwayclass <- KEGGNodes[KEGGNodes$Name%in%pathway_check, "Class"]
    if(any(pathwayclass%in%"Metabolism associated")){class_gene <- append(class_gene, "Metabolism associated")}else{
      class_gene <- append(class_gene, pathwayclass[1])
    }
    
  }else{
    class_gene <- append(class_gene,
                         KEGGNodes[KEGGNodes$Name == pathway_check, "Class"])
  }
}
  
a$Class <- class_gene
colnames(a) <- c("Name","Size","Class")

a <- as.data.frame(a)
KEGGNodes2 <- rbind(KEGGNodes2,a)

KEGGlinks2$Gene=as.numeric(factor(KEGGlinks2$Gene,levels = gene))+57


## Specify colours for specific edges
metaph <- KEGGNodes2[which(KEGGNodes2$Class=="Metabolism associated" & str_detect(KEGGNodes2$Name," ")),]
metabolicindex <- which(KEGGNodes2$Name%in%metaph$Name , arr = TRUE) - 1

#metabolicindex <- which(KEGGNodes2 == "Metabolic pathways" , arr = TRUE)[1] - 1

metabolicindex = which(KEGGlinks2$pathway %in%  metabolicindex)

# Create a colour vector
metabolicCols = ifelse(1:nrow(KEGGlinks2) %in% metabolicindex, "#F2BD46", "#666")


p <- forceNetwork(Links = KEGGlinks2,#线性质数据框
             Nodes = KEGGNodes2,#节点性质数据框
             Source = "pathway",#连线的源变量
             Target = "Gene",#连线的目标变量
             Value = "p-value",#连线的粗细值
             NodeID = "Name",#节点名称
             Group = "Class",#节点的分组
             Nodesize = "Size" ,#节点大小，节点数据框中 
             charge = -0.1,
             linkColour =metabolicCols,
             linkWidth=3,
             radiusCalculation = JS("Math.sqrt(d.nodesize)+6"),
             #radiusCalculation = "d.nodesize",
             fontFamily="serif",#字体设置 
             fontSize = 12, #节点文本标签的数字字体大小（以像素为单位）。
             opacity = 0.8, 
             legend = T,opacityNoHover = TRUE,
             bounded = F,
             zoom = T)
p

forceNetwork(Links = KEGGlinks2,#线性质数据框
             Nodes = KEGGNodes2,#节点性质数据框
             Source = "pathway",#连线的源变量
             Target = "Gene",#连线的目标变量
             Value = "p-value",#连线的粗细值
             NodeID = "Name",#节点名称
             Group = "Class",#节点的分组
             Nodesize = "Size" ,#节点大小，节点数据框中
             ###美化部分
             fontFamily="serif",#字体设置 
             fontSize = 10, #节点文本标签的数字字体大小（以像素为单位）。
             linkColour="black",#连线颜色,black,red,blue,
             #colourScale ,linkWidth,#节点颜色,red，蓝色blue,cyan,yellow等
             charge = -10,#数值表示节点排斥强度（负值）或吸引力（正值）  
             opacity = 0.9,
             legend=T,#显示节点分组的颜色标签
             arrows=F,#是否带方向
             bounded=F,#是否启用限制图像的边框
             opacityNoHover=1.0,#当鼠标悬停在其上时，节点标签文本的不透明度比例的数值
             zoom = T)#允许放缩，双击放大


require(htmlwidgets)
saveWidget(p, file="pathway.html")

require(webshot)
webshot::install_phantomjs("file:/Users/huachunyin/Desktop/Xinqiao/Panx1/2022年8月/pathway.html","pathwaynetwork.pdf")


##############

#data <- matrix(0,ncol=length(unique(KEGGlinks$Gene)),
#               nrow=length(unique(KEGGlinks$pathway)))
#
#colnames(data) <- as.character(unique(KEGGlinks$Gene))
#
#rownames(data) <- unique(KEGGlinks$pathway)
#
#head(data)

for(i in 1:nrow(KEGGlinks)){
  
  data[KEGGlinks[i,1],KEGGlinks[i,2]] <- 1
  
}

KEGGNodes2 <- KEGGNodes2 %>% 
  mutate(
    Class.type=case_when(
      Class=="Metabolism associated" ~ 1,
      Class=="Human diseases" ~ 2,
      Class=="Signal transduction" ~ 3,
      Class=="Cellular processes" ~ 4,
      Class=="Organismal Systems" ~ 5,
      Class=="Genetic Information Processing" ~ 6)
  )

library(igraph)
net <- graph_from_data_frame(d=KEGGlinks, vertices=KEGGNodes2, directed=F) 

class(net)

#E(net)       # The edges of the "net" object

#V(net)       # The vertices of the "net" object

#E(net)$"p-value"  # Edge attribute "p-value"

#V(net)$Size # Vertex attribute "Size"


#plot(net, edge.arrow.size=.4, edge.curved=.1)

# Set edge color to gray, and the node color to orange. 

# Replace the vertex label with the node names 

#plot(net, edge.arrow.size=.2, edge.curved=0,
#     
#     vertex.color="orange", vertex.frame.color=NA,
#     
#     vertex.label=V(net)$Name, vertex.label.color="black",
#     
#     vertex.label.cex=.7) 
#

# Generate colors based on class type:
colrs <- c("#FFCC33", "#32507C", "#97A8C1","#6B9B97","#906C93","#9AB780",alpha=0.7)

V(net)$color <- colrs[V(net)$Class.type]
V(net)$label.color <- "black" #the label color
V(net)$frame.color <- NA #the broder color

# The labels are currently node IDs.
#show the nodes with high degree 
# Setting them to NA will render no labels:
KEGGNodes3 <- KEGGNodes2
KEGGNodes3$Name[which(KEGGNodes3$Size[59:270]<10)+58] <- NA
met.pos <- which(KEGGNodes3$Class[1:58]%in%"Metabolism associated")
pathwy.index <- 1:58
KEGGNodes3$Name[setdiff(pathwy.index,met.pos)] <- NA
V(net)$label <- KEGGNodes3$Name


############scale the node size
a=5;
b=25;
y <- V(net)$Size
Ymax=max(V(net)$Size);#%计算最大值
Ymin=min(V(net)$Size);#%计算最小值
k=(b-a)/(Ymax-Ymin);
norY=a+k*(y-Ymin)

#the edge parameter

#E(net)$width <-E(net)$"p-value"*3 #the edge width
E(net)$edge.color <- "#9E9D9D" #the edge color

#layout
l <- layout_on_sphere(net)
#layout_on_sphere, layout_in_circle, layout_randomly,layout_with_fr,


plot(net, layout=l,
     vertex.label.family="Times",
     vertex.size=norY,#the node size
     edge.arrow.size=.2, 
     edge.curved=0,
     vertex.label=V(net)$label, 
     vertex.label.color="black",
     vertex.label.cex=.7)

legend(x=-1.5, y=-1.1, c("Metabolism associated",
                         "Human diseases",
                         "Signal transduction",
                         "Cellular processes",
                         "Organismal Systems",
                         "Genetic Information Processing"), 
       pch=21,col=NA,
       pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

#cut.off <- mean(links$weight) 
#net.sp <- delete_edges(net, E(net)[weight<cut.off])
#plot(net.sp, layout=layout_with_kk) 



#############GO

upGO <- openxlsx::read.xlsx("upGO.xlsx")

upGO <- upGO[order(upGO$Category,upGO$Hit.Count.in.Query.List,decreasing=T),]

upBP <- upGO[which(upGO$Category%in%"GO: Biological Process"),][1:10,]

downGO <- openxlsx::read.xlsx("downGO.xlsx")

downGO <- downGO[order(downGO$Category,downGO$Hit.Count.in.Query.List,decreasing=T),]

downBP <- downGO[which(downGO$Category%in%"GO: Biological Process"),][1:10,]

colnames(upBP)[8] <- "count"
colnames(downBP)[8] <- "count"

upGOBar<-ggplot(upBP,aes(reorder(Name,count),count))+
  geom_bar(stat="identity",color="white",fill="#7FABCB", width=0.4,size=0.25)+
  coord_flip()+
  labs(x=NULL)+
  theme_classic()+
  theme(text = element_text(family="Times New Roman",size = 10),
        axis.line.y =element_blank())

downGOBar<-ggplot(downBP,aes(reorder(Name,count),count))+
  geom_bar(stat="identity",color="white",fill="#7FABCB", width=0.4,size=0.25)+
  coord_flip()+
  labs(x=NULL)+
  theme_classic()+
  theme(text = element_text(family="Times New Roman",size = 10),
        axis.line.y =element_blank()) 

library(ggpubr)

ggarrange(upGOBar,downGOBar,ncol=2)

#########venn

KEGGmetagene

KEGGmetagene1 <- unlist(KEGGmetagene)

#### metabolites assiated genes
KEGGmetagene1 <- unique(KEGGmetagene1)
#save(KEGGmetagene,file="KEGGmetagene.RData")

library(VennDiagram)

updeg <- com.gene1[as.numeric(as.character(com.gene1$SMD))>0,]
downdeg <- com.gene1[as.numeric(as.character(com.gene1$SMD))<0,]

venn.diagram(list(
  "Metabolism-associated genes" = KEGGmetagene1,
  "DEGs-up" = as.character(updeg$FeatureName),
  "DEGs-down" = as.character(downdeg$FeatureName)),
  filename = "Metabolism-associated.tiff",
  height = 1600,
  width = 1600,
  resolution = 600,imagetype ="tiff",
  col = "transparent",margin = 0.1,
  fill = c("cornflowerblue","orchid3","gray"),
  alpha = 0.70,
  cex = 0.45,cat.cex = 0.45,
  #cat.default.pos = "outer",
  #cat.pos = c(-20,20,180),# ????????????????????????
  cat.dist = c(-0.055,-0.055,-0.055))# ????????????????????????)

metabDEGs <- c(intersect(KEGGmetagene1,as.character(updeg$geneID)),
               intersect(KEGGmetagene1,as.character(downdeg$geneID)))
com.gene$geneID <- as.character(com.gene$geneID)

metabDEGs <- com.gene[which(com.gene$geneID%in%metabDEGs),]


hlgene<-homologene::homologene(metabDEGs$geneID, inTax=10116 , outTax=9606)

metabDEGs1 <- merge(metabDEGs,hlgene,by.x ="geneID",by.y = "10116" )
  
openxlsx::write.xlsx(metabDEGs1,"metabDEGs-294.xlsx")

######################Metablism#############

metamatrix <- openxlsx::read.xlsx("Vector_vs_Panx-OE_info.xlsx",rowNames=TRUE)

metamatrix1 <- metamatrix[,5:16]
 

metamatrix2<-t(metamatrix1)

meta.lab<-group[2:13]

iris.pca <- FactoMineR::PCA(log2(metamatrix2), graph = FALSE)

ind.p <-FactoMineR::fviz_pca_ind(iris.pca,
                     # show points only (nbut not "text") ????????????????????????????????????????????????
                     geom.ind = "point",
                     # ??????????????????
                     col.ind = meta.lab,
                     # ????????????
                     palette = "simpsons",
                     # ???????????? Concentration ellipses
                     addEllipses = TRUE,
                     ellipse.level = 0.80,
                     legend.title = "Groups"
                     
)
#palette = c("#00AFBB","#E7B800")
#color from ggsci (???npg?????????aaas?????????lancet?????????jco?????????ucscgb?????????uchicago?????????simpsons???,???rickandmorty???)
#ellipse.type = "convex"


ggpubr::ggpar(ind.p,
              legend.title = "Group",palette = c("#00468B","#F89B9B"))+
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.x= element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.2), # legend???????????????
        legend.background = element_rect(size = 1, colour = "white"))


pca1 <- as.data.frame(pca$rotation)
 
library(scatterplot3d)
scatterplot3d(pca1$PC1, pca1$PC2, pca1$PC3, 
          groups = as.factor(meta.lab),
          ellipsoid.alpha=0.2,
              color = rep(c('#E9BE3A', '#95A5BD'), c(6, 6)),
              angle=40,cex.symbols=1.5,cex.axis=0.8, 
              xlab = paste('PCA1: 13.2%'), ylab = paste('PCA2: 4.2%'), 
              zlab = paste('PCA3: 3.1%'), pch = 16,grid = F)

#########

###############OPLSDA

set.seed(1)

oplsda = opls(log2(metamatrix2), as.factor(meta.lab), predI = 1, orthoI = NA, crossvalI=7)

sample.score <- oplsda@scoreMN %>%as.data.frame() %>% mutate(label = as.factor(meta.lab),
                                                             o1 = oplsda@orthoScoreMN[,1])

vip <- getVipVn(oplsda)

vip <- as.data.frame(vip)

vip$FeatureName<-rownames(vip)

metamatrix$FeatureName<-rownames(metamatrix)

metamfinal <- merge(metamatrix,vip,by="FeatureName")


metamfinal <- metamfinal%>% 
  dplyr::mutate(Type=case_when(VIP >1 & Fold_Change > 1.5 ~ "up",
                               VIP >1 & Fold_Change < 0.5 ~ "down",
                                    TRUE ~ 'insig')
                
  )


########################

ggplot(metamfinal,
       aes(x=Log2FC,y=VIP,
           color=Type))+
  geom_point( size = 4)+
  scale_color_manual(values = c("insig"="grey",
                                "up"="#EDD3D6",
                                "down"="#95A5BD"))+
  
  theme_classic()+
  ylab('Variable Importance in Projection (VIP)')+
  xlab('Log2FC')+
  ylim(0,1.5)+
  xlim(-2,2)+
  geom_vline(xintercept=c(-1,1),
             lty=3,
             col="black",
             lwd=0.5)+
  geom_hline(yintercept = 1,
             lty=3,
             col="black",
             lwd=0.5)+
  theme(text = element_text(family="serif",size = 12))

###################Joint pathway

Jointpa <- openxlsx::read.xlsx("Joint pathway.xlsx")

colnames(Jointpa)[5] <- "P.value"
ggplot(Jointpa,aes(Impact,reorder(Name,X6)))+
  geom_point(aes(size=Hits,color=X6))+
  scale_color_gradient(low="#95A5BD",high = "#E9BE3A")+
  labs(color=expression(P.value),
       size="Hit number",
       y=NULL,
       x="Impact",
       fill="-log10(P value)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=8),
        text=element_text(family="serif",size=8),
        legend.position = "bottom")


####################

meta.node <- openxlsx::read.xlsx("node-joint netwrok.xlsx")

meta.edge_list <- read.csv("orig_edge_list.csv")

meta.edge_anno <- read.csv("orig_node_list.csv")


colnames(meta.edge_anno)
meta.edge_list <- merge(meta.edge_list,meta.edge_anno[,1:2],by.x = "Source",by.y = "Id")

meta.edge_list <- merge(meta.edge_list[,2:3],meta.edge_anno[,1:2],by.x = "Target",by.y = "Id")

meta.edge_list <-meta.edge_list[,2:3]

colnames(meta.edge_list) <- c("Source","Target")

metabDEGs1 <- openxlsx::read.xlsx("metabDEGs-294.xlsx",rowNames=T)

metadown <- metabDEGs1[as.numeric(as.character(metabDEGs1$SMD))<0,] 


meta.node <- meta.node %>% 
  mutate(Condition=dplyr::case_when(
    meta.node$Label%in%metadown$`9606` ~ 2,
    meta.node$Label%in%"2-Phospho-D-glyceric acid" ~ 2,
    meta.node$Label%in%"Adenosine triphosphate"~ 2,
    meta.node$Label%in%"Uracil"~ 2,
    TRUE ~ 1),
    Type=dplyr::case_when( 
    meta.node$Label%in%metabDEGs1$`9606` ~ 1,
    TRUE ~ 2)
  )


library(igraph)
meta.net <- graph_from_data_frame(d=meta.edge_list, vertices=meta.node, directed=F) 

class(meta.net)
  
V(meta.net)$color <- c("#BF8E9B","#778BA8")[meta.node$Condition]
V(meta.net)$shape <- c("circle","square")[meta.node$Type]


V(meta.net)$label.color <- "black" #the label color
V(meta.net)$frame.color <- NA #the broder color

V(meta.net)$label <- meta.node$Label


############scale the node size
a=9;
b=15;
y <- V(meta.net)$Degree
Ymax=max(V(meta.net)$Degree);#%计算最大值
Ymin=min(V(meta.net)$Degree);#%计算最小值
k=(b-a)/(Ymax-Ymin);
norY=a+k*(y-Ymin)

#the edge parameter

#E(net)$width <-E(net)$"p-value"*3 #the edge width
E(meta.net)$edge.color <- "#9E9D9D" #the edge color

#layout
l <- layout_with_lgl(meta.net)
#layout_on_sphere, layout_in_circle, layout_randomly,layout_with_fr,layout_with_kk,layout_with_sugiyama
#layout_as_tree,

plot(meta.net,layout=l,
     vertex.label.family="Times",
     vertex.size=norY,#the node size
     edge.arrow.size=.2, 
     edge.curved=0,
     vertex.label=V(meta.net)$label, 
     vertex.label.color="black",
     vertex.label.cex=.7)


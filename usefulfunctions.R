
cat("Self-defining functions Collection.\nV3.6   [Autumn Rain]\n")
cat("Author:Zhiming Ye @ Guangzhou Medical University, Nanshan School\n")
cat("load necessary packages...\n")

try(dev.off())
library(tidyverse)
cat("=====Function List=====\n")
FunctionList<-c("Arrange_Table","CalcMad","CorrlatinEstimate0","do.ssGSEA","fpkmToTpm","ggCibersort","MapSymbol","Merge_Table","NMFEstimate","NormalizeTest","PlotPie","rmEnsemblDot","run_DEG_RNAseq","UniCox","WGCNACalcSFT","WGCNAplot2","Write.csvs","Write.tables")
print(FunctionList)
rm(FunctionList)
PlotPie<-function(dt,itemAndfreq=c(1,2),pretty=T,Tag="B",Order="no",col_type="nejm"){
  dt<-as.data.frame(dt)
  colnames(dt)[itemAndfreq[2]]<-"A"
  colnames(dt)[itemAndfreq[1]]<-"B"
  library(ggpubr)
  library(ggsci)
  dt = dt[order(dt$A, decreasing = TRUE),]
  myLabel = as.vector(dt$B)   
  myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")
  colnames(dt)[itemAndfreq[1]]<-"Tag"
  if(Order!="no"){
    dt$Tag<-factor(dt$Tag,levels=Order)
  }
  colnames(dt)[itemAndfreq[1]]<-Tag
  if(pretty){
    if(col_type=="npg"){
      print(ggpie(dt,"A",fill = Tag,label=myLabel)+scale_fill_npg())
    }
    if(col_type=="lancet"){
      print(ggpie(dt,"A",fill = Tag,label=myLabel)+scale_fill_lancet())
    }
    if(col_type=="nejm"){
      print(ggpie(dt,"A",fill = Tag,label=myLabel)+scale_fill_nejm())
    }
  }
  else{
    ggpie(dt,"A",fill = Tag,label=myLabel)
  }
}
Arrange_Table<-function(Mat,FilterList){
  cat("Should have rowname.\n")
  Mat<-as.data.frame(Mat)%>%dplyr::filter(rownames(Mat)%in%FilterList)
  rn1<-rownames(Mat)
  rn2<-FilterList
  Mat<-Mat[rn1[match(rn2, rn1)],]
  return(Mat)
}
do.ssGSEA<-function(Mat,GMTfilePath){
  library(GSVA)
  library(GSEABase)
  genesetHypo <- getGmt(file.path(GMTfilePath))
  CGGA_gsvaHypo<- gsva(as.matrix(Mat),genesetHypo,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE,min.sz=10)
  return(CGGA_gsvaHypo)
}
rmEnsemblDot<-function(Mat,ColNum){
  cat("New Col: SimSymbol.\n")
  Mat<-Mat%>%dplyr::mutate(SimSymbol=sapply(strsplit(Mat[,ColNum],'[.]'),function(x) x[1]))
  return(Mat)
}
MapSymbol<-function(Mat,ColNum){
  cat("New Col: Symbol. \n")
  library(org.Hs.eg.db)
  Mat$Symbol <- mapIds(x = org.Hs.eg.db,
                           keys = Mat[,ColNum],
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
  return(Mat)
}

NMFEstimate<-function(mat,pdfname,pngname,filename){
  library(NMF)
  estimate <- nmf(mat,rank=2:10, method="brunet", nrun=50, seed=486)
  save(estimate,file=filename)
  pdf(pdfname,width = 8,height = 6)
  try(print(plot(estimate)))
  dev.off()
  png(pngname,width = 4400,height = 4000,res=300)
  try(consensusmap(estimate,
               annRow = NA,
               annCol = NA,
               main = "Consensus matrix",
               info = FALSE))
  dev.off()
}

WGCNACalcSFT<-function(Mat,OutPutFile){
  library(WGCNA)
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft = pickSoftThreshold(Mat, powerVector = powers, verbose = 5)
  cat(paste0("sft PowerEstimate is:",sft$powerEstimate,"\n"))
  PowerEst<-sft$powerEstimate
  pdf(OutPutFile,width = 8,height = 5)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
  )
  
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  abline(h = 0.90, col = "red")
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  dev.off()
  return(PowerEst)
}

WGCNAplot2<-function(net,FileName1,FileName2="NA",FileName3){
  cat("FileName1 and 3 is pdf, 2 is png.\n")
  table(net$colors)
  mergedColors = labels2colors(net$colors)
  pdf(FileName1,width = 6,height = 5)
  plotDendroAndColors(
    net$dendrograms[[1]],
    mergedColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  if(FileName2!="NA"){
  geneTree = net$dendrograms[[1]]
  moduleColors = labels2colors(net$colors)
  dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = 14)
  
  plotTOM = dissTOM ^ 7
  diag(plotTOM) = NA
  png("Threeplot.png",width=4000,height = 4000,res=300)
  TOMplot(plotTOM,
          geneTree,
          moduleColors,
          main = "Network heatmap plot, all genes")
  dev.off()}
  MEs = net$MEs
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(
    colnames(MEs), "ME", ""
  ))))
  MEs_col = orderMEs(MEs_col)
  pdf(FileName3,width = 6,height = 7)
  plotEigengeneNetworks(
    MEs_col,
    "Eigengene adjacency heatmap",
    marDendro = c(3, 3, 2, 4),
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = T,
    xLabelsAngle = 90
  )
  dev.off()
}

CalcMad<-function(mat,num){
  mads<-apply(mat,1,mad)
  mat<-mat[rev(order(mads)),]
  return(mat[1:num,])
}

UniCox<-function(StatusAndTimeCol,DATA){
  cat("STATUS=OS,TIME=OStime. Argument 1 is a collection.\n")
  covariates <- colnames(DATA)[-StatusAndTimeCol]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(OStime, OS)~', x)))
  library(survival)
  library(survminer)
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = DATA)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           #获取p值
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           #获取HR
                           HR <-signif(x$coef[2], digits=2);
                           #获取95%置信区间
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(p.value,HR)
                           names(res)<-c("p.value","HR (95% CI for HR)")
                           return(res)
                         })
  #转换成数据框，并转置
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  return(res)
}

NormalizeTest<-function(Table){
  MitoAndHypo2<-as.matrix(Table)
  for(i in 1:ncol(MitoAndHypo2)){
    print(shapiro.test(MitoAndHypo2[,i]))
  }
  
}

Write.tables<-function(table,withRowname=T){
  table<-as.character(table)
  cat("will save to the working dictionary.\n")
  PCAmatnolog<-get(table)
  if(withRowname){
  PCAmatnolog<-PCAmatnolog%>%as.data.frame()%>%rownames_to_column(var="temp3wr4ts")
  colnames(PCAmatnolog)[1]<-""
  }
  write.table(PCAmatnolog, file=paste0(table,".txt"), sep="\t", row.names=F, col.names=T, quote=F)
}
Write.csvs<-function(table){
  table<-as.character(table)
  cat("will save to the working dictionary.\n")
  write.csv(get(table),file=paste0(table,".csv"))
}


Merge_Table<-function(table,cohort,group,col_type,col_id,returntype){
  cat("the cohort colname must be [Group]\n")
  CountsInput<-as.data.frame(table)
  selectvar<-group
  col_geneid<-col_type
  Cohort<-as.data.frame(cohort)%>%dplyr::filter(Group%in%selectvar)
  CountsInput<-as.data.frame(t(CountsInput%>%column_to_rownames(var=col_geneid)))%>%rownames_to_column(var=col_id)
  # return(as.data.frame(merge(Cohort,CountsInput,by=col_id,all.x=T,sort=F))%>%dplyr::select(-Group))
  CountsInputM<-as.data.frame(t(as.data.frame(merge(Cohort,CountsInput,by=col_id,all.x=T,sort=F))%>%dplyr::select(-Group)%>%column_to_rownames(var=col_id)))
  # return(CountsInputM)
  CountsInput<-CountsInputM
  CountsInput2 <- matrix(as.numeric(unlist(CountsInput)),nrow=nrow(CountsInput), dimnames = list(rownames(CountsInput),colnames(CountsInput)))
  data = CountsInput2
  condition<-factor(Cohort$Group)
  coldata <- data.frame(row.names = colnames(CountsInput2), condition)
  if(returntype=="Cohort"){
    return(coldata)
  }
  if(returntype=="WithoutRownames"){
    data<-as.data.frame(data)%>%rownames_to_column(var=col_geneid)
    return(data)
  }
  else{
    return(data)
  }
}

fpkmToTpm <- function(expMatrix)
{
  cat("Might be not correct\n")
  fpkmToTpm2 <- function(fpkm)
  {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  tpms <- apply(expMatrix,2,fpkmToTpm2)
  return(tpms)
}


run_DEG_RNAseq<-function(exprSet,group_list,g1="control",g2="case", pro = "test"){
  ### ---------------
  ###
  ### Firstly run DEseq2 
  ###
  ### ---------------
  cat("Method-1 : DEseq2\n")
  suppressMessages(library(DESeq2))  
  colData <- data.frame(row.names=colnames(exprSet), 
                        group_list=group_list)
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds) 
  res <- results(dds, 
                 contrast=c("group_list",g2,g1))
  resOrdered <- res[order(res$padj),]
  DEG = as.data.frame(resOrdered)
  DEG_DEseq2 = na.omit(DEG)
  
  cat("Method-2 : edgeR\n")
  ### ---------------
  ###
  ### Then run edgeR 
  ###
  ### ---------------
  library(edgeR)
  g=factor(group_list)
  d <- DGEList(counts=exprSet,group=g)
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  DEG_edgeR=as.data.frame(nrDEG)
  
  
  
  # nrDEG=DEG_edgeR[,c(1,5)]
  # colnames(nrDEG)=c('log2FoldChange','pvalue')
  
  ### ---------------
  ###
  ### Then run limma 
  ###
  ### --------------- 
  cat("Method-3 : limma\n")
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  con=paste0(g2,'-',g1)
  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef=con, n=Inf)
  # write.csv(tempOutput,file="limmadegoutput.csv")
  DEG_limma_voom = na.omit(tempOutput)
  # nrDEG=DEG_limma_voom[,c(1,4)]
  # colnames(nrDEG)=c('log2FoldChange','pvalue')
  
  save(DEG_limma_voom,DEG_DEseq2,DEG_edgeR,
       dds,exprSet,group_list,
       file = paste0(pro,'_DEG_results.Rdata')) 
  
  cat("Merge 3 methods' result\n")
  allg=intersect(rownames(DEG_limma_voom),rownames(DEG_edgeR))
  allg=intersect(allg,rownames(DEG_DEseq2))
  nrDEG=cbind(DEG_DEseq2[allg,c(2,6)],
              DEG_edgeR[allg,c(1,5)],
              DEG_limma_voom[allg,c(1,4)]) 
  colnames(nrDEG) = c("logFC_DEseq2","P_DEseq2",
                      "logFC_edgeR","P_edgeR",
                      "logFC_limma","P_limma")
  write.csv(nrDEG,file =  paste0(pro,'_DEG_results.csv'))
}



CorrlatinEstimate0<-function(genegsvamatrix,TargetSet="HALLMARK_HYPOXIA",UseMethod="spearman"){
  data<-genegsvamatrix
  gene_name1<-c()
  gene_name2<-c()
  cor_r<-c()
  pvalue<-c()
  i=which(rownames(data)==TargetSet)
  for (r in 1:nrow(data)){
    g1=rownames(data)[i]
    g2=rownames(data)[r]
    c_r=cor(as.numeric(data[i,]),as.numeric(data[r,]),method=UseMethod)
    p=cor.test(as.numeric(data[i,]),as.numeric(data[r,]),method =UseMethod)[[3]]
    gene_name1=c(gene_name1,g1)
    gene_name2=c(gene_name2,g2)
    cor_r=c(cor_r,c_r)
    pvalue=c(pvalue,p)
    cat(paste0(r,"   done!\n"))
    
  }
  data_cor<-data.frame(gene_name1,gene_name2,cor_r,pvalue)
  return(data_cor)
}

ggCibersort<-function(ResultTable,GroupItems,comparisons,colorList){
  # cat("Sample ID column MUST PLACE in the FIRST column, and its name MUST BE[ID].\n")
  # cat("Sample Group Column Name MUST BE[Group],and place it into the SECOND column.\n")
  # cat("The example of arugument[comparisons](a list):list(c(“Group3”, “Group2”), c(“Group1”, “Group2”), c(“Group3”, “Group1”))\n")
  # cat("The example of arugument[GroupItems](a collection):c(“Group3”,“Group1”,“Group2”)\n")
  library(tidyverse)
  library(reshape2)
  if(ncol(ResultTable)==27){
    k<-melt(ResultTable[,-c(1,25,26,27)],id.vars = "Group")
  }
  if(ncol(ResultTable)==24){
    k<-melt(ResultTable[,-c(1)],id.vars = "Group")
  }
  if(ncol(ResultTable)==23){
    k<-melt(ResultTable,id.vars = "Group")
  }
  library(ggpubr)
  k$Group <- factor(k$Group,levels=GroupItems)
  my_comparisons <- comparisons
  cells<-c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
  p1<-ggboxplot(k%>%dplyr::filter(variable==cells[1]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[1],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p2<-ggboxplot(k%>%dplyr::filter(variable==cells[2]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[2],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p3<-ggboxplot(k%>%dplyr::filter(variable==cells[3]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[3],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p4<-ggboxplot(k%>%dplyr::filter(variable==cells[4]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[4],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p5<-ggboxplot(k%>%dplyr::filter(variable==cells[5]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[5],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p6<-ggboxplot(k%>%dplyr::filter(variable==cells[6]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[6],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p7<-ggboxplot(k%>%dplyr::filter(variable==cells[7]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[7],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p8<-ggboxplot(k%>%dplyr::filter(variable==cells[8]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[8],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p9<-ggboxplot(k%>%dplyr::filter(variable==cells[9]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[9],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p10<-ggboxplot(k%>%dplyr::filter(variable==cells[10]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[10],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p11<-ggboxplot(k%>%dplyr::filter(variable==cells[11]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[11],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p12<-ggboxplot(k%>%dplyr::filter(variable==cells[12]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[12],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p13<-ggboxplot(k%>%dplyr::filter(variable==cells[13]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[13],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p14<-ggboxplot(k%>%dplyr::filter(variable==cells[14]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[14],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p15<-ggboxplot(k%>%dplyr::filter(variable==cells[15]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[15],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p16<-ggboxplot(k%>%dplyr::filter(variable==cells[16]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[16],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p17<-ggboxplot(k%>%dplyr::filter(variable==cells[17]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[17],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p18<-ggboxplot(k%>%dplyr::filter(variable==cells[18]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[18],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p19<-ggboxplot(k%>%dplyr::filter(variable==cells[19]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[19],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p20<-ggboxplot(k%>%dplyr::filter(variable==cells[20]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[20],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p21<-ggboxplot(k%>%dplyr::filter(variable==cells[21]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[21],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  p22<-ggboxplot(k%>%dplyr::filter(variable==cells[22]),x = 'Group',y = 'value',fill='Group',bxp.errorbar = T,ylab=cells[22],palette = colorList)+guides(fill="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test", label = "p.signif")+rotate_x_text()
  pfinal<-cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22, ncol=11,nrow=2,labels=NULL)
  return(pfinal)
}



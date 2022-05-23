# The Mettl3 epitranscriptomic writer amplifies p53 stress responses
# Nitin Raj et al.
# author: Jose A. Seoane joseaseoane@vhio.net

# TCGA PANCAN DE analyses
# Figures S6D, S6E




library(maftools)

#maftools object with mutations (version MC3) and deletions from TCGA
load("/data/datasets/TCGA/PANCAN_33/mut/PANCAN_mut_amp_del_MAF_MC3_v3.RData")

muts = mutCountMatrix(pancan_MAF_ampdel,countOnly = c("Frame_Shift_Del","Frame_Shift_Ins",
                                                      "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                                                      "Nonsense_Mutation","Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))>0
muts.ms = mutCountMatrix(pancan_MAF_ampdel,countOnly = c( "Missense_Mutation" ))>0

muts.tr = mutCountMatrix(pancan_MAF_ampdel,countOnly = c("Frame_Shift_Del","Frame_Shift_Ins",
                                                         "In_Frame_Del", "In_Frame_Ins", 
                                                         "Nonsense_Mutation","Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))>0

dels = mutCountMatrix(pancan_MAF_ampdel,countOnly = "Del")>0

tmp_ms=matrix(F,nrow=2,ncol=ncol(muts),dimnames = list(c("TP53.ms","TP53.tr"),colnames(muts)))
tmp_ms[1,colnames(muts.ms)]=muts.ms["TP53",]
tmp_ms[2,colnames(muts.tr)]=muts.tr["TP53",]
muts = rbind(muts,tmp_ms)


N6_MT_complex=apply(muts[c("METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","ZC3H13"),],2,any)
N6_MT_complex_core=apply(muts[c("METTL3","METTL14","WTAP"),],2,any)

muts2 = rbind(muts,N6_MT_complex,N6_MT_complex_core )

N6_MT_complex=apply(dels[c("METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","ZC3H13"),],2,any)
N6_MT_complex_core=apply(dels[c("METTL3","METTL14","WTAP"),],2,any)

del2 = rbind(dels,N6_MT_complex,N6_MT_complex_core )

allsamples = unique(c(colnames(muts2),colnames(amps2),colnames(del2)))
allgenes = unique(c(rownames(muts2),rownames(amps2),rownames(del2)))

mut3 = matrix(data=F,nrow=length(allgenes),ncol = length(allsamples),dimnames = list(allgenes,allsamples))
del3 = matrix(data=F,nrow=length(allgenes),ncol = length(allsamples),dimnames = list(allgenes,allsamples))
mut3[rownames(muts2),colnames(muts2)]=muts2
del3[rownames(del2),colnames(del2)]=del2

rm(muts2, del2,dels,muts)
gc()

save(mut3,del3,file="mutDistMETTL3.RData")
load("mutDistMETTL3.RData")


library(MultiAssayExperiment)
library(limma)
# expression from TCGA 
load("/data/datasets/TCGA/PANCAN_33/GeneExp_gene_counts.RData")
gc()

coad_exp = voom(all_counts)$E
#coad_exp = t(apply(all_counts,1,function(x) log2(x+1)))
ts = substr(colnames(coad_exp),14,15)
coad_exp.t = coad_exp[,ts !="11"]
colnames(coad_exp.t)=substr(colnames(coad_exp.t),1,12)



cancerTypes = unique(pancan_MAF_ampdel@clinical.data$subtype)
allResults=NULL
for(i in 1:length(cancerTypes)){
  subtype=cancerTypes[i]
  print(subtype)
  try({
    #subtype=cancerTypes[i]
    coad_samples = as.character(pancan_MAF_ampdel@clinical.data$bcr_patient_barcode[pancan_MAF_ampdel@clinical.data$subtype==subtype])
    
    z = intersect(coad_samples,colnames(mut3))
    z = intersect(z,colnames(coad_exp.t))
    
    
    Ks = ifelse(mut3["N6_MT_complex",z]   |  del3["N6_MT_complex",z],"N6_MT_complex",
                ifelse(mut3["TP53.tr",z] |  del3["TP53.tr",z],"TP53 tr",
                       ifelse(mut3["TP53.ms",z] | del3["TP53.ms",z],"TP53 ms" , "WT")))
    print(table(Ks))
    design1 <- model.matrix(~0+Ks)
    
    colnames(design1) = c("N6_MT_complex","TP53_ms","TP53_tr","WT")
    contr.matrix <- makeContrasts(
      N6vsWT = N6_MT_complex-WT, 
      TP53trvsWT = TP53_tr - WT,
      TP53msvsWT = TP53_ms - WT,
      levels = colnames(design1))
    
    
    
    vfit <- lmFit(coad_exp.t[,z], design1)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    
    N6.vs.WT <- topTable(efit, coef=1, n=Inf,confint = T)
    N6.vs.WT$model="N6vsWT"
    N6.vs.WT$subtype=subtype
    N6.vs.WT$genes=rownames(N6.vs.WT)
    allResults=rbind(allResults,N6.vs.WT)
    
    
    TP53tr.vs.WT <- topTable(efit, coef=2, n=Inf,confint = T)
    TP53tr.vs.WT$model="TP53trvsWT"
    TP53tr.vs.WT$subtype=subtype
    TP53tr.vs.WT$genes=rownames(TP53tr.vs.WT)
    allResults=rbind(allResults,TP53tr.vs.WT)
    
    TP53ms.vs.WT <- topTable(efit, coef=3, n=Inf,confint = T)
    TP53ms.vs.WT$model="TP53msvsWT"
    TP53ms.vs.WT$subtype=subtype
    TP53ms.vs.WT$genes=rownames(TP53ms.vs.WT)
    allResults=rbind(allResults,TP53ms.vs.WT)
    

  })
  
}
save(allResults,file="DE_TP53_pancanResults.RData")



library(ggplot2)
library(viridis)

#summary
library(plyr)
library(dplyr)
allResults$se = allResults$logFC-allResults$CI.L



load("/data/datasets/TCGA/TCGA-HNSC/RNAseq/GeneExp_gene_counts.RData")


targets2 = c( "ENSG00000124762","ENSG00000105327","ENSG00000026103",
              "ENSG00000135679","ENSG00000164938")

tumortypes = c("TCGA-LUAD","TCGA-BRCA","TCGA-OV","TCGA-HNSC","TCGA-UCEC","summary")

allResults3$symbol = data.exp@rowRanges$external_gene_name[match(allResults3$genes,data.exp@rowRanges$ensembl_gene_id)]

toPLot = allResults3 %>% filter(genes %in% targets2,model %in% c("N6vsWT","TP53trvsWT"),subtype %in% tumortypes) %>% group_by(subtype,model) %>% select(genes,logFC,CI.L , CI.R,symbol)


# Sup fig 6E
ggplot(toPLot,aes(x=1,y=logFC,ymin=CI.L,ymax=CI.R,col=model))+xlab("TP53 targets")+
  facet_grid(symbol~subtype,switch = "y")+
  geom_abline(slope=0,intercept = 0)+
  geom_pointrange(position=position_dodge2(width=1))+ 
  coord_flip(ylim=c(-1.5,0.3))+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom",axis.text.y = element_blank(),axis.ticks.y=element_blank())+
  scale_color_viridis("Comparison",labels=c("MMettl3 complex vs WT","TP53 trunc vs WT"),discrete = T)

ggsave("figS6E.pdf",width = 9,height = 5)


toPLot3 = allResults3 %>% filter(genes %in% targets,model %in% c("N6vsWT","TP53trvsWT"),subtype =="summary") %>% group_by(subtype,model) %>% select(genes,logFC,CI.L , CI.R,symbol)

# Sup fig 6D
ggplot(toPLot3,aes(x=forcats::fct_rev(symbol),y=logFC,ymin=CI.L,ymax=CI.R,col=model))+xlab("TP53 targets")+
  #  facet_grid(~subtype,switch = "y")+
  geom_abline(slope=0,intercept = 0)+
  geom_pointrange(position=position_dodge2(width=1))+ 
  coord_flip(ylim=c(-0.75,0.3))+theme_classic()+
  theme(legend.position="bottom",)+
  scale_color_viridis("Comparison",labels=c("MMettl3 complex vs WT","TP53 trunc vs WT"),discrete = T)
ggsave("figS6D.pdf",width = 9,height = 5)

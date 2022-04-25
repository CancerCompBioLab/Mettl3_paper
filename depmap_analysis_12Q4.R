# The Mettl3 epitranscriptomic writer amplifies p53 stress responses
# Nitin Raj et al.
# author: Jose A. Seoane joseaseoane@vhio.net

# depmap analyses
# Figures 7F, S7A, S7B, S7C


library(ggpubr)
library(ggrepel)

fsamples= "/data/datasets/CCLE/21Q4/sample_info.csv"
fach2 = "/data/datasets/CCLE/21Q4/CRISPR_gene_effect.csv"


fmut = "/data/datasets/CCLE/21Q4/CCLE_mutations.csv"


ach2 = read.table(fach2,header = T,sep = ",",quote = "",fill=T)


ach2.m = as.matrix(ach2[,-c(1)])
rownames(ach2.m)=samples$CCLE_Name[match(ach2$DepMap_ID,samples$DepMap_ID)]



mutData = read.table(fmut,header = T,sep = ",",stringsAsFactors = F,quote = "",fill=T)
mutData.nonSil = mutData[which(mutData$Variant_Classification!="Silent"),]

mutData.nonSil$isHotSpot = mutData.nonSil$isCOSMIChotspot=="True"|mutData.nonSil$isTCGAhotspot=="True"

fhotspot = function(x){
  any(x)
}


mudDataDam.mat = reshape2::dcast(mutData.nonSil,Hugo_Symbol~DepMap_ID,fun.aggregate = fhotspot,value.var="isHotSpot")
mudDataDam.mat2 =mudDataDam.mat[,-c(1)]
rownames(mudDataDam.mat2) = mudDataDam.mat$Hugo_Symbol
mudDataDam.mat3 = mudDataDam.mat2
colnames(mudDataDam.mat3)=samples$CCLE_Name[match(colnames(mudDataDam.mat2),samples$DepMap_ID)]


z7 = intersect(rownames(ach2.m),colnames(mudDataDam.mat3))
z7 = setdiff(z7,"")


mettl3_scaled = scale(as.numeric(ach2.m[z7,"METTL3..56339."]))
colorsViridis = viridis::viridis(3)

# figS7C
ggboxplot(data.frame(mettl3=mettl3_scaled,TP53=t(mudDataDam.mat3["TP53",z7]),stringsAsFactors = F),
          x="TP53",y="mettl3",color="TP53",add = "jitter",palette = colorsViridis[1:2])+stat_compare_means(label.x.npc = 0.5)+ylab("METTL3 achiles score")


ggsave("figS7C_depmap_achiles_mettl3_byTP53_21Q4.pdf")


wts_3 = apply(mudDataDam.mat3[,z7],1,function(x) tryCatch({wilcox.test(mettl3_scaled~x  )$p.value},error=function(e) 1) )

toplot_wtest = data.frame(pval=wts_3,symbol=names(wts_3),counts=rowSums(mudDataDam.mat3[,z7]),stringsAsFactors = F)
toplot_wtest$symbol2 = ifelse(toplot_wtest$symbol=="TP53","TP53","")
toplot_wtest$logpval = -log10(toplot_wtest$pval)
toplot_wtest$order = order(toplot_wtest$pval) 
toplot_wtest$order = factor(toplot_wtest$order,levels=toplot_wtest$order,ordered = T)



#figS7a

p1 = ggplot(toplot_wtest[which(toplot_wtest$pval<1),],aes(x=reorder(symbol,-logpval),y=logpval,label=symbol2))+
  geom_point()+geom_text_repel()+coord_cartesian(xlim=c(-20,4844))+theme_linedraw()+theme(panel.background = element_blank(),axis.title.x=element_blank(),
                                                                                          axis.text.x=element_blank(),
                                                                                          axis.ticks.x=element_blank(),
                                                                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("-log10 pval")

p1
ggsave("figS7A_1_wilcox_METTL3_cripr_gw_mut_21Q4.pdf",height = 5,width = 8)

p2 = 
  ggplot(toplot_wtest[which(toplot_wtest$pval<0.0262),],aes(x=reorder(symbol,-logpval),y=logpval,label=symbol2))+
  geom_bar(stat = "identity" )+coord_cartesian(xlim=c(0,30))+xlab("symbol")+theme(panel.background = element_blank(),
                                                                                  axis.text.x=element_text(angle = 45, hjust = 1))+ylab("-log10 pval")

p2
ggsave("figS7A_1_wilcox_METTL3_cripr_gw_mut_inner_21Q4.pdf",height = 3,width = 6)




cor_achil.sp = cor(ach2.m,use = "complete.obs",method = "sp")

cor_achil.sp.m3 = cor_achil.sp[,"METTL3..56339."] # 60
cor_achil.sp.m14 = cor_achil.sp[,"METTL14..57721."] # 60
cor_achil.sp.tp = cor_achil.sp[,"TP53..7157."] # 382
cor_achil.sp.tp.rbm15 = cor_achil.sp[,"RBM15..64783."] # 
cor_achil.sp.tp.wtap = cor_achil.sp[,"WTAP..9589."] # 
cor_achil.sp.tp.virma = cor_achil.sp[,"VIRMA..25962."] # 
cor_achil.sp.tp.zc3h13 = cor_achil.sp[,"ZC3H13..23091."] # 

toPlot.corr = data.frame(symbol=gsub("\\.\\..*$","",colnames(cor_achil.sp),fixed=F,perl=T),
                         mettl3_cor=cor_achil.sp.m3,
                         mettl14_cor=cor_achil.sp.m14,
                         rbm15_cor=cor_achil.sp.tp.rbm15,
                         wtap_cor=cor_achil.sp.tp.wtap,
                         virma_cor=cor_achil.sp.tp.virma,
                         zc3h13_cor=cor_achil.sp.tp.zc3h13,
                         tp53_cor = cor_achil.sp.tp,
                         stringsAsFactors = F)

toPlot.corr$symbol2 = ifelse(toPlot.corr$symbol %in% c("MDM2","METTL3","METTL14","WTAP","ZC3H13","RBM15","ATM",
                                                       "TP53","TP53BP1","CBLL1","RBM15B"),toPlot.corr$symbol,"")
toPlot.corr$symbol3 = ifelse(toPlot.corr$symbol %in% c("MDM2","MDM4","TP53","TP53BP1","ATM","RBM15","RBM15B","METTL3",
                                                       "METTL14","RB1","ZMAT3","WTAP","SP1",
                                                       "PTPN14","ZC3H13","CBLL1"),toPlot.corr$symbol,"")

# Fig 7F
qt = quantile(cor_achil.sp.m3,c(0.05,0.95))
ggplot(toPlot.corr,aes(x=mettl3_cor,y=..density..))+geom_density(color=colorsViridis[2],fill=colorsViridis[2],alpha=0.5)+xlim(c(-0.3,0.6))+
  theme_classic()+
  geom_vline(data=toPlot.corr[which(toPlot.corr$symbol2!=""),],mapping=aes(xintercept=mettl3_cor),linetype="dashed",color=colorsViridis[1])+
  geom_vline(xintercept = qt,color="red",linetype="dotted")+
  geom_text_repel(aes(x=mettl3_cor+0.00,y=4,label=symbol2),angle=0,position = position_jitter(width=0,height = 4) ,color=colorsViridis[1],max.overlaps = Inf)+
  xlab("Correlation with Achiles METTL3 effect size")

ggsave("fig7F_achiles_correlation_mettl3_v2_21Q4.pdf",height = 4,width = 7)

# Fig S7B
qt = quantile(cor_achil.sp.tp.wtap,c(0.05,0.95))
ggplot(toPlot.corr,aes(x=wtap_cor,y=..density..))+geom_density(color=colorsViridis[2],fill=colorsViridis[2],alpha=0.5)+xlim(c(-0.4,0.6))+
  theme_classic()+
  geom_vline(data=toPlot.corr[which(toPlot.corr$symbol2!=""),],mapping=aes(xintercept=wtap_cor),linetype="dashed",color=colorsViridis[1])+
  geom_vline(xintercept = qt,color="red",linetype="dotted")+
  geom_text_repel(aes(x=wtap_cor+0.00,y=4,label=symbol2),angle=0,position = position_jitter(width=0,height = 4) ,color=colorsViridis[1],max.overlaps = Inf)+
  xlab("Correlation with Achiles WTAP effect size")

ggsave("figS7B_achiles_correlation_WTAP_v2_21Q4.pdf",height = 4,width = 7)


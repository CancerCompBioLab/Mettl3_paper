# The Mettl3 epitranscriptomic writer amplifies p53 stress responses
# Nitin Raj et al.
# author: Jose A. Seoane joseaseoane@vhio.net

# TCGA PANCAN DE analyses
# Figures 7E, S6B, S6C



library(ComplexHeatmap)

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

rownames(muts)[which(rownames(muts)=="KIAA1429")]<-"VIRMA"


N6_MT_complex=apply(muts[c("METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13"),],2,any)
N6_MT_complex_core=apply(muts[c("METTL3","METTL14","WTAP"),],2,any)

muts2 = rbind(muts,N6_MT_complex,N6_MT_complex_core )


N6_MT_complex=apply(dels[c("METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13"),],2,any)
N6_MT_complex_core=apply(dels[c("METTL3","METTL14","WTAP"),],2,any)

del2 = rbind(dels,N6_MT_complex,N6_MT_complex_core )


allsamples = unique(c(colnames(muts2),colnames(del2)))
allgenes = unique(c(rownames(muts2),rownames(del2)))


mut3 = matrix(data=F,nrow=length(allgenes),ncol = length(allsamples),dimnames = list(allgenes,allsamples))
del3 = matrix(data=F,nrow=length(allgenes),ncol = length(allsamples),dimnames = list(allgenes,allsamples))
mut3[rownames(muts2),colnames(muts2)]=muts2
del3[rownames(del2),colnames(del2)]=del2




plot_cn_mut_v2 = function(subtype,label){
  
  coad_samples = as.character(pancan_MAF_ampdel@clinical.data$bcr_patient_barcode[pancan_MAF_ampdel@clinical.data$subtype==subtype])
  z = intersect(coad_samples,colnames(mut3))
  mat_list = list(mut = mut3[c("TP53.tr","TP53.ms","N6_MT_complex","METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13"),z],
                  del=del3[c("TP53.tr","TP53.ms","N6_MT_complex","METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13"),z])
  
  Ks = ifelse(mut3["N6_MT_complex",z] |
                #amp3["N6_MT_complex",colnames(mut3) %in% coad_samples]  |
                del3["N6_MT_complex",z],"N6_MT_complex",
              ifelse(mut3["TP53.tr",z]  ,"TP53 tr",
                     ifelse(mut3["TP53.ms",z]   ,"TP53 ms","WT")))
  col = c(mut = "#008000", amp = "red",del = "blue")
  v1=250/ncol(mat_list$mut)
  op1=oncoPrint(mat_list,alter_fun = list(
    #    background = function(x, y, w, h) {
    #      grid.rect(x, y, w-unit(v1, "mm"), h-unit(v1, "mm"), 
    #                gp = gpar(fill = "#CCCCCC", col = NA))
    #    },
    amp = function(x, y, w, h) grid.rect(x, y, w-unit(0.5,"mm"), h-unit(0.5,"mm"), 
                                         gp = gpar(fill = col["amp"], col = NA)),
    del = function(x, y, w, h) grid.rect(x, y, w-unit(0.5,"mm"), h-unit(0.5,"mm"), 
                                         gp = gpar(fill = col["del"], col = NA)),
    mut = function(x, y, w, h) grid.rect(x, y, w-unit(0.5,"mm"), h*0.33, 
                                         gp = gpar(fill = col["mut"], col = NA))
  ),col = col,top_annotation = NULL,right_annotation = NULL,column_split = Ks,
  row_order = c("N6_MT_complex","TP53.tr","TP53.ms","ZC3H13","CBLL1","METTL3","RBM15","RBM15B","VIRMA","METTL14","WTAP"))
  
  draw(op1,column_title=label)
  
}



# Figure 7E
pdf(file="fig7E.pdf",height = 4,width = 8)
# LUAD
plot_cn_mut_v2(subtype = "TCGA-LUAD","LUAD ms p=6.01E-6")
dev.off()


# Figure S6B
pdf(file="figS6B.pdf",height = 4,width = 8)
# BRCA
plot_cn_mut_v2(subtype = "TCGA-BRCA","BRCA tr p=1.86E-16, ms p=1.46E-21")
dev.off()

# Figure S6C
pdf(file="figS6C.pdf",height = 4,width = 8)
# HNSC
plot_cn_mut_v2(subtype = "TCGA-HNSC","HNSC tr p=5.26E-11, ms p=1.007E-8")
dev.off()







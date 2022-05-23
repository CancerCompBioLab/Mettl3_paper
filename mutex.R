# The Mettl3 epitranscriptomic writer amplifies p53 stress responses
# Nitin Raj et al.
# author: Jose A. Seoane joseaseoane@vhio.net

# TCGA PANCAN DE analyses
# Table S6A


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
rownames(muts)[which(rownames(muts)=="KIAA1429")]<-"VIRMA"

allsamples = unique(c(colnames(muts),colnames(dels)))
allgenes = unique(c(rownames(muts),rownames(dels)))

mut_del = matrix(data=F,nrow=length(allgenes),ncol = length(allsamples),dimnames = list(allgenes,allsamples))

mut_del[rownames(muts),colnames(muts)]=muts
mut_del[rownames(dels),colnames(dels)]=mut_del[rownames(dels),colnames(dels)]|dels

N6_MT_complex=apply(mut_del[c("METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13"),],2,any)
N6_MT_complex_core=apply(mut_del[c("METTL3","METTL14","WTAP"),],2,any)

mut_del = rbind(mut_del,N6_MT_complex,N6_MT_complex_core )


N6_MT_complex=apply(muts[c("METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13"),],2,any)
N6_MT_complex_core=apply(muts[c("METTL3","METTL14","WTAP"),],2,any)
muts = rbind(muts,N6_MT_complex,N6_MT_complex_core )
save(muts,mut_del,file="boolean_mut_amp_del.RData")

library(discover)
tissue = as.character(pancan_MAF_ampdel@clinical.data$subtype)
names(tissue)=as.character(pancan_MAF_ampdel@clinical.data$bcr_patient_barcode)

events_del = discover.matrix(mut_del,strata = tissue[match(colnames(mut_del),names(tissue))] )
save(events_del,tissue,file="events_mut_del_v4.RData")

events_mut = discover.matrix(muts,strata = tissue[match(colnames(muts),names(tissue))] )
save(events_mut,file="events_mut_v4.RData")

load("events_mut_del_v4.RData")

load("events_mut_v4.RData")


mutex_11 = pairwise.discover.test(events_mut[which(rownames(events_mut) %in% c("N6_MT_complex","N6_MT_complex_core",
                                                                               "METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13",
                                                                               "TP53","TP53.tr")), ],alternative = "less")




# mutex analysis mut only


mutex_12 = pairwise.discover.test(events_mut[which(rownames(events_mut) %in% c("N6_MT_complex","N6_MT_complex_core",
                                                                               "METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13",
                                                                               "TP53","TP53.tr")), ],alternative = "less")



gois = c("N6_MT_complex","N6_MT_complex_core",
         "METTL3","METTL14","WTAP","CBLL1","RBM15","RBM15B","VIRMA","ZC3H13",
         "TP53","TP53.ms","TP53.tr")
set.seed(1980)
gois_sample = unique(c(gois,sample(rownames(events_mut),1000)))

res=NULL
for(i in 1:length(unique(tissue))){
  tstmp = unique(tissue)[i]
  try({
    
    mutex_12_t = pairwise.discover.test(events_mut[which(rownames(events_mut) %in% gois), which(tissue[match(colnames(events_mut),names(tissue))]==tstmp)  ],alternative = "less")
    print(tstmp)
    print(as.data.frame(mutex_12_t,q.threshold = 1))
    restmp = as.data.frame(mutex_12_t,q.threshold = 1)
    restmp$type = tstmp
    res = rbind(res,restmp)
  })
}


library(xlsx)

## fig S6A
write.xlsx2(res[(res$q.value<0.1 & res$gene1 %in% gois & res$gene2 %in% gois),],file="TP53_N6_discover_v5.xlsx",sheetName = "mut_only",row.names = F,append = F)

write.xlsx2(res,file="TP53_N6_discover_v5_allComb.xlsx",sheetName = "mut_only",row.names = F,append = F)


res_cn=NULL
for(i in 1:length(unique(tissue))){
  tstmp = unique(tissue)[i]
  try({

    mutex_12_t = pairwise.discover.test(events_del[which(rownames(events_del) %in% gois), which(tissue[match(colnames(events_del),names(tissue))]==tstmp)  ],alternative = "less")
    print(tstmp)
    print(as.data.frame(mutex_12_t,q.threshold = 1))
    restmp = as.data.frame(mutex_12_t,q.threshold = 1)
    restmp$type = tstmp
    res_cn = rbind(res_cn,restmp)
  })
}
write.xlsx2(res_cn[(res_cn$q.value<0.1 & res_cn$gene1 %in% gois & res_cn$gene2 %in% gois),],file="TP53_N6_discover_v5.xlsx",sheetName = "mut_del",row.names = F,append = T)

write.xlsx2(res_cn,file="TP53_N6_discover_v5_allComb.xlsx",sheetName = "mut_del",row.names = F,append = T)


write.xlsx2(as.data.frame(mutex_12),file="TP53_N6_discover_v5.xlsx",sheetName = "all_types_mut",row.names = F,append = T)

write.xlsx2(as.data.frame(mutex_11),file="TP53_N6_discover_v5.xlsx",sheetName = "all_types_mut_cn",row.names = F,append = T)




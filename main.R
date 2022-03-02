library(Seurat)

# CD45pTumor
fileNames <- list.files("HNSCC/data/CD45pTumor/")
seurat.list <- list()
for(i in 1:length(fileNames)){
  seurat.list[[i]] <- Read10X(paste0("HNSCC/data/CD45pTumor/",fileNames[i],"/"))
  seurat.list[[i]] <- CreateSeuratObject(counts = seurat.list[[i]],min.cells = 3)
}
names(seurat.list) <- fileNames

seurat.list <- lapply(seurat.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
seurat.list$HN01$Sample <- "HN01"
seurat.list$HN02$Sample <- "HN02"
seurat.list$HN03$Sample <- "HN03"
seurat.list$HN04$Sample <- "HN04"
seurat.list$HN05$Sample <- "HN05"
seurat.list$HN06$Sample <- "HN06"
seurat.list$HN07$Sample <- "HN07"
seurat.list$HN08$Sample <- "HN08"
seurat.list$HN09$Sample <- "HN09"
seurat.list$HN10$Sample <- "HN10"
seurat.list$HN11$Sample <- "HN11"
seurat.list$HN12$Sample <- "HN12"
seurat.list$HN13$Sample <- "HN13"
seurat.list$HN14$Sample <- "HN14"
seurat.list$HN15$Sample <- "HN15"
seurat.list$HN16$Sample <- "HN16"
seurat.list$HN17$Sample <- "HN17"
seurat.list$HN18$Sample <- "HN18"

# add clinical information
# HN01
seurat.list$HN01$Gender <- "M"
seurat.list$HN01$AgeGroup <- "70-79"
seurat.list$HN01$Smoking <- "Yes"
seurat.list$HN01$Alcohol <- "No"
seurat.list$HN01$DiseaseSite <- "Oral cavity"
seurat.list$HN01$T_Stage <- "T4A"
seurat.list$HN01$N_Stage <- "N2B"
seurat.list$HN01$M_Stage <- "M0"
seurat.list$HN01$HPV <- "Neg"
seurat.list$HN01$Inflam.status <- "High"
# HN02
seurat.list$HN02$Gender <- "F"
seurat.list$HN02$AgeGroup <- "60-69"
seurat.list$HN02$Smoking <- "No"
seurat.list$HN02$Alcohol <- "No"
seurat.list$HN02$DiseaseSite <- "Oral cavity"
seurat.list$HN02$T_Stage <- "T3"
seurat.list$HN02$N_Stage <- "N2a"
seurat.list$HN02$M_Stage <- "M0"
seurat.list$HN02$HPV <- "Neg"
seurat.list$HN02$Inflam.status <- "Low"
# HN03
seurat.list$HN03$Gender <- "M"
seurat.list$HN03$AgeGroup <- "80-89"
seurat.list$HN03$Smoking <- "No"
seurat.list$HN03$Alcohol <- "No"
seurat.list$HN03$DiseaseSite <- "Oral cavity"
seurat.list$HN03$T_Stage <- "T4a"
seurat.list$HN03$N_Stage <- "N0"
seurat.list$HN03$M_Stage <- "M0"
seurat.list$HN03$HPV <- "Neg"
seurat.list$HN03$Inflam.status <- "NA"
# HN04
seurat.list$HN04$Gender <- "M"
seurat.list$HN04$AgeGroup <- "50-59"
seurat.list$HN04$Smoking <- "Yes"
seurat.list$HN04$Alcohol <- "Yes"
seurat.list$HN04$DiseaseSite <- "Oral cavity"
seurat.list$HN04$T_Stage <- "T3"
seurat.list$HN04$N_Stage <- "N1"
seurat.list$HN04$M_Stage <- "M0"
seurat.list$HN04$HPV <- "Neg"
seurat.list$HN04$Inflam.status <- "Low"
# HN05
seurat.list$HN05$Gender <- "F"
seurat.list$HN05$AgeGroup <- "50-59"
seurat.list$HN05$Smoking <- "Yes"
seurat.list$HN05$Alcohol <- "Yes"
seurat.list$HN05$DiseaseSite <- "Oral cavity"
seurat.list$HN05$T_Stage <- "†T3"
seurat.list$HN05$N_Stage <- "†N3b"
seurat.list$HN05$M_Stage <- "†M0"
seurat.list$HN05$HPV <- "Neg"
seurat.list$HN05$Inflam.status <- "Med"
# HN06
seurat.list$HN06$Gender <- "M"
seurat.list$HN06$AgeGroup <- "30-39"
seurat.list$HN06$Smoking <- "Yes"
seurat.list$HN06$Alcohol <- "Yes"
seurat.list$HN06$DiseaseSite <- "Oral cavity"
seurat.list$HN06$T_Stage <- "T3"
seurat.list$HN06$N_Stage <- "N0"
seurat.list$HN06$M_Stage <- "M0"
seurat.list$HN06$HPV <- "Neg"
seurat.list$HN06$Inflam.status <- "High"
# HN07
seurat.list$HN07$Gender <- "F"
seurat.list$HN07$AgeGroup <- "60-69"
seurat.list$HN07$Smoking <- "Yes"
seurat.list$HN07$Alcohol <- "Yes"
seurat.list$HN07$DiseaseSite <- "Larynx"
seurat.list$HN07$T_Stage <- "T3"
seurat.list$HN07$N_Stage <- "N0"
seurat.list$HN07$M_Stage <- "M0"
seurat.list$HN07$HPV <- "Neg"
seurat.list$HN07$Inflam.status <- "Low"
# HN08
seurat.list$HN08$Gender <- "F"
seurat.list$HN08$AgeGroup <- "70-79"
seurat.list$HN08$Smoking <- "Yes"
seurat.list$HN08$Alcohol <- "Yes"
seurat.list$HN08$DiseaseSite <- "Oral cavity"
seurat.list$HN08$T_Stage <- "†T1"
seurat.list$HN08$N_Stage <- "†N0"
seurat.list$HN08$M_Stage <- "M0"
seurat.list$HN08$HPV <- "Neg"
seurat.list$HN08$Inflam.status <- "Med"
# HN09
seurat.list$HN09$Gender <- "F"
seurat.list$HN09$AgeGroup <- "70-79"
seurat.list$HN09$Smoking <- "Yes"
seurat.list$HN09$Alcohol <- "Yes"
seurat.list$HN09$DiseaseSite <- "Oral cavity"
seurat.list$HN09$T_Stage <- "T3"
seurat.list$HN09$N_Stage <- "N2B"
seurat.list$HN09$M_Stage <- "M0"
seurat.list$HN09$HPV <- "Neg"
seurat.list$HN09$Inflam.status <- "Med"
# HN10
seurat.list$HN10$Gender <- "F"
seurat.list$HN10$AgeGroup <- "50-59"
seurat.list$HN10$Smoking <- "No"
seurat.list$HN10$Alcohol <- "Yes"
seurat.list$HN10$DiseaseSite <- "Oral cavity"
seurat.list$HN10$T_Stage <- "T3"
seurat.list$HN10$N_Stage <- "N0"
seurat.list$HN10$M_Stage <- "M0"
seurat.list$HN10$HPV <- "Neg"
seurat.list$HN10$Inflam.status <- "High"
# HN11
seurat.list$HN11$Gender <- "M"
seurat.list$HN11$AgeGroup <- "80-89"
seurat.list$HN11$Smoking <- "No"
seurat.list$HN11$Alcohol <- "No"
seurat.list$HN11$DiseaseSite <- "Oral cavity"
seurat.list$HN11$T_Stage <- "T2"
seurat.list$HN11$N_Stage <- "N0"
seurat.list$HN11$M_Stage <- "M0"
seurat.list$HN11$HPV <- "Neg"
seurat.list$HN11$Inflam.status <- "Low"
# HN12
seurat.list$HN12$Gender <- "M"
seurat.list$HN12$AgeGroup <- "50-59"
seurat.list$HN12$Smoking <- "Yes"
seurat.list$HN12$Alcohol <- "Yes"
seurat.list$HN12$DiseaseSite <- "Oropharynx"
seurat.list$HN12$T_Stage <- "T2"
seurat.list$HN12$N_Stage <- "N1"
seurat.list$HN12$M_Stage <- "M0"
seurat.list$HN12$HPV <- "Pos"
seurat.list$HN12$Inflam.status <- "Med"
# HN13
seurat.list$HN13$Gender <- "M"
seurat.list$HN13$AgeGroup <- "70-79"
seurat.list$HN13$Smoking <- "No"
seurat.list$HN13$Alcohol <- "No"
seurat.list$HN13$DiseaseSite <- "Oropharynx"
seurat.list$HN13$T_Stage <- "T2"
seurat.list$HN13$N_Stage <- "N0"
seurat.list$HN13$M_Stage <- "M0"
seurat.list$HN13$HPV <- "Pos"
seurat.list$HN13$Inflam.status <- "High"
# HN14
seurat.list$HN14$Gender <- "M"
seurat.list$HN14$AgeGroup <- "50-59"
seurat.list$HN14$Smoking <- "Yes"
seurat.list$HN14$Alcohol <- "Yes"
seurat.list$HN14$DiseaseSite <- "Oropharynx"
seurat.list$HN14$T_Stage <- "T1"
seurat.list$HN14$N_Stage <- "N1"
seurat.list$HN14$M_Stage <- "M0"
seurat.list$HN14$HPV <- "Pos"
seurat.list$HN14$Inflam.status <- "High"
# HN15
seurat.list$HN15$Gender <- "F"
seurat.list$HN15$AgeGroup <- "60-69"
seurat.list$HN15$Smoking <- "Yes"
seurat.list$HN15$Alcohol <- "NA"
seurat.list$HN15$DiseaseSite <- "Oral cavity"
seurat.list$HN15$T_Stage <- "T2"
seurat.list$HN15$N_Stage <- "N0"
seurat.list$HN15$M_Stage <- "M0"
seurat.list$HN15$HPV <- "Neg"
seurat.list$HN15$Inflam.status <- "Low"
# HN16
seurat.list$HN16$Gender <- "M"
seurat.list$HN16$AgeGroup <- "40-49"
seurat.list$HN16$Smoking <- "Yes"
seurat.list$HN16$Alcohol <- "Yes"
seurat.list$HN16$DiseaseSite <- "Oropharynx"
seurat.list$HN16$T_Stage <- "T2"
seurat.list$HN16$N_Stage <- "N1"
seurat.list$HN16$M_Stage <- "M0"
seurat.list$HN16$HPV <- "Pos"
seurat.list$HN16$Inflam.status <- "Med"
# HN17
seurat.list$HN17$Gender <- "M"
seurat.list$HN17$AgeGroup <- "50-59"
seurat.list$HN17$Smoking <- "Yes"
seurat.list$HN17$Alcohol <- "Yes"
seurat.list$HN17$DiseaseSite <- "Oropharynx"
seurat.list$HN17$T_Stage <- "T1"
seurat.list$HN17$N_Stage <- "N1"
seurat.list$HN17$M_Stage <- "M0"
seurat.list$HN17$HPV <- "Pos"
seurat.list$HN17$Inflam.status <- "High"
# HN18
seurat.list$HN18$Gender <- "M"
seurat.list$HN18$AgeGroup <- "50-59"
seurat.list$HN18$Smoking <- "Yes"
seurat.list$HN18$Alcohol <- "Yes"
seurat.list$HN18$DiseaseSite <- "Oropharynx"
seurat.list$HN18$T_Stage <- "T2"
seurat.list$HN18$N_Stage <- "N2"
seurat.list$HN18$M_Stage <- "M0"
seurat.list$HN18$HPV <- "Pos"
seurat.list$HN18$Inflam.status <- "Low"

anchors <- FindIntegrationAnchors(object.list = seurat.list)
saveRDS(anchors,file = "HNSCC/rds/CD45pTumor.anchors.rds")
integrated <- IntegrateData(anchorset = anchors)
integrated <- ScaleData(integrated, verbose = FALSE)
saveRDS(integrated, file = "HNSCC/rds/CD45pTumor.integrated.rds")

rm(list = ls());gc()

# CD45nTumor
fileNames <- list.files("HNSCC/data/CD45nTumor/")
seurat.list <- list()
for(i in 1:length(fileNames)){
  seurat.list[[i]] <- Read10X(paste0("HNSCC/data/CD45nTumor/",fileNames[i],"/"))
  seurat.list[[i]] <- CreateSeuratObject(counts = seurat.list[[i]],min.cells = 3)
}
names(seurat.list) <- fileNames

seurat.list <- lapply(seurat.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


seurat.list$HN01$Sample <- "HN01"
seurat.list$HN05$Sample <- "HN05"
seurat.list$HN06$Sample <- "HN06"
seurat.list$HN07$Sample <- "HN07"
seurat.list$HN08$Sample <- "HN08"
seurat.list$HN09$Sample <- "HN09"
seurat.list$HN10$Sample <- "HN10"
seurat.list$HN11$Sample <- "HN11"
seurat.list$HN12$Sample <- "HN12"
seurat.list$HN13$Sample <- "HN13"
seurat.list$HN14$Sample <- "HN14"
seurat.list$HN15$Sample <- "HN15"
seurat.list$HN16$Sample <- "HN16"
seurat.list$HN17$Sample <- "HN17"
seurat.list$HN18$Sample <- "HN18"

# add clinical information
# HN01
seurat.list$HN01$Gender <- "M"
seurat.list$HN01$AgeGroup <- "70-79"
seurat.list$HN01$Smoking <- "Yes"
seurat.list$HN01$Alcohol <- "No"
seurat.list$HN01$DiseaseSite <- "Oral cavity"
seurat.list$HN01$T_Stage <- "T4A"
seurat.list$HN01$N_Stage <- "N2B"
seurat.list$HN01$M_Stage <- "M0"
seurat.list$HN01$HPV <- "Neg"
seurat.list$HN01$Inflam.status <- "High"
# HN05
seurat.list$HN05$Gender <- "F"
seurat.list$HN05$AgeGroup <- "50-59"
seurat.list$HN05$Smoking <- "Yes"
seurat.list$HN05$Alcohol <- "Yes"
seurat.list$HN05$DiseaseSite <- "Oral cavity"
seurat.list$HN05$T_Stage <- "†T3"
seurat.list$HN05$N_Stage <- "†N3b"
seurat.list$HN05$M_Stage <- "†M0"
seurat.list$HN05$HPV <- "Neg"
seurat.list$HN05$Inflam.status <- "Med"
# HN06
seurat.list$HN06$Gender <- "M"
seurat.list$HN06$AgeGroup <- "30-39"
seurat.list$HN06$Smoking <- "Yes"
seurat.list$HN06$Alcohol <- "Yes"
seurat.list$HN06$DiseaseSite <- "Oral cavity"
seurat.list$HN06$T_Stage <- "T3"
seurat.list$HN06$N_Stage <- "N0"
seurat.list$HN06$M_Stage <- "M0"
seurat.list$HN06$HPV <- "Neg"
seurat.list$HN06$Inflam.status <- "High"
# HN07
seurat.list$HN07$Gender <- "F"
seurat.list$HN07$AgeGroup <- "60-69"
seurat.list$HN07$Smoking <- "Yes"
seurat.list$HN07$Alcohol <- "Yes"
seurat.list$HN07$DiseaseSite <- "Larynx"
seurat.list$HN07$T_Stage <- "T3"
seurat.list$HN07$N_Stage <- "N0"
seurat.list$HN07$M_Stage <- "M0"
seurat.list$HN07$HPV <- "Neg"
seurat.list$HN07$Inflam.status <- "Low"
# HN08
seurat.list$HN08$Gender <- "F"
seurat.list$HN08$AgeGroup <- "70-79"
seurat.list$HN08$Smoking <- "Yes"
seurat.list$HN08$Alcohol <- "Yes"
seurat.list$HN08$DiseaseSite <- "Oral cavity"
seurat.list$HN08$T_Stage <- "†T1"
seurat.list$HN08$N_Stage <- "†N0"
seurat.list$HN08$M_Stage <- "M0"
seurat.list$HN08$HPV <- "Neg"
seurat.list$HN08$Inflam.status <- "Med"
# HN09
seurat.list$HN09$Gender <- "F"
seurat.list$HN09$AgeGroup <- "70-79"
seurat.list$HN09$Smoking <- "Yes"
seurat.list$HN09$Alcohol <- "Yes"
seurat.list$HN09$DiseaseSite <- "Oral cavity"
seurat.list$HN09$T_Stage <- "T3"
seurat.list$HN09$N_Stage <- "N2B"
seurat.list$HN09$M_Stage <- "M0"
seurat.list$HN09$HPV <- "Neg"
seurat.list$HN09$Inflam.status <- "Med"
# HN10
seurat.list$HN10$Gender <- "F"
seurat.list$HN10$AgeGroup <- "50-59"
seurat.list$HN10$Smoking <- "No"
seurat.list$HN10$Alcohol <- "Yes"
seurat.list$HN10$DiseaseSite <- "Oral cavity"
seurat.list$HN10$T_Stage <- "T3"
seurat.list$HN10$N_Stage <- "N0"
seurat.list$HN10$M_Stage <- "M0"
seurat.list$HN10$HPV <- "Neg"
seurat.list$HN10$Inflam.status <- "High"
# HN11
seurat.list$HN11$Gender <- "M"
seurat.list$HN11$AgeGroup <- "80-89"
seurat.list$HN11$Smoking <- "No"
seurat.list$HN11$Alcohol <- "No"
seurat.list$HN11$DiseaseSite <- "Oral cavity"
seurat.list$HN11$T_Stage <- "T2"
seurat.list$HN11$N_Stage <- "N0"
seurat.list$HN11$M_Stage <- "M0"
seurat.list$HN11$HPV <- "Neg"
seurat.list$HN11$Inflam.status <- "Low"
# HN12
seurat.list$HN12$Gender <- "M"
seurat.list$HN12$AgeGroup <- "50-59"
seurat.list$HN12$Smoking <- "Yes"
seurat.list$HN12$Alcohol <- "Yes"
seurat.list$HN12$DiseaseSite <- "Oropharynx"
seurat.list$HN12$T_Stage <- "T2"
seurat.list$HN12$N_Stage <- "N1"
seurat.list$HN12$M_Stage <- "M0"
seurat.list$HN12$HPV <- "Pos"
seurat.list$HN12$Inflam.status <- "Med"
# HN13
seurat.list$HN13$Gender <- "M"
seurat.list$HN13$AgeGroup <- "70-79"
seurat.list$HN13$Smoking <- "No"
seurat.list$HN13$Alcohol <- "No"
seurat.list$HN13$DiseaseSite <- "Oropharynx"
seurat.list$HN13$T_Stage <- "T2"
seurat.list$HN13$N_Stage <- "N0"
seurat.list$HN13$M_Stage <- "M0"
seurat.list$HN13$HPV <- "Pos"
seurat.list$HN13$Inflam.status <- "High"
# HN14
seurat.list$HN14$Gender <- "M"
seurat.list$HN14$AgeGroup <- "50-59"
seurat.list$HN14$Smoking <- "Yes"
seurat.list$HN14$Alcohol <- "Yes"
seurat.list$HN14$DiseaseSite <- "Oropharynx"
seurat.list$HN14$T_Stage <- "T1"
seurat.list$HN14$N_Stage <- "N1"
seurat.list$HN14$M_Stage <- "M0"
seurat.list$HN14$HPV <- "Pos"
seurat.list$HN14$Inflam.status <- "High"
# HN15
seurat.list$HN15$Gender <- "F"
seurat.list$HN15$AgeGroup <- "60-69"
seurat.list$HN15$Smoking <- "Yes"
seurat.list$HN15$Alcohol <- "NA"
seurat.list$HN15$DiseaseSite <- "Oral cavity"
seurat.list$HN15$T_Stage <- "T2"
seurat.list$HN15$N_Stage <- "N0"
seurat.list$HN15$M_Stage <- "M0"
seurat.list$HN15$HPV <- "Neg"
seurat.list$HN15$Inflam.status <- "Low"
# HN16
seurat.list$HN16$Gender <- "M"
seurat.list$HN16$AgeGroup <- "40-49"
seurat.list$HN16$Smoking <- "Yes"
seurat.list$HN16$Alcohol <- "Yes"
seurat.list$HN16$DiseaseSite <- "Oropharynx"
seurat.list$HN16$T_Stage <- "T2"
seurat.list$HN16$N_Stage <- "N1"
seurat.list$HN16$M_Stage <- "M0"
seurat.list$HN16$HPV <- "Pos"
seurat.list$HN16$Inflam.status <- "Med"
# HN17
seurat.list$HN17$Gender <- "M"
seurat.list$HN17$AgeGroup <- "50-59"
seurat.list$HN17$Smoking <- "Yes"
seurat.list$HN17$Alcohol <- "Yes"
seurat.list$HN17$DiseaseSite <- "Oropharynx"
seurat.list$HN17$T_Stage <- "T1"
seurat.list$HN17$N_Stage <- "N1"
seurat.list$HN17$M_Stage <- "M0"
seurat.list$HN17$HPV <- "Pos"
seurat.list$HN17$Inflam.status <- "High"
# HN18
seurat.list$HN18$Gender <- "M"
seurat.list$HN18$AgeGroup <- "50-59"
seurat.list$HN18$Smoking <- "Yes"
seurat.list$HN18$Alcohol <- "Yes"
seurat.list$HN18$DiseaseSite <- "Oropharynx"
seurat.list$HN18$T_Stage <- "T2"
seurat.list$HN18$N_Stage <- "N2"
seurat.list$HN18$M_Stage <- "M0"
seurat.list$HN18$HPV <- "Pos"
seurat.list$HN18$Inflam.status <- "Low"

anchors <- FindIntegrationAnchors(object.list = seurat.list)
saveRDS(anchors,file = "HNSCC/rds/CD45nTumor.anchors.rds")
integrated <- IntegrateData(anchorset = anchors)
integrated <- ScaleData(integrated, verbose = FALSE)
saveRDS(integrated, file = "HNSCC/rds/CD45nTumor.integrated.rds")


# CD45pPBMC
fileNames <- list.files("HNSCC/data/CD45pPBMC/")
seurat.list <- list()
for(i in 1:length(fileNames)){
  seurat.list[[i]] <- Read10X(paste0("HNSCC/data/CD45pPBMC/",fileNames[i],"/"))
  seurat.list[[i]] <- CreateSeuratObject(counts = seurat.list[[i]],min.cells = 3)
}
names(seurat.list) <- fileNames

seurat.list <- lapply(seurat.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
seurat.list$HN01$Sample <- "HN01"
seurat.list$HN02$Sample <- "HN02"
seurat.list$HN03$Sample <- "HN03"
seurat.list$HN04$Sample <- "HN04"
seurat.list$HN05$Sample <- "HN05"
seurat.list$HN06$Sample <- "HN06"
seurat.list$HN07$Sample <- "HN07"
seurat.list$HN08$Sample <- "HN08"
seurat.list$HN09$Sample <- "HN09"
seurat.list$HN11$Sample <- "HN11"
seurat.list$HN12$Sample <- "HN12"
seurat.list$HN13$Sample <- "HN13"
seurat.list$HN14$Sample <- "HN14"
seurat.list$HN15$Sample <- "HN15"
seurat.list$HN16$Sample <- "HN16"
seurat.list$HN17$Sample <- "HN17"
seurat.list$HN18$Sample <- "HN18"

# add clinical information
# HN01
seurat.list$HN01$Gender <- "M"
seurat.list$HN01$AgeGroup <- "70-79"
seurat.list$HN01$Smoking <- "Yes"
seurat.list$HN01$Alcohol <- "No"
seurat.list$HN01$DiseaseSite <- "Oral cavity"
seurat.list$HN01$T_Stage <- "T4A"
seurat.list$HN01$N_Stage <- "N2B"
seurat.list$HN01$M_Stage <- "M0"
seurat.list$HN01$HPV <- "Neg"
seurat.list$HN01$Inflam.status <- "High"
# HN02
seurat.list$HN02$Gender <- "F"
seurat.list$HN02$AgeGroup <- "60-69"
seurat.list$HN02$Smoking <- "No"
seurat.list$HN02$Alcohol <- "No"
seurat.list$HN02$DiseaseSite <- "Oral cavity"
seurat.list$HN02$T_Stage <- "T3"
seurat.list$HN02$N_Stage <- "N2a"
seurat.list$HN02$M_Stage <- "M0"
seurat.list$HN02$HPV <- "Neg"
seurat.list$HN02$Inflam.status <- "Low"
# HN03
seurat.list$HN03$Gender <- "M"
seurat.list$HN03$AgeGroup <- "80-89"
seurat.list$HN03$Smoking <- "No"
seurat.list$HN03$Alcohol <- "No"
seurat.list$HN03$DiseaseSite <- "Oral cavity"
seurat.list$HN03$T_Stage <- "T4a"
seurat.list$HN03$N_Stage <- "N0"
seurat.list$HN03$M_Stage <- "M0"
seurat.list$HN03$HPV <- "Neg"
seurat.list$HN03$Inflam.status <- "NA"
# HN04
seurat.list$HN04$Gender <- "M"
seurat.list$HN04$AgeGroup <- "50-59"
seurat.list$HN04$Smoking <- "Yes"
seurat.list$HN04$Alcohol <- "Yes"
seurat.list$HN04$DiseaseSite <- "Oral cavity"
seurat.list$HN04$T_Stage <- "T3"
seurat.list$HN04$N_Stage <- "N1"
seurat.list$HN04$M_Stage <- "M0"
seurat.list$HN04$HPV <- "Neg"
seurat.list$HN04$Inflam.status <- "Low"
# HN05
seurat.list$HN05$Gender <- "F"
seurat.list$HN05$AgeGroup <- "50-59"
seurat.list$HN05$Smoking <- "Yes"
seurat.list$HN05$Alcohol <- "Yes"
seurat.list$HN05$DiseaseSite <- "Oral cavity"
seurat.list$HN05$T_Stage <- "†T3"
seurat.list$HN05$N_Stage <- "†N3b"
seurat.list$HN05$M_Stage <- "†M0"
seurat.list$HN05$HPV <- "Neg"
seurat.list$HN05$Inflam.status <- "Med"
# HN06
seurat.list$HN06$Gender <- "M"
seurat.list$HN06$AgeGroup <- "30-39"
seurat.list$HN06$Smoking <- "Yes"
seurat.list$HN06$Alcohol <- "Yes"
seurat.list$HN06$DiseaseSite <- "Oral cavity"
seurat.list$HN06$T_Stage <- "T3"
seurat.list$HN06$N_Stage <- "N0"
seurat.list$HN06$M_Stage <- "M0"
seurat.list$HN06$HPV <- "Neg"
seurat.list$HN06$Inflam.status <- "High"
# HN07
seurat.list$HN07$Gender <- "F"
seurat.list$HN07$AgeGroup <- "60-69"
seurat.list$HN07$Smoking <- "Yes"
seurat.list$HN07$Alcohol <- "Yes"
seurat.list$HN07$DiseaseSite <- "Larynx"
seurat.list$HN07$T_Stage <- "T3"
seurat.list$HN07$N_Stage <- "N0"
seurat.list$HN07$M_Stage <- "M0"
seurat.list$HN07$HPV <- "Neg"
seurat.list$HN07$Inflam.status <- "Low"
# HN08
seurat.list$HN08$Gender <- "F"
seurat.list$HN08$AgeGroup <- "70-79"
seurat.list$HN08$Smoking <- "Yes"
seurat.list$HN08$Alcohol <- "Yes"
seurat.list$HN08$DiseaseSite <- "Oral cavity"
seurat.list$HN08$T_Stage <- "†T1"
seurat.list$HN08$N_Stage <- "†N0"
seurat.list$HN08$M_Stage <- "M0"
seurat.list$HN08$HPV <- "Neg"
seurat.list$HN08$Inflam.status <- "Med"
# HN09
seurat.list$HN09$Gender <- "F"
seurat.list$HN09$AgeGroup <- "70-79"
seurat.list$HN09$Smoking <- "Yes"
seurat.list$HN09$Alcohol <- "Yes"
seurat.list$HN09$DiseaseSite <- "Oral cavity"
seurat.list$HN09$T_Stage <- "T3"
seurat.list$HN09$N_Stage <- "N2B"
seurat.list$HN09$M_Stage <- "M0"
seurat.list$HN09$HPV <- "Neg"
seurat.list$HN09$Inflam.status <- "Med"
# HN11
seurat.list$HN11$Gender <- "M"
seurat.list$HN11$AgeGroup <- "80-89"
seurat.list$HN11$Smoking <- "No"
seurat.list$HN11$Alcohol <- "No"
seurat.list$HN11$DiseaseSite <- "Oral cavity"
seurat.list$HN11$T_Stage <- "T2"
seurat.list$HN11$N_Stage <- "N0"
seurat.list$HN11$M_Stage <- "M0"
seurat.list$HN11$HPV <- "Neg"
seurat.list$HN11$Inflam.status <- "Low"
# HN12
seurat.list$HN12$Gender <- "M"
seurat.list$HN12$AgeGroup <- "50-59"
seurat.list$HN12$Smoking <- "Yes"
seurat.list$HN12$Alcohol <- "Yes"
seurat.list$HN12$DiseaseSite <- "Oropharynx"
seurat.list$HN12$T_Stage <- "T2"
seurat.list$HN12$N_Stage <- "N1"
seurat.list$HN12$M_Stage <- "M0"
seurat.list$HN12$HPV <- "Pos"
seurat.list$HN12$Inflam.status <- "Med"
# HN13
seurat.list$HN13$Gender <- "M"
seurat.list$HN13$AgeGroup <- "70-79"
seurat.list$HN13$Smoking <- "No"
seurat.list$HN13$Alcohol <- "No"
seurat.list$HN13$DiseaseSite <- "Oropharynx"
seurat.list$HN13$T_Stage <- "T2"
seurat.list$HN13$N_Stage <- "N0"
seurat.list$HN13$M_Stage <- "M0"
seurat.list$HN13$HPV <- "Pos"
seurat.list$HN13$Inflam.status <- "High"
# HN14
seurat.list$HN14$Gender <- "M"
seurat.list$HN14$AgeGroup <- "50-59"
seurat.list$HN14$Smoking <- "Yes"
seurat.list$HN14$Alcohol <- "Yes"
seurat.list$HN14$DiseaseSite <- "Oropharynx"
seurat.list$HN14$T_Stage <- "T1"
seurat.list$HN14$N_Stage <- "N1"
seurat.list$HN14$M_Stage <- "M0"
seurat.list$HN14$HPV <- "Pos"
seurat.list$HN14$Inflam.status <- "High"
# HN15
seurat.list$HN15$Gender <- "F"
seurat.list$HN15$AgeGroup <- "60-69"
seurat.list$HN15$Smoking <- "Yes"
seurat.list$HN15$Alcohol <- "NA"
seurat.list$HN15$DiseaseSite <- "Oral cavity"
seurat.list$HN15$T_Stage <- "T2"
seurat.list$HN15$N_Stage <- "N0"
seurat.list$HN15$M_Stage <- "M0"
seurat.list$HN15$HPV <- "Neg"
seurat.list$HN15$Inflam.status <- "Low"
# HN16
seurat.list$HN16$Gender <- "M"
seurat.list$HN16$AgeGroup <- "40-49"
seurat.list$HN16$Smoking <- "Yes"
seurat.list$HN16$Alcohol <- "Yes"
seurat.list$HN16$DiseaseSite <- "Oropharynx"
seurat.list$HN16$T_Stage <- "T2"
seurat.list$HN16$N_Stage <- "N1"
seurat.list$HN16$M_Stage <- "M0"
seurat.list$HN16$HPV <- "Pos"
seurat.list$HN16$Inflam.status <- "Med"
# HN17
seurat.list$HN17$Gender <- "M"
seurat.list$HN17$AgeGroup <- "50-59"
seurat.list$HN17$Smoking <- "Yes"
seurat.list$HN17$Alcohol <- "Yes"
seurat.list$HN17$DiseaseSite <- "Oropharynx"
seurat.list$HN17$T_Stage <- "T1"
seurat.list$HN17$N_Stage <- "N1"
seurat.list$HN17$M_Stage <- "M0"
seurat.list$HN17$HPV <- "Pos"
seurat.list$HN17$Inflam.status <- "High"
# HN18
seurat.list$HN18$Gender <- "M"
seurat.list$HN18$AgeGroup <- "50-59"
seurat.list$HN18$Smoking <- "Yes"
seurat.list$HN18$Alcohol <- "Yes"
seurat.list$HN18$DiseaseSite <- "Oropharynx"
seurat.list$HN18$T_Stage <- "T2"
seurat.list$HN18$N_Stage <- "N2"
seurat.list$HN18$M_Stage <- "M0"
seurat.list$HN18$HPV <- "Pos"
seurat.list$HN18$Inflam.status <- "Low"

anchors <- FindIntegrationAnchors(object.list = seurat.list)
saveRDS(anchors,file = "HNSCC/rds/CD45pPBMC.anchors.rds")
integrated <- IntegrateData(anchorset = anchors)
integrated <- ScaleData(integrated, verbose = FALSE)
saveRDS(integrated, file = "HNSCC/rds/CD45pPBMC.integrated.rds")

rm(list = ls());gc()

#-------------------- CD45nTumor -----------------
if(F){
  library(Seurat)
  integrated <- readRDS("HNSCC/rds/CD45nTumor.integrated.rds")
  DefaultAssay(integrated) <- "RNA"
  integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(integrated, features = VariableFeatures(integrated))
  # Run PCA and Determine Dimensions for 90% Variance
  PCDeterminators <- function(object){
    stdev <- object@reductions$pca@stdev
    var <- stdev^2
    endVar <- 0
    for(i in 1:length(var)){
      total <- sum(var)
      numerator <- sum(var[1:i])
      exp.var <- numerator/total
      if(endVar == 0){
        if(exp.var > 0.9){
          endVar <- endVar + 1
          PC.num <- i
        }
      }
    }
    sum(var[1:PC.num])/sum(var)
    return(PC.num)
  }
  PC.num <- PCDeterminators(integrated)
  integrated <- FindNeighbors(integrated, dims = 1:PC.num, reduction = "pca")
  integrated <- FindClusters(integrated, resolution = seq(0.2, 1.0, 0.2), reduction = "pca")
  integrated <- RunUMAP(integrated, dims = 1:PC.num, reduction = "pca")
  integrated <- RunTSNE(integrated, dims = 1:PC.num)
  saveRDS(integrated,file = "HNSCC/rds/CD45nTumor.integrated.dr.rds")
  
  DefaultAssay(integrated) <- "RNA"
  # DotPlot(integrated, features = c("CD3D","CD3E","CD3G","CD4","CD8A","NCAM1","TNF","IFNG","NKG7"),group.by = "seurat_clusters")
  
  # CD3-,NKG7+, GNLY+, KLRD1+
  # DotPlot(integrated, features = c("CD3D","NKG7","GNLY","KLRD1","NCAM1","IFNG"))
  
  rm(list = ls());gc()
}




#-------------------- CD45pTumor -----------------
library(Seurat)
integrated <- readRDS("HNSCC/rds/CD45pTumor.integrated.rds")
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)
integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, features = VariableFeatures(integrated))
# Run PCA and Determine Dimensions for 90% Variance
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(integrated)
integrated <- FindNeighbors(integrated, dims = 1:PC.num, reduction = "pca")
integrated <- FindClusters(integrated, resolution = seq(0.2, 1.0, 0.2), reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:PC.num, reduction = "pca")
integrated <- RunTSNE(integrated, dims = 1:PC.num)
saveRDS(integrated,file = "HNSCC/rds/CD45pTumor.integrated.dr.rds")

DefaultAssay(integrated) <- "RNA"
# DotPlot(integrated, features = c("CD3D","CD3E","CD3G","CD4","CD8A","NCAM1","TNF","IFNG","NKG7"),group.by = "seurat_clusters")

# CD3-,NKG7+, GNLY+, KLRD1+, NCAM1-
DotPlot(integrated, features = c("CD3D","NKG7","GNLY","KLRD1","NCAM1"))

nk.cells <- rownames(subset(integrated@meta.data, seurat_clusters == 14 | seurat_clusters == 15))
CD45pTumor.NK <- subset(integrated, cells = nk.cells)
saveRDS(CD45pTumor.NK, file = "HNSCC/rds/CD45pTumor.NK.rds")

rm(list = ls());gc()






#-------------------- CD45pPBMC -----------------
library(Seurat)
integrated <- readRDS("HNSCC/rds/CD45pPBMC.integrated.rds")
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)
integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, features = VariableFeatures(integrated))
# Run PCA and Determine Dimensions for 90% Variance
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(integrated)
integrated <- FindNeighbors(integrated, dims = 1:PC.num, reduction = "pca")
integrated <- FindClusters(integrated, resolution = seq(0.2, 1.0, 0.2), reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:PC.num, reduction = "pca")
integrated <- RunTSNE(integrated, dims = 1:PC.num)
saveRDS(integrated,file = "HNSCC/rds/CD45pPBMC.integrated.dr.rds")

# visualize
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey25"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "none",
                 plot.title = element_blank(),
                 legend.text = element_text(size = 14))
p.ct <- DimPlot(integrated, reduction = "tsne", label = T,label.size = 5) + mytheme
ggsave("HNSCC/figures/CD45pPBMC.cluster.pdf",p.ct, width = 5,height = 5)

cell.highlight <- rownames(subset(integrated@meta.data, seurat_clusters == 2 | seurat_clusters == 11 | seurat_clusters == 17))
p.NK <- DimPlot(integrated,reduction = "tsne", label = T, cells.highlight = cell.highlight,label.size = 4, 
                label.color = rainbow(length(unique(as.character(integrated$seurat_clusters))))) + mytheme
ggsave("CD45pPBMC.nk.highlight.pdf",p.NK,width = 5,height = 5)

DefaultAssay(integrated) <- "RNA"
# DotPlot(integrated, features = c("CD3D","CD3E","CD3G","CD4","CD8A","NCAM1","TNFA","IFNG"))

# CD3-,NKG7+, GNLY+, KLRD1+, NCAM1-
fea <- c("CD3D","NKG7","GNLY","KLRD1","NCAM1")
library(MySeuratWrappers)
p.nk <- VlnPlot(integrated,features = c("CD3D","NKG7","GNLY","KLRD1","NCAM1"),stacked = T, pt.size = 0,
        direction = "horizontal",combine = T, x.lab = '', y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave("HNSCC/figures/violin.nk.markers.pdf",p.nk,width = 6,height = 10)

nk.cells <- rownames(subset(integrated@meta.data, seurat_clusters == 2 | seurat_clusters == 11 | seurat_clusters == 17))
CD45pPBMC.NK <- subset(integrated, cells = nk.cells)
saveRDS(CD45pPBMC.NK, file = "HNSCC/rds/CD45pPBMC.NK.rds")

rm(list = ls());gc()

# ------------------ NK subsets --------------------
CD45pPBMC.NK <- readRDS("HNSCC/rds/CD45pPBMC.NK.rds")
CD45pPBMC.NK <- NormalizeData(CD45pPBMC.NK, normalization.method = "LogNormalize", scale.factor = 10000)
CD45pPBMC.NK <- FindVariableFeatures(CD45pPBMC.NK, selection.method = "vst", nfeatures = 2000)
DefaultAssay(CD45pPBMC.NK) <- "integrated"
CD45pPBMC.NK <- ScaleData(CD45pPBMC.NK, features = rownames(CD45pPBMC.NK))
CD45pPBMC.NK <- RunPCA(CD45pPBMC.NK, features = VariableFeatures(CD45pPBMC.NK))
# Run PCA and Determine Dimensions for 90% Variance
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(CD45pPBMC.NK)
CD45pPBMC.NK <- FindNeighbors(CD45pPBMC.NK, dims = 1:PC.num, reduction = "pca")
CD45pPBMC.NK <- FindClusters(CD45pPBMC.NK, resolution = seq(0.2, 0.6, 0.1), reduction = "pca")
CD45pPBMC.NK <- RunUMAP(CD45pPBMC.NK, dims = 1:PC.num, reduction = "pca")
CD45pPBMC.NK <- RunTSNE(CD45pPBMC.NK, dims = 1:PC.num)
CD45pPBMC.NK$seurat_clusters <- CD45pPBMC.NK$integrated_snn_res.0.2
Idents(CD45pPBMC.NK) <- CD45pPBMC.NK$seurat_clusters
p.ct <- DimPlot(CD45pPBMC.NK, reduction = "tsne",label = T,label.size = 5, cols = c("#4E9BD3","#DB4840","#DF78C1")) + mytheme
ggsave("HNSCC/figures/CD45pPBMC.NK.cluster.pdf",p.ct, width = 4,height = 4)

p.ct2 <- DimPlot(CD45pPBMC.NK, reduction = "tsne",label = T,label.size = 5,group.by = "classI", 
                 cols = c("#4E9BD3","#DB4840","#DF78C1")) + mytheme
ggsave("HNSCC/figures/CD45pPBMC.NK.Blood.pdf",p.ct2, width = 4,height = 4)

DefaultAssay(CD45pPBMC.NK) <- "RNA"
CD45pPBMC.NK <- ScaleData(CD45pPBMC.NK, features = rownames(CD45pPBMC.NK))

# CD56dim, CD16, CD57
# CD56brightCD16-NK
# CD56dimCD16+CD57-NK
# CD56dimCD16+CD57+NK
library(MySeuratWrappers)
DefaultAssay(CD45pPBMC.NK) <- "RNA"
p.nk.sub <- VlnPlot(CD45pPBMC.NK, features = c("NCAM1","FCGR3A","B3GAT1"),
                    cols = c("#4E9BD3","#DB4840","#DF78C1"),
                    pt.size = 0,direction = "horizontal",stacked = T,combine = T,x.lab = '', y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave("HNSCC/figures/violin.nk.sub.pdf",p.nk.sub,width = 2.5,height = 4)

# GZMB
p.GZMB <- VlnPlot(CD45pPBMC.NK, features = "GZMB",pt.size = 0,direction = "horizontal",
                  cols = c("#4E9BD3","#DB4840","#DF78C1"),
                  stacked = T,combine = T,x.lab = '', y.lab = '') +
  NoLegend() + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_rect(size = 1.5,colour = "grey45"))
ggsave("HNSCC/figures/violin.nk.gzmb.pdf",p.GZMB,width = 2,height = 4)

ILs <- c("IL12A","IL15","IL18")
p.ILs <- VlnPlot(CD45pPBMC.NK, features = ILs,pt.size = 0,direction = "horizontal",
                 cols = c("#4E9BD3","#DB4840","#DF78C1"),
                 stacked = T,combine = T,x.lab = '', y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_rect(size = 1.5,colour = "grey45"))
ggsave("HNSCC/figures/violin.nk.ILs.pdf",p.ILs,width = 4,height = 6)

# KIRs: Killer Cell Ig-Like Receptors
# KIR3DX1,KIR2DL3,KIR2DL1,KIR2DL4,KIR3DL1,KIR3DL2,KIR3DL3
# KIRs <-c("KIR3DX1",'KIR2DL3','KIR2DL1','KIR2DL4','KIR3DL1','KIR3DL2','KIR3DL3')
# chemokines <- c("CCL4","CCL3")

Idents(CD45pPBMC.NK) <- CD45pPBMC.NK$integrated_snn_res.0.2
all.markers <- FindAllMarkers(CD45pPBMC.NK,assay = "RNA",only.pos = T)
save(all.markers, file = "HNSCC/rds/CD45pPBMC.nk.all.markers.rda")
curr.ids <- c(0:2)
new.ids <- c("Bl1","Bl2","Bl3")
CD45pPBMC.NK$classI <- plyr::mapvalues(CD45pPBMC.NK$seurat_clusters, from = curr.ids, to = new.ids)
saveRDS(CD45pPBMC.NK,file = "HNSCC/rds/CD45pPBMC.integrated.nk.rds")

## Cell proportion
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey25"),
                 panel.grid = element_blank(),
                 axis.text = element_text(size = 12,face = "bold",colour = "grey20"),
                 axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
                 legend.position = "top", legend.direction = "horizontal",
                 legend.text = element_text(size = 14),legend.title = element_text(size = 15))

p.NK.prop <- ggplot(CD45pPBMC.NK@meta.data, aes(x = Sample, fill = classI)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = c("#4E9BD3","#DB4840","#DF78C1")) +
  mytheme +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction", x = '')
ggsave("HNSCC/figures/CD45pPBMC.NK.cell.prop.pdf",p.NK.prop,width = 4,height = 8)

# cluster 0
genes.c0 <- all.markers[all.markers$cluster == 0,]$gene
write(genes.c0, file = "HNSCC/marker.c0.txt")
# cluster 1
genes.c1 <- all.markers[all.markers$cluster == 1,]$gene
write(genes.c1, file = "HNSCC/marker.c1.txt")
# cluster 2
genes.c2 <- all.markers[all.markers$cluster == 2,]$gene
write(genes.c2, file = "HNSCC/marker.c2.txt")

mat <- GetAssayData(CD45pPBMC.NK, assay = "RNA", slot = "scale.data")
mat[mat > 2] = 2
mat[mat < -2] = -2
cls.info <- sort(CD45pPBMC.NK$classI)
mat <- as.matrix(mat[all.markers$gene, names(cls.info)])

gene <- c("SPON2","FGFBP2","PRF1","PNF1","ACTB","FCGR3A","GZMB","S100A4",
          "CD52","CCL5","GZMH","KLRC1",
          "GZMK","XCL1","XCL2","SELL","COTL1")
gene <- unique(gene)
gene.pos <- which(rownames(mat) %in% gene)
row.anno <- ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = gene.pos, labels = rownames(mat)[which(rownames(mat) %in% gene)]))
cols <- c("#4E9BD3","#DB4840","#DF78C1")
names(cols) <- levels(cls.info)
top.anno <- ComplexHeatmap::HeatmapAnnotation(
  cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill=cols), labels = levels(cls.info), labels_gp = grid::gpar(cex=0.5,col='white'))
)
col.fun <- circlize::colorRamp2(seq(min(mat),max(mat),length=3), c('#377EB8','white','#E41A1C'))
pdf("HNSCC/diff.cls.heatmap.pdf",width = 10,height = 10)
ComplexHeatmap::Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F,
                        show_column_names = F, show_row_names = F,
                        column_split = cls.info, right_annotation = row.anno,
                        column_title = NULL, top_annotation = top.anno,
                        heatmap_legend_param = list(title = 'Expression', title_position = 'leftcenter-rot'),col = col.fun)
dev.off()

# Biocarta
pathway <- xlsx::read.xlsx("HNSCC/CD45pPBMC.NK.diff.Biocarta.xlsx",sheetIndex = 1)
pathway <- reshape2::melt(pathway)
colnames(pathway) <- c("Biocarta","Sample","-log10(p-value)")
p.bioc <- ggplot(data = pathway, aes(x=`-log10(p-value)`, y=Biocarta, fill=Sample)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#4E9BD3","#DB4840","#DF78C1"))

ggsave("HNSCC/figures/CD45pPBMC.NK.diff.pathway.Biocarta.pdf",p.bioc,width = 5,height = 4)


#---------- Monocle -------------
library(Seurat)
library(monocle)
CD45pPBMC.NK <- readRDS("HNSCC/rds/CD45pPBMC.integrated.nk.rds")
dat <- as(as.matrix(CD45pPBMC.NK@assays$RNA@counts), "sparseMatrix")
pd <- new('AnnotatedDataFrame', data = CD45pPBMC.NK@meta.data)
fData <- data.frame(gene_short_name = rownames(dat), row.names = rownames(dat))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(dat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds,relative_expr = T)
# 使用monocle选择的高变基因
disp.tbl <- dispersionTable(mycds)
disp.genes <- subset(disp.tbl, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p1 <- plot_ordering_genes(mycds)
# 使用disp.genes展开后续分析
mycds <- reduceDimension(mycds, max_components = 2, method = "DDRTree")
mycds <- orderCells(mycds)
saveRDS(mycds, file = "HNSCC/rds/CD45pPBMC.monocle.cds.rds")
source("HNSCC/pseudotime_heatmap.R")
p1 <- plot_cell_trajectory(mycds, color_by = "classI") + 
  scale_color_manual(values = c("#4E9BD3","#DB4840","#DF78C1")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("HNSCC/figures/CD45pPBMC.NK.trajectory.pdf",p1,width = 6,height = 4)
  
p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime") +
  scale_color_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256))) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("HNSCC/figures/CD45pPBMC.NK.pseudo.pdf",p2,width = 6,height = 4)

p3 <- plot_cell_trajectory(mycds, color_by = "HPV") + facet_wrap(~HPV, nrow = 2) +
  scale_color_manual(values = c("#4E9BD3","#DB4840")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("HNSCC/figures/CD45pPBMC.NK.splitByHPV.pdf",p3,width = 6,height = 8)

# Mature NK signatures
# CD160,FCGR3A, GZMB, CCL3
mature <- c("CD160","FCGR3A","GZMB","CCL3")

cols <- c("grey88","DarkCyan")

for(i in 1:length(mature)){
  p <- FeaturePlot(CD45pPBMC.NK, features = mature[i], reduction = "tsne", cols = cols) +
    theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
          axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank())
  ggsave(paste0("HNSCC/figures/CD45pPBMC.",mature[i],".pdf"),p,width = 5,height = 5)
}

# gene <- c("SPON2","FGFBP2","PRF1","PNF1","ACTB","FCGR3A","GZMB","S100A4",
#           "CD52","CCL5","GZMH","KLRC1",
#           "GZMK","XCL1","XCL2","SELL","COTL1")
gene <- unique(all.markers$gene)

# cluster-based differential genes
disp.tbl <- dispersionTable(mycds)
disp.genes <- subset(disp.tbl, mean_expression >= 0.5 & dispersion_empirical >= 1 & dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff.test <- differentialGeneTest(mycds[disp.genes,], cores = 1)
# diff.test <- differentialGeneTest(mycds[gene,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig.genename <- rownames(subset(diff.test, qval < 0.01))

p1 <- pseudotime_heatmap(mycds[sig.genename,],show_rownames = T,return_heatmap = T,num_clusters = 1)

save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p1$ph_res, filename = "HNSCC/figures/CD45pPBMC.NK.pseudo.heatmap2.pdf",width = 5,height = 6)

rm(list = ls());gc()


#--------- Compare between HPV+ and HPV-
CD45pPBMC.NK <- readRDS("HNSCC/rds/CD45pPBMC.integrated.nk.rds")
Idents(CD45pPBMC.NK) <- CD45pPBMC.NK$classI
Bl1 <- FindMarkers(CD45pPBMC.NK,assay = "RNA",only.pos = F,group.by = "HPV", ident.1 = "Pos",subset.ident = "Bl1")
Bl1$NKsub <- "Bl1"
Bl1$compare <- "HPV+_vs_HPV-"
Bl2 <- FindMarkers(CD45pPBMC.NK,assay = "RNA",only.pos = F,group.by = "HPV", ident.1 = "Pos",subset.ident = "Bl2")
Bl2$NKsub <- "Bl2"
Bl2$compare <- "HPV+_vs_HPV-"
Bl3 <- FindMarkers(CD45pPBMC.NK,assay = "RNA",only.pos = F,group.by = "HPV", ident.1 = "Pos",subset.ident = "Bl3")
Bl3$NKsub <- "Bl3"
Bl3$compare <- "HPV+_vs_HPV-"
all.markers <- rbind(Bl1,Bl2,Bl3)
rm(Bl1,Bl2,Bl3);gc()

write(rownames(Bl1[Bl1$avg_logFC > 0,]), file = "HNSCC/out/CD45pPBMC.NK.Bl1.HPV+.genes.txt")
write(rownames(Bl1[Bl1$avg_logFC < 0,]), file = "HNSCC/out/CD45pPBMC.NK.Bl1.HPV-.genes.txt")
write(rownames(Bl2[Bl2$avg_logFC > 0,]), file = "HNSCC/out/CD45pPBMC.NK.Bl2.HPV+.genes.txt")
write(rownames(Bl2[Bl2$avg_logFC < 0,]), file = "HNSCC/out/CD45pPBMC.NK.Bl2.HPV-.genes.txt")
write(rownames(Bl3[Bl3$avg_logFC > 0,]), file = "HNSCC/out/CD45pPBMC.NK.Bl3.HPV+.genes.txt")
write(rownames(Bl3[Bl3$avg_logFC < 0,]), file = "HNSCC/out/CD45pPBMC.NK.Bl3.HPV-.genes.txt")

library(ggplot2)
library(ggrepel)
mytheme<- theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 12))
# Bl1
Bl1 <- all.markers[all.markers$NKsub=="Bl1",]
Bl1$Significance <- ifelse(Bl1$p_val < 0.01, TRUE, FALSE)
Bl1 <- Bl1[order(Bl1[,2],decreasing = T),]
top10 <- rbind(head(Bl1,10),tail(Bl1,10))
p.Bl1 <- ggplot(Bl1, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#C72719","#BEBBBD")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave("HNSCC/figures/CD45pPBMC.NK.Bl1.DEgenes.pdf",p.Bl1,width = 6,height = 5)
# Bl2
Bl2 <- all.markers[all.markers$NKsub=="Bl2",]
Bl2$Significance <- ifelse(Bl2$p_val < 0.01, TRUE, FALSE)
Bl2 <- Bl2[order(Bl2[,2],decreasing = T),]
top10 <- rbind(head(Bl2,10),tail(Bl2,10))
p.Bl2 <- ggplot(Bl2, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave("HNSCC/figures/CD45pPBMC.NK.Bl2.DEgenes.pdf",p.Bl2,width = 6,height = 5)
# Bl3
Bl3 <- all.markers[all.markers$NKsub=="Bl3",]
Bl3$Significance <- ifelse(Bl3$p_val < 0.01, TRUE, FALSE)
Bl3 <- Bl3[order(Bl3[,2],decreasing = T),]
top10 <- rbind(head(Bl3,10),tail(Bl3,10))
p.Bl3 <- ggplot(Bl3, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave("HNSCC/figures/CD45pPBMC.NK.Bl3.DEgenes.pdf",p.Bl3,width = 6,height = 5)

save(all.markers, file = "HNSCC/out/CD45pPBMC.NK.sub.HPV+_vs_PHV-.rda")


#----- Extract Diff genes for visualization ----
# Bl1
Bl1 <- subset(CD45pPBMC.NK, cells = rownames(subset(CD45pPBMC.NK@meta.data, classI == "Bl1")))
Idents(Bl1) <- Bl1$HPV
Bl1.diff <- all.markers[all.markers$NKsub == "Bl1",]
Bl1.diff <- Bl1.diff[order(Bl1.diff[,2],decreasing = T),]
Bl1.gene <- unique(rownames(Bl1.diff))
Bl1.m <- match(Bl1.gene, rownames(Bl1))
Bl1.gene <- rownames(Bl1)[Bl1.m]
Bl1.gene <- Bl1.gene[!is.na(Bl2.gene)]
diff.exp <- AverageExpression(Bl1, assays = "RNA", features = Bl1.gene,slot = "scale.data")
diff.exp <- diff.exp$RNA 
diff.exp$`HPV+` <- diff.exp$Pos; diff.exp$`HPV-` <- diff.exp$Neg
diff.exp$Pos <- NULL; diff.exp$Neg <- NULL
diff.exp <- t(diff.exp)
pheatmap(diff.exp, cluster_rows = F,cluster_cols = F,border_color = NA,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),
         gaps_row = 1, filename = "Bl1.tmp.pdf",width = 14,height = 3)
pheatmap(diff.exp, cluster_rows = F,cluster_cols = F,border_color = NA,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),
         filename = "Bl1.pdf",width = 8,height = 1.5,legend = F)
# Bl2
Bl2 <- subset(CD45pPBMC.NK, cells = rownames(subset(CD45pPBMC.NK@meta.data, classI == "Bl2")))
Idents(Bl2) <- Bl2$HPV
Bl2.diff <- all.markers[all.markers$NKsub == "Bl2",]
Bl2.diff <- Bl2.diff[order(Bl2.diff[,2],decreasing = T),]
Bl2.gene <- unique(rownames(Bl2.diff))
Bl2.m <- match(Bl2.gene, rownames(Bl2))
Bl2.gene <- rownames(Bl2)[Bl2.m]
Bl2.gene <- Bl2.gene[!is.na(Bl2.gene)]
diff.exp <- AverageExpression(Bl2, assays = "RNA", features = Bl2.gene,slot = "scale.data")
diff.exp <- diff.exp$RNA 
diff.exp$`HPV+` <- diff.exp$Pos; diff.exp$`HPV-` <- diff.exp$Neg
diff.exp$Pos <- NULL; diff.exp$Neg <- NULL
diff.exp <- t(diff.exp)
tmp <- diff.exp
tmp[tmp > 0.4] <- 0.4
tmp[tmp < -0.2] <- -0.2
pheatmap(tmp, cluster_rows = F,cluster_cols = F,border_color = NA,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),
         gaps_row = 1, filename = "Bl2.tmp.pdf",width = 14,height = 3)
pheatmap(tmp, cluster_rows = F,cluster_cols = F,border_color = NA,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),
         filename = "Bl2.pdf",width = 8,height = 1.5,legend = F)
# Bl3
Bl3 <- subset(CD45pPBMC.NK, cells = rownames(subset(CD45pPBMC.NK@meta.data, classI == "Bl3")))
Idents(Bl3) <- Bl3$HPV
Bl3.diff <- all.markers[all.markers$NKsub == "Bl3",]
Bl3.diff <- Bl3.diff[order(Bl3.diff[,2],decreasing = T),]
Bl3.gene <- unique(rownames(Bl3.diff))
Bl3.m <- match(Bl3.gene, rownames(Bl3))
Bl3.gene <- rownames(Bl3)[Bl3.m]
Bl3.gene <- Bl3.gene[!is.na(Bl3.gene)]
diff.exp <- AverageExpression(Bl3, assays = "RNA", features = Bl3.gene,slot = "scale.data")
diff.exp <- diff.exp$RNA 
diff.exp$`HPV+` <- diff.exp$Pos; diff.exp$`HPV-` <- diff.exp$Neg
diff.exp$Pos <- NULL; diff.exp$Neg <- NULL
diff.exp <- t(diff.exp)
tmp <- diff.exp
tmp[tmp > 0.4] <- 0.4
tmp[tmp < -0.2] <- -0.2
pheatmap(tmp, cluster_rows = F,cluster_cols = F,border_color = NA,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),
         gaps_row = 1, filename = "Bl3.tmp.pdf",width = 14,height = 3)
pheatmap(tmp, cluster_rows = F,cluster_cols = F,border_color = NA,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),
         filename = "Bl3.pdf",width = 8,height = 1.5,legend = F)

rm(list = ls());gc()

#------- survival analysis ------------

mytheme <- theme(legend.title = element_text(size = 14), 
                 legend.text = element_text(size = 14),
                 axis.text.x = element_text(size = 14),
                 axis.text.y = element_text(size = 14),
                 axis.line.x = element_line(size = 1),
                 axis.line.y = element_line(size = 1),
                 axis.title.x = element_text(size = 16),
                 axis.title.y = element_text(size = 16))

# make TCGA expression data frame
tcga.exp.df <- read.table("HNSCC/data/cBioPortal/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",header = T)
tcga.exp.df <- tcga.exp.df[!duplicated(tcga.exp.df$Hugo_Symbol),]
rownames(tcga.exp.df) <- tcga.exp.df$Hugo_Symbol
tcga.exp.df <- tcga.exp.df[,-c(1,2)]
colnames(tcga.exp.df) <- gsub("\\.","-",colnames(tcga.exp.df))
tcga.exp.df <- data.frame(t(tcga.exp.df))
colnames(tcga.exp.df) <- gsub("\\.","-",colnames(tcga.exp.df))
# make clinical info
clin <- read.table("HNSCC/data/cBioPortal/data_clinical_patient.txt",skip = 4,header = T,sep = "\t")
clin <- data.frame(clin$PATIENT_ID, clin$OS_MONTHS, clin$OS_STATUS)
clin <- clin[clin$clin.OS_STATUS != "[Not Avaibale]",]
Dead <- rep(NA, dim(clin)[1])
for(i in 1:length(clin$clin.OS_STATUS)){
  if(stringr::str_split(clin$clin.OS_STATUS, ":", simplify = T)[,2][i] == "DECEASED"){
    Dead[i] <- 1
  }else{
    Dead[i] <- 0
  }
}
clin$Dead <- Dead
clin$clin.OS_MONTHS <- as.integer(clin$clin.OS_MONTHS)

load("HNSCC/out/CD45pPBMC.NK.sub.HPV+_vs_PHV-.rda")
# Bl1
Bl1.diff <- all.markers[all.markers$NKsub=="Bl1",]
Bl1.diff <- Bl1.diff[order(Bl1.diff[,2],decreasing = T),]
Bl1.gene <- rownames(Bl1.diff)
Bl1.gene <- colnames(tcga.exp.df)[match(Bl1.gene,colnames(tcga.exp.df))]
Bl1.gene <- Bl1.gene[!is.na(Bl1.gene)]
tcga.exp.df.2 <- data.frame(tcga.exp.df[,Bl1.gene])
removeColsAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
tcga.exp.df.2 <- removeColsAllNA(tcga.exp.df.2)
tcga.exp.df.2 <- data.frame(scale(tcga.exp.df.2))
tcga.exp.df.2$PATIENT_ID <- gsub("-01","",rownames(tcga.exp.df.2))
avg.exp <- aggregate(.~PATIENT_ID, FUN = mean, data = tcga.exp.df.2)
merged.data <- merge(clin,avg.exp, by.x = "clin.PATIENT_ID", by.y = "PATIENT_ID")
# survival analysis for signatures
# HPV+
merged.data.1 <- merged.data[,c(1:19)]
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data.1[,c(2,4:ncol(merged.data.1))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data.1[,5:ncol(merged.data.1)] - colMeans(merged.data.1[,5:ncol(merged.data.1)]))) %*% coefs
group <- rep(NA, dim(merged.data.1)[1])
for(i in 1:dim(merged.data.1)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data.1$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data.1)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "HPV+")) +
  mytheme +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave("Bl1.HPV+.surv.pdf",p.surv,width = 6,height = 5)

# # HPV-
merged.data.2 <- merged.data[,c(1:4,20:ncol(merged.data))]
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data.2[,c(2,4:ncol(merged.data.2))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data.2[,5:ncol(merged.data.2)] - colMeans(merged.data.2[,5:ncol(merged.data.2)]))) %*% coefs
group <- rep(NA, dim(merged.data.2)[1])
for(i in 1:dim(merged.data.2)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data.2$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data.2)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "HPV-")) +
  mytheme +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave("Bl1.HPV-.surv.pdf",p.surv,width = 6,height = 5)

Bl1 <- list(Bl1.HPVpos = merged.data.1, Bl1.HPVneg = merged.data.2)
save(Bl1, file = "HNSCC/out/CD45pPBMC.NK.Bl1.survdat.rda")
rm(Bl1)

# Bl2
Bl2.diff <- all.markers[all.markers$NKsub=="Bl2",]
Bl2.diff <- Bl2.diff[order(Bl2.diff[,2],decreasing = T),]
Bl2.gene <- rownames(Bl2.diff)
Bl2.gene <- colnames(tcga.exp.df)[match(Bl2.gene,colnames(tcga.exp.df))]
Bl2.gene <- Bl2.gene[!is.na(Bl2.gene)]
tcga.exp.df.2 <- data.frame(tcga.exp.df[,Bl2.gene])
removeColsAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
tcga.exp.df.2 <- removeColsAllNA(tcga.exp.df.2)
tcga.exp.df.2 <- data.frame(scale(tcga.exp.df.2))
tcga.exp.df.2$PATIENT_ID <- gsub("-01","",rownames(tcga.exp.df.2))
avg.exp <- aggregate(.~PATIENT_ID, FUN = mean, data = tcga.exp.df.2)
merged.data <- merge(clin,avg.exp, by.x = "clin.PATIENT_ID", by.y = "PATIENT_ID")
# survival analysis for signatures
# HPV+
merged.data.1 <- merged.data[,c(1:47)]
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data.1[,c(2,4:ncol(merged.data.1))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data.1[,5:ncol(merged.data.1)] - colMeans(merged.data.1[,5:ncol(merged.data.1)]))) %*% coefs
group <- rep(NA, dim(merged.data.1)[1])
for(i in 1:dim(merged.data.1)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data.1$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data.1)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "HPV+")) +
  mytheme +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave("Bl2.HPV+.surv.pdf",p.surv,width = 6,height = 5)

# # HPV-
merged.data.2 <- merged.data[,c(1:4,48:ncol(merged.data))]
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data.2[,c(2,4:ncol(merged.data.2))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data.2[,5:ncol(merged.data.2)] - colMeans(merged.data.2[,5:ncol(merged.data.2)]))) %*% coefs
group <- rep(NA, dim(merged.data.2)[1])
for(i in 1:dim(merged.data.2)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data.2$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data.2)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "HPV-")) +
  mytheme +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave("Bl2.HPV-.surv.pdf",p.surv,width = 6,height = 5)

Bl2 <- list(Bl2.HPVpos = merged.data.1, Bl2.HPVneg = merged.data.2)
save(Bl2, file = "HNSCC/out/CD45pPBMC.NK.Bl2.survdat.rda")
rm(Bl2)

# Bl3
Bl3.diff <- all.markers[all.markers$NKsub=="Bl3",]
Bl3.diff <- Bl3.diff[order(Bl3.diff[,2],decreasing = T),]
Bl3.gene <- rownames(Bl3.diff)
Bl3.gene <- colnames(tcga.exp.df)[match(Bl3.gene,colnames(tcga.exp.df))]
Bl3.gene <- Bl3.gene[!is.na(Bl3.gene)]
tcga.exp.df.2 <- data.frame(tcga.exp.df[,Bl3.gene])
removeColsAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
tcga.exp.df.2 <- removeColsAllNA(tcga.exp.df.2)
tcga.exp.df.2 <- data.frame(scale(tcga.exp.df.2))
tcga.exp.df.2$PATIENT_ID <- gsub("-01","",rownames(tcga.exp.df.2))
avg.exp <- aggregate(.~PATIENT_ID, FUN = mean, data = tcga.exp.df.2)
merged.data <- merge(clin,avg.exp, by.x = "clin.PATIENT_ID", by.y = "PATIENT_ID")
# survival analysis for signatures
# HPV+
merged.data.1 <- merged.data[,c(1:15)]
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data.1[,c(2,4:ncol(merged.data.1))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data.1[,5:ncol(merged.data.1)] - colMeans(merged.data.1[,5:ncol(merged.data.1)]))) %*% coefs
group <- rep(NA, dim(merged.data.1)[1])
for(i in 1:dim(merged.data.1)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data.1$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data.1)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "HPV+")) +
  mytheme +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave("HNSCC/figures/Bl3.HPV+.surv.pdf",p.surv,width = 6,height = 5)

# # HPV-
merged.data.2 <- merged.data[,c(1:4,16:ncol(merged.data))]
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data.2[,c(2,4:ncol(merged.data.2))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data.2[,5:ncol(merged.data.2)] - colMeans(merged.data.2[,5:ncol(merged.data.2)]))) %*% coefs
group <- rep(NA, dim(merged.data.2)[1])
for(i in 1:dim(merged.data.2)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data.2$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data.2)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "HPV-")) +
  mytheme +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave("HNSCC/figures/Bl3.HPV-.surv.pdf",p.surv,width = 6,height = 5)

Bl3 <- list(Bl3.HPVpos = merged.data.1, Bl3.HPVneg = merged.data.2)
save(Bl3, file = "HNSCC/out/CD45pPBMC.NK.Bl3.survdat.rda")
rm(Bl3)

## survival analysis for each gene
load("HNSCC/out/CD45pPBMC.NK.Bl1.survdat.rda")
# Bl1.HPV+
Bl1.pos <- Bl1$Bl1.HPVpos
mysurv <- survival::Surv(Bl1.pos$clin.OS_MONTHS,Bl1.pos$Dead)
Group <- ifelse(Bl1.pos[,5] > median(na.omit(Bl1.pos[,5])), 'high', 'low')
sfit <- survival::survfit(mysurv~Group,data = Bl1.pos)
survp <- survminer::ggsurvplot(sfit)
if(!dir.exists("HNSCC/figures/surv/Bl1/pos")) dir.create("HNSCC/figures/surv/Bl1/pos",recursive = T)
for(i in 5:(ncol(Bl1.pos)-1)){
  Group <- ifelse(Bl1.pos[,i] > median(na.omit(Bl1.pos[,i])), 'high', 'low')
  sfit <- survival::survfit(mysurv~Group,data = Bl1.pos)
  survp <- survminer::ggsurvplot(sfit, conf.int = F, pval = F)
  p1 <- survp$plot + annotate("text", label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                              x = 40, y = 0.15, size = 5) +
    annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
    guides(color = guide_legend(title = "HPV+")) +
    mytheme +
    scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
  ggsave(paste0("HNSCC/figures/surv/Bl1/pos/",colnames(Bl1.pos)[i],".pdf"),p1,width = 6,height = 5)
}
# Bl1.HPV-
Bl1.neg <- Bl1$Bl1.HPVneg
mysurv <- survival::Surv(Bl1.neg$clin.OS_MONTHS,Bl1.neg$Dead)
Group <- ifelse(Bl1.neg[,5] > median(na.omit(Bl1.neg[,5])), 'high', 'low')
sfit <- survival::survfit(mysurv~Group,data = Bl1.neg)
survp <- survminer::ggsurvplot(sfit)
if(!dir.exists("HNSCC/figures/surv/Bl1/neg")) dir.create("HNSCC/figures/surv/Bl1/neg",recursive = T)
for(i in 5:(ncol(Bl1.neg)-1)){
  Group <- ifelse(Bl1.neg[,i] > median(na.omit(Bl1.neg[,i])), 'high', 'low')
  sfit <- survival::survfit(mysurv~Group,data = Bl1.neg)
  survp <- survminer::ggsurvplot(sfit, conf.int = F, pval = F)
  p1 <- survp$plot + annotate("text", label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                              x = 40, y = 0.15, size = 5) +
    annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
    guides(color = guide_legend(title = "HPV-")) +
    mytheme +
    scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
  ggsave(paste0("HNSCC/figures/surv/Bl1/neg/",colnames(Bl1.neg)[i],".pdf"),p1,width = 6,height = 5)
}

load("HNSCC/out/CD45pPBMC.NK.Bl2.survdat.rda")
# Bl2.HPV+
Bl2.pos <- Bl2$Bl2.HPVpos
mysurv <- survival::Surv(Bl2.pos$clin.OS_MONTHS,Bl2.pos$Dead)
Group <- ifelse(Bl2.pos[,5] > median(na.omit(Bl2.pos[,5])), 'high', 'low')
sfit <- survival::survfit(mysurv~Group,data = Bl2.pos)
survp <- survminer::ggsurvplot(sfit)
if(!dir.exists("HNSCC/figures/surv/Bl2/pos")) dir.create("HNSCC/figures/surv/Bl2/pos",recursive = T)
for(i in 5:(ncol(Bl2.pos)-1)){
  Group <- ifelse(Bl2.pos[,i] > median(na.omit(Bl2.pos[,i])), 'high', 'low')
  sfit <- survival::survfit(mysurv~Group,data = Bl2.pos)
  survp <- survminer::ggsurvplot(sfit, conf.int = F, pval = F)
  p1 <- survp$plot + annotate("text", label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                              x = 40, y = 0.15, size = 5) +
    annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
    guides(color = guide_legend(title = "HPV+")) +
    mytheme +
    scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
  ggsave(paste0("HNSCC/figures/surv/Bl2/pos/",colnames(Bl2.pos)[i],".pdf"),p1,width = 6,height = 5)
}
# Bl2.HPV-
Bl2.neg <- Bl2$Bl2.HPVneg
mysurv <- survival::Surv(Bl2.neg$clin.OS_MONTHS,Bl2.neg$Dead)
Group <- ifelse(Bl2.neg[,5] > median(na.omit(Bl2.neg[,5])), 'high', 'low')
sfit <- survival::survfit(mysurv~Group,data = Bl2.neg)
survp <- survminer::ggsurvplot(sfit)
if(!dir.exists("HNSCC/figures/surv/Bl2/neg")) dir.create("HNSCC/figures/surv/Bl2/neg",recursive = T)
for(i in 5:(ncol(Bl2.neg)-1)){
  Group <- ifelse(Bl2.neg[,i] > median(na.omit(Bl2.neg[,i])), 'high', 'low')
  sfit <- survival::survfit(mysurv~Group,data = Bl2.neg)
  survp <- survminer::ggsurvplot(sfit, conf.int = F, pval = F)
  p1 <- survp$plot + annotate("text", label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                              x = 40, y = 0.15, size = 5) +
    annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
    guides(color = guide_legend(title = "HPV-")) +
    mytheme +
    scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
  ggsave(paste0("HNSCC/figures/surv/Bl2/neg/",colnames(Bl2.neg)[i],".pdf"),p1,width = 6,height = 5)
}

load("HNSCC/out/CD45pPBMC.NK.Bl3.survdat.rda")
# Bl3.HPV+
Bl3.pos <- Bl3$Bl3.HPVpos
mysurv <- survival::Surv(Bl3.pos$clin.OS_MONTHS,Bl3.pos$Dead)
Group <- ifelse(Bl3.pos[,5] > median(na.omit(Bl3.pos[,5])), 'high', 'low')
sfit <- survival::survfit(mysurv~Group,data = Bl3.pos)
survp <- survminer::ggsurvplot(sfit)
if(!dir.exists("HNSCC/figures/surv/Bl3/pos")) dir.create("HNSCC/figures/surv/Bl3/pos",recursive = T)
for(i in 5:(ncol(Bl3.pos)-1)){
  Group <- ifelse(Bl3.pos[,i] > median(na.omit(Bl3.pos[,i])), 'high', 'low')
  sfit <- survival::survfit(mysurv~Group,data = Bl3.pos)
  survp <- survminer::ggsurvplot(sfit, conf.int = F, pval = F)
  p1 <- survp$plot + annotate("text", label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                              x = 40, y = 0.15, size = 5) +
    annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
    guides(color = guide_legend(title = "HPV+")) +
    mytheme +
    scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
  ggsave(paste0("HNSCC/figures/surv/Bl3/pos/",colnames(Bl3.pos)[i],".pdf"),p1,width = 6,height = 5)
}
# Bl2.HPV-
Bl3.neg <- Bl3$Bl3.HPVneg
mysurv <- survival::Surv(Bl3.neg$clin.OS_MONTHS,Bl3.neg$Dead)
Group <- ifelse(Bl3.neg[,5] > median(na.omit(Bl3.neg[,5])), 'high', 'low')
sfit <- survival::survfit(mysurv~Group,data = Bl3.neg)
survp <- survminer::ggsurvplot(sfit)
if(!dir.exists("HNSCC/figures/surv/Bl3/neg")) dir.create("HNSCC/figures/surv/Bl3/neg",recursive = T)
for(i in 5:(ncol(Bl3.neg)-1)){
  Group <- ifelse(Bl3.neg[,i] > median(na.omit(Bl3.neg[,i])), 'high', 'low')
  sfit <- survival::survfit(mysurv~Group,data = Bl3.neg)
  survp <- survminer::ggsurvplot(sfit, conf.int = F, pval = F)
  p1 <- survp$plot + annotate("text", label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                              x = 40, y = 0.15, size = 5) +
    annotate("text", label = paste0("n = 178"), x = 20, y = 0.22, size = 5) +
    guides(color = guide_legend(title = "HPV-")) +
    mytheme +
    scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
  ggsave(paste0("HNSCC/figures/surv/Bl3/neg/",colnames(Bl3.neg)[i],".pdf"),p1,width = 6,height = 5)
}

# ---------------------------------------------------------------------------
# seuratToLoom(obj, dir) - write a saved Seurat .rds out as a loom directory.
#
# A helper, not a pipeline step: step 02 sources this file and calls the
# function once per timepoint subset. The loom exports exist so the object can
# be read by python tools (scanpy/scVelo) without going through SeuratDisk.
#
# Unchanged from the cluster original except for this header.
# ---------------------------------------------------------------------------

seuratToLoom <- function(obj, dir){
  library(loomR)
  library(hdf5r)
  library(Seurat)
  
  seur <- readRDS(obj)
  dir.create(dir)
  
  if(!grepl('^seurat',class(seur)[1],ignore.case = T)){
    showNotification('Error: The selected object is of class ', class(seur)[1], ' but must be of class Seurat', type = 'error')
  }
  
  current_version <- 3
  if(strsplit(as.character(seur@version), split = '\\.')[[1]][1]<current_version){
    showNotification(paste0("Warning: The selected seurat object is out of date (version: ", seur@version, "). Updating to object now..."), type = 'error')
    seur <- UpdateSeuratObject(seur)
    showNotification(paste0("Update complete"), type = 'message')
  }
  
  # make sure percent.mito is added to the object
  if(any(names(seur@assays)=='RNA')){
    DefaultAssay(seur) <- 'RNA'
  }
  if(any(colnames(seur@meta.data) == 'percent.mito') == FALSE){
    if(any(grepl('^MT-', rownames(seur)))){
      seur$percent.mito <- PercentageFeatureSet(seur, pattern = '^MT-')
    }else if(any(grepl('^mt-', rownames(seur)))){
      seur$percent.mito <- PercentageFeatureSet(seur, pattern = '^mt-')
    }else if(any(grepl("^Mt-", rownames(seur)))){
      seur$percent.mito <- PercentageFeatureSet(seur, pattern = '^Mt-')
    }else{
      seur$percent.mito <- 0
      showNotification(paste0("Mitochondrial genes were not identified."), type = 'warning')
    }
  }
  
  project_dir <- paste0(dir,'/')
  
  #down_sample <- round(6000/length(unique(seur$seurat_clusters)))
  #seur <- subset(seur, downsample = down_sample, idents = 'seurat_clusters')
  # if(dim(seur)[2]>6000){
  #   seur <- subset(seur, cells = sample(Cells(seur), 6000))
  # }
  
  assays <- names(seur@assays)
  
  if(length(assays)>1){
    dims <- lapply(seur@assays, function(x){
      return(dim(x@data))
    })
    names(dims) <- assays
    assay_index <- NULL
    if(any(names(assays)=="RNA" & dims[["RNA"]][2] > 0)){
      for(i in 1:length(dims)){
        if(dims[[i]][2]!=dims[["RNA"]][2]){
          showNotification(paste0("Warning: the ", assays[i], "assay does not have the same number of cells (",dims[[i]][2],") as the RNA assay (", dims[["RNA"]][2],") and will be removed."), type = 'warning')
          assay_index <- c(assay_index,i)
        }
      }
    }else{
      first_assay <- dims[[1]][2]
      for(i in 2:length(dims)){
        if(first_assay!=dims[[i]][2]){
          showNotification(paste0("Warning: the ", assays[i], "assay does not have the same number of cells (",dims[[i]][2],") as the reference assay (", assay[1],": ",dims[["RNA"]][2],") and will be removed."), type = 'warning')
          assay_index <- c(assay_index,i)
        }
      }
    }
    if(!is.null(assay_index)){
      assays <- assays[-assay_index] 
    }
  }
  
  for(assay in assays){
    DefaultAssay(seur) <- assay
    
    filename <- paste0(project_dir,assay,".loom")
    loomR::create(filename = filename, data = seur[[assay]]@data, calc.count = F, overwrite = T)
    data <- loomR::connect(filename = filename, mode = "r+")
    data$link_delete('row_attrs/Gene')
    
    # add metadata
    meta.data <- seur@meta.data
    colnames(meta.data) <- paste0(colnames(meta.data),"_meta_data")
    
    data$add.col.attribute(as.list(meta.data))
    data$add.row.attribute(list(features = rownames(seur[[assay]])))
    
    # add reduction embeddings
    reduction_names <- names(seur@reductions)
    for(i in 1:length(reduction_names)){
      reductions <- as.data.frame(seur@reductions[[i]]@cell.embeddings)
      if(nrow(reductions)==0){
        showNotification(paste0("Warning: There was no data found for the ", reduction_names[i], " reduction. Skipping this reduction."), type = 'warning')
        next
      }
      assay_used <- tolower(seur@reductions[[i]]@assay.used)
      #if(!grepl(paste0('(?![a-z])(?<![a-z])(',assay_used,')'),reduction_names[i],perl = T,ignore.case = T)){
      if(!grepl(assay_used,reduction_names[i], ignore.case = T)){
        reduction_names[i] <- paste0(reduction_names[i],'_',assay_used)
      }
      n <- ncol(reductions)
      if(n>3){
        for(j in 2:3){
          tmp <- reductions[,1:j]
          tmp_name <- reduction_names[i]
          if(!grepl(paste0('(?![a-z])(?<![a-z])(',j,'d)'),tmp_name,perl=T,ignore.case = T)){
            tmp_name <- paste0(tmp_name,'_',j,'d')
          }
          if(grepl(' ', tmp_name)){
            tmp_name <- gsub(' ','_',tmp_name)
          }
          colnames(tmp) <- paste0(tmp_name,'_',1:j,'_reduction')
          data$add.col.attribute(as.list(tmp))
        }
      }else{
        if(!grepl(paste0('(?![a-z])(?<![a-z])(',n,'d)'),reduction_names[i],perl=T,ignore.case = T)){
          reduction_names[i] <- paste0(reduction_names[i],'_',n,'d')
        }
        if(grepl(' ', reduction_names[i])){
          reduction_names[i] <- gsub(' ','_',reduction_names[i])
        }
        colnames(reductions) <- paste0(reduction_names[i],'_',1:n,'_reduction')
        data$add.col.attribute(as.list(reductions))
      }
    }
    data$close_all()
  }
  return(1)
}


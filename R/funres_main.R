#'@description FunRes function is the main function of this code it uses tissue RNA-seq data and celltype annotation
#' to predict ligand-receptor interactions that maintain cellular phenotypes in different celltypes of the tissue.
#' @param data TPM count data for a tissue with celltype information as column names and rownames in form of Ensemblid_Genename.
#' @param anno.tbl Annotation dataframe for the TPM data with cell ids in first column and cell types in second.
#' @param LR Dataframe contatining ligand-receptor interactions with first column in form of ligand_receptor.
#' @param tfs Character vector with names of transcription factors.
#' @param species Species information (HUMAN or MOUSE)
#' @param sighot.cutoff cutoff parameter for SigHOtSpotter
#' @param sighot.percentile percentile parameter for SigHotSpotter
#' @param consv.thrs Conservation threshold
#' @param n Number of iterations durinh bootstrapping
#' @param ncores Number of cores 
#' @param cutoff1 Significance threshold for Cell-Cell interaction during bootstrapping (Not being used currently)
#' @param cutoff2 Sognificance threshold for LIgand-Receptor interaction during bootstrapping
#' @author Kartikeya Singh

  suppressPackageStartupMessages({
    require(textshape)
    require(ggplot2)
    require(selma)
    require(dplyr)
    require(doParallel)
    require(stringr)
    require(plyr)
    require(igraph)
    require(Matrix)
    require(reshape2)
    require(RSpectra)
    require(snow)
    require(taRifx)
    require(stats)
    require(gtools)
    require(data.table)
    require(rlist)
    require(pheatmap)
    require(clusterCrit)
    require(org.Mm.eg.db)
    require(org.Hs.eg.db)
    require(clusterProfiler)
    require(Seurat)
  })
  

FunRes <- function(data,anno.tbl,species,sighot.cutoff=0.1,sighot.percentile=70,consv.thrs=0.1,n=1000,ncores=4,z.score.cutoff=2,tissue.name,temp.folder.name,out.path,gen.markers=TRUE){
  
  #source("funres_library.R")
  #source("for_plotting_network_functions.R")
  

  system(paste0("mkdir ",out.path))
  system(paste0("mkdir ",out.path,"/",temp.folder.name))
  
  
  cat("Creating input parameters file\n\n")
  
  parms <- c("tissue.name","out.path","temp.folder.name","species","sighot.cutoff","sighot.percentile","consv.thrs","n","ncores")
  
  parms.file <- do.call("rbind",lapply(X = parms, function(x){
    return(paste0(x," = ",get(x)))
  }))
  
  write.table(x = parms.file,file = paste0(out.path,"/","input_parameters_",Sys.Date(),".txt"), sep = "\t", row.names = F, quote = F,col.names = F)
  
  cat("Tissue : ",tissue.name,"\n")
  cat(" Preparing data\n")
  
  # Prepare data
  
  colnames(anno.tbl) <- c("cell.name","cell.type")
  
  celltype.freq <- as.data.frame(table(anno.tbl$cell.type))
  rm.celltype <- as.character(celltype.freq[which(celltype.freq$Freq <= 10),][,1])
  
  
  new.colnames <- sapply(colnames(data),function(x) {
    make.names(as.character(anno.tbl[["cell.type"]][x == as.character(anno.tbl[["cell.name"]])]),
               unique=FALSE)
  })
  colnames(data) <- new.colnames
  
  # Dynamic loading of maxsub function shared object required for running the get.max.cluster function.
  # 
  # dyn.load("maxsubf.so")
  
  if(species == "MOUSE"){
    load("MOUSE_Background_data.RData",envir = .GlobalEnv)
    mmu_tab <- read.table("uniprot-filtered-organism-MMU.tab", sep = '\t', stringsAsFactors = FALSE, quote = "", comment.char = "")
    mmu_tab <- mmu_tab[mmu_tab$V4 != "",]
    secreted <- mmu_tab[base::grepl("Secreted",mmu_tab$V4),]
    secreted <- do.call("c",lapply(secreted$V3,function(x){
      out <- strsplit(x," ")[[1]]
      return(out)
    }))
    secreted <- secreted[!duplicated(secreted)]
    secreted <- sort(secreted)
    secreted <- setdiff(secreted,c("App", "Gnai2"))
    
    LR <- LR[which(LR$Ligand %in% secreted),]
    Ligands <- intersect(Ligands,secreted)
  } else {
    if(species == "HUMAN"){
      load("HUMAN_Background_data.RData",envir = .GlobalEnv)
      hsa_tab <- read.table("uniprot-filtered-organism-HSA.tab", sep = '\t', stringsAsFactors = FALSE, quote = "", comment.char = "")
      hsa_tab <- hsa_tab[hsa_tab$V3 != "",]
      secreted <- hsa_tab[base::grepl("Secreted",hsa_tab$V3),]
      secreted <- do.call("c",lapply(secreted$V2,function(x){
        out <- strsplit(x," ")[[1]]
        return(out)
      }))
      secreted <- secreted[!duplicated(secreted)]
      secreted <- sort(secreted)
      
      LR <- LR[which(LR$Ligand %in% secreted),]
      Ligands <- intersect(Ligands,secreted)
    } else {
      stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
    }
  }
  
  row.names(Background_signaling_interactome) <- 1:nrow(Background_signaling_interactome)
  non.rec.id <- which(Background_signaling_interactome$Source == dummy.var & !(Background_signaling_interactome$Target %in% Receptors))
  keep.id <- setdiff(row.names(Background_signaling_interactome),non.rec.id)
  Background_signaling_interactome <- Background_signaling_interactome[keep.id,]
  
  TF_TF_interactions <- remove.factors(TF_TF_interactions)
  
  all.pops <- setdiff(unique(anno.tbl$cell.type),rm.celltype)
  anno.tbl <- anno.tbl[which(anno.tbl$cell.type %in% all.pops),]
  
  data.lig.exp <- get.gene.expr(exp.tbl = data,genes = intersect(Ligands,row.names(data)),cell.type = all.pops)
  colnames(data.lig.exp) <- paste0("Ligand.",colnames(data.lig.exp))
  
  L.frame <- dplyr::inner_join(x = data.lig.exp,y = LR[,-1], by = c("Ligand.gene" = "Ligand"))
  L.frame <- L.frame[which(L.frame$Ligand.exp.perc > consv.thrs),]
  
  save(list = c("data","anno.tbl","data.lig.exp","L.frame"),file = paste0(out.path,"/",temp.folder.name,"/temp_data.RData"))
  
  rm(data.lig.exp)
  
  all.pops <- setdiff(all.pops,c("Unknown","Uknown"))
  
  invisible(lapply(all.pops, function(celltype){
    cell.exp.tbl <- data[,which(colnames(data) == celltype)]
    saveRDS(object = cell.exp.tbl,file = paste0(out.path,"/",temp.folder.name,"/temp_",celltype,".Rds"))
    rm(cell.exp.tbl)
  }))
  
  rm(data)
  
  invisible(gc())
  
  lapply(X = all.pops,FUN = function(celltype){
    
    if(length(list.files(path = paste0(out.path,"/",temp.folder.name),pattern = paste0(celltype,"_results.RData"))) == 0){
      
      cat(paste0(" Celltype : ",celltype,"\n"))
      
      cell.exp.tbl <- readRDS(file = paste0(out.path,"/",temp.folder.name,"/temp_",celltype,".Rds"))
      
      cat("   Finding maximum sum subcluster in expression space\n")
      
      #max.cluster.info <- get.max.cluster(exp.tbl = cell.exp.tbl)
      max.cluster.info <- get.cons.tfs(exp.tbl = cell.exp.tbl)
      if(class(max.cluster.info) == "list"){
        
        sig.input <- cell.exp.tbl[,max.cluster.info$tf.max.mat.cell]
        sig.input <- cbind.data.frame(row.names(sig.input),sig.input,stringsAsFactors = F)
        
        cons.tfs <- as.data.frame(unique(max.cluster.info$tf.count$Gene),stringsAsFactors = F)
        
        colnames(cons.tfs)[1] <- "Gene"
        
        cons.tfs$bool <- as.numeric(1)
        
        cat("   Starting SigHotSpotter analysis for the sub-cluster identified\n")
        
        hotspot.out <- SigHotSpotter_pipeline(idata = sig.input,species = species,cutoff = sighot.cutoff,DE_Genes = cons.tfs,percentile = sighot.percentile,ncores = ncores)
        
        save(list = c("max.cluster.info","sig.input","cons.tfs","hotspot.out"),file = paste0(out.path,"/",temp.folder.name,"/temp_",celltype,"_results.RData"))
        cat("   Saving results\n")
        
        rm(list = c("cell.exp.tbl","max.cluster.info","sig.input","cons.tfs","hotspot.out"))
        
      }else{
        cat(max.cluster.info,"\n")
        save(list = "max.cluster.info",file = paste0(out.path,"/",temp.folder.name,"/temp_",celltype,"_results.RData"))
        rm(max.cluster.info)
      }
    }
    

  })
  
  cat(" Collating Results\n")
  
  result.files <- list.files(path = paste0(out.path,"/",temp.folder.name),pattern = "_results.RData")
  
  all.pops <- as.character(sapply(result.files, function(file){
    # unlist(str_split(string = file,pattern = "_"))[2]
    out <- gsub(x = file,pattern = "temp_",replacement = "")
    out <- gsub(x = out,pattern = "_results.RData",replacement = "")
    return(out)
  }))
  
  collate <- lapply(X = all.pops,FUN = function(celltype){
    
    load(paste0(out.path,"/",temp.folder.name,"/temp_",celltype,"_results.RData"))
    
    if(class(max.cluster.info) == "list"){
      cat(" Collating data for ... ",celltype,"\n")
      # R.perturbs <- remove.factors(hotspot.out$perturbations)
      
      out <- list()
      
      hotspot.recs <- unique(c(hotspot.out$active,hotspot.out$inactive))
      
      iTF.targets <- hotspot.out$iTF.targets
      
      path.sums <- hotspot.out$path.sums
      
      if(class(path.sums) == "data.frame"){
        path.sums <- path.sums[which(path.sums$z.score > 0),]
      }
      
      if(!is.null(path.sums) & (class(path.sums) == "data.frame")){
        R.TF.info <- inner_join(x = path.sums, y = iTF.targets,by = "Gene")
        
        R.TF.info <- R.TF.info[which(R.TF.info$Receptor %in% hotspot.recs),]
        
        if(dim(R.TF.info)[1] > 0){
          
          cell.exp.tbl <- readRDS(file = paste0(out.path,"/",temp.folder.name,"/temp_",celltype,".Rds"))
          
          coexp.info <- gene.coexp(exp.tbl = sig.input[-1],gene.frame = R.TF.info[,c(1,2,5)],ncores = ncores)
          coexp.info$submat.coexp.perc <- coexp.info$coexp.count/(dim(sig.input)[2]-1)
          coexp.info$coexp.perc <- coexp.info$coexp.count/dim(cell.exp.tbl)[2]
          
          R.TF.info <- inner_join(x = R.TF.info, y = coexp.info, by = c("Receptor","Gene","Target"))
          R.TF.info <- R.TF.info[which(R.TF.info$coexp.perc > consv.thrs),]
          
          out$maxmat <- max.cluster.info
          out$hotspot <- hotspot.out
          
          if(dim(R.TF.info)[1] > 0){
            
            out$R.TF.info <- R.TF.info
            
            rec.expr <- get.gene.expr(exp.tbl = cell.exp.tbl,genes = unique(R.TF.info$Receptor),cell.type = celltype)
            
            colnames(rec.expr) <- paste0("Receptor.",colnames(rec.expr))
            
            LR.frame <- dplyr::inner_join(x = L.frame,y = rec.expr, by = c("Receptor" = "Receptor.gene"))
            
            # satisf.test <- check.R.satisfaction(LR.frame = LR.frame)
            
            # LR.frame <- dplyr::inner_join(x = LR.frame, y = satisf.test, by = c("Receptor" = "x"))
            
            LR.frame <- inner_join(x = LR.frame,y = R.TF.info,by = "Receptor")
            LR.frame <- unique(LR.frame[,c(1,2,4,6,7,9)])
            colnames(LR.frame) <- c("Lig.pop","Ligand","Lig.exp.perc","Receptor","Rec.pop","Rec.exp.perc")
            out$LR.frame <- LR.frame
            return(out)
            
          } else {
            return(NULL)
          }
        }else {
          return(NULL)
        }
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  })
  collate <- setNames(object = collate,nm = all.pops)
  collate <- list.clean(.data = collate,fun = is.null)
  
  tissue.LR <- data.frame()
  tissue.R.TF <- data.frame()
  
  all.pops <- names(collate)
  
  for (Celltype in all.pops) {
    if(exists(x = "LR.frame",where = collate[[Celltype]])){
      tissue.LR <- rbind.data.frame(tissue.LR,collate[[Celltype]][["LR.frame"]],stringsAsFactors = F)
    }
    if(exists(x = "R.TF.info",where = collate[[Celltype]])){
      tissue.R.TF <- rbind.data.frame(tissue.R.TF,cbind.data.frame(Celltype,collate[[Celltype]][["R.TF.info"]],stringsAsFactors = F),stringsAsFactors = F)
    }
  }
  
  collate$tissue.LR <- tissue.LR
  collate$tissue.R.TF <- tissue.R.TF
  
  
  saveRDS(object = tissue.LR,file = paste0(out.path,"/tissue_LR_no_bootstrap.Rds"))
  write.table(x = tissue.LR,file = paste0(out.path,"/tissue_LR_no_bootstrap.txt"), sep = "\t", row.names = F, quote = F)
  
  cat(" Running bootstrap for finding significant LR interactions\n")
  # cat(paste0(" Number of iterations : ",n,"\n"))

  load(paste0(out.path,"/",temp.folder.name,"/temp_data.RData"))

  #bootstrap.LR <- Bootstrap.interactions(data = data,anno.tbl = anno.tbl,LR = LR,ncores = ncores,n = n,consv.thrs = consv.thrs)
  # bootstrap.LR <- Bootstrap.interactions.minimal(data = data,anno.tbl = anno.tbl,tissue.LR = tissue.LR,ncores = ncores,n = n,consv.thrs = consv.thrs)
  cat("Performing bootstrapping ... \n")
  final <- Bootstrap.NewScoring(data = data,LR = tissue.LR,R.TF = tissue.R.TF,significance.cutoff = 0.7)

  saveRDS(object = final,file = paste0(out.path,"/tissue_LR_with_bootstrap.Rds"))

  # final <- final[which(final$Strength.z.score > z.score.cutoff),]

  write.table(x = final,file = paste0(out.path,"/tissue_LR_with_bootstrap.txt"), sep = "\t", row.names = F, quote = F)

  collate$final <- final

  output <- collate

  save(list = c("output"),file = paste0(out.path,"/output_",tissue.name,".RData"))

  partitions <- funres_partitions(output = output,out.path = out.path,tissue.name = tissue.name)

  if(gen.markers == TRUE){

    gen.markers(partitions = partitions,out.path = out.path,tissue.name = tissue.name)

  }
}

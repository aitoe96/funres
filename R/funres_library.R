set.seed(100)

############################
#  FunRes Functions
############################


#### prepare.data function prepares the data fro downstream analysis by using expression matrix and annotation table. Same as in commap.

#' @description The function prepares the input data for downstream analysis.
#' @param expr A dataframe with TPM normalized gene expression data, gene names as rownames in form Ensemblid_Genename and
#' cell ids as column names.
#' @param ph A dataframe with two columns. First column contains cell ids corresponding to column names in expr file and
#' second column containing cell type information.
#' @param ph.cell.id.col Name of column in ph with cell ids
#' @param ph.cell.type.col Nmae of column in ph with cell types.
#' @return A dataframe with gene names as rownames and cell types as column names.
#' 

prepare.data <- function(expr,ph,placeholder=NULL,ph.cell.id.col="cell.id",ph.cell.type.col="cell.type"){
  ## Prepares the data of a certain format (see below) to be supplied to the commap.
  ## Format:
  ## rownames(expr) have a format ENSMBLID_GENESYMBOL
  ## colnames(expr) are cell id's
  ## phenotype (ph) data have cell id column (ph.cell.id.col)
  ##  and cell type column (ph.cell.type.col)
  ## Remove ENSMBLID part
  ## print(grep("HLA",rownames(expr),value=T))    
  new.rownames <- sapply(strsplit(rownames(expr),"_"),function(x) x[2])
  
  if(!is.null(placeholder)){
    ## Just remove rows with un-identified gene symbols marked with the placeholder
    idx.to.remove <- new.rownames == placeholder
    expr <- expr[!idx.to.remove,]
    rownames(expr) <- new.rownames[!idx.to.remove]
    ## it will have a warning here automatically about the duplicated rownames
  }
  else {
    ## Just remove duplicates in the rownames, leaving potentially a single unknown placeholder
    idx.to.remove <- duplicated(new.rownames)
    cat(" ",length(which(idx.to.remove)),"duplicated entries found in input data and were removed\n")
    expr <- expr[!idx.to.remove,]
    rownames(expr) <- new.rownames[!idx.to.remove]
  }
  ## Substitue cell.id with cell.type information
  if(!(ph.cell.id.col %in% colnames(ph)) || !(ph.cell.type.col %in% colnames(ph)))
    stop("Indicate the proper names of cell.id and cell.type columns")
  new.colnames <- sapply(colnames(expr),function(x) {
    make.names(as.character(ph[[ph.cell.type.col]][x == as.character(ph[[ph.cell.id.col]])]),
               unique=FALSE)
  })
  colnames(expr) <- new.colnames
  return(expr)
}


get.gene.expr <- function(exp.tbl,anno.tbl,genes,cell.type=NULL){
  gene.exp.tbl <- exp.tbl[genes,,drop=FALSE]
  
  all.pops <- cell.type
  
  if(length(all.pops) == 1){
    cell.gene.exp <- gene.exp.tbl[,which(grepl(x = colnames(gene.exp.tbl),pattern = paste0("^",cell.type,"[\\.0-9]*$"),ignore.case = F)),drop=FALSE]
    cell.gene.abs.exp <- rowSums(cell.gene.exp)
    cell.gene.abs.exp <- cbind.data.frame(row.names(cell.gene.exp),cell.gene.abs.exp,stringsAsFactors = F)
    colnames(cell.gene.abs.exp) <- c("gene","abs.expr")
    
    cell.gene.bool <- as.data.frame(bool.data(exp.tbl = cell.gene.exp))
    
    cell.gene.cons <- rowSums(cell.gene.bool)
    cell.gene.cons <- cbind.data.frame(row.names(cell.gene.bool),cell.gene.cons,stringsAsFactors = F)
    colnames(cell.gene.cons) <- c("gene","cell.count")
    cell.gene.cons$exp.perc <- cell.gene.cons$cell.count/dim(cell.gene.bool)[2]
    cell.gene.out <- dplyr::inner_join(x = cell.gene.cons,y = cell.gene.abs.exp, by = "gene")
    cell.gene.out$celltype <- cell.type
    out <- cell.gene.out[,c(5,1:4)]
  }else{
    
    out <- do.call("rbind",lapply(X = all.pops,FUN = function(celltype1){
      # print(celltype1)
      # celltype.clms <- which(unlist(lapply(colnames(gene.exp.tbl),FUN = function(cl) paste0(unlist(strsplit(x = cl,split = "\\."))[1:2],collapse = "."))) == celltype1)
      # cell.gene.exp <- gene.exp.tbl[,celltype.clms]
      # cell.gene.exp <- gene.exp.tbl[,which(grepl(x = colnames(gene.exp.tbl),pattern = paste0("^",celltype1,"[\\.0-9]*$"),ignore.case = F))]
      cell.gene.exp <- gene.exp.tbl[,which(colnames(gene.exp.tbl) == celltype1)]
      cell.gene.abs.exp <- rowSums(cell.gene.exp) 
      cell.gene.abs.exp <- cbind.data.frame(row.names(cell.gene.exp),cell.gene.abs.exp,stringsAsFactors = F)
      colnames(cell.gene.abs.exp) <- c("gene","abs.expr")
      
      cell.gene.bool <- bool.data(exp.tbl = cell.gene.exp)
      
      cell.gene.cons <- rowSums(cell.gene.bool)
      cell.gene.cons <- cbind.data.frame(row.names(cell.gene.bool),cell.gene.cons,stringsAsFactors = F)
      colnames(cell.gene.cons) <- c("gene","cell.count")
      cell.gene.cons$exp.perc <- cell.gene.cons$cell.count/dim(cell.gene.bool)[2]
      
      cell.gene.out <- dplyr::inner_join(x = cell.gene.cons,y = cell.gene.abs.exp, by = "gene")
      cell.gene.out$celltype <- celltype1
      return(cell.gene.out[,c(5,1:4)])
      
    }))
  }
  invisible(gc())
  return(out)
}


#### bool.data function booleanizes the expression matrix based on expression threshold. It also removes genes not expression in any cell.
## Input : 1. exp.tbl : Gene expression count data matrix
##         2. expr.thrs : Expression threshold for booleanization (default 0).
## Output : Booleanized matrix

#' @description The function booleanizes the expression matrix based on expression threshold. 
#' It also removes genes not expression in any cell.
#' @param exp.tbl Gene expression dataframe.
#' @param expr.thrs Expression threshold.
#' @return Gene expression dataframe with booleanized expression.


bool.data <- function(exp.tbl,expr.thrs=0){
  bool.tbl <- matrix(data = as.numeric(exp.tbl > expr.thrs),nrow = nrow(exp.tbl),ncol = ncol(exp.tbl))
  colnames(bool.tbl) <- colnames(exp.tbl)
  row.names(bool.tbl) <- row.names(exp.tbl)
  bool.tbl <- bool.tbl[which(rowSums(bool.tbl) > 0),,drop=FALSE]
  return(bool.tbl)
}


#### maxsub2d function finds the rectangle in the matrix with maximum sum. The function and its fortran implementation has been from the adagio R package (after correction).
## Input : 1. A : Booleanized matrix (the matrix must have booleanization in 1/-1 form and not 1/0 form.)
## Output : Sum, index and the maximum sum sub-matrix.

#' @description The function finds the maximum sum sub-matrix in a matrix.
#' The function and its fortran implementation has been adapted from the adagio R package (after correction).
#' @param A Matrix with booleanization in 1/-1 form and not 1/0 form.
#' @return List containing Sum, index of the sub-matrix and the sub-matrix.
#'


maxsub2d <- function(A) {
  stopifnot(is.numeric(A), is.matrix(A))
  n <- nrow(A); m <- ncol(A)
  
  if (all(A <= 0))
    stop("At least on element of matrix 'A' must be positive.")
  if (all(A >= 0))
    return(list(sum = sum(A), inds = c(1, n, 1, m), submat = A))
  
  mi <- vector("integer", 4)
  S <- matrix(0, nrow = n+1, ncol = m)
  aa <- numeric(m)
  b <- 0.0
  
  fm <- 0.0
  R <- .Fortran("maxsub2f", as.numeric(A), as.numeric(S),
                as.integer(n), as.integer(m),
                fmax = as.numeric(fm), mind = as.integer(mi),
                as.numeric(aa), as.numeric(b))
  
  fm <- R$fmax
  mi <- R$mind
  
  invisible(gc())
  
  return(list(sum = fm, inds = mi,
              submat = A[mi[1]:mi[2], mi[3]:mi[4], drop = FALSE]))
}

get.cons.tfs <- function(exp.tbl,quantile = 0.95){
  tf.tbl <- exp.tbl[which(row.names(exp.tbl) %in% tfs),]
  tf.bool <- bool.data(exp.tbl = tf.tbl)
  freq_df <- apply(tf.bool,1,sum)
  freq_df <- data.frame(Gene = names(freq_df), Freq = freq_df,stringsAsFactors = FALSE)
  cutoff <- unname(quantile(freq_df$Freq[freq_df$Freq > 0],0.95))
  return(list(tf.max.mat.cell = colnames(exp.tbl), tf.count = freq_df[freq_df$Freq >= cutoff,]))
}

#### get.max.cluster function is used to find the sub-matrix with conserved TFs and receptors in the gene expression matrix.
## Input : 1. exp.tbl : Gene expression matrix

#' @description The function computes the sub-matrix with conserved TFs and receptors in the gene expression matrix
#' @param exp.tbl Gene expression matrix.
#' @return list containing TFs and receptors in the sub-matrix and their cell counts.

get.max.cluster <- function(exp.tbl){
  
  tf.tbl <- exp.tbl[which(row.names(exp.tbl) %in% tfs),]
  
  tf.bool <- bool.data(exp.tbl = tf.tbl)
  
  tf.bool[tf.bool == 0] <- -1
  
  tf.clust <- cluster_matrix(x = tf.bool,dim = "both")
  
  maxsum.tf <- maxsub2d(A = tf.clust)
  
  maxsum.tf.idx <- maxsum.tf$inds
  
  maxsum.matrix.tf <- tf.clust[maxsum.tf.idx[1]:maxsum.tf.idx[2],maxsum.tf.idx[3]:maxsum.tf.idx[4]]
  
  maxsum.matrix.tf[maxsum.matrix.tf == -1] <- 0
  
  if(abs(maxsum.tf.idx[3] - maxsum.tf.idx[4]) < 0.1 * as.numeric(dim(exp.tbl)[2])){
    return("No conserved TFs found.")
  }
  
  maxsum.tf.count <- as.data.frame(rowSums(maxsum.matrix.tf),stringsAsFactors = F)
  
  maxsum.tf.count <- cbind.data.frame(row.names(maxsum.tf.count),maxsum.tf.count,stringsAsFactors = F) 
  
  colnames(maxsum.tf.count) <- c("Gene","Cell.count")
  
  maxsum.tf.count$Cell.perc <- maxsum.tf.count$Cell.count/dim(exp.tbl)[2]
  maxsum.tf.count$TF.matrix.frac <- maxsum.tf.count$Cell.count/dim(maxsum.matrix.tf)[2]
  
  # rec.tbl <- exp.tbl[which(row.names(exp.tbl) %in% Receptors),]
  # 
  # rec.tbl <- rec.tbl[,which(colnames(rec.tbl) %in% colnames(maxsum.matrix.tf))]
  # 
  # rec.bool <- bool.data(rec.tbl)
  # 
  # rec.bool[rec.bool == 0] <- -1
  # 
  # rec.clust <- cluster_matrix(x = rec.bool,dim = "both")
  # 
  # maxsum.rec <- maxsub2d(A = rec.clust)
  # 
  # maxsum.rec.idx <- maxsum.rec$inds
  # 
  # maxsum.matrix.rec <- rec.clust[maxsum.rec.idx[1]:maxsum.rec.idx[2],maxsum.rec.idx[3]:maxsum.rec.idx[4]]
  # 
  # if(abs(maxsum.rec.idx[3] - maxsum.rec.idx[4]) < 0.1 * as.numeric(dim(exp.tbl)[2])){
  #   return("No conserved Receptors found.")
  # }
  # 
  # maxsum.matrix.rec[maxsum.matrix.rec == -1] <- 0
  # 
  # if(class(maxsum.matrix.rec) != "matrix"){
  #   maxsum.rec.count <- sum(maxsum.matrix.rec)
  #   maxsum.rec.count <- cbind.data.frame(row.names(rec.clust)[maxsum.rec.idx[1]:maxsum.rec.idx[2]],maxsum.rec.count,stringsAsFactors = F)
  #   
  # }else{
  #   maxsum.rec.count <- as.data.frame(rowSums(maxsum.matrix.rec),stringsAsFactors = F)
  #   maxsum.rec.count <- cbind.data.frame(row.names(maxsum.rec.count),maxsum.rec.count,stringsAsFactors = F) 
  #   
  # }
  # 
  # colnames(maxsum.rec.count) <- c("Gene","Cell.count")
  # 
  # maxsum.rec.count$Cell.perc <- maxsum.rec.count$Cell.count/dim(exp.tbl)[2]
  # maxsum.rec.count$TF.matrix.perc <- maxsum.rec.count$Cell.count/dim(maxsum.matrix.tf)[2]
  
  invisible(gc())
  
  return(list("tf.max.mat.cell" = colnames(maxsum.matrix.tf),"tf.count" = maxsum.tf.count))
}


gene.coexp <- function(exp.tbl,gene.set.1,gene.set.2,ncores=4){
  out.frame <- do.call("rbind",lapply(X = gene.set.1, FUN = function(x){
    out.row <- do.call("rbind",parallel::mclapply(mc.cores = ncores,X = gene.set.2, FUN = function(y){
      # print(x)
      # print(y)
      exp.bool <- bool.data(exp.tbl = exp.tbl)
      Combo <- paste0(x,"_",y)
      # print(head(row.names(exp.bool)))
      Set1.exp <- names(which(exp.bool[x,] == 1 ))
      Set2.exp <- names(which(exp.bool[y,] == 1))
      Set1.not.exp <- names(which(exp.bool[x,] == 0))
      Set2.not.exp <- names(which(exp.bool[y,] == 0))
      n11 <- length(intersect(Set1.exp,Set2.exp))
      n01 <- length(intersect(Set1.not.exp,Set2.exp))
      n10 <- length(intersect(Set1.exp,Set2.not.exp))
      n00 <- length(intersect(Set1.not.exp,Set2.not.exp))
      # asc.tbl[[R_TF]] <- matrix(c(n11,n01,n10,n00),nrow=2)
      # asc.test <- fisheSet1.test(x = asc.tbl[[R_TF]])
      asc.test <- fisheSet1.test(x = matrix(c(n11,n01,n10,n00),nrow=2))
      p.val <- asc.test$p.value
      coexp.count <- n11
      Set1.count <- length(Set1.exp)
      Set2.count <- length(Set2.exp)
      temp.row <- cbind.data.frame(Combo,Set1.count,Set2.count,coexp.count,p.val,stringsAsFactors = F)      
      return(temp.row)
    }))
    return(out.row)
  }))
  invisible(gc())
}

gene.frame.coexp <- function(exp.tbl,gene.set.frame,ncores=4){
  exp.bool <- bool.data (exp.tbl = exp.tbl)
  out.frame <- do.call("rbind",mclapply(X = 1:nrow(gene.set.frame), FUN = function(x){
    
    # print(head(row.names(exp.bool)))
    Set1.exp <- names(which(exp.bool[as.character(gene.set.frame[x,1]),] == 1 ))
    Set2.exp <- names(which(exp.bool[as.character(gene.set.frame[x,2]),] == 1))
    Set1.not.exp <- names(which(exp.bool[as.character(gene.set.frame[x,1]),] == 0))
    Set2.not.exp <- names(which(exp.bool[as.character(gene.set.frame[x,2]),] == 0))
    n11 <- length(intersect(Set1.exp,Set2.exp))
    n01 <- length(intersect(Set1.not.exp,Set2.exp))
    n10 <- length(intersect(Set1.exp,Set2.not.exp))
    n00 <- length(intersect(Set1.not.exp,Set2.not.exp))
    # asc.tbl[[R_TF]] <- matrix(c(n11,n01,n10,n00),nrow=2)
    # asc.test <- fisheSet1.test(x = asc.tbl[[R_TF]])
    asc.test <- fisher.test(x = matrix(c(n11,n01,n10,n00),nrow=2))
    p.val <- asc.test$p.value
    coexp.count <- n11
    Set1.count <- length(Set1.exp)
    Set2.count <- length(Set2.exp)
    temp.row <- cbind.data.frame(gene.set.frame[x,],Set1.count,Set2.count,coexp.count,p.val,stringsAsFactors = F)
    return(temp.row)
  },mc.cores = ncores))
  invisible(gc())
  return(out.frame)
}

#' @description The function calculates the coexpression of genes in each row of a given dataframe using the gene expression table.
#' @param exp.tbl Gene expression matrix to calculate coexpression of genes.
#' @param gene.frame A dataframe with genes. The coexpression is calculated for genes in each row.
#' @param ncores Number of cores to be used in the function. (Default : 4) 
#' @return gene.frame dataframe with an additional column corresponding to the number of cells coexpressing the genes. 

gene.coexp <- function(exp.tbl,gene.frame,ncores=4){
  gene.frame <- taRifx::remove.factors(gene.frame)
  row.names(gene.frame) <- 1:nrow(gene.frame)
  bool.exp.tbl <- bool.data(exp.tbl = exp.tbl)
  out <- do.call("rbind",lapply(X = 1:nrow(gene.frame), FUN = function(i){
    genes <- as.character(gene.frame[i,])
    len <- length(intersect(genes,row.names(bool.exp.tbl)))
    if (len == length(genes)) {
      bool.gene.tbl <- bool.exp.tbl[genes,]
      bool.sum <- colSums(bool.gene.tbl)
      coexp.count <- length(bool.sum[bool.sum == ncol(gene.frame)])
      out.row <- cbind.data.frame(gene.frame[i,],coexp.count,stringsAsFactors = F)
      return(out.row)
    }else{
      return(NULL)
    }
  }))
  return(out)
}

#'@description The function provides the satisfaction ratio for receptors. 
#'@details Satisfaction ratio is defined as the ratio of sum of absolute expression of all cognate ligands to the absolute expression of receptor. 
#'@param LR.frame A dataframe containing absolute gene exprssion data for all receptors and their cognate ligands.
#'@return A dataframe containing receptors and satisfaction ratios. 
#' 

check.R.satisfaction <- function(LR.frame){
  out <- data.frame()
  for (x in unique(LR.frame$Receptor)) {
    temp.LR.frame <- LR.frame[which(LR.frame$Receptor == x),]
    temp <- cbind.data.frame(x,sum(temp.LR.frame$Ligand.abs.expr),unique(temp.LR.frame$Receptor.abs.expr),stringsAsFactors = F)
    temp$satisf.ratio <- temp[,2]/temp[,3]
    out <- rbind.data.frame(out,temp[,c(1,4)],stringsAsFactors = F)
  }
  return(out)
}


get.top.lig.pop <- function(LR.frame,percent=0.9){
  new.frame <- data.frame()
  receptor <- unique(LR.frame$Receptor)
  for (x in receptor) {
    temp.l.frame <- LR.frame[which(LR.frame$Receptor == x),]
    if(dim(temp.l.frame)[1] > 0){
      temp.l.frame <- temp.l.frame[order(-temp.l.frame$Ligand.abs.expr),]
      l.sum <- sum(temp.l.frame$Ligand.abs.expr)
      temp.sum <- 0
      for (i in c(1:nrow(temp.l.frame))) {
        temp.sum <- temp.sum + as.numeric(temp.l.frame$Ligand.abs.expr[i])
        logic <- temp.sum > percent * l.sum
        if(logic){
          break()
        }
      }
      new.frame <- rbind.data.frame(new.frame,temp.l.frame[1:i,])
    }
  }
  return(new.frame)
}

#' The function computes the Median Absolute Deviation (MAD) values for the gene expression of input genes.
#' @param exp.tbl Gene expressoin table.
#' @param genes Set of genes for which MAD values are to be calculated.
#' @return Dataframe containing genes and corresponding MAD values. 

MAD.values <- function(exp.tbl,genes){
  out <- do.call("rbind",lapply(genes, function(y){
    temp <- cbind.data.frame(y,mad(x = as.numeric(exp.tbl[y,])))
    colnames(temp) <- c("Gene","MAD.value")
    return(temp)
  }))
}


############################
#  SighotSpotter Functions
############################


#' General pipeline for SigHotSpotter
#'
#' The function computes compatibility scores for signaling intermediates
#'
#' @param species Currently supported species: "HUMAN", "MOUSE"
#' @param input_data File name for input gene expression data
#' @param cutoff Maximum number of zero-value genes, above this cutoff the genes are excluded
#' @param DE_Genes_data Differential expression dataset (1 for up-regulated, -1 for down-regulated genes)
#' @param percentile Predicted intermediates are taken into account above this threshold
#' @param invert_DE If the differential expression should be inverted, default = FALSE
#' @param showprogress shows progress bar in shiny app if set to TRUE, set it to FALSE in batch mode without GUI, default = TRUE
#' @return Compatibility scores
#' @export
SigHotSpotter_pipeline <- function(species, idata, cutoff, DE_Genes, percentile, invert_DE = FALSE, showprogress = TRUE,ncores=4){
  
  subg=Data_preprocessing(input_data = idata,cutoff = cutoff,species = species)
  
  ## Calculate stationary distribution of the MC
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  
  prob.matrix <- summary(Steady_state_true$prob.matrix)
  
  edge.id.name <- cbind.data.frame(as_edgelist(graph = subg,names = T), as_edgelist(graph = subg,names = F),stringsAsFactors = F)
  colnames(edge.id.name) <- c("a","b","c","d")
  
  prob.matrix <- inner_join(x = edge.id.name,y = prob.matrix, by = c("c" = "i" , "d" = "j"))[,c(1,2,5)] 
  
  Steady_state_true <- Steady_state_true$SD
  
  ## Retrieves high probability intermediates
  int=high_probability_intermediates(x = Steady_state_true, intermediates = intermediates,percentile =  percentile)
  gintg=integrate_sig_TF(g = subg,x = Steady_state_true,deg = DE_Genes, non_interface_TFs = non_interface_TFs,TF_TF_interactions = TF_TF_interactions )
  
  # nTF=nonterminal_DE_TFs(g = gintg,deg = DE_Genes,non_interface_TFs = non_interface_TFs)
  
  target.TFs <- inner_join(x = DE_Genes, y = TF_TF_interactions, by = c("Gene" = "Source"))
  target.TFs <- target.TFs[which(target.TFs$Target %in% V(gintg)$name),]
  if(nrow(target.TFs) == 0){
    return(NULL)
  }
  target.coexp <- gene.coexp(exp.tbl = idata,gene.frame = target.TFs[,c(1,3)],ncores = ncores)
  target.coexp <- inner_join(x = target.coexp, y = target.TFs[,-2], by = c("Gene", "Target"))
  target.coexp$perc <- target.coexp$coexp.count/dim(idata)[2]
  target.coexp <- target.coexp[which(target.coexp$perc > 0.1 & target.coexp$Effect == 1),]
  iTF.target.info <- target.coexp[,c(1,2,4)]
  nTF <- iTF.target.info$Target
  
  nTF_scoring <- unique(nTF)
  names(nTF_scoring) <- nTF_scoring
  
  if (length(int) == 0){
    cat("No intermediates found. You may decrease the percentile in order to find intermediates.")
    return(NULL)
  }else if (length(nTF) == 0){
    cat("No non-terminal conserved TFs found at lower conservation also.")
    nTF <- DE_Genes
    return(NULL)
  }
  ## Computing compatibility scores
  #comp_score for each TF individually
  score <- lapply(nTF_scoring,function(x){return(comp_score_tf(x,int,gintg))})
  score <- lapply(nTF,function(x){score[[x]]})
  #score=lapply(nTF,comp_score_tf, int, gintg)
  if(is.null(score)){
    cat("No shortest path found. You may decrease the cutoff in order to find shortest path.")
    return(NULL)
  }else{
    #converting the nested list into a matrix whose row sum will give probability of each intermediate
    score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
    score_m_means=as.list(rowMeans(score_m))
    final_score=compatability_score(score_m_means,Steady_state_true,int)
    
    ## Computing networks for visualization
    # trimmed_score_I=.trimResults(final_score,F)
    # trimmed_score_A=.trimResults(final_score,T)
    # toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
    # toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
    # twoi=c(toiintI,toiintA)
    toiintA <- as.character(final_score[which(final_score$Activation_probability > 0.5),]$Gene)
    toiintI <- as.character(final_score[which(final_score$Activation_probability < 0.5),]$Gene)
    
    
    #pruning the integrated networks
    gintg.p=prun.int.g(gintg)
    
    #building networks for all intermediates for active signaling hotspots
    sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
    
    #building networks for inactive signaling hotspots
    sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
    
    #retrieve receptors for active intermediates & inactive:
    u.gr <- Reduce(graph.union,sp_int_A)
    if(class(u.gr) == "igraph"){
      aa <- incident(u.gr,"NICHE",mode="out")
      bb <- igraph::ends(u.gr,aa)
      active_receptors <- bb[,2]
      active_receptors <- intersect(x = active_receptors,y = LR$Receptor)
    }else{
      active_receptors <- NULL
    }
    
    u.gr <- Reduce(graph.union,sp_int_I)
    if(class(u.gr) == "igraph"){
      aa <- incident(u.gr,"NICHE",mode="out")
      bb <- igraph::ends(u.gr,aa)
      inactive_receptors <- bb[,2]
      inactive_receptors <- intersect(x = inactive_receptors, y = LR$Receptor)
    }else{
      inactive_receptors <- NULL
    }
    
    
    ### Link signaling molecules to receptors
    del <- incident(gintg.p,"NICHE",mode = "out")
    gintg.p.noNiche <- delete.edges(gintg.p,del)
    dists <- distances(gintg.p.noNiche,to = as.character(final_score$Gene), v = as.character(unique(c(active_receptors,inactive_receptors))), mode = "out", weights = NA, algorithm = "johnson")
    dists <- melt(dists)
    dists <- dists[which(!is.infinite(dists$value)),]
    
    recs.to.perturb <- unique(c(active_receptors,inactive_receptors))
    
    #idata_bool <- idata
    #idata_bool[idata_bool > 0] <- 1
    #save(idata,DE_Genes, file = paste0(gsub("\\.[0-9]*$","",colnames(idata)[2]),".RData"))
    #for(rec in recs.to.perturb){
    #  cep <- coexp.pairs(idata_bool,rec,Ligands,"T.cell")
    #  cps <- get.conserved.paths(idata_bool,g_tf_tf_lr, rec, unique(iTF.target.info$Gene), intersect(do.call(rbind,strsplit(cep$Pair,"_"))[,2],V(g_tf_tf_lr)$name), TF_TF_LR_interactions, cnsrv=0.8)
    #}
    
    #cat("   Starting receptor perturbation simulations\n")
    #perturb.tests <- do.call("rbind",mclapply(X = recs.to.perturb,FUN = function(rec){
    #  temp.test <- perturbCompare(r = rec,g = subg,SS = Steady_state_true,percentile = percentile,DE_Genes = DE_Genes
    #                              ,intermediates = intermediates,TF_TF_interactions = TF_TF_interactions,non_interface_TFs = non_interface_TFs)
    #  return(temp.test)
    #},mc.cores = ncores))
    
    cat("   Calculating shortest path weights\n")
    
    path.sums <- path.prob.sum(subg = subg,iTF = unique(iTF.target.info$Gene), Receptors = Receptors,prob.matrix = prob.matrix,ncores = ncores)
    
    # "perturbations" = as.data.frame(perturb.tests,stringsAsFactors = F)
    OUTPUT <- list("active"= active_receptors, "inactive"= inactive_receptors,"iTF.targets" = iTF.target.info,"final.score" = as.data.frame(final_score,stringsAsFactors = F),"path.sums" = path.sums, "rec.hotspot" = dists)
    return(OUTPUT)
  }
}


Data_preprocessing <- function(input_data,cutoff,species){
  
  b=input_data
  ##COMVERT CELLS WITH LESS THAT 1 FPKM TO 0
  b[b < 1] <- 0
  ##Add a new gene Dummy with expression value 1, Only works if Dummy is present in the initial network
  b[nrow(b)+1, ] <- c(dummy.var, rep(1,(ncol(b)-1)))
  #This is to convert chr into numericversion
  b[2:ncol(b)]<-as.data.frame(lapply(b[2:ncol(b)],as.numeric,b[2:ncol(b)]))
  ##Renaming the first column for finding the union
  colnames(b)[1] <- "Source"
  a=Background_signaling_interactome
  
  ## Removing vertices connected to Dummy which are not receptors
  non.recs <- which(a$Source == dummy.var & !(a$Target %in% Receptors))
  recs <- setdiff(c(1:nrow(a)),non.recs)
  a <- a[recs,]
  
  
  ab=join(a,b,by=c("Source"),type="left",match="first")
  colnames(ab)[3:ncol(ab)]="source"
  colnames(b)[1] <- "Target"
  ab1=join(a,b,by=c("Target"),type="left",match="first")
  names(ab1) <- NULL
  names(ab) <- NULL
  ab=ab[,4:ncol(ab)]
  ab1=ab1[,4:ncol(ab1)]
  ab=as.matrix(ab)
  ab1=as.matrix(ab1)
  ########Elementwise product
  g=ab * ab1
  ########Sum of elementwise product
  sum_product_expression=rowSums(g, na.rm = FALSE, dims = 1)
  g3=cbind(a,sum_product_expression)
  ########Calculation of precentage of cells expressed
  h=rowSums(g != 0)
  percent_expressed=(h*100)/ncol(ab)
  g3=cbind(g3,percent_expressed)
  g3[is.na(g3)]<-0
  #names(g3)=c("V1","V2","Effect","V3","V4")
  #g3[g3$percent_expressed <as.numeric(cutoff) ,3]=0
  #write.table(g3,out, sep="\t",row.names=FALSE)
  ######NETWORK preprocessing
  g <- graph.data.frame(as.data.frame(g3))
  #g=simplify(g,edge.attr.comb=list("first"))
  #DELETE those edges whoese sum is zero
  #del=E(g)[V3==0]
  #DELETE those edges whoese sum is zero (V3) OR whose average expression is lessthan args(3) (V4)
  #convert cutoff to numeric ##NOT DELETING NODES ANYMORE
  del=E(g)[sum_product_expression==0|percent_expressed<as.numeric(cutoff)]
  g <- delete.edges(g,del)
  #SINCE THE TFs AND RECEPTORS ARE ALREADY CONNECTED TO DUMMY, REMOVE ANY NODE THAT HAS ZERO in degree or zero out degree
  #To ensure reachability for the Markov chain
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
  #Select Nodes to be deleted
  del=V(g)[degree==0]
  #delete vertices from graph
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
    del=V(g)[degree==0]
  }
  #Same as above but remove nodes with with zero out degree
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
  #Select Nodes to be deleted
  del=V(g)[degree==0]
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
    del=V(g)[degree==0]
  }
  #####TO EXTRACT THE LARGEST STRONGLY CONNECTED COMPONENT
  members <- membership(clusters(g, mode="strong"))
  #l=lapply(unique(members), function (x) induced.subgraph(g, which(members == x)))
  #g=l[[1]]
  SCC <- clusters(g, mode="strong")
  subg <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
  subg=igraph::simplify(subg,edge.attr.comb=list("first"))
  subg
}

##STATIONARY DISTRIBUTION IN R ITSELF
##make a sparce-matrix from the edgelist
Markov_chain_stationary_distribution <- function(subg){ #The function takes the edgelist with probabilitys and computes the SS probability
  ####Write this subgraph as edgelist to make normal graph with ids ie.e names=F
  out <- list()
  transition_probability=as.data.frame(as_edgelist(subg,names = F),stringsAsFactors = FALSE)
  transition_probability$probability=paste(E(subg)$sum_product_expression)
  transition_probability[3]=as.numeric(transition_probability[[3]])
  myMatrix = sparseMatrix(i = transition_probability[1:nrow(transition_probability),1], j = transition_probability[1:nrow(transition_probability),2],x = transition_probability[1:nrow(transition_probability),3])
  #Making a stochastic matrix
  myMatrix = (myMatrix)/Matrix::rowSums((myMatrix))
  #write.table(myMatrix,"matrix.txt",sep="\t")
  #print(myMatrix)
  ##Eigen value of the sparse matrix
  #ev=eigen(Matrix::t(myMatrix))
  ##eigen value by Rspectra package
  el=eigs(Matrix::t(myMatrix),1,which="LR")
  ##eigen value that matches 1
  #match(1.0000,Re(round(ev$values, digits = 5)))
  #col=which.min(abs(ev$values - 1.00000000))
  ##STATIONARY DISTRIBUTION
  #SD=(abs(ev$vectors[,col]))/sum(abs(ev$vectors[,col]))
  #Stationary distribution from Rspectra eigen vectors
  SD=(abs(el$vectors))/sum(abs(el$vectors))
  SD=as.data.frame(SD,stringsAsFactors=FALSE)
  SD
  SD=cbind((as.data.frame(V(subg)$name,stringsAsFactors = FALSE)),SD)
  SD=as.data.frame(SD,stringsAsFactors=FALSE)
  colnames(SD)[1] <- "Gene"
  out_SS=paste("Steady_state",sep="")
  colnames(SD)[2] <- out_SS
  out$SD <- SD
  out$prob.matrix <- myMatrix
  return(out)
}

#----------------------------------------------------------------------------------
# TO CALCULATE THE COMPATABILITY SCORE FOR THE HIGH PROBABILITY INTERMEDIATES
#----------------------------------------------------------------------------------

high_probability_intermediates <- function(x, intermediates, percentile)
{
  intermediates=unique(join(intermediates,x,by=c("Gene"),type="inner"))
  if(nrow(intermediates) == 0){
    return(c())
  }
  #Selecting top 90 percentile of intermediates with steady state probabilities
  percentile=as.numeric(percentile)/100
  SS_90percentile=as.vector(quantile(intermediates[,2], as.numeric(percentile)))
  ##Shortlisting high-probability intermediates > 90 percentile of SS probability
  int=as.vector((subset(intermediates, intermediates[2] > SS_90percentile , select=c("Gene")))$Gene)
  int
}


#Function for integrating signaling and TF networks, g=signaling graph, x steady state vector, deg=differentially expressed genes
integrate_sig_TF <- function(g,x,deg, non_interface_TFs, TF_TF_interactions ){
  
  el=as_edgelist(g)
  graph_markov=as.data.frame(cbind(el,E(g)$Effect))
  colnames(graph_markov)=c("Source","Target","Effect")
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner")
  #Get the phenotype from the non_interface_DE_TFs and map it to the TF-TF interaction network
  colnames(non_interface_DE_TFs)[1]="Target"
  DE_TF_TF_interactions_target=join(TF_TF_interactions,non_interface_DE_TFs,by=c("Target"),type="left")
  DE_TF_TF_interactions_target=na.omit(DE_TF_TF_interactions_target)
  ab=DE_TF_TF_interactions_target
  #ab$Effect=ab$Effect*ab$Phenotype1 #second Phenotype will be opposite of this #must input second data for second phenotype
  names(ab)<-NULL
  ab=as.matrix(ab)
  ab[,3]=as.numeric(ab[,3])*as.numeric(ab[,4])
  ab=as.data.frame(ab)
  names(ab)=c("Source","Target","Effect","DEG")
  graph_markov$Effect=as.numeric(as.character(graph_markov$Effect))
  graph_markov=rbind(graph_markov,ab[1:3]) #merging the nTF interaction with appropriate sign Effect with the original graph
  colnames(x)[1] <- "Source"
  ab=join(graph_markov,x,by=c("Source"),type="left",match="first")
  colnames(ab)[3:ncol(ab)]="source"
  colnames(x)[1] <- "Target"
  ab1=join(graph_markov,x,by=c("Target"),type="left",match="first")
  names(ab1) <- NULL
  names(ab) <- NULL
  ab=ab[,4:ncol(ab)]
  ab1=ab1[,4:ncol(ab1)]
  #creating node SS as the edge property
  weight=as.numeric(as.matrix(ab))
  #edge_P=as.data.frame(weight)
  graph_markov=(cbind(graph_markov,weight))
  graph_markov[is.na(graph_markov)] <- 1  #Making TF-TF interactions dependent only on the expression status
  graph_markov$Effect=as.numeric(as.matrix((graph_markov$Effect)))
  g3 <- graph.data.frame(as.data.frame(graph_markov))
  #updating the graph attribute for the adjacency matrix i.e. product SS (weight) and effect
  E(g3)$weight=E(g3)$weight*E(g3)$Effect
  #deleting TF nodes with no indegree
  V(g3)$degree=igraph::degree(g3, v=V(g3), mode = c("in"))
  #Select Nodes to be deleted
  del=V(g3)[degree==0]
  #delete vertices from graph
  while(length(del)!=0)
  {
    g3 <- delete.vertices(g3,del)
    V(g3)$degree=igraph::degree(g3, v=V(g3), mode = c("in"))
    del=V(g3)[degree==0]
  }
  g3
}



#Function of shortlisting non-terminal differentially expressed genes
nonterminal_DE_TFs <- function(g,deg,non_interface_TFs){
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  #load non-terminal TFs
  nTF=non_interface_DE_TFs[1]
  names(nTF)=NULL
  nTF=as.vector(t(nTF))
  nTF<-intersect(nTF,V(g)$name) #Some TFs must still be missing in the final g3
  # if (length(nTF) == 0) stop ('No downstream TF found for the cutoff employed. You may decrease the cutoff in order to find TFs.') # decrease percentile cutoff
  # nTF
  if (length(nTF) == 0){
    cat("No downstream TF found for the cutoff employed. You may decrease the cutoff in order to find TFs.")
    return(NULL)
  }else{
    return(nTF)
  } 
}

#Function of shortlisting non-terminal differentially expressed genes But without stop command for building networks
#Function for classifying nTF as up or down regulated
up_down_tfs <- function(g,deg,non_interface_TFs)
{
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  #load non-terminal TFs
  nTF=non_interface_DE_TFs[1]
  names(nTF)=NULL
  nTF=as.vector(t(nTF))
  nTF<-intersect(nTF,V(g)$name) #Some TFs must still be missing in the final g3
  nTF
}

##Path_weight
#product_path_weight <- function(path, graph) sum(E(graph, path=path)$weight)/length(path)
product_path_weight <- function(path, graph){
  edge_weights=E(graph, path=path)$weight
  v_Edge_weight=sum(edge_weights)
  if (v_Edge_weight == 0) {
    return(0)
  } else{
    return(prod(edge_weights[edge_weights!=0]))
  }
}

path_weights <- function(path, graph) (E(graph, path=path)$weight)

##--Function for spliting and taking product of shortest path res file---
##This function takes ONE shortest path (x) and the adjacency matrix (l) as as input, and returns the product of SS probability of the intermediates in that shortest path.


##Function for Compatability score
##This function takes s="source" (one source gene), t="target" (a vector of target gene),g="graph", l="adjacency matrix" as input and finds the shortest paths and passes the argument to spsplit.
##Then it gets the product of intermediates of each path for a source and all its targets and returs its product as final output.

spcal_path_weight <- function(s,t,g){
  # if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
  # if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
  # paths=(get.all.shortest.paths(g, s, t, mode = c("out"), weights=NA)$res)
  # if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
  # if (length(paths) == 0){
  #   return (0)
  # }
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"), weights=NA)$res)
  if (length(paths) == 0){
    return(0)
  } 
  
  weight=lapply(paths,product_path_weight,g)
  #s=skewness(unlist(weight))
  s=weight_probability(unlist(weight))
  #s=sum(unlist(weight))
  return(s)
}

#Function to parallelize the Compatability Score calculations
parallel.function.comp <- function(i){
  mclapply(i,spcal,nTF,g3,l)
}

#FUNCTION FOR SKEWNESS # FROM MOMENTS PACKAGE
skewness <- function (x, na.rm = FALSE)
{
  if (is.matrix(x))
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x))
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

weight_probability <- function(x)
{
  x=unlist(x)
  x_pos=x[x>0]
  x_neg=x[x<0]
  x_tot=sum(abs(x_pos),abs(x_neg))
  #probability of the intermediate to be compatible: closer to 1 more compatible it is and closer t zero more incompatible it is
  #p_neg=sum(x_neg)/x_tot
  p_pos=sum(x_pos)/x_tot
}

comp_score_tf <- function(t,s,g) #x is a list of comp score
{
  comp_score=lapply(s,spcal_path_weight,t,g)
}

compatability_score <- function(x,y,int) #where x is compatability score as list and y is steady state
{
  x=unlist(x)
  x=cbind(as.data.frame(int),as.data.frame(unlist(x)))
  colnames(x)=c("Gene", "Activation_probability")
  x=join(x,y,by=c("Gene"),type="inner")
  x=x[order(x$Activation_probability,decreasing = TRUE),]
}


perturb_receptor <- function(r,g,intermediates,percentile,DE_Genes){
  
  #delete all receptors except 1
  
  g=delete.vertices(g,r)
  
  #To ensure reachability for the Markov chain
  
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
  
  #Select Nodes to be deleted
  
  del=V(g)[degree==0]
  
  #delete vertices from graph
  
  while(length(del)!=0)
    
  {
    
    g <- delete.vertices(g,del)
    
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
    
    del=V(g)[degree==0]
    
  }
  
  #Same as above but remove nodes with with zero out degree
  
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
  
  #Select Nodes to be deleted
  
  del=V(g)[degree==0]
  
  while(length(del)!=0)
    
  {
    
    g <- delete.vertices(g,del)
    
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
    
    del=V(g)[degree==0]
    
  }
  
  #####TO EXTRACT THE LARGEST STRONGLY CONNECTED COMPONENT
  
  members <- membership(clusters(g, mode="strong"))
  
  #l=lapply(unique(members), function (x) induced.subgraph(g, which(members == x)))
  
  #g=l[[1]]
  
  SCC <- clusters(g, mode="strong")
  
  subg <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
  
  #subg=simplify(subg,edge.attr.comb=list("first"))
  
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  Steady_state_true <- Steady_state_true$SD
  
  
  ## Retrieves high probability intermediates
  
  int=high_probability_intermediates(Steady_state_true, intermediates = intermediates, percentile = percentile)
  
  gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )
  
  nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)
  
  #comp_score for each TF individually
  
  score=lapply(nTF,comp_score_tf, int, gintg)
  
  if(is.null(score)){
    cat("Skipping rest of the analysis")
    return(NULL)
  }else{
    #converting the nested list into a matrix whose row sum will give probability of each intermediate
    
    score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
    
    score_m_means=as.list(rowMeans(score_m))
    
    final_score=compatability_score(score_m_means,Steady_state_true,int)
    
    
    
    ## Computing networks for visualization
    
    # trimmed_score_I=.trimResults(final_score,F)
    # 
    # trimmed_score_A=.trimResults(final_score,T)
    # 
    # toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
    # 
    # toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
    # 
    # twoi=c(toiintI,toiintA)
    
    
    
    #pruning the integrated networks
    
    # gintg.p=prun.int.g(gintg)
    
    
    
    #building networks for all intermediates for active signaling hotspots
    
    # sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
    
    
    
    #building networks for inactive signaling hotspots
    
    # sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
    
    
    
    #converting to visNetwork
    
    # vis_net_A <- lapply(sp_int_A,toVisNetworkData)
    # 
    # vis_net_I <- lapply(sp_int_I,toVisNetworkData)
    
    
    
    #for edge color
    
    # vis_net_A <- lapply(vis_net_A,vis.edge.color)
    # 
    # vis_net_I <- lapply(vis_net_I,vis.edge.color)
    
    ret_value = list(Steady_state_true=Steady_state_true,final_score=final_score)
    
    return(ret_value)
  }
}


perturbCompare<-function(r,g,SS,TF_TF_interactions,non_interface_TFs,percentile,intermediates,DE_Genes) {
  
  interface_TFs<-as.data.frame(setdiff(V(graph.data.frame(TF_TF_interactions))$name,as.character(as.vector(non_interface_TFs)$Gene)))
  names(interface_TFs)<-"Gene"
  
  rP=perturb_receptor(r,g,intermediates = intermediates,percentile = percentile,DE_Genes = DE_Genes)
  
  ssReal=merge(interface_TFs,SS)
  
  ssRp=merge(interface_TFs,rP$Steady_state_true)
  
  merge_ss=merge(ssReal,ssRp,by="Gene")
  
  merge_ss <- merge_ss[which(merge_ss$Gene %in% as.character(DE_Genes[,1])),]
  
  merge_ss=cbind(merge_ss,round(merge_ss[2]/merge_ss[3],digits = 3))
  colnames(merge_ss)[2:4] <- c("Steady.state.unpert","Steady.state.perturb","Steady.state.change") 
  
  merge_ss$Receptor <- r
  
  return(merge_ss)
  
}

path.prob.sum <- function(subg,iTF,Receptors,prob.matrix,ncores=4){
  # edges <- as.data.frame(as_edgelist(subg))
  # edge.weights <- as.data.frame(edge_attr(graph = subg)$sum_product_expression)
  # edge.weights <- cbind.data.frame(edges,edge.weights)
  
  Receptors <- intersect(Receptors,V(subg)$name)
  iTF <- intersect(iTF,V(subg)$name)
  
  if(length(iTF) == 0){
    cat("  No interface TFs in graph\n")
    return(NULL)
  }else{
    
    path.sum.frame <- do.call("rbind",mclapply(X = Receptors,FUN = function(rec){
      #print(rec)
      paths <- all_shortest_paths(graph = subg,from = rec, to = iTF,mode = "out",weights = NA)$res
      path.weights <- lapply(X = paths,FUN = function(path){
        path <- names(path)
        prod = 1
        if(length(path) > 1){
          for (j in 1:(length(path)-1)) {
            weight <- as.numeric(prob.matrix[which(prob.matrix$a == path[j] & prob.matrix$b == path[j+1]),][,3])
            prod = prod * weight
          }
        }
        names(prod) <- path[length(path)]
        return(prod)
      })
      names(path.weights) <- sapply(path.weights,function(x){names(x)})
      path.weights.sum <- aggregate(unlist(path.weights),by = list(Gene = names(path.weights)), FUN = function(x){sum(x,na.rm = TRUE)}, simplify = FALSE, drop = FALSE)
      path.weights.sum$Receptor <- rep(rec,nrow(path.weights.sum))
      path.weights.sum <- path.weights.sum[,c("Receptor","Gene","x")]
      colnames(path.weights.sum)[3] <- "path.weight.sum"
      path.weights.sum$path.weight.sum <- as.numeric(path.weights.sum$path.weight.sum)
      return(path.weights.sum)
    },mc.cores = ncores))
    vals  <- path.sum.frame$path.weight.sum
    path.sum.frame$z.score <- (vals-mean(vals))/sd(vals)
    # path.sum.frame <- path.sum.frame[which(path.sum.frame$z.score > 0),]
    return(path.sum.frame)
  }
}



############################
#  Bootstrap Functions
############################

Bootstrap.NewScoring <- function(data,LR,R.TF,significance.cutoff = 0.9){
  
  # LR <- LR[which(LR$Lig.exp.perc > 0.1 & LR$Rec.exp.perc > 0.1),]
  
  pops <- union(LR$Lig.pop,LR$Rec.pop)
  out.tissue.lr.scored <- list()
  for(lig.pop in pops){
    for(rec.pop in pops){
      # print(paste0(lig.pop,"_",rec.pop))
      lrsub <- LR[LR$Lig.pop == lig.pop & LR$Rec.pop == rec.pop,]
      if(nrow(lrsub) == 0){
        next
      }
      dat_lig <- data[lrsub$Ligand,which(colnames(data) == lig.pop)]
      dat_rec <- data[lrsub$Receptor,which(colnames(data) == rec.pop)]
      lig_means <- apply(dat_lig,1,function(x){mean(x[x>0])})
      rec_means <- apply(dat_rec,1,function(x){mean(x[x>0])})
      #lig_means <- apply(dat_lig,1,function(x){mean(x)})
      #rec_means <- apply(dat_rec,1,function(x){mean(x)})
      lrsub$score <- apply(lrsub,1,function(x){as.numeric(lig_means[x[2]])*as.numeric(rec_means[x[4]])})
      lrsub$score[is.na(lrsub$score)] <- 0
      test_lig <- data[LR$Ligand,which(colnames(data) == lig.pop)]
      test_lig_means <- apply(test_lig,1,function(x){mean(as.numeric(x[x>0]))})
      #test_lig_means <- apply(test_lig,1,function(x){mean(as.numeric(x))})
      test_rec <- data[LR$Receptor,which(colnames(data) == rec.pop)]
      test_rec_means <- apply(test_rec,1,function(x){mean(as.numeric(x[x>0]))})
      #test_rec_means <- apply(test_rec,1,function(x){mean(as.numeric(x))})
      scores <- unname(test_lig_means)*unname(test_rec_means)
      scores[is.na(scores)] <- 0
      #scores <- scores[!is.na(scores)]
      e <- ecdf(scores)
      lrsub$significance <- e(lrsub$score)
      z2 = (mean(scores)* 2)/sd(scores)
      lrsub$z2.diff <- lrsub$score-z2
      out.tissue.lr.scored[[paste0(lig.pop,"_",rec.pop)]] <- lrsub
    }
  }
  out.tissue.lr.scored.joined <- do.call("rbind",out.tissue.lr.scored)
  
  m <- R.TF[R.TF$z.score >= 0,c(1,2)]
  m <- m[!duplicated(m),]
  mm <- merge(out.tissue.lr.scored.joined,m,by.x = c("Rec.pop","Receptor"), by.y = c("Celltype","Receptor"))
  # mm <- inner_join(x = out.tissue.lr.scored.joined,y = m,by = c("Rec.pop" = "Celltype","Receptor"))
  mm <- mm[mm$significance >= significance.cutoff,]
  mm <- mm[,c(3,4,5,2,1,6,7,8)]
  # output$final <- mm
  # save(list = "output",file = paste0(path,"/output_",sample,"_NewScoring.RData"))
  return(mm)
}

############################
#  Feedback Loop Functions
############################


get.feedback.loops <- function(data,anno.tbl,cnsrv_rec = 0.1,cnsrv_lig = 0.1,zscore = 2,cnsrv_rl_pos = 0.4,cnsrv_rl_neg = 0.1,cond_prob_thr = 0.5,data.path,sample.name){
  
  cat(" Finding feedback loops\n")
  
  assign("TF_TF_LR_interactions",rbind.data.frame(TF_TF_interactions,Background_signaling_interactome[Background_signaling_interactome$Source %in% tfs & Background_signaling_interactome$Target %in% Ligands,]),envir = globalenv())
  assign("g_tf_tf_lr",graph_from_data_frame(TF_TF_LR_interactions, directed = TRUE, vertices = NULL),envir = globalenv())
  
  TF_TF_LR_interactions <- TF_TF_LR_interactions[base::intersect(which(as.character(TF_TF_LR_interactions$Source) %in% rownames(data)),which(as.character(TF_TF_LR_interactions$Target) %in% rownames(data))),]
  
  load(paste0(data.path,"/output_",sample.name,".RData"),envir = .GlobalEnv)
  
  subpops <- unique(intersect(output$final$Lig.pop,output$final$Rec.pop))
  
  cps_total <- list()
  # ligs <- unique(output$final$Ligand)
  
  for (pop in subpops) {
    
    recs <- unique(output$final$Receptor[which(output$final$Rec.pop == pop)])
    ligs <- unique(output$final$Ligand[which(output$final$Lig.pop == pop)])
    
    # print("Preparing gene expression data")
    cell.exp.tbl <- data[,which(colnames(data) == pop)]
    sig.input.bool <- cell.exp.tbl
    sig.input.bool[sig.input.bool > 0] <- 1
    
    cat("  Computing correlated receptors and ligands in ",pop,"\n")
    #ceps1 <- coexp.pairs.condprob(exp.tbl = sig.input.bool,rec.list = recs,lig.list = ligs,pop = pop,all_ligs = Ligands,all_recs = Receptors,num_samples = 10000)
    #ceps2 <- coexp.pairs.condprob.single(exp.tbl = sig.input.bool,rec.list = recs,lig.list = ligs,pop = pop,all_ligs = Ligands,all_recs = Receptors,num_samples = 10000)
    ceps <- coexp.pairs.coexpr(exp.tbl = sig.input.bool,rec.list = recs,lig.list = ligs,pop = pop,all_ligs = Ligands,all_recs = Receptors,num_samples = 10000)
    
    if(nrow(ceps) == 0){
      next
    }
    ceps$Receptor <- do.call(rbind,strsplit(ceps$Pair,"_"))[,1]
    ceps$Ligand <- do.call(rbind,strsplit(ceps$Pair,"_"))[,2]
    
    rl_pairs <- paste(ceps$Receptor,ceps$Ligand,sep = "_")
    
    recToTF <- output[[pop]]$R.TF.info[,c(1,2)]
    recToTF <- recToTF[!duplicated(recToTF),]
    recToTF <- recToTF[which(recToTF$Receptor %in% ceps$Receptor),]
    
    cat("  Computing conserved paths in ",pop,"\n")
    cps_total[[pop]] <- get.conserved.paths(expr_bool = sig.input.bool,g_tf_tf_lr = g_tf_tf_lr, rec = unique(ceps$Receptor),iTF.vec =  unique(output[[pop]]$R.TF.info$Gene[output[[pop]]$R.TF.info$Receptor %in% unique(ceps$Receptor)]), lig.vec = intersect(unique(ceps$Ligand),V(g_tf_tf_lr)$name),TF_TF_LR_interactions =  TF_TF_LR_interactions,recToTF = recToTF, cnsrv=0)
    cps_total[[pop]]$Pairs <- cps_total[[pop]]$Pairs[which(paste(cps_total[[pop]]$Pairs$Receptor,cps_total[[pop]]$Pairs$Ligand,sep = "_") %in% paste(ceps$Receptor,ceps$Ligand,sep = "_")),]
    cps_total[[pop]]$Neg_Pairs <- cps_total[[pop]]$Neg_Pairs[which(paste(cps_total[[pop]]$Neg_Pairs$Receptor,cps_total[[pop]]$Neg_Pairs$Ligand,sep = "_") %in% paste(ceps$Receptor,ceps$Ligand,sep = "_")),]
    cps_total[[pop]]$cerl <- ceps
  }
  
  popsizes <- table(anno.tbl$cell.type)
  for(n in names(cps_total)){
    colnames(cps_total[[n]]$Pairs) <- c("Receptor", "TF", "Ligand", "Frac.Pos")
    colnames(cps_total[[n]]$Neg_Pairs) <- c("Receptor", "TF", "Ligand", "Frac.Neg")
    cps_total[[n]]$Joined <- join(cps_total[[n]]$Pairs,cps_total[[n]]$Neg_Pairs)
    cps_total[[n]]$Joined[is.na(cps_total[[n]]$Joined)] <- 0
    #nsize <- as.numeric(popsizes[n])
    #pvals <- apply(cps_total[[n]]$Joined,1,function(x){prop.test(c(as.numeric(x[4])*nsize,as.numeric(x[5])*nsize),n = c(nsize,nsize),alternative = "greater")$p.value})
    #qvals <- p.adjust(pvals,method = "fdr")
    #cps_total[[n]]$Joined$pvals <- pvals
    #cps_total[[n]]$Joined$qvals <- qvals
  }
  
  saveRDS(object = cps_total,file = paste0(data.path,"/",sample.name,"_Cell_Internal_Ligand_Exp_Through_Recs.Rds"))
  
  # load(paste0(data.path,"/",sample.name,"_Cell_Internal_Ligand_Exp_Through_Recs.RData"))
  # 
  # ## Subset LR interactions ##
  # lr_ints_df <- output[["final"]][output[["final"]]$Lig.exp.perc >= cnsrv_lig & output[["final"]]$Rec.exp.perc >= cnsrv_rec & output[["final"]]$Strength.z.score >= zscore,]
  # lr_ints <- data.frame(Receptor=paste(lr_ints_df$Receptor,lr_ints_df$Rec.pop,sep="_"),Ligand=paste(lr_ints_df$Ligand,lr_ints_df$Lig.pop,sep="_"), stringsAsFactors = FALSE)
  # 
  # ## Subset RL interactions ##
  # cps_total_tmp <- cps_total
  # for(n in names(cps_total_tmp)){
  #   # idx <- which(paste(cps_total_tmp[[n]]$Joined$Receptor,cps_total_tmp[[n]]$Joined$Ligand,sep="_") %in% cps_total_tmp[[n]]$cerl$Pair[cps_total_tmp[[n]]$cerl$Cond.Prob >= cond_prob_thr])
  #   idx <- which(paste(cps_total_tmp[[n]]$Joined$Receptor,cps_total_tmp[[n]]$Joined$Ligand,sep="_") %in% cps_total_tmp[[n]]$cerl$Pair)
  #   if(length(idx) == 0){
  #     print(n)
  #     next
  #   }
  #   cps_total_tmp[[n]]$Joined <- cps_total_tmp[[n]]$Joined[idx,]
  #   cps_total_tmp[[n]]$Joined$Receptor <- paste(cps_total_tmp[[n]]$Joined$Receptor,n,sep="_")
  #   cps_total_tmp[[n]]$Joined$Ligand <- paste(cps_total_tmp[[n]]$Joined$Ligand,n,sep="_")
  # }
  # rl_ints_df <- do.call("rbind",lapply(cps_total_tmp,function(x){x$Joined}))
  # rl_ints_df <- rl_ints_df[rl_ints_df$Frac.Pos >= cnsrv_rl_pos & (is.na(rl_ints_df$Frac.Neg) | (rl_ints_df$Frac.Neg <= cnsrv_rl_neg)),]
  # #rl_ints_df <- rl_ints_df[rl_ints_df$qvals <= 0.00001 & rl_ints_df$Frac.Pos >= cnsrv_rl_pos & rl_ints_df$Frac.Neg <= cnsrv_rl_neg,]
  # #rl_ints_df <- rl_ints_df[which((rl_ints_df$Frac.Pos > rl_ints_df$Frac.Neg*2) | (rl_ints_df$Frac.Pos >=0.1 & is.na(rl_ints_df$Frac.Neg))),]
  # 
  # rl_ints <- rl_ints_df[,c("Receptor","Ligand")]
  # rl_ints <- rl_ints[!duplicated(rl_ints),]
  # 
  # total_ints <- data.frame(From=c(lr_ints$Ligand,rl_ints$Receptor),To=c(lr_ints$Receptor,rl_ints$Ligand),stringsAsFactors = FALSE)
  # tmp <- apply(total_ints,1,function(x){strsplit(x[1],"_")[[1]]})
  # total_ints$Gene1 <- tmp[1,]
  # tmp <- apply(total_ints,1,function(x){strsplit(x[2],"_")[[1]]})
  # total_ints$Gene2 <- tmp[1,]
  # all.lr <- unique(c(LR$Ligand,LR$Receptor))
  # 
  # final.lr <- unique(c(output$final$Ligand,output$final$Receptor))
  # 
  # Mbotc <- read.table("MBotCOnt.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = "",quote = "")
  # Mbotc$Symbol <- str_to_sentence(Mbotc$Symbol)
  # Mbotc <- Mbotc[Mbotc$Symbol %in% final.lr,]
  # all.lr.Mbotc <- Mbotc
  # all.lr.Mbotc <- all.lr.Mbotc[,c(5,3)]
  # #colnames(all.lr.Mbotc) <- c("SYMBOL","GO")
  # system(paste0("mkdir ",data.path,"/",sample.name,"_Mbotc"))
  # categories <- unique(Mbotc$ProcessID)
  # for (x in 1:length(categories)) {
  #   print(categories[x])
  #   temp.Mbotc <- Mbotc[Mbotc$ProcessID == categories[x],c(2,5),drop = FALSE]
  #   genes <- unique(temp.Mbotc$Symbol)
  #   temp.ints <- total_ints[which(total_ints$Gene1 %in% genes & total_ints$Gene2 %in% genes),]
  #   
  #   Mbotc.name <- gsub(temp.Mbotc$ProcessName,pattern = " ",replacement = "_")
  #   Mbotc.name <- gsub(Mbotc.name,pattern = "/",replacement = "_")
  #   Mbotc.name <- unique(Mbotc.name)
  #   
  #   if (nrow(temp.ints) > 0) {
  #     g_sub <- graph_from_data_frame(temp.ints[,c(1,2)],directed = TRUE)
  #     assign(x = "graph_ids",value = V(g_sub)$name,envir = .GlobalEnv)
  #     g_sub_new <- set.vertex.attribute(g_sub, "name", value=1:length(graph_ids))
  #     sink(paste0(data.path,"/",sample.name,"_Mbotc/",Mbotc.name,"_FeedbackLoops.txt"))
  #     test <- get.elementary.circuits(g_sub_new)
  #     sink()
  #   }
  # }
}

#######################################################################################################
# This is an implementation of Johnson's algorithm for finding all the elementary circuits of a       #
# directed graph.                                                                                     #
# An elementary circuit is a circuit where no vertex but the first and last vertex appears twice.     #
# Two elementary circuits are distinct if one is not a cyclic permutation of the other.               #
# The details on this algorithm may be found in;                                                      #
#                                                                                                     #
# Finding all the elementary circuits of a directed graph - Donald B. Johnson, 1974.                  #
#                                                                                                     #  
# Author: Corrie Jacobien Carstens                                                                    #
#######################################################################################################

# library(igraph)

# Function to obtain the strong component K with least vertex in subgraph of g induced by {s,s+1, ..., n}
# So the strong component of the subgraph on vertex s,s+1,...,n that contains the least vertex. 
# Note that this function may return NULL if no such component exists. 
get.induced.strong <- function(g, s){
  # Create the induced subgraph on {s, s+1, ..., n} and compute the strong components
  sg <- induced.subgraph(g, vids=s:vcount(g))
  sc <- clusters(sg, mode="strong")
  
  # Obtain the names for the remaining nodes - this has to be done to make sure that we use
  # the right order of nodes, we want to find the strong component with the least vertex.
  # Igraph always uses ids 1:n' for a graph, so we need to use the names. 
  ids <- as.numeric(get.vertex.attribute(sg, "name", 1:vcount(sg)))
  order <- sort(ids, index.return=T)$ix
  
  # Obtain the vertices of the strong component with the least vertex
  others <- c() 
  for(v in order)
    if(length(others) <= 1)
      others <- which(sc$membership == sc$membership[v])
  
  # If there is a strong component with more than 1 vertex, return this component
  if(length(others) > 1)
    return(induced.subgraph(sg, others))
  # Else return NULL
  else 
    return(NULL)
}

# This function returns a list where u is unblocked and all vertices in B(u)
# are unblocked (recursively)
unblock <- function(u, b, B){
  b[u] <- F
  for(w in B[[u]]){
    B[[u]] <- B[[u]][-which(B[[u]]==w)]
    if(b[w]){
      bB <- unblock(w,b,B)
      b <- bB$b
      B <- bB$B
    }
  }
  return(list(b=b, B=B))
}

# Prints out the circuits starting at s 
circuit <- function(s, v, Ak, B, b, f, stack, ids){
  stack <- c(stack, v)
  b[v] <- T
  for(w in neighbors(Ak, v, mode="out")){
    if(w==s){
      cat(graph_ids[sapply(c(stack,s), FUN=function(i){return(ids[i])})])
      cat("\n")
      f = T
    }else if (!b[w]){
      updated <- circuit(s,w,Ak,B,b,f,stack, ids)
      B <- updated$B
      b <- updated$b
      stack <- updated$stack      
      if(updated$f)
        f = T
    }
  }
  if(f){
    updated <- unblock(v, b, B)
    b <- updated$b
    B <- updated$B
  }else{for(w in neighbors(Ak, v, mode="out"))
    if (! v %in% B[[w]])
      B[[w]] <- c(B[[w]], v)
  }
  stack <- stack[-length(stack)]
  return(list(B=B, b=b, f=f, stack=stack))
}

# Main loop (L3)
get.elementary.circuits <- function(g){
  b <- rep(F, vcount(g))
  B <- vector("list", vcount(g))
  s = 1
  while(s < vcount(g)){
    Ak <- get.induced.strong(g,s)
    if(!is.null(Ak)){
      ids <- as.numeric(get.vertex.attribute(Ak, "name", 1:vcount(Ak)))
      s <- min(ids)
      for(i in ids){
        b[i] <- F
        B[[i]] <- numeric(0)
      }
      s_indx <- which(ids == s)
      circuit(s_indx, s_indx, Ak, B, b, F, numeric(0), ids)
      s <- s + 1
    }else
      s <- vcount(g) 
  }
}


#' Compute incompatible interactions.
#' @param ints Interaction table - 1st col: Source, 2nd col: Target, 3rd col: Effect.
#' @param expr data frame of expression values (only one column).
#' @return vector containing all incompatible interactions as data frame. 
get.incompatible.interactions <- function(ints,expr){
  incomp_idx <- c()
  incomp_idx <- c(incomp_idx,which(expr[as.character(ints$Source),1] == 0))
  incomp_idx <- c(incomp_idx,which((expr[as.character(ints$Source),1] == 1 & expr[as.character(ints$Target),1] == 0 & ints$Effect == 1)))
  incomp_idx <- c(incomp_idx,which((expr[as.character(ints$Source),1] == 1 & expr[as.character(ints$Target),1] == 1 & ints$Effect == -1)))
  incomp_df <- ints[incomp_idx,]
  if(nrow(incomp_df) > 0){
    return(unlist(as.list(t(incomp_df[,1:2]))))
  }
  return(c())
}

get.conserved.paths <- function(expr_bool, g_tf_tf_lr, rec, iTF.vec, lig.vec, TF_TF_LR_interactions, recToTF, cnsrv=0.1){
  pairs <- data.frame(Receptor=character(0),Ligand=character(0),cnsrv=numeric(0),stringsAsFactors = FALSE)
  pairs_non <- data.frame(Receptor=character(0),Ligand=character(0),cnsrv=numeric(0),stringsAsFactors = FALSE)
  iTF.vec <- iTF.vec[iTF.vec %in% V(g_tf_tf_lr)$name]
  c_total <- c()
  numCells <- c()
  c_total_non <- c()
  numCells_non <- c()
  for(c in 1:ncol(expr_bool)){
    idx_non <- rec[which(expr_bool[rec,c] == 0)]
    idx <- rec[which(expr_bool[rec,c] == 1)]
    todel <- get.edge.ids(g_tf_tf_lr,get.incompatible.interactions(TF_TF_LR_interactions,expr_bool[,c,drop=FALSE]))
    g_test <- delete.edges(g_tf_tf_lr,todel)
    paths <- do.call("c",sapply(iTF.vec[which(expr_bool[iTF.vec,c] == 1)],function(x){all_shortest_paths(g_test,x,intersect(V(g_test)$name,lig.vec[which(expr_bool[lig.vec,c] == 1)]),mode = "out")$res}))
    if(length(idx_non) > 0){
      if(length(paths) > 0){
        tmp_df <- as.data.table(table(unlist(lapply(paths,function(x){paste0(as.character(names(x)),collapse = " ")}))))
        tmp_df <- tmp_df[rep(1:nrow(tmp_df),length(idx_non))]
        tmp_df <- tmp_df[, "Rec" := do.call("c",lapply(idx_non,function(x){rep(x,nrow(tmp_df)/length(idx_non))}))]
        setcolorder(tmp_df,c(1,3,2))
        c_df <- rbindlist(list(c_total_non,tmp_df),use.names = FALSE)
        c_total_non <- c_df[,lapply(.SD,sum), by = .(V1,Rec)]
        #tmp_df <- do.call("rbind", replicate(length(idx_non), as.data.frame(table(unlist(lapply(paths,function(x){paste0(as.character(names(x)),collapse = " ")})))), simplify = FALSE))
        #tmp_df$Rec <- do.call("c",lapply(idx_non,function(x){rep(x,nrow(tmp_df)/length(idx_non))}))
        #tmp_df <- tmp_df[,c(1,3,2)]
        #c_df <- rbind(c_total_non,tmp_df)
        #c_total_non <- aggregate(c_df$Freq,by = list(Var1 = c_df$Var1, Rec = c_df$Rec),FUN="sum")
        #colnames(c_total_non) <- c("Var1","Rec","Freq")
      }
      if(length(which(expr_bool[lig.vec,c] == 1)) > 0){
        numCells_tmp <- rbind(numCells_non,as.data.frame(table(lig.vec[which(expr_bool[lig.vec,c] == 1)])))
        numCells_non <- aggregate(numCells_tmp$Freq,by = list(Var1 = numCells_tmp$Var1),FUN="sum")
        colnames(numCells_non) <- c("Var1","Freq")
      }
    }
    if(length(idx) > 0){
      if(length(paths) > 0){
        tmp_df <- as.data.table(table(unlist(lapply(paths,function(x){paste0(as.character(names(x)),collapse = " ")}))))
        tmp_df <- tmp_df[rep(1:nrow(tmp_df),length(idx))]
        tmp_df <- tmp_df[, "Rec" := do.call("c",lapply(idx,function(x){rep(x,nrow(tmp_df)/length(idx))}))]
        setcolorder(tmp_df,c(1,3,2))
        c_df <- rbindlist(list(c_total,tmp_df),use.names = FALSE)
        c_total <- c_df[,lapply(.SD,sum), by = .(V1,Rec)]
        #tmp_df <- do.call("rbind", replicate(length(idx), as.data.frame(table(unlist(lapply(paths,function(x){paste0(as.character(names(x)),collapse = " ")})))), simplify = FALSE))
        #tmp_df$Rec <- do.call("c",lapply(idx,function(x){rep(x,nrow(tmp_df)/length(idx))}))
        #tmp_df <- tmp_df[,c(1,3,2)]
        #c_df <- rbind(c_total,tmp_df)
        #c_total <- aggregate(c_df$Freq,by = list(Var1 = c_df$Var1, Rec = c_df$Rec),FUN="sum")
        #colnames(c_total) <- c("Var1","Rec","Freq")
      }
      if(length(which(expr_bool[lig.vec,c] == 1)) > 0){
        numCells_tmp <- rbind(numCells,as.data.frame(table(lig.vec[which(expr_bool[lig.vec,c] == 1)])))
        numCells <- aggregate(numCells_tmp$Freq,by = list(Var1 = numCells_tmp$Var1),FUN="sum")
        colnames(numCells) <- c("Var1","Freq")
      }
    }
    # if(expr_bool[rec,c] == 0){
    #   todel <- get.edge.ids(g_tf_tf_lr,get.incompatible.interactions(TF_TF_LR_interactions,expr_bool[,c,drop=FALSE]))
    #   g_test <- delete.edges(g_tf_tf_lr,todel)
    #   paths <- do.call("c",sapply(iTF.vec[which(expr_bool[iTF.vec,c] == 1)],function(x){all_shortest_paths(g_test,x,intersect(V(g_test)$name,lig.vec[which(expr_bool[lig.vec,c] == 1)]),mode = "out")$res}))
    #   tmp_df <- do.call("rbind", replicate(length(idx_non), as.data.frame(table(unlist(lapply(paths,function(x){paste0(as.character(names(x)),collapse = " ")})))), simplify = FALSE))
    #   tmp_df$Rec <- do.call("c",lapply(idx_non,function(x){rep(x,nrow(tmp_df)/length(idx_non))}))
    #   tmp_df <- tmp_df[,c(1,3,2)]
    #   c_df <- rbind(c_total_non,tmp_df)
    #   c_total_non <- aggregate(c_df$Freq,by = list(Var1 = c_df$Var1, Rec = c_df$Rec),FUN="sum")
    #   colnames(c_total_non) <- c("Var1","Rec","Freq")
    #   numCells_tmp <- rbind(numCells_non,as.data.frame(table(lig.vec[which(expr_bool[lig.vec,c] == 1)])))
    #   numCells_non <- aggregate(numCells_tmp$Freq,by = list(Var1 = numCells_tmp$Var1),FUN="sum")
    #   colnames(numCells_non) <- c("Var1","Freq")
    # }else{
    #   todel <- get.edge.ids(g_tf_tf_lr,get.incompatible.interactions(TF_TF_LR_interactions,expr_bool[,c,drop=FALSE]))
    #   g_test <- delete.edges(g_tf_tf_lr,todel)
    #   paths <- do.call("c",sapply(iTF.vec[which(expr_bool[iTF.vec,c] == 1)],function(x){all_shortest_paths(g_test,x,intersect(V(g_test)$name,lig.vec[which(expr_bool[lig.vec,c] == 1)]),mode = "out")$res}))
    #   c_df <- rbind(c_total,as.data.frame(table(unlist(lapply(paths,function(x){paste0(as.character(names(x)),collapse = " ")})))))
    #   c_total <- aggregate(c_df$Freq,by = list(Var1 = c_df$Var1),FUN="sum")
    #   colnames(c_total) <- c("Var1","Freq")
    #   numCells_tmp <- rbind(numCells,as.data.frame(table(lig.vec[which(expr_bool[lig.vec,c] == 1)])))
    #   numCells <- aggregate(numCells_tmp$Freq,by = list(Var1 = numCells_tmp$Var1),FUN="sum")
    #   colnames(numCells) <- c("Var1","Freq")
    #   #allpaths <- c()
    # }
  }
  c_total <- as.data.frame(c_total)
  
  
  if(nrow(c_total) != 0){
    c_total <- c_total[c_total[,1] != "",]
  }
  if(nrow(c_total) == 0){
    pairs <- data.frame(Receptor=character(),TF=character(),Ligand=character(),Frac=numeric(),stringsAsFactors = FALSE)
  }else{
    colnames(c_total) <- c("Var1","Rec","Freq")
    c_total$Ligand <- unlist(lapply(strsplit(as.character(c_total$Var1)," "),function(x){x[length(x)]}))
    c_total$TF <- unlist(lapply(strsplit(as.character(c_total$Var1)," "),function(x){x[1]}))
    c_total$Var1 <- c_total$Ligand
    c_total <- c_total[!duplicated(c_total),]
    c_merged <- merge(c_total,numCells,by = c("Var1"))
    c_merged$Frac <- c_merged$Freq.x/c_merged$Freq.y
    c_merged <- c_merged[which(paste(c_merged$Rec,c_merged$TF,sep = "_") %in% paste(recToTF$Receptor,recToTF$Gene,sep = "_")),]
    pairs <- c_merged[c_merged$Frac >= cnsrv,]
    pairs <- pairs[,c(2,5,4,7)]
    pairs <- aggregate(pairs$Frac,by=list(Rec = pairs$Rec, TF = pairs$TF, Ligand = pairs$Ligand),FUN = max)
    colnames(pairs) <- c("Receptor", "TF", "Ligand","Frac")
  }
  
  
  if(is.null(c_total_non) == FALSE){
    c_total_non <- as.data.frame(c_total_non)
    
    if(nrow(c_total_non) != 0){
      c_total_non <- c_total_non[c_total_non[,1] != "",]
    }
    if(nrow(c_total_non) == 0){
      pairs_non <- data.frame(Receptor=character(),TF=character(),Ligand=character(),Frac=numeric(),stringsAsFactors = FALSE)
    }else{
      colnames(c_total_non) <- c("Var1","Rec","Freq")
      c_total_non$Ligand <- unlist(lapply(strsplit(as.character(c_total_non$Var1)," "),function(x){x[length(x)]}))
      c_total_non$TF <- unlist(lapply(strsplit(as.character(c_total_non$Var1)," "),function(x){x[1]}))
      c_total_non$Var1 <- c_total_non$Ligand
      c_total_non <- c_total_non[!duplicated(c_total_non),]
      c_merged_non <- merge(c_total_non,numCells_non,by = c("Var1"))
      c_merged_non$Frac <- c_merged_non$Freq.x/c_merged_non$Freq.y
      c_merged_non <- c_merged_non[which(paste(c_merged_non$Rec,c_merged_non$TF,sep = "_") %in% paste(recToTF$Receptor,recToTF$Gene,sep = "_")),]
      pairs_non <- c_merged_non[c_merged_non$Frac >= cnsrv,]
      pairs_non <- pairs_non[,c(2,5,4,7)]
      pairs_non <- aggregate(pairs_non$Frac,by=list(Rec = pairs_non$Rec, TF = pairs_non$TF, Ligand = pairs_non$Ligand),FUN = max)
      colnames(pairs_non) <- c("Receptor", "TF", "Ligand","Frac")
    }
  }else{
    pairs_non <- data.frame(Receptor=character(),TF=character(),Ligand=character(),Frac=numeric(),stringsAsFactors = FALSE)
  }
  
  
  return(list("Pairs" = pairs,"Neg_Pairs" = pairs_non))
  #return(pairs[!(paste(pairs$Receptor,pairs$Ligand) %in% paste(pairs_non$Receptor,pairs_non$Ligand)),])
}

asc.tbl.perm <- function(dat.b.p,gene,lookup.gene,norm=TRUE,cnsrv=0.1,n.perm=1000){
  ## Make co-expression (association) tables between a pair
  ## of "gene" and "lookup.gene"
  ## Input: dat.b.p - expression table, bool and single pop
  asc.tbls <- vector("list",n.perm)
  for(i in 1:n.perm){
    names(asc.tbls)[i] <- paste(paste(gene,lookup.gene,sep="_"),as.character(i),sep=".")
    gON <- which(dat.b.p[gene,])
    gOFF <- which(!dat.b.p[gene,])
    s <- sample(dat.b.p[lookup.gene,])
    lgON <- which(s)
    lgOFF <- which(!s)
    n11 <- length(intersect(gON,lgON))
    n10 <- length(intersect(gON,lgOFF))
    n01 <- length(intersect(gOFF,lgON))
    n00 <- length(intersect(gOFF,lgOFF))
    asc.tbls[[i]] <- matrix(c(n11,n01,n10,n00),nrow=2)
    if(norm)
      asc.tbls[[i]] <- asc.tbls[[i]] / sum(asc.tbls[[i]])
  }
  return(asc.tbls)
}


asc.tbl <- function(dat.b.p,genes,lookup.genes,norm=TRUE,cnsrv=0.1){
  ## Make co-expression (association) tables between each pair
  ## of "genes" and "lookup.genes"
  ## Input: dat.b.p - expression table, bool and single pop
  asc.tbls <- vector("list",length(genes)*length(lookup.genes))
  N <- dim(dat.b.p)[2]
  for(i in seq_along(genes)){
    for(j in seq_along(lookup.genes)){
      if(!all(c(genes[i],lookup.genes[j]) %in% rownames(dat.b.p)))
        next ## two genes must be present in data
      if(any((apply(dat.b.p[c(genes[i],lookup.genes[j]),],1,sum)/N) < cnsrv))
        next ## two genes must be expressed in more then cnsrv% cells
      idx <- (i-1)*length(lookup.genes) + j
      names(asc.tbls)[idx] <- paste(genes[i],lookup.genes[j],sep="_")
      gON <- which(dat.b.p[genes[i],])
      gOFF <- which(!dat.b.p[genes[i],])
      lgON <- which(dat.b.p[lookup.genes[j],])
      lgOFF <- which(!dat.b.p[lookup.genes[j],])
      n11 <- length(intersect(gON,lgON))
      n10 <- length(intersect(gON,lgOFF))
      n01 <- length(intersect(gOFF,lgON))
      n00 <- length(intersect(gOFF,lgOFF))
      asc.tbls[[idx]] <- matrix(c(n11,n01,n10,n00),nrow=2)
      if(norm)
        asc.tbls[[idx]] <- asc.tbls[[idx]] / sum(asc.tbls[[idx]])
    }
  }
  ## if(sum(sapply(asc.tbls,is.null)))
  ##     cat("Warning: some genes were not found.\n")
  asc.tbls <- asc.tbls[!sapply(asc.tbls,is.null)]
  return(asc.tbls)
}

# cond.prob.test <- function(dat.p,r.genes,l.genes,sign.lev=0.05,cnsrv=0.1){
## Make a table of the conditional probability test of co-expression of a pair
## of genes.
## Input:
## dat.p -- population data matrix (1 population only). May not be booleanized.
## r.genes -- one set of genes ("receptors")
## l.genes -- another set of genes ("ligands")
## sign.lev -- significance level to keep the associations
## cnsrv -- minimum conservation level to expect from both genes in a pair
## ==========================================================================
#   dat.b.p <- dat.p == 1#make.exp.bool(dat.p,0.0)
#   n <- length(r.genes)*length(l.genes)
#   cond.p.tbl <- data.frame(Pair=rep(NA,n),
#                            Cond.Prob1=rep(NA,n),
#                            Cond.Prob2=rep(NA,n),
#                            P.Val1=rep(NA,n),
#                            P.Val2=rep(NA,n),
#                            Sign=rep(NA,n),
#                            Cnsrv1=rep(NA,n),
#                            Cnsrv2=rep(NA,n))
#   N <- dim(dat.b.p)[2]
#   for(i in seq_along(r.genes)){
#     for(j in seq_along(l.genes)){
#       if(!all(c(r.genes[i],l.genes[j]) %in% rownames(dat.b.p)))
#         next ## two genes must be present in data
#       if(any((apply(dat.b.p[c(r.genes[i],l.genes[j]),],1,sum)/N) < cnsrv))
#         next ## two genes must be expressed in more then cnsrv% cells
#       
#       idx <- (i-1)*length(l.genes) + j # order index
#       pair <- paste(r.genes[i],l.genes[j],sep="_") # pair name
#       a.tbl <- asc.tbl(dat.b.p,r.genes[i],l.genes[j],norm=FALSE,cnsrv=cnsrv)[[1]] # calc asc tbl
#       ## cond prob L=ON | R=ON
#       cond.p1 <- a.tbl[1,1] / (a.tbl[1,1] + a.tbl[1,2])
#       cond.p2 <- a.tbl[1,1] / (a.tbl[1,1] + a.tbl[2,1])
# 
#       a.tbl.p <- asc.tbl.perm(dat.b.p,r.genes[i],l.genes[j],norm=FALSE,cnsrv=cnsrv,n.perm=1000) # permuted asc tbl
#       
#       cond.p.dist1 <- unlist(lapply(a.tbl.p,function(t) t[1,1]/(t[1,1]+t[1,2])))
#       cond.p.den1 <- density(cond.p.dist1,n=1024)
#       p.val1 <- sum(cond.p.den1$y[cond.p.den1$x > cond.p1])/sum(cond.p.den1$y)
#       
#       cond.p.dist2 <- unlist(lapply(a.tbl.p,function(t) t[1,1]/(t[1,1]+t[2,1])))
#       cond.p.den2 <- density(cond.p.dist2,n=1024)
#       p.val2 <- sum(cond.p.den2$y[cond.p.den2$x > cond.p2])/sum(cond.p.den2$y)
#       
#       cnsrv1 <- sum(dat.b.p[r.genes[i],])/N # conservation 1
#       cnsrv2 <- sum(dat.b.p[l.genes[j],])/N # conservation 2
#       cond.p.tbl$Pair[idx] <- pair
#       cond.p.tbl$Cond.Prob1[idx] <- cond.p1
#       cond.p.tbl$Cond.Prob2[idx] <- cond.p2
#       cond.p.tbl$P.Val1[idx] <- p.val1
#       cond.p.tbl$P.Val2[idx] <- p.val2
#       cond.p.tbl$Sign[idx] <- (p.val1 < sign.lev & p.val2 < sign.lev) 
#       cond.p.tbl$Cnsrv1[idx] <- cnsrv1
#       cond.p.tbl$Cnsrv2[idx] <- cnsrv2
#     }
#   }    
#   cond.p.tbl <- cond.p.tbl[!is.na(cond.p.tbl$Pair),]
#   return(cond.p.tbl[order(-cond.p.tbl$Cond.Prob),])
# }


#' Compute receptor-ligand co-occurrance for a given population.
#' @param dat Full data matrix (may not be bool)
#' @param rec.list vector of (important) receptors for population pop.
#' @param lig.list vector of (important) ligands for population pop.
#' @param pop The population under consideration
#' @return vector containing all incompatible interactions in character format "Gene1->Gene2". 

coexp.pairs <- function(dat,rec.list,lig.list=Ligands,pop,pop.exclude=NULL){
  ## Calculate coexpression tables between receptors and ligands,
  ## based on the conditional probability test.
  ## Input:
  ## dat -- full data matrix (may not be bool)
  ## rec.list -- list of (important) receptors for each population
  ## lig.list -- list of (important) ligands per population to select them further
  ##             (if NULL, cognate ligands to all niche receptors taken)
  ## pop.exclude -- some populations to exclude (like "Unknown")
  ## ======================================================================
  
  # cat("Population under consideration:",pop,"\n")
  all.recs <- unique(unlist(rec.list)) # all receptors REDUNDANT
  ## Take all ligands that could be cognate to any rec in the niche
  all.ligs <- unique(unlist(lig.list)) # all ligands REDUNDANT
  ## cat("All receptors:",all.recs,"\n")
  cat("*** Population:",pop,"\n")
  cat("Computing correlated receptors and ligands\n")
  ## niche.recs <- setdiff(all.recs,rec.list[[pops[ip]]]) # niche receptors
  ## cat("*** Niche receptors:",niche.recs,"\n")
  rec.p <- rec.list
  if(!is.null(lig.list))
    lig.p <- lig.list
  else
    lig.p <- all.ligs
  
  st <- system.time({
    cp.tbl <- cond.prob.test(dat.p = dat,r.genes = rec.p,l.genes = lig.p)
  })
  cat("Test table time:\n")
  print(st)
  cp.tbl <- cp.tbl[cp.tbl$Sign,]
  return(cp.tbl)
}


coexp.pairs.pearson <- function(exp.tbl,rec.list,lig.list,pop,all_ligs,all_recs,num_samples = 1000){
  r_ligs <- sample(all_ligs,num_samples,replace = TRUE)
  r_recs <- sample(all_recs,num_samples,replace = TRUE)
  #Take random pairs
  tmp <- 2*exp.tbl[r_ligs,]-exp.tbl[r_recs,]
  tmp[tmp != 1] <- 0
  backg <- rowSums(tmp)/rowSums(exp.tbl[r_recs,])
  backg[is.nan(backg)] <- 0
  backg <- backg[!is.na(backg)]
  q <- ecdf(backg)
  out.frame <- data.frame(Pair=character(0),
                          Cond.prob=numeric(0),
                          P.Val=numeric(0),
                          Sign=logical(0),
                          stringsAsFactors = FALSE)
  i <- 1
  for(r in rec.list){
    tmp <- 2*exp.tbl[lig.list,]-exp.tbl[rep(r,length(lig.list)),]
    tmp[tmp != 1] <- 0
    tmp <- rowSums(tmp)/rowSums(exp.tbl[r,])
    tmp[is.nan(tmp)] <- 0
    tmp <- tmp[!is.na(tmp)]
    pairs <- paste(r,names(tmp),sep = "_")
    pval <- 1-q(unname(tmp))
    sign <- pval <= 0.05
    out.frame <- rbind.data.frame(out.frame,data.frame(Pair = pairs,
                                                       Cond.prob = unname(tmp),
                                                       P.Val = pval,
                                                       Sign = sign,
                                                       stringsAsFactors = FALSE),stringsAsFactors = FALSE)
  }
  return(out.frame)
}

coexp.pairs.condprob <- function(exp.tbl,rec.list,lig.list,pop,all_ligs,all_recs,num_samples = 1000){
  r_ligs <- sample(all_ligs,num_samples,replace = TRUE)
  r_recs <- sample(all_recs,num_samples,replace = TRUE)
  bg <- do.call("rbind",lapply(seq(1,num_samples),function(i){
    tab <- table(x = factor(simplify2array(exp.tbl[c(r_recs[i]),]),levels = c(0,1)), y=factor(simplify2array(exp.tbl[r_ligs[i],]),levels = c(0,1)))
    return(c(prop.table(tab,margin = 1)[2,2],prop.table(tab,margin = 2)[2,2]))
  }))
  bg <- cbind(bg,(bg[,1]+bg[,2])/2)
  q_1 <- ecdf(bg[,1])
  q_2 <- ecdf(bg[,2])
  q_b <- ecdf(bg[,3])
  n <- length(rec.list)*length(lig.list)
  out.frame <- data.frame(Pair=rep(NA,n),
                          Cond.Prob=rep(-1,n),
                          Prob1=rep(-1,n),
                          Prob2=rep(-1,n),
                          P.Val.1=rep(-1,n),
                          P.Val.2=rep(-1,n),
                          P.Val.b=rep(-1,n),
                          Sign.1=rep(FALSE,n),
                          Sign.2=rep(FALSE,n),
                          Sign.b=rep(FALSE,n),
                          Sign=rep(FALSE,n),
                          Cnsrv1=rep(-1,n),
                          Cnsrv2=rep(-1,n))
  i <- 1
  for(l in lig.list){
    for(r in rec.list){
      pair <- paste(r,l,sep = "_")
      #correlation <- cor.test(simplify2array(exp.tbl[l,]),simplify2array(exp.tbl[r,]))
      #est <- prop.table(table(x = factor(simplify2array(sig.input.bool[r,]),levels = c(0,1)), y=factor(simplify2array(sig.input.bool[l,]),levels = c(0,1))),margin = 1)[2,2]
      tab <- table(x = factor(simplify2array(exp.tbl[r,]),levels = c(0,1)), y=factor(simplify2array(exp.tbl[l,]),levels = c(0,1)))
      prob1 <- prop.table(tab,margin = 1)[2,2]
      prob2 <- prop.table(tab,margin = 2)[2,2]
      est <- (prob1 + prob2)/2
      #est <- correlation$estimate
      #pval <- correlation$p.value
      pval.1 <- 1-q_1(prob1)
      pval.2 <- 1-q_2(prob2)
      pval.b <- 1-q_b(est)
      sign.1 <- pval.1 <= 0.05
      sign.2 <- pval.2 <= 0.05
      sign.b <- pval.b <= 0.05
      sign <- sign.1 & sign.2 & sign.b
      Cnsrv1 <- unname(sum(exp.tbl[l,]))
      Cnsrv2 <- unname(sum(exp.tbl[r,]))
      out.frame[i,] <- c(pair,est,prob1,prob2,pval.1,pval.2,pval.b,sign.1,sign.2,sign.b,sign,Cnsrv1,Cnsrv2)
      i <- i+1
    }
  }
  return(out.frame[which(out.frame$Sign == TRUE),])
}

coexp.pairs.condprob.single <- function(exp.tbl,rec.list,lig.list,pop,all_ligs,all_recs,num_samples = 1000){
  r_ligs <- sample(all_ligs,num_samples,replace = TRUE)
  r_recs <- sample(all_recs,num_samples,replace = TRUE)
  bg <- sapply(seq(1,num_samples),function(i){
    tab <- table(x = factor(simplify2array(exp.tbl[c(r_recs[i]),]),levels = c(0,1)), y=factor(simplify2array(exp.tbl[r_ligs[i],]),levels = c(0,1)))
    return(prop.table(tab,margin = 1)[2,2])
  })
  q <- ecdf(bg)
  n <- length(rec.list)*length(lig.list)
  out.frame <- data.frame(Pair=rep(NA,n),
                          Cond.Prob=rep(-1,n),
                          P.Val=rep(-1,n),
                          Sign=rep(FALSE,n),
                          Cnsrv1=rep(-1,n),
                          Cnsrv2=rep(-1,n))
  i <- 1
  for(l in lig.list){
    for(r in rec.list){
      pair <- paste(r,l,sep = "_")
      #correlation <- cor.test(simplify2array(exp.tbl[l,]),simplify2array(exp.tbl[r,]))
      #est <- prop.table(table(x = factor(simplify2array(sig.input.bool[r,]),levels = c(0,1)), y=factor(simplify2array(sig.input.bool[l,]),levels = c(0,1))),margin = 1)[2,2]
      tab <- table(x = factor(simplify2array(exp.tbl[r,]),levels = c(0,1)), y=factor(simplify2array(exp.tbl[l,]),levels = c(0,1)))
      prob1 <- prop.table(tab,margin = 1)[2,2]
      est <- prob1
      #est <- correlation$estimate
      #pval <- correlation$p.value
      pval <- 1-q(est)
      sign <- pval <= 0.05
      Cnsrv1 <- unname(sum(exp.tbl[l,]))
      Cnsrv2 <- unname(sum(exp.tbl[r,]))
      out.frame[i,] <- c(pair,est,pval,sign,Cnsrv1,Cnsrv2)
      i <- i+1
    }
  }
  return(out.frame[which(out.frame$Sign == TRUE),])
}

coexp.pairs.coexpr<- function(exp.tbl,rec.list,lig.list,pop,all_ligs,all_recs,num_samples = 1000){
  r_ligs <- sample(all_ligs,num_samples,replace = TRUE)
  r_recs <- sample(all_recs,num_samples,replace = TRUE)
  #bg <- sapply(seq(1,num_samples),function(i){
  #  prop.table(table(x = factor(simplify2array(exp.tbl[c(r_recs[i]),]),levels = c(0,1)), y=factor(simplify2array(exp.tbl[r_ligs[i],]),levels = c(0,1))),margin = 1)[2,2]
  #})
  #q <- ecdf(bg)
  n <- length(rec.list)*length(lig.list)
  out.frame <- data.frame(Pair=rep(NA,n),
                          Cond.Prob=rep(-1,n),
                          P.Val=rep(-1,n),
                          Sign=rep(FALSE,n),
                          Cnsrv1=rep(-1,n),
                          Cnsrv2=rep(-1,n))
  i <- 1
  for(l in lig.list){
    for(r in rec.list){
      pair <- paste(r,l,sep = "_")
      #correlation <- cor.test(simplify2array(exp.tbl[l,]),simplify2array(exp.tbl[r,]))
      #est <- prop.table(table(x = factor(simplify2array(sig.input.bool[r,]),levels = c(0,1)), y=factor(simplify2array(sig.input.bool[l,]),levels = c(0,1))),margin = 1)[2,2]
      est <- table(x = factor(simplify2array(exp.tbl[r,]),levels = c(0,1)), y=factor(simplify2array(exp.tbl[l,]),levels = c(0,1)))
      est <- est[1,1]+est[2,2]
      #est <- correlation$estimate
      #pval <- correlation$p.value
      pval <- 1-pbinom(est,ncol(exp.tbl),0.5)
      #pval <- 1-q(est)
      sign <- pval <= 0.05
      Cnsrv1 <- unname(sum(exp.tbl[l,]))
      Cnsrv2 <- unname(sum(exp.tbl[r,]))
      out.frame[i,] <- c(pair,est/ncol(exp.tbl),pval,sign,Cnsrv1,Cnsrv2)
      i <- i+1
    }
  }
  return(out.frame[which(out.frame$Sign == TRUE),])
}

########################################
# Functional state partition functions
########################################



clust.h.pop <- function(mat,method = "ward.D",n = 5,min.cells=5){
  
  map <- pheatmap(mat = mat,clustering_method = method,silent = T)
  # max.h <- floor(max(map$tree_col$height))
  heights <- sort(map$tree_col$height,decreasing = T)
  heights <- heights[1:n]
  heights <- heights[!is.na(heights)]
  map$tree_col$height <- sort(map$tree_col$height,decreasing = F)
  
  cutree.info <- unlist(lapply(heights, function(h){
    h <- ceiling(h)
    h.tmp <- cutree(tree = map$tree_col,h = h)
    # h.tmp[,1] <- as.integer(h.tmp[,1])
    crit <- intCriteria(traj = as.matrix(t(mat)),part = h.tmp,crit = "dunn")
    out <- as.data.frame(unlist(crit))
    return(out)
  }))
  cutree.info <- cutree.info[!is.infinite(cutree.info)]
  max.h <- as.numeric(which(cutree.info == max(cutree.info[-1]))[1])
  
  checks <- unlist(lapply(1:max.h, function(h){
    h.tmp <- as.numeric(table(cutree(tree = map$tree_col,k = h)))
    # print(h.tmp)
    if(length(h.tmp) == length(h.tmp[h.tmp >= min.cells])){
      # print(h)
      return(h)
    }else{
      return(NULL)
    }
  }))
  return(max(checks))
}

clust.h.gene <- function(mat,method = "ward.D",n = 5,min.geneset=5){
  map <- pheatmap(mat = mat,clustering_method = method,silent = T)
  # max.h <- floor(max(map$tree_col$height))
  heights <- sort(map$tree_row$height,decreasing = T)
  heights <- heights[1:n]
  heights <- heights[!is.na(heights)]
  map$tree_row$height <- sort(map$tree_row$height,decreasing = F)
  
  cutree.info <- unlist(lapply(heights, function(h){
    # cat(h,"\n")
    h <- ceiling(h)
    h.tmp <- cutree(tree = map$tree_row,h = h)
    # h.tmp[,1] <- as.integer(h.tmp[,1])
    crit <- intCriteria(traj = as.matrix(mat),part = h.tmp,crit = "dunn")
    out <- as.data.frame(unlist(crit))
    return(out)
  }))
  cutree.info <- cutree.info[!is.infinite(cutree.info)]
  max.h <- as.numeric(which(cutree.info == max(cutree.info[-1]))[1])
  
  checks <- unlist(lapply(1:max.h, function(h){
    h.tmp <- as.numeric(table(cutree(tree = map$tree_row,k = h)))
    # print(h.tmp)
    if(length(h.tmp) == length(h.tmp[h.tmp >= min.geneset])){
      # print(h)
      return(h)
    }else{
      return(NULL)
    }
  }))
  return(max(checks))
}


pop.R.lig <- function(cps){
  cnsrv_rl_pos = 0.1
  cnsrv_rl_neg = 1
  
  R.lig <- do.call("rbind",lapply(names(cps), function(pop){
    temp <- cps[[pop]][["Joined"]]
    if(nrow(temp) > 0){
      temp$celltype <- pop
      return(temp)
    }else{
      return(NULL)
    }
  }))
  colnames(R.lig)[3] <- "Lig.out"
  R.lig <- R.lig[R.lig$Frac.Pos >= cnsrv_rl_pos & (is.na(R.lig$Frac.Neg) | (R.lig$Frac.Neg <= cnsrv_rl_neg)),c(1:3,6)]
  return(R.lig)
}  

pop.R.lig.inf <- function(cps){
  cnsrv_rl_pos = 0.4
  cnsrv_rl_neg = 0.1
  
  R.lig <- do.call("rbind",lapply(names(cps), function(pop){
    temp <- cps[[pop]][["Joined"]]
    if(nrow(temp) > 0){
      temp$celltype <- pop
      return(temp)
    }else{
      return(NULL)
    }
  }))
  colnames(R.lig)[3] <- "Lig.out"
  R.lig <- R.lig[R.lig$Frac.Pos >= cnsrv_rl_pos & (is.na(R.lig$Frac.Neg) | (R.lig$Frac.Neg <= cnsrv_rl_neg)) & R.lig$pvals < 0.05,c(1:3,8)]
  return(R.lig)
}  

tree.go <- function(tree.frame){
  tree.frame <- cbind.data.frame(row.names(tree.frame),tree.frame,stringsAsFactors = F)
  go <- lapply(unique(tree.frame[,2]), function(clust){
    geneset <- tree.frame[which(tree.frame[,2] == clust),1]
    genes <- unique(unlist(lapply(geneset, function(set){
      return(unlist(strsplit(x = set,split = "_",fixed = 3)))
    })))
    return(simplify(enrichGO(gene = genes,OrgDb = Org.db,keyType = "SYMBOL",ont = "BP",universe = background.genes)))
  })
  go <- setNames(object = go,nm = paste0("Cluster_",unique(tree.frame[,2])))
  return(go)
}

diff.check <- function(mat,k.col,k.row,sig.diff = 0.3){
  
  map <-  pheatmap(mat = mat,show_colnames = F,clustering_method = "ward.D",cutree_rows = k.row,cutree_cols = k.col,silent = T)
  
  col.data <- as.data.frame(cutree(tree = map$tree_col,k = k.col))
  col.data <- cbind.data.frame(row.names(col.data),col.data,stringsAsFactors = F)
  colnames(col.data) <- c("cell","cluster")
  
  
  row.data <- as.data.frame(cutree(tree = map$tree_row,k = k.row))
  row.data <- cbind.data.frame(row.names(row.data),row.data,stringsAsFactors = F)
  colnames(row.data) <- c("gene","cluster")
  
  
  data.means <- as.data.frame(do.call("rbind",lapply(1:k.row, function(row){
    rows <- row.data[which(row.data$cluster == row),1]
    row.means <- do.call("cbind",lapply(1:k.col, function(col,genes = rows){
      cols <- col.data[which(col.data$cluster == col),1]
      temp <- as.data.frame(mat[genes,cols])
      mean <- sum(rowSums(temp))/(ncol(temp) * nrow(temp))
      return(mean)
    }))
    return(row.means)
  })))
  
  check.diff <- do.call("cbind",lapply(1:ncol(data.means), function(i){
    
    temp <- do.call("rbind",lapply(1:ncol(data.means), function(j){
      
      col.diff <- data.means[,i]-data.means[,j]
      return(as.numeric(any(abs(col.diff) > sig.diff)))
    }))
    return(temp)
  }))
  
  return(sum(colSums(check.diff) < (ncol(check.diff) - 1))-1)
}


gen.map <- function(pop.data,pop.name,out.path,folder.name = "RTF_heatmaps"){
  
  if(dim(pop.data)[1] < 5){
    return(NULL)
  }else{
    k.row <- clust.h.gene(mat = pop.data,n = 10)
    k.col <- clust.h.pop(mat = pop.data,n = 10)
    
    col.diff <- diff.check(mat = pop.data,k.col = k.col,k.row = k.row)
    
    if(col.diff >= 0){
      k.col <- k.col - col.diff
    }
    
    if(k.col < 1){
      k.col = 1
    }
    
    map <-  pheatmap(mat = pop.data,show_colnames = F,clustering_method = "ward.D",cutree_rows = k.row,
                     cutree_cols = k.col,silent = T)
    
    tree.frame <- as.data.frame(cutree(tree = map$tree_col,k = k.col))
    # tree.frame <- rbind.data.frame(row.names(tree.frame),tree.frame,stringsAsFactors = F)
    
    col.anno <- as.data.frame(as.factor(cutree(tree = map$tree_col,k = k.col)))
    colnames(col.anno) <- colnames(tree.frame) <- "State"
    
    # map <-  pheatmap(mat = pop.data,annotation_col = col.anno,show_colnames = F,clustering_method = "ward.D",
    #                  cutree_rows = k.row,cutree_cols = k.col,
    #                  filename = paste0(out.path,"/",folder.name,"/",pop.name,".tiff"),silent = T)
    map <-  pheatmap(mat = pop.data,annotation_col = col.anno,show_colnames = F,clustering_method = "ward.D",
                     cutree_rows = k.row,cutree_cols = k.col,silent = T)
    
    return(list("pop.data" = as.data.frame(pop.data),"map" = map,"tree.frame" = tree.frame,"k.row" = k.row,"k.col" = k.col))
    
    # return(list("map" = map))
  }
}

funres_partitions <- function(output,out.path,tissue.name,folder.name = "RTF_heatmaps"){
  
  system(paste0("mkdir ",out.path,"/",folder.name))
  
  final <- output$final
  R.TF <- output$tissue.R.TF
  
  sig.R.TF <- unique(inner_join(x = final[,c(4,5)], y = R.TF[,c(1:3)], by = c("Rec.pop" = "Celltype","Receptor")))
  
  R.TF.pop <- lapply(unique(sig.R.TF$Rec.pop), function(pop){
    temp.pop <- sig.R.TF[which(sig.R.TF$Rec.pop == pop),c(1,3)]
    ints <- unlist(lapply(1:nrow(temp.pop),function(i){
      return(paste(temp.pop[i,],collapse = "_"))
    }))
    return(ints)
  })
  
  R.TF.pop <- setNames(object = R.TF.pop,unique(sig.R.TF$Rec.pop))
  
  ### Get coexpression matrix for each population based on the gene sets.
  
  cat(" Creating genesets coexpression matrices\n")
  
  matrix.data <- lapply(X = names(R.TF.pop), function(pop){
    
    cell.exp <- readRDS(file = paste0(out.path,"/temp/temp_",pop,".Rds"))
    cell.exp[cell.exp > 0] <- 1
    
    data <- do.call("cbind",lapply(X = 1:ncol(cell.exp),FUN = function(cell.id){
      col <- unlist(lapply(X = R.TF.pop[[pop]],FUN = function(int,cell = cell.id){
        int <- unlist(strsplit(x = int,split = "_",fixed = 2))
        if(sum(cell.exp[int,cell]) == 2){
          return(1)
        }else{
          return(0)
        }
      }))
      col <- as.data.frame(col)
      row.names(col) <- R.TF.pop[[pop]]
      return(col)
    }))
    colnames(data) <- colnames(cell.exp)
    return(data)
  })
  matrix.data <- setNames(object = matrix.data,nm = names(R.TF.pop))
  
  cat(" Generating partitioned gene sets\n")
  
  pop.map <- lapply(names(matrix.data), function(pop){
    cat(" Celltype : ",pop,"\n")
    return(gen.map(pop.data = matrix.data[[pop]],pop.name = pop,out.path = out.path,folder.name = folder.name))
  })
  pop.map <- setNames(object = pop.map,nm = names(matrix.data))
  
  save(list = c("matrix.data","pop.map"),file = paste0(out.path,"/",tissue.name,"_heatmap_data.RData"))
  return(list("matrix.data" = matrix.data, "pop.map" = pop.map))
}


gen.markers <- function(partitions,out.path,tissue.name){
  
  cat("Identifying cluster markers (Bimod test)\n")
  
  bimod.markers <- lapply(X = names(partitions$pop.map),FUN = function(pop){
    
    if(!is.null(partitions$pop.map[[pop]])){
      
      cat(" Celltype : ",pop,"\n")
      
      tree.frame <- as.data.frame(as.factor(cutree(tree = partitions$pop.map[[pop]]$map$tree_col,
                                                   k = partitions$pop.map[[pop]]$k.col)))
      colnames(tree.frame)[1] <- "State"
      
      data <- readRDS(paste0(out.path,"/temp/temp_",pop,".Rds"))
      
      data <- data[intersect(row.names(data),c(Receptors,tf.db$Symbol)),]
      
      data <- CreateSeuratObject(counts = data,meta.data = tree.frame)
      Idents(data) <- data$State
      
      markers <- suppressWarnings(FindAllMarkers(object = data,test.use = "bimod"))
      markers <- markers[which(markers$p_val_adj < 0.05),]
      
      return(markers)
    }
  })
  bimod.markers <- setNames(object = bimod.markers,nm = names(partitions$pop.map))
  
  
  cat("Performing proportion of success test\n")
  
  gene.prop.test <- lapply(X = names(partitions$pop.map),FUN = function(pop){
    
    cat(" Celltype : ",pop,"\n")
    if(!is.null(partitions$pop.map[[pop]])){
      matrix <- readRDS(paste0(out.path,"/temp/temp_",pop,".Rds"))
      matrix <- matrix[intersect(row.names(matrix),c(Receptors,tf.db$Symbol)),]
      matrix[matrix > 0 ] <- 1
      tree <- partitions$pop.map[[pop]][["tree.frame"]]
      
      # tree <- tree[order(tree$State,decreasing = F),]
      
      state.lengths <- as.numeric(table(tree$State))
      
      state.test <- lapply(X = 1:nrow(matrix), function(i){
        row <- matrix[i,]
        exp <- do.call("c",lapply(X = unique(tree$State), function(state){
          return(sum(row[,row.names(tree)[which(tree$State == state)]]))
        }))
        # cat(exp,"\n")
        test <- suppressWarnings(prop.test(x = exp,n = state.lengths,correct = F))
        return(test)
      })
      
      assoc.test <- as.data.frame(do.call("rbind",lapply(1:length(state.test), function(j){
        return(c(row.names(matrix)[j],state.test[[j]][["p.value"]],state.test[[j]][["estimate"]]))
      })))
      colnames(assoc.test) <- c("gene","p.value",paste0("Cluster_",unique(tree$State)))
      
      assoc.test$p.value <- as.numeric(as.character(assoc.test$p.value))
      
      assoc.test <- assoc.test[which(assoc.test$p.value < 0.05),]
      
      return(assoc.test)
    }else{
      return(NULL)
    }
  })
  gene.prop.test <- setNames(object = gene.prop.test,nm = names(partitions$pop.map))
  
  save(list = c("bimod.markers","gene.prop.test"),file = paste0(out.path,"/",tissue.name,"_markers_test.RData"))
}







#' @title FunRes: Functional Resolution Pipeline
#' @description
#' FunRes is the main function of the package. It uses tissue scRNA-seq data
#' and cell-type annotations to predict ligand–receptor interactions that
#' maintain cellular phenotypes across cell types in a tissue.
#'
#' @param data A numeric matrix (genes × cells), rownames = EnsemblID_GeneName.
#' @param anno.tbl A data.frame with two columns: first = cell IDs, second = cell types.
#' @param species Character; either "HUMAN" or "MOUSE".
#' @param sighot.cutoff Numeric; cutoff for SigHotSpotter (default 0.1).
#' @param sighot.percentile Numeric; percentile for SigHotSpotter (default 70).
#' @param consv.thrs Numeric; conservation threshold (default 0.1).
#' @param n Integer; number of bootstrapping iterations (default 1000).
#' @param ncores Integer; number of cores to use (default 4).
#' @param z.score.cutoff Numeric; z-score cutoff after bootstrap (default 2).
#' @param tissue.name Character; label for your tissue (used in output filenames).
#' @param temp.folder.name Character; name of temporary subfolder (e.g. "tmp").
#' @param out.path Character; path where results will be written.
#' @param gen.markers Logical; whether to generate marker heatmaps (default TRUE).
#'
#' @return Invisibly returns a list with:
#'   - `tissue.LR_no_bootstrap`
#'   - `tissue.LR_with_bootstrap`
#'   - per-cell-type `R.TF.info`, `LR.frame`, etc.
#'
#' @export
FunRes <- function(
  data,
  anno.tbl,
  species,
  sighot.cutoff     = 0.1,
  sighot.percentile = 70,
  consv.thrs        = 0.1,
  n                 = 1000,
  ncores            = 4,
  z.score.cutoff    = 2,
  tissue.name,
  temp.folder.name,
  out.path,
  gen.markers       = TRUE
) {
  # Create output directories
  dir.create(out.path, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out.path, temp.folder.name), showWarnings = FALSE)

  cat("Creating input parameters file\n\n")
  parms <- c("tissue.name","out.path","temp.folder.name","species",
             "sighot.cutoff","sighot.percentile","consv.thrs","n","ncores","z.score.cutoff")
  parms.file <- do.call(rbind, lapply(parms, function(x) paste0(x, " = ", get(x))))
  write.table(parms.file, file = file.path(out.path, paste0("input_parameters_", Sys.Date(), ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat("Tissue : ", tissue.name, "\n")
  cat(" Preparing data\n")

  # Standardize annotation
  colnames(anno.tbl) <- c("cell.name", "cell.type")
  celltype.freq <- as.data.frame(table(anno.tbl$cell.type))
  rm.celltype   <- as.character(celltype.freq[celltype.freq$Freq <= 10, 1])
  
  # Rename data columns by cell type
  new.colnames <- sapply(colnames(data), function(x) {
    make.names(
      anno.tbl$cell.type[anno.tbl$cell.name == x],
      unique = FALSE
    )
  })
  colnames(data) <- new.colnames

  # Load background data and filter for species
  if (toupper(species) == "HUMAN") {
    data("HUMAN_Background_data", package = "funres", envir = environment())
    hsa_tab <- read.table(
      system.file("extdata", "uniprot-filtered-organism-HSA.tab", package = "funres"),
      sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    hsa_tab <- hsa_tab[hsa_tab$V3 != "", ]
    secreted <- unique(unlist(strsplit(hsa_tab$V2, " ")))
    LR       <- LR[LR$Ligand %in% secreted, ]
    Ligands  <- intersect(Ligands, secreted)
  } else if (toupper(species) == "MOUSE") {
    load(system.file("data", "MOUSE_Background_data.rda", package = "funres"), envir = environment())
    mmu_tab <- read.table(
      system.file("extdata", "uniprot-filtered-organism-MMU.tab", package = "funres"),
      sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    mmu_tab <- mmu_tab[mmu_tab$V4 != "", ]
    secreted <- unique(unlist(strsplit(mmu_tab$V3, " ")))
    LR       <- LR[LR$Ligand %in% secreted, ]
    Ligands  <- intersect(Ligands, secreted)
  } else {
    stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
  }

  # Filter the background signaling interactome
  rownames(Background_signaling_interactome) <- seq_len(nrow(Background_signaling_interactome))
  non.rec.id <- which(
    Background_signaling_interactome$Source == dummy.var &
      !(Background_signaling_interactome$Target %in% Receptors)
  )
  keep.id <- setdiff(rownames(Background_signaling_interactome), non.rec.id)
  Background_signaling_interactome <- Background_signaling_interactome[keep.id, ]

  TF_TF_interactions <- remove.factors(TF_TF_interactions)

  # Subset populations
  all.pops <- setdiff(unique(anno.tbl$cell.type), rm.celltype)
  anno.tbl <- anno.tbl[anno.tbl$cell.type %in% all.pops, ]

  # Compute ligand expression per cell type
    # Compute ligand expression per cell type
  data.lig.exp <- get.gene.expr(data, intersect(Ligands, rownames(data)), all.pops)
  # Ensure matrix format even if single population
  if (is.null(dim(data.lig.exp)) || length(dim(data.lig.exp)) < 2) {
    tmp <- data.lig.exp
    data.lig.exp <- matrix(tmp, nrow = length(tmp), dimnames = list(names(tmp), NULL))
  }
  colnames(data.lig.exp) <- paste0("Ligand.", colnames(data.lig.exp))
  L.frame <- dplyr::inner_join(
    data.lig.exp, LR[-1], by = c("Ligand.gene" = "Ligand")
  )::inner_join(
    data.lig.exp, LR[-1], by = c("Ligand.gene" = "Ligand")
  )
  L.frame <- L.frame[L.frame$Ligand.exp.perc > consv.thrs, ]

  # Save temporary data
  save(data, anno.tbl, data.lig.exp, L.frame,
       file = file.path(out.path, temp.folder.name, "temp_data.RData"))
  rm(data.lig.exp)

  # Per-population analysis
  invisible(lapply(all.pops, function(celltype) {
    cell.exp.tbl <- readRDS(file.path(out.path, temp.folder.name,
                                       paste0("temp_", celltype, ".Rds")))
    cat(" Celltype : ", celltype, "\n")
    max.cluster.info <- get.cons.tfs(cell.exp.tbl)
    if (is.list(max.cluster.info)) {
      sig.input <- cell.exp.tbl[, max.cluster.info$tf.max.mat.cell]
      sig.input <- cbind.data.frame(row.names(sig.input), sig.input,
                                     stringsAsFactors = FALSE)
      cons.tfs <- data.frame(Gene = unique(max.cluster.info$tf.count$Gene),
                             bool = 1, stringsAsFactors = FALSE)

      cat("   Starting SigHotSpotter analysis\n")
      hotspot.out <- SigHotSpotter_pipeline(
        idata      = sig.input,
        species    = species,
        cutoff     = sighot.cutoff,
        DE_Genes   = cons.tfs,
        percentile = sighot.percentile,
        ncores     = ncores
      )
      save(max.cluster.info, sig.input, cons.tfs, hotspot.out,
           file = file.path(out.path, temp.folder.name,
                            paste0("temp_", celltype, "_results.RData")))
      rm(cell.exp.tbl, max.cluster.info, sig.input, cons.tfs, hotspot.out)
    } else {
      save(max.cluster.info,
           file = file.path(out.path, temp.folder.name,
                            paste0("temp_", celltype, "_results.RData")))
      rm(cell.exp.tbl, max.cluster.info)
    }
  }))

  cat(" Collating Results\n")
  result.files <- list.files(file.path(out.path, temp.folder.name),
                             pattern = "_results.RData")
  pops <- gsub("(temp_|_results.RData)", "", result.files)
  collate <- setNames(lapply(pops, function(celltype) {
    load(file.path(out.path, temp.folder.name,
                   paste0("temp_", celltype, "_results.RData")))
    if (is.list(max.cluster.info)) {
      out <- list(
        maxmat  = max.cluster.info,
        hotspot = hotspot.out
      )
      if (exists("R.TF.info")) out$R.TF.info <- R.TF.info
      if (exists("LR.frame"))  out$LR.frame  <- LR.frame
      return(out)
    } else {
      return(NULL)
    }
  }), pops)
  collate <- collate[!vapply(collate, is.null, FALSE)]

  # Combine tissue-wide results
  tissue.LR <- do.call(rbind, lapply(collate, `[[`, "LR.frame"))
  tissue.R.TF <- do.call(rbind, lapply(names(collate), function(ct) {
    df <- collate[[ct]][["R.TF.info"]]
    cbind(Celltype = ct, df)
  }))

  saveRDS(tissue.LR, file = file.path(out.path, "tissue_LR_no_bootstrap.Rds"))
  write.table(tissue.LR, file = file.path(out.path, "tissue_LR_no_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat("Performing bootstrapping ... \n")
  load(file.path(out.path, temp.folder.name, "temp_data.RData"))
  final <- Bootstrap.NewScoring(data, LR = tissue.LR, R.TF = tissue.R.TF,
                                significance.cutoff = 0.7)
  saveRDS(final, file = file.path(out.path, "tissue_LR_with_bootstrap.Rds"))
  write.table(final, file = file.path(out.path, "tissue_LR_with_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  collate$final <- final
  save(collate, file = file.path(out.path, paste0("output_", tissue.name, ".RData")))

  partitions <- funres_partitions(output = collate,
                                  out.path = out.path,
                                  tissue.name = tissue.name)
  if (gen.markers) gen.markers(partitions, out.path = out.path,
                                 tissue.name = tissue.name)

  invisible(collate)
}

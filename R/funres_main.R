#' @title FunRes: Functional Resolution Pipeline
#' @description
#' FunRes is the main function of the package. It uses tissue scRNA-seq data
#' and cell-type annotations to predict ligand–receptor interactions that
#' maintain cellular phenotypes across cell types in a tissue.
#'
#' @param data A numeric matrix (genes × cells), rownames = EnsemblID_GeneName.
#' @param anno.tbl A data.frame with two columns: cell IDs and cell types.
#' @param species Character; either "HUMAN" or "MOUSE" (case-insensitive).
#' @param sighot.cutoff Numeric; cutoff for SigHotSpotter (default 0.1).
#' @param sighot.percentile Numeric; percentile for SigHotSpotter (default 70).
#' @param consv.thrs Numeric; conservation threshold (default 0.1).
#' @param n Integer; number of bootstrapping iterations (default 1000).
#' @param ncores Integer; number of cores to use (default 4).
#' @param z.score.cutoff Numeric; z-score cutoff after bootstrap (default 2).
#' @param tissue.name Character; label for your tissue (used in output filenames).
#' @param temp.folder.name Character; name of temporary subfolder (default "tmp").
#' @param out.path Character; path where results will be written.
#' @param gen.markers Logical; whether to generate marker heatmaps (default TRUE).
#'
#' @return Invisibly returns a list with collated results, including
#'   tissue-level and per-cell-type outputs.
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
  temp.folder.name  = "tmp",
  out.path,
  gen.markers       = TRUE
) {
  # === 1. Create output directories ===
  dir.create(out.path, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out.path, temp.folder.name), showWarnings = FALSE)

  # === 2. Save input parameters ===
  params <- data.frame(
    parameter = c("tissue.name", "out.path", "temp.folder.name", "species",
                  "sighot.cutoff", "sighot.percentile", "consv.thrs",
                  "n", "ncores", "z.score.cutoff"),
    value = c(tissue.name, out.path, temp.folder.name, species,
              sighot.cutoff, sighot.percentile, consv.thrs,
              n, ncores, z.score.cutoff),
    stringsAsFactors = FALSE
  )
  write.table(
    params,
    file = file.path(out.path, paste0("input_parameters_", Sys.Date(), ".txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # === 3. Standardize annotation ===
  colnames(anno.tbl)[1:2] <- c("cell.name", "cell.type")
  freq <- table(anno.tbl$cell.type)
  rm.types <- names(freq)[freq <= 10]

  # === 4. Rename data columns by cell type ===
  colnames(data) <- sapply(colnames(data), function(cid) {
    ct <- anno.tbl$cell.type[anno.tbl$cell.name == cid]
    make.names(ct)
  })

  # === 5. Load background data ===
  species_up <- toupper(species)
  if (species_up == "HUMAN") {
    data("HUMAN_Background_data", package = "funres", envir = environment())
    tab <- read.table(system.file("extdata", "uniprot-filtered-organism-HSA.tab", package = "funres"),
                      sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    lig_col <- 2
  } else if (species_up == "MOUSE") {
    data("MOUSE_Background_data", package = "funres", envir = environment())
    tab <- read.table(system.file("extdata", "uniprot-filtered-organism-MMU.tab", package = "funres"),
                      sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    lig_col <- 3
  } else {
    stop("species must be either 'HUMAN' or 'MOUSE'")
  }
  secreted <- unique(unlist(strsplit(tab[[lig_col]], " ")))
  LR       <- LR[LR$Ligand %in% secreted, ]
  Ligands  <- intersect(Ligands, secreted)

  # === 6. Filter interactome and TF interactions ===
  Background_signaling_interactome <-
    Background_signaling_interactome[!(Background_signaling_interactome$Source == dummy.var &
      !(Background_signaling_interactome$Target %in% Receptors)), ]
  TF_TF_interactions <- remove.factors(TF_TF_interactions)

  # === 7. Subset populations ===
  pops <- setdiff(unique(anno.tbl$cell.type), rm.types)
  anno.tbl <- anno.tbl[anno.tbl$cell.type %in% pops, ]

  # === 8. Compute ligand expression per population ===
  data_lig_exp <- get.gene.expr(
    exp.tbl   = data,
    anno.tbl  = anno.tbl,
    genes     = intersect(Ligands, rownames(data)),
    cell.type = pops
  )
  if (is.null(dim(data_lig_exp))) {
    data_lig_exp <- matrix(
      data_lig_exp,
      ncol = 1,
      dimnames = list(names(data_lig_exp), NULL)
    )
  }
  colnames(data_lig_exp) <- paste0("Ligand.", colnames(data_lig_exp))

  L_frame <- dplyr::inner_join(
    data.frame(Ligand.gene = rownames(data_lig_exp), data_lig_exp, check.names = FALSE),
    LR[, -1],
    by = c("Ligand.gene" = "Ligand")
  )

  # === 9. Save temp data ===
  save(data, anno.tbl, data_lig_exp, L_frame,
       file = file.path(out.path, temp.folder.name, "temp_data.RData"))
  rm(data_lig_exp)

  # === 10. Per-population analysis ===
  results_list <- vector("list", length(pops))
  names(results_list) <- pops
  for (ct in pops) {
    data_ct <- data[, colnames(data) == ct, drop = FALSE]
    saveRDS(data_ct, file = file.path(out.path, temp.folder.name, paste0("temp_", ct, ".Rds")))

    max_info <- get.cons.tfs(data_ct)
    if (is.list(max_info)) {
      sig_in <- data_ct[, max_info$tf.max.mat.cell, drop = FALSE]
      sig_in <- data.frame(Gene = rownames(sig_in), sig_in, check.names = FALSE)
      tf_df <- data.frame(Gene = unique(max_info$tf.count$Gene), stringsAsFactors = FALSE)

      hs_out <- SigHotSpotter_pipeline(
        idata      = sig_in,
        species    = species_up,
        cutoff     = sighot.cutoff,
        DE_Genes   = tf_df,
        percentile = sighot.percentile,
        ncores     = ncores
      )

      results_list[[ct]] <- list(
        max_info = max_info,
        hotspot  = hs_out
      )
    }
  }

  # === 11. Collate results ===
  LR_all   <- do.call(rbind, lapply(results_list, function(x) x$hotspot$path.sums))
  RTF_info <- do.call(rbind, lapply(names(results_list), function(ct) {
    df <- results_list[[ct]]$hotspot$iTF.targets
    cbind(celltype = ct, df)
  }))

  # Save collated outputs
  saveRDS(LR_all,   file = file.path(out.path, "tissue_LR_no_bootstrap.Rds"))
  write.table(LR_all, file = file.path(out.path, "tissue_LR_no_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # === 12. Bootstrapping ===
  load(file.path(out.path, temp.folder.name, "temp_data.RData"))
  final <- Bootstrap.NewScoring(data, LR = LR_all, R.TF = RTF_info,
                                significance.cutoff = 0.7)
  saveRDS(final, file = file.path(out.path, "tissue_LR_with_bootstrap.Rds"))
  write.table(final, file = file.path(out.path, "tissue_LR_with_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  results_list$final <- final
  saveRDS(results_list, file = file.path(out.path, paste0("output_", tissue.name, ".Rds")))

  # === 13. Generate markers if requested ===
  if (gen.markers) {
    partitions <- funres_partitions(
      output      = results_list,
      out.path    = out.path,
      tissue.name = tissue.name
    )
    gen.markers(partitions, out.path = out.path, tissue.name = tissue.name)
  }

  invisible(results_list)
}

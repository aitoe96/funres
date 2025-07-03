#' @title FunRes: Functional Resolution Pipeline
#' @description
#' Main function of the funres package: performs ligand–receptor analysis on scRNA-seq data
#' to identify interactions maintaining cellular phenotypes.
#'
#' @param data Numeric matrix (genes × cells). Rownames: "EnsemblID_GeneName".
#' @param anno.tbl Data.frame with two columns: cell IDs and cell types.
#' @param species Character; "HUMAN" or "MOUSE" (case-insensitive).
#' @param sighot.cutoff Numeric; SigHotSpotter cutoff (default 0.1).
#' @param sighot.percentile Numeric; SigHotSpotter percentile (default 70).
#' @param consv.thrs Numeric; conservation threshold (default 0.1).
#' @param n Integer; bootstrapping iterations (default 1000).
#' @param ncores Integer; number of cores (default 4).
#' @param z.score.cutoff Numeric; z-score cutoff post-bootstrap (default 2).
#' @param tissue.name Character; tissue label for outputs.
#' @param temp.folder.name Character; temporary folder name (default "tmp").
#' @param out.path Character; output directory path.
#' @param gen.markers Logical; generate marker heatmaps (default TRUE).
#'
#' @return Invisibly returns a list of results (per-cell-type and tissue-wide).
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
  # 1) Create output dirs
  dir.create(out.path, showWarnings = FALSE, recursive = TRUE)
  tmpdir <- file.path(out.path, temp.folder.name)
  dir.create(tmpdir, showWarnings = FALSE)

  # 2) Save input parameters
  params <- data.frame(
    parameter = c(
      "tissue.name", "out.path", "temp.folder.name", "species",
      "sighot.cutoff", "sighot.percentile", "consv.thrs",
      "n", "ncores", "z.score.cutoff"
    ),
    value = as.character(c(
      tissue.name, out.path, temp.folder.name, species,
      sighot.cutoff, sighot.percentile, consv.thrs,
      n, ncores, z.score.cutoff
    )),
    stringsAsFactors = FALSE
  )
  write.table(
    params,
    file = file.path(out.path, paste0("input_parameters_", Sys.Date(), ".txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # 3) Standardize annotation
  colnames(anno.tbl)[1:2] <- c("cell.name", "cell.type")
  freq <- table(anno.tbl$cell.type)
  rm.types <- names(freq)[freq <= 10]
  anno.tbl <- anno.tbl[!anno.tbl$cell.type %in% rm.types, , drop = FALSE]

  # 4) Filter and rename columns of data by cell type
  common_ids <- intersect(colnames(data), anno.tbl$cell.name)
  if (length(common_ids) == 0) {
    stop("No matching cell IDs between expression data and annotation table.")
  }
  data <- data[, common_ids, drop = FALSE]
  colnames(data) <- sapply(colnames(data), function(cid) {
    idx <- which(anno.tbl$cell.name == cid)
    make.names(anno.tbl$cell.type[idx])
  })

  # 5) Load background data and filter ligands
  species_up <- toupper(species)
  if (species_up == "HUMAN") {
    data("HUMAN_Background_data", package = "funres", envir = environment())
    tab <- data.table::fread(
      file = system.file(
        "extdata", "uniprot-filtered-organism-HSA.tab", package = "funres"
      ), sep = "\t", header = FALSE, fill = TRUE
    )
    lig_col <- 2
  } else if (species_up == "MOUSE") {
    data("MOUSE_Background_data", package = "funres", envir = environment())
    tab <- data.table::fread(
      file = system.file(
        "extdata", "uniprot-filtered-organism-MMU.tab", package = "funres"
      ), sep = "\t", header = FALSE, fill = TRUE
    )
    lig_col <- 3
  } else {
    stop("species must be 'HUMAN' or 'MOUSE'")
  }
  secreted <- unique(unlist(strsplit(tab[[lig_col]], " ")))
  LR <- LR[LR$Ligand %in% secreted, , drop = FALSE]
  Ligands <- intersect(Ligands, secreted)

  # 6) Filter interactome
  Background_signaling_interactome <-
    Background_signaling_interactome[!(
      Background_signaling_interactome$Source == dummy.var &
      !(Background_signaling_interactome$Target %in% Receptors)
    ), ]
  TF_TF_interactions <- remove.factors(TF_TF_interactions)

  # 7) Populations to analyze
  pops <- unique(colnames(data))

  # 8) Compute ligand expression
  data_lig_exp <- get.gene.expr(
    exp.tbl   = data,
    anno.tbl  = anno.tbl,
    genes     = intersect(Ligands, rownames(data)),
    cell.type = pops
  )
  if (is.null(dim(data_lig_exp))) {
    data_lig_exp <- matrix(
      data_lig_exp, ncol = 1,
      dimnames = list(names(data_lig_exp), NULL)
    )
  }
  colnames(data_lig_exp) <- paste0("Ligand.", colnames(data_lig_exp))

  # Prepare for join and remove duplicates
  x_df <- data.frame(
    Ligand.gene = rownames(data_lig_exp),
    data_lig_exp,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  dups <- duplicated(x_df$Ligand.gene)
  if (any(dups)) {
    warning(
      "Detected duplicated genes in data_lig_exp: ",
      paste(unique(x_df$Ligand.gene[dups]), collapse = ", "),
      ". Keeping first occurrence."
    )
    x_df <- x_df[!dups, , drop = FALSE]
  }

  # 9) Join ligands to expression
  L_frame <- dplyr::inner_join(
    x_df,
    LR[, -1, drop = FALSE],
    by = c("Ligand.gene" = "Ligand")
  )

  # 10) Save temp data
  save(data, anno.tbl, data_lig_exp, L_frame,
       file = file.path(tmpdir, "temp_data.RData"))
  rm(data_lig_exp)

  # 11) Per-population processing
  res_list <- setNames(vector("list", length(pops)), pops)
  for (ct in pops) {
    data_ct <- data[, colnames(data) == ct, drop = FALSE]
    saveRDS(data_ct, file = file.path(tmpdir, paste0("temp_", ct, ".Rds")))
    max_info <- get.cons.tfs(data_ct)
    if (is.list(max_info)) {
      sig_in <- data.frame(
        Gene = rownames(data_ct)[max_info$tf.max.mat.cell],
        data_ct[, max_info$tf.max.mat.cell, drop = FALSE],
        check.names = FALSE
      )
      hs_out <- SigHotSpotter_pipeline(
        idata = sig_in,
        species = species_up,
        cutoff = sighot.cutoff,
        DE_Genes = data.frame(Gene = unique(max_info$tf.count$Gene)),
        percentile = sighot.percentile,
        ncores = ncores
      )
      res_list[[ct]] <- list(max_info = max_info, hotspot = hs_out)
    }
  }

  # 12) Collate results
  LR_all <- do.call(rbind, lapply(res_list, function(x) x$hotspot$path.sums))
  RTF_df <- do.call(rbind, lapply(names(res_list), function(ct) {
    df <- res_list[[ct]]$hotspot$iTF.targets
    cbind(celltype = ct, df)
  }))
  saveRDS(LR_all, file = file.path(out.path, "tissue_LR_no_bootstrap.Rds"))
  write.table(LR_all, file = file.path(out.path, "tissue_LR_no_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # 13) Bootstrapping
  load(file.path(tmpdir, "temp_data.RData"))
  final <- Bootstrap.NewScoring(
    data, LR = LR_all, R.TF = RTF_df, significance.cutoff = 0.7
  )
  saveRDS(final, file = file.path(out.path, "tissue_LR_with_bootstrap.Rds"))
  write.table(final, file = file.path(out.path, "tissue_LR_with_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # 14) Final outputs
  res_list$final <- final
  saveRDS(res_list, file = file.path(out.path, paste0("output_", tissue.name, ".Rds")))

  # 15) Generate markers
  if (gen.markers) {
    parts <- funres_partitions(
      output = res_list,
      out.path = out.path,
      tissue.name = tissue.name
    )
    gen.markers(parts, out.path = out.path, tissue.name = tissue.name)
  }

  invisible(res_list)
}

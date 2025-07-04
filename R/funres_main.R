#' @title FunRes: Functional Resolution Pipeline
#' @description
#' Main function of the funres package: performs ligand–receptor analysis on scRNA-seq data
#' to identify interactions maintaining cellular phenotypes.
#'
#' @param data Numeric matrix (genes × cells). Rownames: "EnsemblID_GeneName".
#' @param anno.tbl Data.frame with two columns: cell.name (ID) and cell.type.
#' @param species Character; "HUMAN" or "MOUSE" (case-insensitive).
#' @param sighot.cutoff Numeric; SigHotSpotter cutoff (default 0.1).
#' @param sighot.percentile Numeric; SigHotSpotter percentile (default 70).
#' @param consv.thrs Numeric; conservation threshold (default 0.1).
#' @param n Integer; number of bootstrap iterations (default 1000).
#' @param ncores Integer; number of cores to use (default 4).
#' @param z.score.cutoff Numeric; z-score cutoff post-bootstrap (default 2).
#' @param tissue.name Character; tissue label for output files.
#' @param temp.folder.name Character; name of temporary subfolder (default "tmp").
#' @param out.path Character; path where results will be written.
#' @param gen.markers Logical; whether to generate marker heatmaps (default TRUE).
#'
#' @return Invisibly returns a list with per-cell-type and tissue-wide results.
#' @importFrom dplyr select inner_join
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
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
  # 1) Create output directories
  dir.create(out.path, recursive = TRUE, showWarnings = FALSE)
  tmpdir <- file.path(out.path, temp.folder.name)
  dir.create(tmpdir, showWarnings = FALSE)

  # 2) Save input parameters
  params <- data.frame(
    parameter = c(
      "tissue.name","out.path","temp.folder.name","species",
      "sighot.cutoff","sighot.percentile","consv.thrs",
      "n","ncores","z.score.cutoff"
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

  # 3) Standardize annotation and filter rare cell types
  colnames(anno.tbl)[1:2] <- c("cell.name","cell.type")
  freq <- table(anno.tbl$cell.type)
  keep.types <- names(freq)[freq > 10]
  anno.tbl <- anno.tbl[anno.tbl$cell.type %in% keep.types, , drop = FALSE]

  # 4) Align and rename data columns
  valid.cells <- intersect(colnames(data), anno.tbl$cell.name)
  if (length(valid.cells) == 0) {
    stop("No matching cell IDs between data and annotation table.")
  }
  data <- data[, valid.cells, drop = FALSE]
  colnames(data) <- vapply(valid.cells, function(cid) {
    ct <- anno.tbl$cell.type[anno.tbl$cell.name == cid]
    make.names(ct)
  }, FUN.VALUE = "")

  # 5) Load background and filter ligands
  sp <- toupper(species)
  if (sp == "HUMAN") {
    data("HUMAN_Background_data", package = "funres", envir = environment())
    tab <- data.table::fread(
      system.file("extdata","uniprot-filtered-organism-HSA.tab", package = "funres"),
      sep = "\t", header = FALSE, fill = TRUE
    )
    lig_col <- 2
  } else if (sp == "MOUSE") {
    data("MOUSE_Background_data", package = "funres", envir = environment())
    tab <- data.table::fread(
      system.file("extdata","uniprot-filtered-organism-MMU.tab", package = "funres"),
      sep = "\t", header = FALSE, fill = TRUE
    )
    lig_col <- 3
  } else {
    stop("species must be 'HUMAN' or 'MOUSE'.")
  }
  secreted <- unique(unlist(strsplit(tab[[lig_col]], " ")))
  LR <- LR[LR$Ligand %in% secreted, , drop = FALSE]
  Ligands <- intersect(Ligands, secreted)

  # 6) Filter interactome and TF-TF interactions
  BSI <- Background_signaling_interactome
  BSI <- BSI[!(BSI$Source == dummy.var & !(BSI$Target %in% Receptors)), ]
  TF_TF <- TF_TF_interactions

  # 7) Define populations
  pops <- unique(colnames(data))

  # 8) Compute ligand expression and pivot to matrix
  lig_df <- get.gene.expr(
    exp.tbl   = data,
    anno.tbl  = anno.tbl,
    genes     = intersect(Ligands, rownames(data)),
    cell.type = pops
  )
  lig_sub <- dplyr::select(lig_df, gene, celltype, exp.perc)
  wide <- tidyr::pivot_wider(
    lig_sub,
    names_from  = celltype,
    values_from = exp.perc,
    values_fill = list(exp.perc = 0)
  )
  data_lig_exp <- tibble::column_to_rownames(wide, var = "gene")
  colnames(data_lig_exp) <- paste0("Ligand.", colnames(data_lig_exp))

  # 9) Prepare x_df and ensure unique genes
  x_df <- data.frame(
    Ligand.gene = rownames(data_lig_exp),
    data_lig_exp,
    check.names     = FALSE,
    stringsAsFactors = FALSE
  )

  # 10) Join with LR
  L_frame <- dplyr::inner_join(
    x_df,
    LR[, -1, drop = FALSE],
    by = c("Ligand.gene" = "Ligand")
  )

  # 11) Save temp data
  save(data, anno.tbl, data_lig_exp, L_frame,
       file = file.path(tmpdir, "temp_data.RData"))
  rm(data_lig_exp)

  # 12) Per-population analysis
  res_list <- setNames(vector("list", length(pops)), pops)
  for (ct in pops) {
    sub <- data[, colnames(data) == ct, drop = FALSE]
    saveRDS(sub, file = file.path(tmpdir, paste0("temp_", ct, ".Rds")))
    mi <- get.cons.tfs(sub)
    if (is.list(mi)) {
      # Prepare sig input: select rows corresponding to tf.max.mat.cell
      sig <- data.frame(
        Gene = rownames(sub)[mi$tf.max.mat.cell],
        sub[mi$tf.max.mat.cell, , drop = FALSE],
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      hs <- SigHotSpotter_pipeline(
        idata      = sig,
        species    = sp,
        cutoff     = sighot.cutoff,
        DE_Genes   = data.frame(Gene = unique(mi$tf.count$Gene)),
        percentile = sighot.percentile,
        ncores     = ncores
      )
      res_list[[ct]] <- list(max_info = mi, hotspot = hs)
    }
  }

# 13) Collate results
  LR_all <- do.call(rbind, lapply(res_list, function(x) x$hotspot$path.sums))
  RTF_df <- do.call(rbind, lapply(names(res_list), function(ct) {
    df <- res_list[[ct]]$hotspot$iTF.targets
    cbind(celltype = ct, df)
  }))
  saveRDS(LR_all, file = file.path(out.path, "tissue_LR_no_bootstrap.Rds"))
  write.table(LR_all,
              file = file.path(out.path, "tissue_LR_no_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # 14) Bootstrap
  load(file.path(tmpdir, "temp_data.RData"))
  final <- Bootstrap.NewScoring(
    data, LR = LR_all, R.TF = RTF_df,
    significance.cutoff = 0.7
  )
  saveRDS(final, file = file.path(out.path, "tissue_LR_with_bootstrap.Rds"))
  write.table(final,
              file = file.path(out.path, "tissue_LR_with_bootstrap.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # 15) Final results and markers
  res_list$final <- final
  saveRDS(res_list,
          file = file.path(out.path, paste0("output_", tissue.name, ".Rds")))
  if (gen.markers) {
    parts <- funres_partitions(
      output     = res_list,
      out.path   = out.path,
      tissue.name= tissue.name
    )
    gen.markers(parts, out.path = out.path, tissue.name = tissue.name)
  }

  invisible(res_list)
}

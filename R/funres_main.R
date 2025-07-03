#' @title FunRes: Functional Resolution Pipeline
#' @description
#' FunRes is the main function of the package. It uses tissue scRNA-seq data
#' and cell-type annotations to predict ligand–receptor interactions that
#' maintain cellular phenotypes across cell types in a tissue.
#'
#' @param data A numeric matrix (genes × cells), rownames = EnsemblID_GeneName.
#' @param anno.tbl A data.frame with two columns: first = cell IDs, second = cell types.
#' @param species Character; either `"HUMAN"` or `"MOUSE"`.
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
  
  # create output directories
  dir.create(out.path,           showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out.path, temp.folder.name), showWarnings = FALSE)
  
  message("Preparing background data for species: ", species)
  
  if (toupper(species) == "HUMAN") {
    # load data/HUMAN_Background_data.rda
    load(system.file("data", "HUMAN_Background_data.rda", package = "funres"),
         envir = environment())
    # read inst/extdata/uniprot-filtered-organism-HSA.tab
    hsa_tab <- read.table(
      system.file("extdata", "uniprot-filtered-organism-HSA.tab", package = "funres"),
      sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    secreted <- hsa_tab[grepl("Secreted", hsa_tab$V3), "V2"]
    secreted <- unique(unlist(strsplit(secreted, " ")))
    
    LR       <- LR[LR$Ligand %in% secreted, ]
    Ligands  <- intersect(Ligands, secreted)
    
  } else if (toupper(species) == "MOUSE") {
    load(system.file("data", "MOUSE_Background_data.rda", package = "funres"),
         envir = environment())
    mmu_tab <- read.table(
      system.file("extdata", "uniprot-filtered-organism-MMU.tab", package = "funres"),
      sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    secreted <- mmu_tab[grepl("Secreted", mmu_tab$V4), "V3"]
    secreted <- unique(unlist(strsplit(secreted, " ")))
    
    LR      <- LR[LR$Ligand %in% secreted, ]
    Ligands <- intersect(Ligands, secreted)
    
  } else {
    stop("`species` must be either 'HUMAN' or 'MOUSE'.")
  }
  
  # standardize annotation table
  colnames(anno.tbl) <- c("cell.name", "cell.type")
  
  # ... (rest of your pipeline follows exactly as before) ...
  # e.g. compute data.lig.exp, L.frame, save temp data, per-celltype loops, 
  # SigHotSpotter_pipeline, collate results, bootstrap, final save.
  
  # at end, return invisibly
  invisible(output)
}

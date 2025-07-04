#' @title FunRes: Functional Resolution Pipeline
#' @description
#' Main function of the funres package: performs ligand–receptor analysis on scRNA-seq data
#' to identify interactions maintaining cellular phenotypes.
#'
#' @param data Numeric matrix (genes × cells). Rownames: "EnsemblID_GeneName".
#' @param anno.tbl Data.frame con dos columnas: cell.name (ID) y cell.type.
#' @param species Character; "HUMAN" o "MOUSE".
#' @param sighot.cutoff Numeric; SigHotSpotter cutoff (default 0.1).
#' @param sighot.percentile Numeric; percentile para SigHotSpotter (default 70).
#' @param consv.thrs Numeric; conservation threshold (default 0.1).
#' @param n Integer; número de iteraciones de bootstrap (default 1000).
#' @param ncores Integer; núm. de cores (default 4).
#' @param z.score.cutoff Numeric; z-score cutoff tras bootstrap (default 2).
#' @param tissue.name Character; etiqueta del tejido para ficheros de salida.
#' @param temp.folder.name Character; nombre de carpeta temporal (default "tmp").
#' @param out.path Character; ruta de salida para resultados.
#' @param gen.markers Logical; generar mapas de marcadores (default TRUE).
#'
#' @return Invisiblemente, una lista con los resultados por célula y globales.
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
  # 1) Crear dirs de salida
  dir.create(out.path, recursive=TRUE, showWarnings=FALSE)
  tmpdir <- file.path(out.path, temp.folder.name)
  dir.create(tmpdir, showWarnings=FALSE)

  # 2) Guardar parámetros
  params <- data.frame(
    parameter = c("tissue.name","out.path","temp.folder.name","species",
                  "sighot.cutoff","sighot.percentile","consv.thrs",
                  "n","ncores","z.score.cutoff"),
    value     = as.character(c(tissue.name, out.path, temp.folder.name,
                  species, sighot.cutoff, sighot.percentile,
                  consv.thrs, n, ncores, z.score.cutoff)),
    stringsAsFactors = FALSE
  )
  write.table(params,
              file = file.path(out.path, paste0("input_parameters_", Sys.Date(), ".txt")),
              sep = "\t", row.names=FALSE, quote=FALSE)

  # 3) Estandarizar anotación y filtrar cell types con pocas células
  colnames(anno.tbl)[1:2] <- c("cell.name","cell.type")
  frec <- table(anno.tbl$cell.type)
  keep.types <- names(frec)[frec > 10]
  anno.tbl <- anno.tbl[anno.tbl$cell.type %in% keep.types, , drop=FALSE]

  # 4) Alinear y renombrar columnas de data
  valid.cells <- intersect(colnames(data), anno.tbl$cell.name)
  if (length(valid.cells)==0) stop("No matching cell IDs entre data y anno.tbl")
  data <- data[, valid.cells, drop=FALSE]
  colnames(data) <- vapply(valid.cells, function(cid) {
    ct <- anno.tbl$cell.type[anno.tbl$cell.name==cid]
    make.names(ct)
  }, FUN.VALUE="")

  # 5) Cargar datos de fondo y filtrar ligandos
  sp_up <- toupper(species)
  if (sp_up=="HUMAN") {
    data("HUMAN_Background_data", package="funres", envir=environment())
    tab <- data.table::fread(
      system.file("extdata","uniprot-filtered-organism-HSA.tab",package="funres"),
      sep="\t", header=FALSE, fill=TRUE
    )
    col.lig <- 2
  } else if (sp_up=="MOUSE") {
    data("MOUSE_Background_data", package="funres", envir=environment())
    tab <- data.table::fread(
      system.file("extdata","uniprot-filtered-organism-MMU.tab",package="funres"),
      sep="\t", header=FALSE, fill=TRUE
    )
    col.lig <- 3
  } else {
    stop("species debe ser 'HUMAN' o 'MOUSE'")
  }
  secreted <- unique(unlist(strsplit(tab[[col.lig]], " ")))
  LR       <- LR[LR$Ligand %in% secreted, , drop=FALSE]
  Ligands  <- intersect(Ligands, secreted)

  # 6) Filtrar interactoma y TF-TF
  BSI <- Background_signaling_interactome
  BSI <- BSI[!(BSI$Source==dummy.var & !(BSI$Target %in% Receptors)), ]
  TFint <- TF_TF_interactions

  # 7) Definir poblaciones
  pops <- unique(colnames(data))

  # 8) Calcular expresión de ligandos
  data_lig_exp <- get.gene.expr(
    exp.tbl   = data,
    anno.tbl  = anno.tbl,
    genes     = intersect(Ligands, rownames(data)),
    cell.type = pops
  )
  if (is.null(dim(data_lig_exp))) {
    data_lig_exp <- matrix(data_lig_exp, ncol=1,
      dimnames=list(names(data_lig_exp),NULL))
  }
  colnames(data_lig_exp) <- paste0("Ligand.", colnames(data_lig_exp))

  # 9) Preparar x_df y deduplicar
  x_df <- data.frame(Ligand.gene=rownames(data_lig_exp),
                     data_lig_exp, check.names=FALSE, stringsAsFactors=FALSE)
  dups <- duplicated(x_df$Ligand.gene)
  if (any(dups)) {
    warning("Genes duplicados en data_lig_exp: ",
            paste(unique(x_df$Ligand.gene[dups]), collapse=", "),
            ". Se mantiene primera ocurrencia.")
    x_df <- x_df[!dups, , drop=FALSE]
  }

  # 10) Unir con LR
  L_frame <- dplyr::inner_join(
    x_df,
    LR[ , -1, drop=FALSE],
    by = c("Ligand.gene"="Ligand")
  )

  # 11) Guardar datos temporales
  save(data, anno.tbl, data_lig_exp, L_frame,
       file=file.path(tmpdir,"temp_data.RData"))
  rm(data_lig_exp)

  # 12) Procesar por población
  res_list <- setNames(vector("list",length(pops)), pops)
  for (ct in pops) {
    sub <- data[, colnames(data)==ct, drop=FALSE]
    saveRDS(sub, file=file.path(tmpdir,paste0("temp_",ct,".Rds")))
    mi <- get.cons.tfs(sub)
    if (is.list(mi)) {
      sig <- data.frame(Gene=rownames(sub)[mi$tf.max.mat.cell],
                        sub[,mi$tf.max.mat.cell,drop=FALSE],
                        check.names=FALSE)
      hs <- SigHotSpotter_pipeline(
        idata   = sig,
        species = sp_up,
        cutoff  = sighot.cutoff,
        DE_Genes   = data.frame(Gene=unique(mi$tf.count$Gene)),
        percentile = sighot.percentile,
        ncores     = ncores
      )
      res_list[[ct]] <- list(max_info=mi, hotspot=hs)
    }
  }

  # 13) Colapsar resultados
  LR_all <- do.call(rbind,lapply(res_list,function(x)x$hotspot$path.sums))
  RTF_df <- do.call(rbind,lapply(names(res_list),function(ct){
    df <- res_list[[ct]]$hotspot$iTF.targets
    cbind(celltype=ct,df)
  }))
  saveRDS(LR_all,file=file.path(out.path,"tissue_LR_no_bootstrap.Rds"))
  write.table(LR_all,file=file.path(out.path,"tissue_LR_no_bootstrap.txt"),
              sep="\t",row.names=FALSE,quote=FALSE)

  # 14) Bootstrap
  load(file.path(tmpdir,"temp_data.RData"))
  final <- Bootstrap.NewScoring(data, LR=LR_all, R.TF=RTF_df,
                                significance.cutoff=0.7)
  saveRDS(final,file=file.path(out.path,"tissue_LR_with_bootstrap.Rds"))
  write.table(final,file=file.path(out.path,"tissue_LR_with_bootstrap.txt"),
              sep="\t",row.names=FALSE,quote=FALSE)

  # 15) Salida final y marcadores
  res_list$final <- final
  saveRDS(res_list,file=file.path(out.path,paste0("output_",tissue.name,".Rds")))
  if (gen.markers) {
    parts <- funres_partitions(output=res_list,
                               out.path=out.path,
                               tissue.name=tissue.name)
    gen.markers(parts, out.path=out.path, tissue.name=tissue.name)
  }
  invisible(res_list)
}

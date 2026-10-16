## DepMap and Sanger CRISPR inputs (spec Section 4.5).

#' depmapSymbols
#' @param cols DepMap gene column names such as "KRAS (3845)".
#' @param hgnc Data frame with symbol, entrez_id.
#' @return HGNC symbols from the Entrez id, falling back to the text before " (".
#' @export
depmapSymbols <- function(cols, hgnc) {
  entrez <- ifelse(grepl("\\([0-9]+\\)$", cols), sub("^.*\\(([0-9]+)\\)$", "\\1", cols), NA_character_)
  sym <- hgnc$symbol[match(entrez, as.character(hgnc$entrez_id), incomparables = NA)]
  ifelse(is.na(sym), trimws(sub("\\s*\\(.*$", "", cols)), sym)
}

#' readDepmapEffect
#' @param path Gene-effect CSV whose first column holds ids: model ids (models in rows, the DepMap
#' layout) or gene names (`genes_in_rows = TRUE`).
#' @param hgnc HGNC table.
#' @param genes_in_rows TRUE when genes are rows and models are columns.
#' @return Genes x models matrix (duplicate symbols: the first is kept).
#' @export
readDepmapEffect <- function(path, hgnc, genes_in_rows = FALSE) {
  x <- data.table::fread(path, sep = ",", data.table = FALSE)
  m <- as.matrix(x[, -1, drop = FALSE])
  rownames(m) <- x[[1]]
  if (genes_in_rows) m <- t(m)
  sym <- depmapSymbols(colnames(m), hgnc)
  keep <- !is.na(sym) & !duplicated(sym)
  out <- t(m[, keep, drop = FALSE])
  rownames(out) <- sym[keep]
  out
}

#' readNonessential
#' @param path DepMap non-essential control list (first column: "SYMBOL (ENTREZ)").
#' @param hgnc HGNC table.
#' @return HGNC symbols.
#' @export
readNonessential <- function(path, hgnc) {
  x <- data.table::fread(path, sep = ",", data.table = FALSE)
  depmapSymbols(as.character(x[[1]]), hgnc)
}

#' hcmiTokenFromRow
#' @param df Data frame of model metadata.
#' @return The first HCMI case id found in each row's text columns, or NA.
#' @export
hcmiTokenFromRow <- function(df) {
  chr <- vapply(df, is.character, logical(1))
  tok <- do.call(cbind, lapply(df[chr], hcmiToken))
  apply(tok, 1, function(r) {
    r <- r[!is.na(r)]
    if (length(r)) r[1] else NA_character_
  })
}

#' combine3dSources
#'
#' Combines DepMap NextGen and Sanger 3D gene effects when they agree: the correlation of
#' per-gene means over shared genes is at least `min_cor`. In that case each model is z-scored
#' within its source. Otherwise the source with more models is used alone.
#' @param E_dm,E_sg Genes x models matrices (E_sg may be NULL).
#' @param min_cor Agreement threshold.
#' @return A list with effect and decision (a sentence recorded in gene_graph_summary.csv).
#' @export
combine3dSources <- function(E_dm, E_sg = NULL, min_cor = 0.8) {
  if (is.null(E_sg) || !ncol(E_sg)) return(list(effect = E_dm, decision = "DepMap 3D only"))
  genes <- intersect(rownames(E_dm), rownames(E_sg))
  r <- stats::cor(rowMeans(E_dm[genes, , drop = FALSE], na.rm = TRUE),
                  rowMeans(E_sg[genes, , drop = FALSE], na.rm = TRUE), use = "complete.obs")
  if (is.finite(r) && r >= min_cor) {
    z <- function(E) {
      out <- scale(E)
      attr(out, "scaled:center") <- NULL
      attr(out, "scaled:scale") <- NULL
      out
    }
    return(list(effect = cbind(z(E_dm[genes, , drop = FALSE]), z(E_sg[genes, , drop = FALSE])),
                decision = sprintf("combined (gene-mean r = %.2f)", r)))
  }
  if (ncol(E_sg) > ncol(E_dm)) {
    list(effect = E_sg, decision = sprintf("Sanger alone (gene-mean r = %.2f < %.2f)", r, min_cor))
  } else {
    list(effect = E_dm, decision = sprintf("DepMap alone (gene-mean r = %.2f < %.2f)", r, min_cor))
  }
}

#' screenClass
#'
#' Arm membership of each gene-effect column (spec Section 4.5): "2D" for traditional screens
#' (ScreenType 2DS), "3D" for NextGen screens run in 3D culture (3DO dome organoids, 3DN CNS
#' spheroids), "exclude" otherwise (including 2D coat screens of NextGen models and screens
#' failing QC). Columns can be screen ids (matched to `screens$ScreenID`) or model ids (typed by
#' `models$CRISPRScreenType`).
#' @param ids Column ids of a gene-effect matrix.
#' @param screens NextGen screen metadata (ScreenID, ModelID, ScreenType, PassesQC).
#' @param models NextGen model metadata (ModelID, CRISPRScreenType).
#' @return Data frame with id, model_id, screen_type, class.
#' @export
screenClass <- function(ids, screens, models) {
  if (mean(ids %in% screens$ScreenID) > 0.5) {
    s <- screens[match(ids, screens$ScreenID, incomparables = NA), , drop = FALSE]
    type <- ifelse(toupper(as.character(s$PassesQC)) %in% "TRUE", s$ScreenType, NA_character_)
    model <- s$ModelID
  } else {
    type <- models$CRISPRScreenType[match(ids, models$ModelID, incomparables = NA)]
    model <- ids
  }
  cls <- ifelse(type %in% "2DS", "2D", ifelse(type %in% c("3DO", "3DN"), "3D", "exclude"))
  data.frame(id = ids, model_id = model, screen_type = type, class = cls, stringsAsFactors = FALSE)
}

#' collapseScreens
#' @param effect Genes x screens matrix.
#' @param model_id Model of each column.
#' @return Genes x models matrix, averaging the screens of each model.
#' @export
collapseScreens <- function(effect, model_id) {
  groups <- split(seq_len(ncol(effect)), model_id)
  out <- vapply(groups, function(i) rowMeans(effect[, i, drop = FALSE], na.rm = TRUE), numeric(nrow(effect)))
  matrix(out, nrow = nrow(effect), dimnames = list(rownames(effect), names(groups)))
}

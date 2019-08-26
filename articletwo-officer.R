characteristics <- function(dset, tableone.names, tableone.factors, extras = list()) {
    title <- "Characteristics"
    overall <- paste0("Cases, n=", dim(dset)[1])
    tableobject <- tableone::CreateTableOne(data = dset, vars = names(tableone.names), factorVars = tableone.factors)
    tablecsv <- print(tableobject,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 1,
                      pDigits = 3,
                      contDigits=1)
    tableone.fullnames <- c(tableone.names, extras)
    tablecsv %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        dplyr::filter(row_number() > 1) %>%
        format(justify = "left", trim = TRUE) %>%
        rowwise() %>%
        mutate(id = gsub("^ *([A-Za-z_0-9]+).*", "\\1", rowname)) %>%
        mutate(present = id %in%  names(tableone.fullnames)) %>%
        mutate(rowname = ifelse(present == TRUE, tableone.fullnames[[id]], rowname)) %>%
        select(rowname, Overall)
}

tableone <- function(dset) {

    tableone.names <- list("BL_AGE" = "Age, y (SD)",
                           "SEX" = "Female, N (%)",
                           "BMI" = "BMI, kg/m² (SD)",
                           "oSYSTM" = "Systolic BP, mmHg (SD)",
                           "oDIASM" = "Diastolic BP, mmHg (SD)",
                           "oPULSEPRESSURE" = "Pulse pressure, mmHg (SD)",
                           "oMAP" = "Mean arterial pressure, mmHg (SD)",
                           "HYPERTENSION" = "Hypertension, N (%)",
                           "CURR_SMOKE" = "Current smoker, N (%)",
                           "PREVAL_DIAB" = "Diabetes mellitus, N (%)",
                           "Q57X" = "Exercise, N (%)",
                           "ANYDRUG" = "Antihypertensive medication, N (%)",
                           "BL_USE_RX_C03" = "  Diuretics, N (%)",
                           "BL_USE_RX_C07" = "  Beta blockers, N (%)",
                           "BL_USE_RX_C08" = "  Calcium channel blockers, N (%)",
                           "BL_USE_RX_C09" = "  RAS blockers, N (%)")
    extras <- list("1"  =  "  Light",
                           "2"  =  "  Moderate",
                           "3"  =  "  Heavy")
    tableone.factors <- c("SEX", "HYPERTENSION",  "ANYDRUG", "CURR_SMOKE", "PREVAL_DIAB",
                          "BL_USE_RX_C03", "BL_USE_RX_C07", "BL_USE_RX_C08", "BL_USE_RX_C09")
    data <- characteristics(dset, tableone.names, tableone.factors, extras)

    flextable(data = data) %>%
        set_header_labels(rowname = "Characteristics",
                          Overall = paste0("Cases, n=", dim(dset)[1])) %>%
        flextable::width(j = 1, width = 3) %>%
        flextable::width(j = 2, width = 1.2) %>%
        flextable::border(border = fp_border(width=0), part="body") %>%
        flextable::border(border = fp_border(width=0), part="header") %>%
        flextable::border(part="header", border.bottom = fp_border(width=1)) %>%
        flextable::border(i = nrow(data), part="body", border.bottom = fp_border(width=1)) %>%
        flextable::bold(bold = FALSE, part = "header") %>%
        flextable::bold(bold = FALSE, part = "body") %>%
        flextable::fontsize(size = 12, part = "header") %>%
        flextable::fontsize(size = 12, part = "body") %>%
        flextable::align(align = "center", part = "all") %>%
        flextable::align(align = "left", part = "header", j = 1) %>%
        flextable::align(align = "left", part = "body", j = 1)
}

typologyformatter <- function(data, font = 12, typology, left = c(1), hleft = c(1)) {
  flex <- flextable(data = data) %>%
    flextable::theme_booktabs() %>%
    flextable::border(border = fp_border(width=0), part="body") %>%
    flextable::border(border = fp_border(width=0), part="header") %>%
    flextable::border(part="header", border.bottom = fp_border(width=1))

  if (!missing(typology)) {
      flex <- flex %>%
          set_header_df(mapping = typology, key = "col_keys") %>%
          merge_h(part = "header") %>%
          flextable::border(part="header", border.bottom = fp_border(width=1))
      if (missing(hleft)) {
          hleft <- c(2)
      }
  }

  flex %>%
      flextable::border(i = nrow(data), part="body", border.bottom = fp_border(width=1)) %>%
      flextable::bold(bold = FALSE, part = "header") %>%
      flextable::bold(bold = FALSE, part = "body") %>%
      flextable::fontsize(size = font, part = "header") %>%
      flextable::fontsize(size = font, part = "body") %>%
      flextable::align(align = "center", part = "all") %>%
      flextable::align(align = "left", part = "header", j = left, i = hleft) %>%
      flextable::align(align = "left", part = "body", j = left)
}

xtableformatter <- function(table, font = 12) {
  xtable_to_flextable(x = xtable(table),
                      NA.string = "",
                      include.rownames = FALSE) %>%
    flextable::border(border.top = fp_border(width=0), part="header", border.bottom = fp_border()) %>%
    flextable::bold(bold = FALSE, part = "header") %>%
    flextable::bold(bold = FALSE, part = "body") %>%
    flextable::fontsize(size = font, part = "header") %>%
    flextable::fontsize(size = font, part = "body") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(align = "left", part = "header", j = 1) %>%
    flextable::align(align = "left", part = "body", j = 1)
}

riskmodel.getn <- function(df, riskclass, inviduals = FALSE, cols = c("htn")) {
  ret <- df %>%
    group_by(.dots = riskclass) %>%
    summarize(n=n(), htn=sum(HT == 1)) %>%
    select(cols)
  return(rbind(c("", ""), c("", ""), c("", ""), ret))
}


riskmodel.onecolumn <- function(riskmodels) {
  lapply(riskmodels, function(result) {
    cbind(
      rbind(result$unadj_persd$mean_ci,
            result$unadj_persd$p,
            "",
            "1.0 (referent)",
            cbind(result$unadj_class$mean_ci)),
      rbind(result$adj_persd$mean_ci,
            result$adj_persd$p,
            "",
            "1.0 (referent)",
            cbind(result$adj_class$mean_ci)))
  })
}


riskmodel.joincolumns <- function(rset, riskmodels, idforn = "fwdbonf_riskclass", order = list()) {
  beginning.list <- cbind(
    c("Per 1-SD", "p", "Quartile", "   Q1", "   Q2", "   Q3", "   Q4"),
    riskmodel.getn(rset, idforn, cols = "n")
  )
  col.list <- lapply(order, function(ord) {
    cbind(riskmodel.getn(rset, paste0(ord, "_riskclass")),
          get(ord, riskmodels))
  })
  col.names <- lapply(order, function(ord) {
    c(paste0(ord, ".event"), paste0(ord, "unadjusted"), paste0(ord, ".adjusted"))
  })
  list <- do.call(cbind, list(beginning.list, col.list)[lapply(list(beginning.list, col.list), length)>0])
  colnames(list) <- c("term", "n", unlist(col.names))
  list
}

writetable <- function(doc, tbl, number = 0, head, foot) {
	text.bold <- fp_text(font.size = 12, bold = TRUE, font.family = "Arial")
	text.normal <- fp_text(font.size = 12, font.family = "Arial")

  doc <- doc %>%
    body_add_fpar(fpar(ftext(paste0("Table ", number, ". "), prop = text.bold), ftext(head, prop = text.normal)), style = "Table Caption") %>%
    body_add_flextable(tbl, align = "left") %>%
    body_add_par(foot, style = "footnote text") %>%
    body_add_break(pos = "after")
}

writeimage <- function(doc, number = 0, filename, header, footer, width = 7) {
text.bold <- fp_text(font.size = 12, bold = TRUE, font.family = "Arial")
text.normal <- fp_text(font.size = 12, font.family = "Arial")

  img <- readPNG(filename)
  doc <- doc %>%
    body_add_fpar(fpar(ftext(paste0("Figure ", number, ". "), prop = text.bold), ftext(header, prop = text.normal)), style = "Image Caption") %>%
    body_add_img(src = filename, width = width, height = width*dim(img)[1]/dim(img)[2]) %>%
    body_add_par(footer, style = "footnote text") %>%
    body_add_break(pos = "after")
}

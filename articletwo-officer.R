characteristics <- function(dset, type = "fr02") {
  if (type == "fr02") {
    factors <- c("HT", "HRX", "sex", "curr_smk", "curr_diab", "plate")
    names <- list("Age, y (SD)" = "AGE",
                  "Female, N (%)" = "sex",
                  "BMI, kg/m² (SD)" = "BMI",
                  "Systolic blood pressure, mmHg (SD)" = "SBP",
                  "Diastolic blood pressure, mmHg (SD)" = "DBP",
                  "Pulse pressure, mmHg (SD)" = "PP",
                  "Mean arterial pressure, mmHg (SD)" = "MAP",
                  "Hypertension, N (%)" = "HT",
                  "Antihypertensive medication, N (%)" = "HRX",
                  "Current smoker, N (%)" = "curr_smk",
                  "Diabetes mellitus, N (%)" = "curr_diab")
  } else {
    factors <- c("HT8", "HRX8", "sex", "CURRSMK8", "curr_diab8", "plate")
    names <- list("Age, y (SD)" = "AGE8",
                  "Female, N (%)" = "sex",
                  "BMI, kg/m² (SD)" = "BMI8",
                  "Systolic blood pressure, mmHg (SD)" = "SBP8",
                  "Diastolic blood pressure, mmHg (SD)" = "DBP8",
                  "Pulse pressure, mmHg (SD)" = "PP8",
                  "Mean arterial pressure, mmHg (SD)" = "MAP8",
                  "Hypertension, N (%)" = "HT8",
                  "Antihypertensive medication, N (%)" = "HRX8",
                  "Current smoker, N (%)" = "CURRSMK8",
                  "Diabetes mellitus, N (%)" = "curr_diab8")
  }
  tableobject <- tableone::CreateTableOne(data = dset, vars = unlist(names), factorVars = factors)
  tablecsv <- print(tableobject,
                    exact = "stage",
                    quote = FALSE,
                    noSpaces = TRUE,
                    printToggle = FALSE,
                    digits = 1,
                    pDigits = 3,
                    contDigits=1)
  
  title <- "Characteristics"
  overall <- paste0("Cases, n=", dim(dset)[1])
  
  tablecsv <- tablecsv %>%
    as.data.frame %>%
    dplyr::filter(row_number() > 1) %>%
    dplyr::mutate(!!title := names(names), !!overall := Overall) %>%
    dplyr::select(-Overall) %>%
    format(justify = "left", trim = TRUE)
  return(tablecsv)
}

typologyformatter <- function(data, font = 12, typology) {
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
  }
  
  flex %>%
    flextable::border(i = nrow(data), part="body", border.bottom = fp_border(width=1)) %>%
    flextable::bold(bold = FALSE, part = "header") %>%
    flextable::bold(bold = FALSE, part = "body") %>%
    flextable::fontsize(size = font, part = "header") %>%
    flextable::fontsize(size = font, part = "body") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(align = "left", part = "header", j = 1) %>%
    flextable::align(align = "left", part = "body", j = 1)
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
  doc <- doc %>% 
    body_add_fpar(fpar(ftext(paste0("Table ", number, "."), prop = text.bold), ftext(head, prop = text.normal)), style = "Table Caption") %>%
    body_add_flextable(tbl, align = "left") %>%
    body_add_par(foot, style = "footnote text") %>%
    body_add_break(pos = "after")
}

writeimage <- function(doc, number = 0, filename, header, footer, width = 7) {
  img <- readPNG(filename)
  doc <- doc %>% 
    body_add_fpar(fpar(ftext(paste0("Figure ", number, "."), prop = text.bold), ftext(header, prop = text.normal)), style = "Image Caption") %>%
    body_add_img(src = filename, width = width, height = width*dim(img)[1]/dim(img)[2]) %>% 
    body_add_par(footer, style = "footnote text") %>%
    body_add_break(pos = "after") 
}

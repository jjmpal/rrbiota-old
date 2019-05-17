# Helper functions for articletwo

replace.brackets <- function (genus) gsub('(.*) \\(+(.+)\\)', '\\1.\\2', genus)

rename.genus  <- function (genus, markword = "Plasmid", mark = "*") {
    name <- gsub('(.*) \\(.*\\)', '\\1', genus)
    star <- ifelse(grepl(markword, genus), mark, "")
    paste0(gsub("_", " ", name), star)
}

myscale <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

filter.phenotype.data <- function(pseq,
                                  included = c("MEN", "BL_AGE", "SYSTM", "DIASM", "BMI",
                                               "CURR_SMOKE", "DIAB", "BL_USE_RX_C09", "Q57X",
                                               "BL_USE_RX_C03",
                                               "BL_USE_RX_C07", "BL_USE_RX_C08"),
                                  excluded) {
    meta(pseq) %>%
        tibble::rownames_to_column(var = "rowname") %>%
        dplyr::select(rowname, included[!included %in% excluded]) %>%
        stats::na.omit() %>%
        dplyr::mutate(ANYDRUG = factor(ifelse(BL_USE_RX_C03  == 1 | BL_USE_RX_C07 == 1 |
                                       BL_USE_RX_C08 == 1  | BL_USE_RX_C09 == 1, 1, 0)),
                      ANYEXERCICE = factor(ifelse(Q57X == 1, 0, 1)),
                      PULSEPRESSURE = SYSTM - DIASM,
                      HYPERTENSION = factor(ifelse(SYSTM >= 140 | DIASM >= 90 | ANYDRUG == 1, 1, 0)),
                      SEX = factor(ifelse(MEN == "Female", 1, 0)),
                      MAP = 2./3.*DIASM + 1./3.*SYSTM) %>%
        dplyr::mutate(oSYSTM = SYSTM,
                      oDIASM = DIASM,
                      oPULSEPRESSURE = PULSEPRESSURE,
                      oMAP = MAP,
                      SYSTM = myscale(SYSTM),
                      DIASM = myscale(DIASM),
                      PULSEPRESSURE = myscale(PULSEPRESSURE),
                      MAP = myscale(MAP),
                      DIAB = ifelse("DIAB" %in% excluded, NA, as.factor(DIAB)),
                      SEX = as.factor(SEX),
                      Q57X = as.factor(Q57X),
                      BL_USE_RX_C03 = as.factor(BL_USE_RX_C03),
                      BL_USE_RX_C07 = as.factor(BL_USE_RX_C07),
                      BL_USE_RX_C08 = as.factor(BL_USE_RX_C08),
                      BL_USE_RX_C09 = as.factor(BL_USE_RX_C09)) %>%
        dplyr::select(-MEN) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = "rowname")
}

meta.merge.alphadiversity <- function(pseq, index = "shannon") {
    alphadiversity  <- microbiome::global(pseq, index = index)
    base::merge(meta(pseq), alphadiversity, by=0, all=TRUE) %>%
        dplyr::rename(Sample_ID = Row.names) %>%
            dplyr::mutate(BL_AGE_GROUP = group_age(BL_AGE))
}

calculate.beta.matrix <- function(pseq) {
    compositional.profile <- microbiome::transform(pseq, 'compositional')
    otu <- microbiome::abundances(compositional.profile)
    bray.dist.m <- vegan::vegdist(t(otu), method="bray")
    dist.matrix <- as.matrix(bray.dist.m)
    attr(dist.matrix, "method") <- "bray"
    dist.matrix
}
calculateglm <- function(dset,
                         responses = list(model_1 = "MAP", model_2 = "SYSTM", model_3 = "DIASM",
                                          model_4 = "PULSEPRESSURE", model_5 = "HYPERTENSION"),
                         min_n_for_continuous = 10,
                         modelstring = "%s ~ BL_AGE + SEX",
                         filterstr = "diversities_") {
    glmlist <- lapply(responses, function(response) {
        fo.family <- ifelse(length(unique(pull(dset, response))) > min_n_for_continuous, stats::gaussian, stats::binomial)
        fo <- sprintf(modelstring, response)
        stats::glm(formula = as.formula(fo), family = fo.family, data = dset) %>%
            broom::tidy() %>%
                dplyr::filter(grepl(filterstr, term)) %>%
                dplyr::mutate(response = response, fo = fo)
    })
    data.table::rbindlist(glmlist, id = "model_name") %>%
        as.data.frame %>%
        dplyr::mutate(conf.low = estimate - qnorm(1- 0.05/2) * std.error,
                      conf.high = estimate + qnorm(1- 0.05/2) * std.error)
}

import_filter_data <- function(file, excluded = c()) {
    pseq.full <- readRDS(file)
    pseq.meta <- filter.phenotype.data(pseq.full, excluded = excluded)
    phyloseq::sample_data(pseq.full) <- phyloseq::sample_data(pseq.meta)
    return(pseq.full)
}

import_salt_pseq <- function(file = "data/phfinrisk_genus_all_drop50k_2018-11-16.RDs") {
    pseq.salt <- readRDS(file)

    pseq.salt.meta <- meta(pseq.salt) %>%
        tibble::rownames_to_column(var = "rowname") %>%
        mutate(SEX = factor(ifelse(MEN == "Female", 1, 0)),
               ANYDRUG = factor(ifelse(BL_USE_RX_C03  == 1 | BL_USE_RX_C07 == 1 |
                                       BL_USE_RX_C08 == 1  | BL_USE_RX_C09 == 1, 1, 0)),
               HYPERTENSION = factor(ifelse(SYSTM >= 140 | DIASM >= 90 | ANYDRUG == 1, 1, 0)),
               PULSEPRESSURE = SYSTM - DIASM,
               MAP = 2./3.*DIASM + 1./3.*SYSTM,
               SYSTM = myscale(SYSTM),
               DIASM = myscale(DIASM),
               PULSEPRESSURE = myscale(PULSEPRESSURE),
               MAP = myscale(MAP),
               dUNA = myscale(NA.)) %>%
        select("rowname", "BL_AGE", "SEX", "dUNA", "HYPERTENSION", "SYSTM",
               "DIASM", "PULSEPRESSURE", "MAP") %>%
        na.omit %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = "rowname")
    
    phyloseq::sample_data(pseq.salt) <- phyloseq::sample_data(pseq.salt.meta)
    return(pseq.salt)
}

calculateadonis <- function(dset,
                            matrix,
                            responses = list(model_1 = "MAP", model_2 = "SYSTM", model_3 = "DIASM",
                                          model_4 = "PULSEPRESSURE", model_5 = "HYPERTENSION"),
                            covariates = c("BL_AGE", "SEX"),
                            npermutations = 99,
                            maxcores = 100) {
    mclapply(responses, function(response) {
        fo <- sprintf("matrix ~ %s + %s", response, paste(covariates, collapse = " + "))
        ad <- adonis(formula = as.formula(fo), data = dset, permutations = npermutations)
        ad$aov.tab %>%
            as.data.frame %>%
            tibble::rownames_to_column(var = "term") %>%
            dplyr::mutate(response = response)
      }, mc.cores = min(maxcores, length(responses)))
}

prepare.maaslin <- function (pseq,
                             tsvfile = "merged_phfinrisk_genus.tsv",
                             conffile = "maaslin_config.config",
                             includedvars = c("BL_AGE", "BMI", "SEX"),
                             readpclrows = "Archaea.Archaea-",
                             core = FALSE,
                             core.detection = 0.001,
                             core.prevalence = 0.01) {
    pseq.comp <- microbiome::transform(pseq, "compositional")
    if (core) {
         pseq.comp  <- microbiome::core(pseq.comp, detection = core.detection, prevalence = core.prevalence)
    }
    pheno <- t(microbiome::meta(phyloseq::sample_data(pseq.comp)))
    sample.ids <- colnames(pheno)
    abud.rel <- abundances(pseq.comp)
    rownames(abud.rel) <- sapply(rownames(abud.rel), replace.brackets)
    write.table(t(rbind(sample.ids, pheno, abud.rel)), file = tsvfile, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
    data <- file(conffile , open='wt')
    pcl.contents <- c("Matrix: Metadata",
                      paste(c("Read_PCL_Rows: ", paste(includedvars, collapse=",")), collapse = ""),
                      "",
                      "Matrix: Adundance",
                      paste("Read_PCL_Rows:", readpclrows))
    writeLines(pcl.contents, con=data)
    close(data)
}

 merge.pheno.abu <- function(pseq,
                             core = FALSE,
                             core.detection = 0.001,
                             core.prevalence = 0.01) {
     pseq.comp <- microbiome::transform(pseq, "compositional")
     if (core) {
         pseq.comp  <- microbiome::core(pseq.comp, detection = core.detection, prevalence = core.prevalence)
    }
     pheno <- microbiome::meta(phyloseq::sample_data(pseq.comp))
     sample.ids <- rownames(pheno)
     abud.rel <- abundances(pseq.comp)
     rownames(abud.rel) <- sapply(rownames(abud.rel), replace.brackets)
     cbind(pheno, t(abud.rel))
}


legend_extractor <- function (tmp) {
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

ylabl_extractor <- function(tmp) {
  yl <- which(grepl('axis.title.y.left', sapply(tmp$grobs, function(x) x$name)))
  tmp$grobs[[yl]]
}

xlabb_extractor <- function(tmp) {
  xl <- which(grepl('axis.title.x.bottom', sapply(tmp$grobs, function(x) x$name)))
  tmp$grobs[[xl]]
}

vectortolist <- function(c) {
  l <- as.list(c)
  names(l) <- l
  l
}

mergediversities <- function(alphadiversity, betadiversity,
                             responses = c("MAP", "SYSTM", "DIASM", "PULSEPRESSURE", "HYPERTENSION")) {
    lapply(vectortolist(names(alphadiversity)), function(x) {
        alpha.select <- alphadiversity[[x]] %>% select(estimate, p.value, response)
        beta.rbind <- betadiversity[[x]] %>%
            data.table::rbindlist(id = "model_name") %>%
            dplyr::filter(term %in% responses) %>%
            dplyr::select(response, R2)
        merge( alpha.select , beta.rbind, by = "response") %>%
            merge(names.dset %>% select(Covariate, Category, Name), by.x ="response", by.y = "Covariate") %>%
            mutate(Qstar = ifelse(p.value < 0.05, '*', ' '))
    })
}

getdescriptions <- function() {
    data("phfinrisk_metadatadesc", package="finriskmetagcommon")
    names.dset <- filter(phfinrisk_metadatadesc,
                         Ignored.Covariate.in.Cross.Sectional.Analysis.Aaro20181115==0)
    bind_rows(names.dset,
              data.frame(Covariate = c("MAP", "PULSEPRESSURE", "HYPERTENSION", "SEX"),
                         Category = rep("Physical", 4),
                         Name =c( "Mean arterial pressure", "Pulse pressure", "Hypertension", "Female"),
                         Desc = c("Mean arterial presure, 2./3.*DIASM + 1./3.*SYSTM)",
                                  "Pulse pressure, SYSTM - DIASM",
                                  "Hypertension, SYSTM >=140 or DIAS >= 90 or BP mediaction",
                                  "Sex is female, True for female, False for male")))
}

maaslinwrapper <- function(pseq, looped, forced, taxa, tempstr = "%s/temp/maaslin-%s")  {
    tempdir <- sprintf(tempstr, Sys.getenv("HOME"), format(Sys.time(), '%s'))
    dir.create(tempdir)
    cwd <- getwd()
    setwd(tempdir)
    prepare.maaslin(pseq,
                    tsvfile = "maaslin.tsv",
                    conffile = "maaslin.config",
                    includedvars = c(looped, forced),
                    readpclrows = taxa)
    Maaslin("maaslin.tsv",
            "maaslin",
            strInputConfig = "maaslin.config",
            strModelSelection="none",
            dSignificanceLevel = 1.0,
            fAllvAll = TRUE,
            strForcedPredictors = forced)
    ret <- read.table("maaslin/maaslin.txt", header=TRUE, sep='\t')
    setwd(cwd)
    ret
}

mygrep <- function(word, list, ignorecase = TRUE) {
    list[grepl(word, list, ignore.case = ignorecase)]
}

cc <- function(...) {
    c2l(c(...))
}


c2l <- function(...) {
    l <- as.list(c(...))
    names(l) <- c(...)
    l
}

pub.p <- function(p, Nfdr = FALSE) {
    p <- as.numeric(p)
    if (Nfdr) p <- p.adjust(p, method="BH", n = Nfdr)
    ifelse(p < 0.01, ifelse(p<0.001, sprintf("%.2e", p), sprintf("%.3f", p)), sprintf("%.2f", p))
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

mykable <- function(x, ...) {
  capture.output(x <- print(x))
  knitr::kable(x, ...)
}

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
                                  included = c("MEN",
                                               "BL_AGE",
                                               "SYSTM",
                                               "DIASM",
                                               "BMI",
                                               "CURR_SMOKE",
                                               "PREVAL_DIAB",
                                               "Q57X",
                                               "BL_USE_RX_C09",
                                               "BL_USE_RX_C03",
                                               "BL_USE_RX_C07",
                                               "BL_USE_RX_C08",
                                               "NA.",
                                               "BP_TREAT",
                                               "PREVAL_CHD",
                                               "PREVAL_CR_ANYCANC"),
                                  allowna = c("NA.",
                                              "BP_TREAT",
                                              "PREVAL_CHD",
                                              "PREVAL_CR_ANYCANC")) {
    meta(pseq) %>%
        tibble::rownames_to_column(var = "rowname") %>%
        dplyr::select(rowname, included) %>%
        tidyr::drop_na(myin(included, allowna, complement = TRUE)) %>%
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
    alphadiversity  <- microbiome::alpha(pseq, index = index)
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
                         filterstr = ".") {
    glmlist <- lapply(responses, function(response) {
        is.logistic <- length(unique(pull(dset, response))) < min_n_for_continuous
        fo.family <- ifelse(is.logistic, stats::binomial, stats::gaussian)
        fo <- sprintf(modelstring, response)
        stats::glm(formula = as.formula(fo), family = fo.family, data = dset) %>%
            broom::tidy(conf.int = TRUE, conf.level = 0.95, exponentiate = is.logistic) %>%
                dplyr::filter(grepl(filterstr, term)) %>%
                dplyr::mutate(response = response, fo = fo)
    })
    data.table::rbindlist(glmlist, id = "model_name") %>%
        as.data.frame
}

import_filter_data <- function(file, included) {
    pseq.full <- readRDS(file)
    if (missing(included)) {
        pseq.meta <- filter.phenotype.data(pseq.full)
    } else {
        pseq.meta <- filter.phenotype.data(pseq.full, included = included)
    }
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
                            covariates = ".",
                            npermutations = 99,
                            maxcores = 100) {
    mclapply(responses, function(response) {
        fo <- sprintf("matrix ~ %s + %s", paste(covariates, collapse = " + "), response)
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
                             readpclrows = "") {
    pheno <- t(microbiome::meta(phyloseq::sample_data(pseq)))
    sample.ids <- colnames(pheno)
    abud.rel <- abundances(pseq)
    rownames(abud.rel) <- sapply(rownames(abud.rel), replace.brackets)
    write.table(t(rbind(sample.ids, pheno, abud.rel)), file = tsvfile, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
    data <- file(conffile , open='wt')
    pcl.contents <- c("Matrix: Metadata",
                      paste(c("Read_PCL_Rows: ", paste(includedvars, collapse=",")), collapse = ""),
                      "",
                      "Matrix: Adundance",
                      paste(c("Read_PCL_Rows: ", paste(readpclrows, collapse=",")), collapse = ""))
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
                             responses = c("MAP", "SYSTM", "DIASM", "PULSEPRESSURE", "HYPERTENSION"),
                             names.dset) {
    lapply(vectortolist(names(alphadiversity)), function(x) {
        alpha.select <- alphadiversity[[x]] %>% select(estimate, p.value, response)
        beta.rbind <- betadiversity[[x]] %>%
            data.table::rbindlist(id = "model_name") %>%
            dplyr::filter(term %in% responses) %>%
            dplyr::rename(R2.p = `Pr(>F)`) %>%
            dplyr::select(response, R2, R2.p)
        merge( alpha.select , beta.rbind, by = "response") %>%
            merge(names.dset %>% select(Covariate, Category, Name), by.x ="response", by.y = "Covariate") %>%
            mutate(Qstar = ifelse(p.value < 0.05, '*', ' '))
    })
}

getdescriptions <- function() {
    names.dset <- readRDS("mydata/phfinrisk_metadatadesc.rds") %>%
        dplyr::filter(Ignored.Covariate.in.Cross.Sectional.Analysis.Aaro20181115==0)
    bind_rows(names.dset,
              data.frame(Covariate = c("MAP", "PULSEPRESSURE", "HYPERTENSION", "SEX"),
                         Category = rep("Physical", 4),
                         Name =c( "Mean arterial pressure", "Pulse pressure", "Hypertension", "Female"),
                         Desc = c("Mean arterial presure, 2./3.*DIASM + 1./3.*SYSTM)",
                                  "Pulse pressure, SYSTM - DIASM",
                                  "Hypertension, SYSTM >=140 or DIAS >= 90 or BP mediaction",
                                  "Sex is female, True for female, False for male")))
}

maaslinwrapper <- function(pseq, looped, forced, taxa, tempstr = "%s/maaslinruns/maaslin-%s")  {
    tempdir <- sprintf(tempstr, Sys.getenv("HOME"), format(Sys.time(), '%s'))
    dir.create(tempdir)
    cwd <- getwd()
    setwd(tempdir)
    prepare.maaslin(pseq,
                    tsvfile = "maaslin.tsv",
                    conffile = "maaslin.config",
                    includedvars = c(looped, forced),
                    readpclrows = taxa)

    fileConn<-file("run.R")
    writeLines(c('.libPaths(c(.libPaths(), "~/maaslinruns/R-3.4"))',
                 'library("Maaslin")',
                 'Maaslin("maaslin.tsv",',
                 '\t"maaslin",',
                 '\tstrInputConfig = "maaslin.config",',
                 '\tstrModelSelection="none",',
                 '\tdSignificanceLevel = 1.0,',
                 '\tfAllvAll = TRUE,',
                 sprintf("\tstrForcedPredictors = \"%s\")", paste0(forced, collapse=","))), fileConn)
    close(fileConn)

    system("/apps/statistics2/R-3.4.3/bin/Rscript run.R > output.log 2> error.log")

    ret <- read.table("maaslin/maaslin.txt", header=TRUE, sep='\t') %>% as.data.frame
    setwd(cwd)
    return(ret)
}

mygrep <- function(..., word, ignorecase = TRUE, complement = FALSE) {
    c(...)[xor(grepl(word, c(...), ignore.case = ignorecase), (complement == TRUE))]
}


myin <- function(x,y, complement = FALSE) {
    x[xor(x %in% y, complement)]
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
    ifelse(p < 0.01, ifelse(p<0.001, "<0.001", sprintf("%.3f", p)), sprintf("%.2f", p))
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

mykable <- function(x, ...) {
  capture.output(x <- print(x))
  knitr::kable(x, ...)
}

myinstall.packages <- function(...) {
    list.of.packages <- c(...)
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if (length(new.packages) == 0) { return(TRUE) }
    for (package in new.packages) {
        message(sprintf("Installing: %s", package))
        myinstall.packages(gtools::getDependencies(package))
        install.packages(package)
    }
}


mydropna <- function(...) {
    c(...)[!is.na(c(...))]
}


calculate.alphadiversity <- function(pseq, vars, modelstr = "%%s ~ %s + diversity_shannon") {
    dset <- meta.merge.alphadiversity(pseq) %>%
        mutate(diversity_shannon = scale(diversity_shannon),
               MAP = oMAP,
               PULSEPRESSURE = oPULSEPRESSURE,
               SYSTM = oSYSTM,
               DIASM = oDIASM)

    lapply(vars, function(var)
        min = calculateglm(dset,
                           modelstring = sprintf(modelstr, paste0(var, collapse = " + ")),
                           filterstr = "shannon"))

}

calculate.betadiversity <- function(pseq, matrix, vars, npermutations = 999) {
    lapply(vars, function(var)
        calculateadonis(dset = meta(pseq),
                        matrix = matrix,
                        covariates = var,
                        npermutations = npermutations))
}

deseqresults <- function(modellist, names.dset) {
    vars <- c2l(names(modellist))
    lapply(c2l(vars), function(x, models = modellist) {
        name <- ifelse(x %in% c("HYPERTENSION"), sprintf("%s_1_vs_0", x), x)
        results(models[[x]], name = name) %>%
            as.data.frame %>%
            tibble::rownames_to_column("Feature") }) %>%    
        map_df(., ~as.data.frame(.x), .id="name") %>%
        dplyr::mutate(qval = p.adjust(pvalue, method="BH"),
                      Feature = gsub("g_", "", Feature)) %>%
        dplyr::filter(qval < 0.05) %>%
        merge(names.dset %>% select(Covariate, Name), by.x ="name", by.y = "Covariate")
}

prune_lactobacillus <- function(pseq, transform = "compositional", drop = TRUE) {
    microbiome::transform(pseq, transform = transform) %>%
        prune_taxa(c("g_Lactobacillus"), .) %>%
        psmelt %>%
        spread(OTU, Abundance) %>%
        dplyr::mutate(lacto.prevalent = ifelse(g_Lactobacillus > 0.001, TRUE, FALSE),
                      duna.present = ifelse(is.na(NA.), 0, 1)) %>%
        { if (drop == TRUE) dplyr::filter(., duna.present == 1) else . }
}

lm_ribbon <- function(dset, model = "g_Lactobacillus ~ NA. + BL_AGE + SEX", term = "NA.", taxon = "g_Lactobacillus") {
    ret <- lm(as.formula(model), data=dset)
    ggpredict(ret, terms = term) %>%
        dplyr::rename(!!term := x, !!taxon := predicted)
}

mytableone <- function(dset, variables) {
    { if (!missing(variables)) dplyr::select(dset, variables) else dset } %>%
        select(-Sample) %>% 
        gather(key,value) %>%
        group_by(key) %>%
        summarise_all(funs(N = n(),
                           Mean = mean(as.numeric(value), na.rm = TRUE),
                           ones = sum(value == 1, na.rm = TRUE),
                           sd = sd(as.numeric(value), na.rm = TRUE),
                           NAs = sum(is.na(value))))  %>%
        mutate_if(is.numeric, round, 3)
}

deseq.formula <- function(..., fo = "~ %s") {
    vars <- c(...)
    if ("HYPERTENSION" %in% vars)
        vars <- vars[!grepl("BL_USE", vars)]
    as.formula(sprintf(fo, paste0(vars, collapse = "+")))
}

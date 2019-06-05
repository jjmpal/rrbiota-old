# some useful functions if importing finriskmetagcommon does not work outside FIMM Atlas

# also I can't compile the package anyway as long Atlas is down


my_na_filter <- function(phf_f, variables) {
  tmp_samples <- sample_data(phf_f)
  tmp_samples$names <- rownames(tmp_samples)
  sample_data(phf_f) <- tmp_samples
  filtered_data <- filter_at(sample_data(phf_f),
                                  variables,
                                  all_vars(!is.na(.)))

  rownames(filtered_data) <- filtered_data$names

  sample_data(phf_f) <- filtered_data
  phf_f
}


#' Remove variables that have only one non-NA value when given fixed covariate is present.
#'
#' @importFrom phyloseq get_variable
#' @importFrom dplyr n_distinct filter
#' @importFrom stats na.omit
#'
#' @export
remove_unnec_vars <- function(variables, phf_f, fixed.variable) {
  if (class(phf_f) %in% c("phyloseq", "sample_data"))
    sample.df <- meta(phf_f)
  if (class(phf_f) == "data.frame")
    sample.df <- phf_f
  sample.df <- filter(sample.df, !is.na(sample.df[,fixed.variable]))

  variables_tmp=c()
  for (vi in 1:length(variables)) {
    vr.content <- sample.df[,variables[vi]]
    if (is.factor(vr.content) && (n_distinct(na.omit(vr.content)) < 2)) {
      cat(paste("dropping factor", variables[vi], ", only one level\n"))
    } else if (is.numeric(vr.content) && n_distinct(na.omit(vr.content)) < 2) {
      cat(paste("dropping numeric", variables[vi], ", only one value\n"))
    } else {
      variables_tmp <- c(variables_tmp, variables[vi])
    }
  }
  variables <- variables_tmp
}



fit_adonis <- function (phf.rel, variables, save.intermediate=FALSE, intermediate.name="phyla",
                        dist.matrix=NULL, ...) {
  # dist.matrix if supplied must be dist object and contain labels

  # filter out rows with na's
  phf.rel <- my_na_filter(phf.rel, variables)

  variables <- as.vector(remove_unnec_vars(variables, phf.rel, last(variables)))

  pheno <- meta(sample_data(phf.rel))

  pheno <- pheno[,variables]

  modelstr <- paste(variables, collapse=" + ")

  if (is.null(dist.matrix)) {
    otu <- abundances(phf.rel)
    modelstr <- paste0("t(otu) ~ ", modelstr)
    # adonis does not need +1 against log(0) issues, we can use OTU table as-is
  } else {
    # use precomputed distance matrix in model formula
    modelstr <- paste0("dist.matrix ~ ", modelstr)
    # however, one must drop the row/col pairs from distance matric that correspond to na-filtered rows
    dist.tmp <- as.matrix(dist.matrix)
    dist.method <- attr(dist.matrix, "method")
    included.samples <- rownames(dist.tmp)[rownames(dist.tmp) %in% rownames(pheno)]
    dist.matrix <- as.dist(dist.tmp[included.samples, included.samples])
    attr(dist.matrix, "method") <- dist.method
  }
  phf.adonis <- adonis(as.formula(modelstr), data=pheno, ...)

  if (save.intermediate) {
    saveRDS(phf.adonis, file=paste0("/csc/fr_metagenome/microbiome_scratch/scratch/",
                                    "data_private/adonis_intermediate/",
                                    paste(variables, collapse="_"),
                                    "_", intermediate.name,
                                    "_adonis.RDs"))
  }
  return(phf.adonis)

}


fit_asine <- function(phf_f,
                      variables,
                      intermediate.name="_asin1_",
                      extra.printouts=FALSE,
                      save.intermediate=FALSE,
                      intermediate.list.file="data_private/processed_asine.txt",
                      return.only.filename=FALSE) {

    phf_f <- my_na_filter(phf_f, variables)
    cat("\nafter na filter, fitting asin(sqrt(OTU))")

    variables <- as.vector(remove_unnec_vars(variables, phf_f, last(variables)))

    modelstr <- paste(variables, collapse=" + ")
    modelstr <- paste0("~ ", modelstr)
    cat(modelstr)
    cat("\n")

    ytaxa <- taxa(phf_f)
    cat(paste("for ", length(ytaxa), " OTUs\n"))

    pheno <- meta(sample_data(phf_f))[,variables]
    pheno[,"sample.ids"] <- rownames(pheno)

    fits <- mclapply(c2l(ytaxa), function (g) {
        x <- as.data.frame(t(otu_table(phf_f)[as.character(g),]))
        tmp <- rownames(x)
        genus.clear <- taxa2underscore(g)
        colnames(x)[[1]] <- genus.clear
        x[,"sample.ids"] <- tmp
        x <- merge(x, pheno, by="sample.ids", all.x=TRUE)
        modelstr <- paste0("asin(sqrt(", genus.clear, ")) ", modelstr)
        lmwrap <- function (modelstr, x) {
            lm(as.formula(modelstr), data=x)
        }
        lmout <- tryCatch(
            lmwrap (modelstr, x),
            error=function(cond) {
                message("lm error, result item will be NA")
                message(cond)
                return(NA)
            }
        )
        return(lmout)
    }, mc.cores = min(8, parallel::detectCores()), mc.preschedule = FALSE)
    return(fits)
}


get_fit <- function (fits.all, vi, read.from.files) {
  if (read.from.files) {
    fit <- readRDS(fits.all[[vi]])
    return(fit)
  } else {
    return(fits.all[[vi]])
  }
}


# collect statistics from a list of lists of lm fits
# (designed for a cross-sectional analysis of several variables and taxa of interest:
# a list of variables of interest, a list of linear models for all OTUs per each v-o-i)
collect_lm_stats <- function (our_vars, otu_names, fits.all, read.from.files=FALSE) {

  varlist <- c()
  otulist <- c()
  pvals <- c()
  betas <- c()
  std.errors <- c()
  model.r2 <- c()
  model.pvals <- c()
  for (vi in 1:length(our_vars)) {
    varname <- our_vars[[vi]]
    fits <- get_fit(fits.all, vi, read.from.files)
    for (gi in 1:length(otu_names)) {
      g <- otu_names[[gi]] # not used: obtain length(otu_names) duplicates of varname
      varlist <- c(varlist, varname)
      otulist <- c(otulist, g)

      # assume the variable of interest is the last coefficient
      tmp_summary <- summary(fits[[gi]])
      nr <-nrow(coef(tmp_summary))

      pvals <- c(pvals, as.numeric(coef(tmp_summary)[nr,4]))

      betas <- c(betas, as.numeric(coef(tmp_summary)[nr,1]))

      std.errors <- c(std.errors, as.numeric(coef(tmp_summary)[nr,2]))

      if ("r.squared" %in% names(summary(fits.all[[1]][[1]]))) {
        model.r2 <- c(model.r2, as.numeric(tmp_summary[["r.squared"]]))
      } else {
        model.r2 <- c(model.r2, NA)  
      }

      # model p val

      f <- tmp_summary$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      model.pvals <- c(model.pvals, p)
    }
  }

  qvals <- p.adjust(pvals, method="BH")
  model.qvals <- p.adjust(model.pvals, method="BH")


  stats.df <- data.frame(Variable = varlist,
                         OTU = otulist,
                         Pvalue = pvals,
                         Qvalue = qvals,
                         Beta = betas,
                         Std.Error = std.errors,
                         Model.R2 = model.r2,
                         Model.Pvalue = model.pvals,
                         Model.Qvalue = model.qvals,
                         Qvalue.stars5 = cut(qvals, c(-Inf, 0.001, 0.01, 0.05, 0.1, 0.25, Inf),
                                            labels=c('***','**','*','..','.', ' ')),
                         Qvalue.stars = cut(qvals, c(-Inf, 0.01, Inf),
                                            labels=c('*',' ')))

}

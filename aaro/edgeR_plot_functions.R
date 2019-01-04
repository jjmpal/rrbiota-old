library(reshape2)
library(scales)
library(ggpubr)
# plot heatmap output of edgeR heatmap
plot_edgeR_fc_heatmap <- function (tables.df, xlab.text = "Phyla", logFC.maxabs=NULL) {

  if (nrow(tables.df) == 0) {
    cat("empty heatmap\n")
    p <- ggplot(data.frame()) + geom_blank()
    return(p)
  }

  p <- ggplot(tables.df, aes(x = OTU, y = Name, fill=logFC)) + geom_tile()

  if (is.null(logFC.maxabs)) {
    p <- p + scale_fill_distiller(palette = "Spectral", name = 'log fold-change')
  } else {
    p <- p + scale_fill_distiller(palette = "Spectral", name = 'log fold-change',
                                  limits=c(-logFC.maxabs, logFC.maxabs))
  }
  p <- p + xlab(xlab.text) + ylab("")
  p <- p + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p <- p + theme(axis.text.y = element_text(size=5))

  p
}

plot_edgeR_fc_heatmap_reps <- function (tables.df, xlab.text = "Phyla", logFC.maxabs=NULL) {

  if (nrow(tables.df) == 0) {
    cat("empty heatmap\n")
    p <- ggplot(data.frame()) + geom_blank()
    return(p)
  }

  p <- ggplot(tables.df, aes(x = OTU, y = Name, fill=logFC)) + geom_tile()

  p <- p + geom_text(aes(label = paste0("P ",
                                        sprintf("%0.4f", round(PValue, digits=4)),
                                        " rep ",
                                        sprintf("%0.4f", round(PValue.reps, digits=4)),
                                        "\nQ ",
                                        sprintf("%0.4f", round(Qval, digits=4)),
                                        " rep ",
                                        sprintf("%0.4f", round(Qval.reps, digits=4)),
                                        "\nlog FC ",
                                        sprintf("%0.4f", round(logFC, digits=4)),
                                        " rep ",
                                        sprintf("%0.4f", round(logFC.reps, digits=4)))),
                         size=1)

  if (is.null(logFC.maxabs)) {
    p <- p + scale_fill_distiller(palette = "Spectral", name = 'log fold-change')
  } else {
    p <- p + scale_fill_distiller(palette = "Spectral", name = 'log fold-change',
                                  limits=c(-logFC.maxabs, logFC.maxabs))
  }
  p <- p + xlab(xlab.text) + ylab("")
  p <- p + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p <- p + theme(axis.text.y = element_text(size=5))

  p
}


plot_edgeR_beta_heatmap <- function (tables.df,  xlab.text = "Phyla") {

  if (nrow(tables.df) == 0) {
    cat("empty heatmap\n")
    p <- ggplot(data.frame()) + geom_blank()
    return(p)
  }

  p <- ggplot(tables.df, aes(x = OTU, y = Name, fill=beta)) + geom_tile()

  p <- p + scale_fill_distiller(palette = "Spectral", name = 'beta')
  p <- p + xlab(xlab.text) + ylab("")
  p <- p + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p <- p + theme(axis.text.y = element_text(size=5))

  p
}

plot_edgeR_pval_heatmap <- function (tables.df,  xlab.text = "Phyla") {
  if (nrow(tables.df) == 0) {
    cat("empty heatmap\n")
    p <- ggplot(data.frame()) + geom_blank()
    return(p)
  }

  p <- ggplot(tables.df, aes(x = OTU, y = Name, fill=Qval)) + geom_tile()

  p <- p + scale_fill_distiller(palette = "Spectral", name = 'BH-adjusted P-value')
  p <- p + xlab(xlab.text) + ylab("")
  p <- p + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p <- p + theme(axis.text.y = element_text(size=5))

  p
}

plot_edgeR_fc_heatmap_final <- function (tables.df, xlab.text = "Phyla", logFC.maxabs=NULL, text.size=10) {

  if (nrow(tables.df) == 0) {
    cat("empty heatmap\n")
    p <- ggplot(data.frame()) + geom_blank()
    return(p)
  }

  # ensure not significant pairs (no rows) are set to NA

  all.otus <- unique(tables.df[,'OTU'])
  all.names <- unique(tables.df[,pname])

  all.pairs <- expand.grid(OTU=all.otus, Name=all.names)

  tables.df <- merge(all.pairs, tables.df, all.x=TRUE, by=c('OTU',pname))

  # plot

  p <- ggplot(tables.df, aes(x = OTU, y = Name, fill=logFC)) + geom_tile(colour="white",size=0.25)

  if (is.null(logFC.maxabs)) {
    p <- p + scale_fill_distiller(palette = "RdBu", name = 'log fold-change', na.value="grey30")
  } else {
    p <- p + scale_fill_distiller(palette = "RdBu", name = 'log fold-change',
                                  limits=c(-logFC.maxabs, logFC.maxabs), na.value="grey30")
  }
  p <- p + coord_fixed()
  p <- p + xlab(xlab.text) + ylab("")
  p <- p + theme_classic(text.size)
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=floor(0.7*text.size)))
  p <- p + theme(axis.text.y = element_text(size=floor(0.7*text.size)))

  p
}

plot_maaslin_heatmap_final <- function (tables.df, xlab.text = "Phyla", beta.maxabs=NULL, text.size=10,
                                        cluster=TRUE, signif=NULL, signif.th=0.05,
                                        fname='Feature', pname='Name',vname='beta',
                                        qvalues=NULL, colours=NULL, breaks=NULL, names.df=NULL,
                                        order.by.category=FALSE, ystr=NULL,
                                        v.title.name="Effect size",
                                        sign.name="Significance", no.sign.legend=FALSE,
                                        discrete.colors=FALSE,
                                        vname_discrete=NULL) {


  # otherwise the same as edgeR heatmap plotter, but relevant df columns are different

  if (is.null(beta.maxabs) && !discrete.colors) {
    beta.maxabs <- max(abs(tables.df[,vname]))
    cat("auto choose bet.maxbs: ")
    cat(beta.maxabs)
    cat("\n")
  }

  if (nrow(tables.df) == 0) {
    cat("empty heatmap\n")
    p <- ggplot(data.frame()) + geom_blank()
    return(p)
  }

  if (!is.null(signif)) {
    tables.df[,'Qlevel'] <- tables.df[,signif] < signif.th
    tables.df[,'Qstar'] <- sapply(tables.df[,'Qlevel'], function (x) ifelse(x,'*',' '))
  }

  # ensure not significant pairs (no rows) are set to NA

  all.otus <- unique(tables.df[,fname])
  all.names <- unique(tables.df[,pname])

  all.pairs <- expand.grid(Feature=all.otus, Name=all.names)

  tables.df <- merge(tables.df, all.pairs, all.y=TRUE, by.y=c('Feature','Name'), by.x=c(fname,pname))
  cat(colnames(tables.df))
  if (!is.null(names.df)) {
    for (pn in as.character(tables.df[[pname]])) {
      tables.df[tables.df[,pname] == pn, 'Category'] <- as.character(names.df[names.df[,pname]==pn, 'Category'])
    }
  }



  # plot
  if (cluster) {
    rnams <- unique(as.character(tables.df[[fname]]))
    cnams <- unique(as.character(tables.df[[pname]]))

    # Rearrange into matrix
    # FIXME: could be better done with cast
    mat <- matrix(0, nrow=length(rnams), ncol=length(cnams))
    rownames(mat) <- rnams
    colnames(mat) <- cnams

    for (i in 1:nrow(tables.df)) {
      mat[as.character(tables.df[i, fname]),
          as.character(tables.df[i, pname])] <- tables.df[i, vname]
    }

    mat <- t(mat)

    # hclust like clustering, idea from microbiome / heat.R
    scor.mat <- cor(mat, use="pairwise.complete.obs")
    scor.mat[is.na(scor.mat)] <-0
    cind <- hclust(as.dist(1 - scor.mat))$order
    order.cols <- colnames(mat)[cind]

    scor.tmat <- cor(t(mat), use="pairwise.complete.obs")
    scor.tmat[is.na(scor.tmat)] <-0
    rind <- hclust(as.dist(1 - scor.tmat))$order
    order.rows <- rownames(mat)[rind]
    tables.df[[pname]] <- factor(tables.df[[pname]], levels=order.rows)
    tables.df[[fname]] <- factor(tables.df[[fname]], levels=order.cols)
  }

  if (!is.null(ystr)) {
    ystr <- ystr
  } else if (order.by.category) {
    # TODO this causes buggy behavior
    ystr <- paste0("reorder(paste0(", pname, ",' (', Category, ')') , as.numeric(as.ordered(Category)))")

  } else {
    ystr <- pname
  }

  vname2plot <- ifelse(discrete.colors, vname_discrete, vname)
  p <- ggplot(tables.df, aes_string(x = paste0("reorder(", fname, " , as.numeric(",fname,"))"),
                                    y = ystr,
                                    fill=vname2plot)) + geom_tile(colour="white", size=0.25)

  if (!is.null(signif)){
    p <- p + geom_point(data=tables.df, aes_string(x = fname, y = pname,  shape='Qstar'),
                                                               show.legend=TRUE, color='black', size=floor(text.size/2))
    p <- p + scale_shape_manual(name=sign.name, values=c('*'='*', ' '=' '),
                                                labels=c("*"=paste("significant at FDR", signif.th),
                                                         ' '=' '),
                                                breaks=c("*", ' '))
  }

  if (no.sign.legend) {
    p <- p + guides(shape=FALSE)
  }

  if (!discrete.colors && !is.null(qvalues) && !is.null(colours)) {
    p <- p + scale_fill_gradientn(name=v.title.name,
                             colours=colours,
                             values=rescale(c(qvalues, c(-beta.maxabs, beta.maxabs)))[1:length(qvalues)],
                             limits=c(-beta.maxabs, beta.maxabs),
                             breaks=breaks)
  } else if (discrete.colors && !is.null(colours)) {
    p <- p + scale_fill_manual(values=colours, name=v.title.name, na.value="gray50")
  
  } else if (discrete.colors) {
    p <- p + scale_fill_brewer(palette="RdBu", drop=FALSE, name=v.title.name, na.value="gray50")
  } else {
    p <- p + scale_fill_distiller(palette = "RdBu", name = v.title.name,
                                  limits=c(-beta.maxabs, beta.maxabs), na.value="grey30")
  }

  p <- p + coord_fixed()
  p <- p + xlab(xlab.text) + ylab("")
  p <- p + theme_classic(text.size)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust=0.5, size=floor(0.7*text.size)))
  p <- p + theme(axis.text.y = element_text(size=floor(0.7*text.size)))

  p
}


plot_maaslin_heatmap_final_categories <- function (categories, tables.df, ...) {

  pl <- lapply(categories, function (ctg) {
                                  tables.df.ctg <- tables.df[as.character(tables.df[,'Category']) == ctg,]
                                  plot_maaslin_heatmap_final(tables.df.ctg, ...)
                                  })

  hs <- sapply(categories, function (ctg) {
                                  tables.df.ctg <- tables.df[as.character(tables.df[,'Category']) == ctg,]
                                  nrow(tables.df.ctg)
                                  })



  do.call(ggarrange, c(pl, ncol=1, nrow=length(categories), align='v'))
}





plot_genus_covariate_crossplot <- function(pseq, genus, covariate, logY=FALSE,
                                           scale_cv=FALSE, norm_count=FALSE) {
  genus.counts <- abundances(pseq)[genus,]
  genus.counts.df <- as.data.frame(genus.counts)
  if (logY) {
    genus.counts.df <- genus.counts.df + 0.1
  }
  if (norm_count) {
    genus.counts.df <- genus.counts.df / colSums(abundances(pseq))
  }
  genus.clear <- gsub("\\)", "", gsub(" \\(", ".", genus))
  colnames(genus.counts.df)[1] <- genus.clear
  genus.counts.df[,'ids'] <- names(genus.counts)
  covariate.values <- sample_data(pseq)[,covariate]
  covariate.values.df <- meta(covariate.values)
  if (scale_cv) {
    covariate.values.df[,covariate] <- scale(as.numeric(covariate.values.df[,covariate]))
  }
  covariate.values.df[,'ids'] <- rownames(covariate.values)

  combined.df <- merge(genus.counts.df, covariate.values.df, by='ids')
  p <- ggplot(data=combined.df, aes_string(x=covariate, y=genus.clear))

  p <- p + geom_jitter(size=0.8, alpha=0.8) + theme_minimal(20)

  if (logY) {
    p <- p + scale_y_log10()
    p <- p + ylab("count (log scale)")
  }
  p
}

plot_genus_covariate_crossplot_arcsin <- function(pseq, genus, covariate, logY=FALSE,
                                           scale_cv=FALSE, norm_count=FALSE) {
  genus.counts <- abundances(pseq)[genus,]
  genus.counts.df <- as.data.frame(genus.counts)
  if (logY) {
    genus.counts.df <- genus.counts.df + 0.1
  }
  if (norm_count) {
    genus.counts.df <- genus.counts.df / colSums(abundances(pseq))
  }
  genus.counts.df <- asin(sqrt(genus.counts.df))
  genus.clear <- gsub("\\)", "", gsub(" \\(", ".", genus))
  colnames(genus.counts.df)[1] <- genus.clear
  genus.counts.df[,'ids'] <- names(genus.counts)
  covariate.values <- sample_data(pseq)[,covariate]
  covariate.values.df <- meta(covariate.values)
  if (scale_cv) {
    covariate.values.df[,covariate] <- scale(as.numeric(covariate.values.df[,covariate]))
  }
  covariate.values.df[,'ids'] <- rownames(covariate.values)

  combined.df <- merge(genus.counts.df, covariate.values.df, by='ids')
  p <- ggplot(data=combined.df, aes_string(x=covariate, y=genus.clear))

  p <- p + geom_jitter(size=0.8, alpha=0.8) + theme_minimal(20)

  if (logY) {
    p <- p + scale_y_log10()
    p <- p + ylab("count (log scale)")
  }
  p
}

plot_genus_covariate_crossplot_with_fitted <- function(pseq, genus, covariate, fitted.values, logY=FALSE,
                                                       scale_cv=FALSE, name="model") {
  genus.counts <- abundances(pseq)[genus,]
  genus.counts.df <- as.data.frame(genus.counts)
  if (logY) {
    genus.counts.df <- genus.counts.df + 0.1
  }
  genus.clear <- gsub("\\)", "", gsub(" \\(", ".", genus))
  colnames(genus.counts.df)[1] <- genus.clear
  genus.counts.df[,'ids'] <- names(genus.counts)

  fitted.counts.df <- as.data.frame(fitted.values[genus,])
  colnames(fitted.counts.df)[1] <- "fitted"
  fitted.counts.df[,'ids'] <- rownames(fitted.counts.df)

  genus.counts.df.m <- merge(genus.counts.df, fitted.counts.df, by="ids")

  covariate.values <- sample_data(pseq)[,covariate]
  covariate.values.df <- meta(covariate.values)
  if (scale_cv) {
    covariate.values.df[,covariate] <- scale(as.numeric(covariate.values.df[,covariate]))
  }
  covariate.values.df[,'ids'] <- rownames(covariate.values)

  combined.df <- merge(genus.counts.df.m, covariate.values.df, by='ids')
  combined.df[,"observed"] <- combined.df[,genus.clear]
  combined.df.m <- melt(combined.df, measure.vars=c("observed", "fitted"))
  p <- ggplot(data=combined.df.m, aes_string(x=covariate, y="value", shape="variable", colour="variable"))

  p <- p + scale_colour_manual(values=c("black", "red"), name=name)
  p <- p + scale_shape_manual(values=c(16, 4), name=name)

  p <- p + geom_jitter(size=1.5, alpha=0.5) + theme_minimal(20)
  p <- p + ylab("count")
  if (logY) {
    p <- p + scale_y_log10()
    p <- p + ylab("count (log scale)")
  }
  p <- p + ggtitle(genus.clear)
  p
}

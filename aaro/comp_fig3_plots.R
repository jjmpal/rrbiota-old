# Plots for Figure 3

library(biomformat)
library(microbiome)
library(phyloseq)
library(data.table)
library(RColorBrewer)
library(dplyr)
packageVersion('phyloseq')
library(finriskmetagcommon)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(grid)
library(gtable)
library(scales)
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/<Paste>

source('code/phyloseq/plot_catcol_palette.R')

pltdir <- 'docs/assets/figure3/'
if (!dir.exists(pltdir)) {
  dir.create(pltdir)
}

### Preliminaries: load and subset data

maindir <- "/csc/fr_metagenome/microbiome_scratch/scratch/"
datadir <- file.path(maindir, "data_private/")


# alpha diversity linear models


shannon.df <- readRDS(file.path(datadir, 'lm_stats_dichodiversities_shannon.spt0.2018-10-22.RDs'))
#isimpsn.df <- readRDS(file.path(datadir, 'lm_stats.diversities_inverse_simpson.Rds'))
#observd.df <- readRDS(file.path(datadir, 'lm_stats.observed.Rds'))


# PCA regression results


## PERMANOVA results
# (preprocessing lifted from comp_adonis_plots.R)
ident <- "20181030_nperm50"

adon.res <- readRDS(paste0("/csc/fr_metagenome/microbiome_scratch/scratch/",
                           "data_private/adonis.result.", ident, ".RDs"))


adon.res.df <- adon.res[[1]] # only the variate of interest
adon.res.all.df <- readRDS(paste0("/csc/fr_metagenome/microbiome_scratch/scratch/",
                           "data_private/adonis.result.", ident, ".all.df.RDs"))# all output of aov.tab
rm(adon.res)


# PERMANOVA: merge covariate category into df

if (!exists('phfinrisk_metadatadesc')) {
  data("phfinrisk_metadatadesc", package="finriskmetagcommon")
}
names.df <- filter(phfinrisk_metadatadesc, Ignored.Covariate.in.Cross.Sectional.Analysis.Aaro20181016==0)

md <- names.df



adon.res.df <- merge(adon.res.df, md[,c('Covariate','Category','Name')],
                       by='Covariate')

# PERMANOVA: extract residual R2 and compute model R2 from them

adon.res.all.df.rds <- filter(adon.res.all.df, Model.Covariates == "Residuals")
adon.res.all.df.rds[,"Model.R2"] <- 1 - adon.res.all.df.rds[,"R2"]

adon.res.df <- merge(adon.res.df, adon.res.all.df.rds[,c("Model.R2", "Covariate")], by="Covariate", all.x=TRUE)
adon.res.df[,"Variable.R2"] <- adon.res.df[,"R2"]
adon.res.df[,"Rest.R2"] <- adon.res.df[,"Model.R2"]

adon.res.df.sig <- filter(adon.res.df, Qval < 0.05)
adon.res.df.sub1 <- adon.res.df.sig[order(adon.res.df.sig$Variable.R2, decreasing=TRUE)[1:40],]


# merge adonis results to colorbar dfs
shannon.df <- merge(shannon.df, adon.res.df.sub1[,c('Covariate', 'Variable.R2')], by='Covariate', all.y=TRUE)




fig.height <- 4.0 + nrow(adon.res.df.sub1) * 0.33
fig.permar2.text <- "R2 of variable \nin Bray-Curtis distance" # TODO latexify expression?
fig.textsize <- 15

### Create ggplot2 objects 1: functions

# separate heatmap object so that we can use different scales if needed (results for 'observed' have different magnitude from others)

# alternative version of "heatmap": hack a barplot code instead
single_col_hm_bar <- function(m_index.df, md, m_index.name, signf.text="", nonames=FALSE, amaxb=NULL, size=20,
                              colours=c("darkblue", "blue", "white", "red", "darkred"),
                              colorvalues=rescale(c(-amaxb,-0.1, 0, 0.1, amaxb)),
                              m_index.shortname=" ") {

  m_index.df <- merge(m_index.df, md, by='Covariate')
  m_index.df.f <- filter(m_index.df, Covariate %in% adon.res.df.sub1[,'Covariate'])
  m_index.df.f[,'Qlevel'] <- m_index.df.f['Qval'] < 0.05
  m_index.df.f[,'Qstar'] <- sapply(m_index.df.f[,'Qlevel'], function (x) ifelse(x,'*',' '))
  p.m_index <- ggplot(m_index.df.f,
                      aes(x=reorder(Name, Variable.R2), y=1, fill=Estimate))
  p.m_index <- p.m_index + geom_bar(stat="identity", color="black") + coord_flip()
  if (is.null(amaxb)){
    amaxb <- signif(1.05*max(abs(m_index.df[,'Estimate'])), digits=2)
  }
  p.m_index <- p.m_index + scale_fill_gradientn(colours=colours,
                                                values=colorvalues,
                                                name=m_index.name, limits=c(-amaxb, amaxb), na.value="black")
  p.m_index <- p.m_index + theme_classic(size)


  p.m_index <- p.m_index + geom_point(data=m_index.df.f, aes(x=reorder(Name, Variable.R2), y=0.5, shape=Qstar),
                                                               show.legend=TRUE, color='black', size=size)
  p.m_index <- p.m_index + scale_shape_manual(name=signf.text, values=c('*'='*', ' '=' '),
                                              labels=c("*"="significant at FDR 0.05", ' '=' '),
                                              breaks=c("*", ' '))

  p.m_index <- p.m_index + ylab("") + xlab("")

  p.m_index <- p.m_index + scale_y_continuous(breaks=c(0.5), labels=c(m_index.shortname))

  p.m_index <- p.m_index + theme(legend.position="right", legend.title=element_text(size=size),
        			 legend.direction="vertical", legend.text=element_text(size=floor(0.75*size)),
                                 legend.box.margin=margin(0,0,0,0,"pt"),
                                 legend.margin=margin(0,0,0,0,"pt"),
                                 legend.justification="left",
                                 #axis.ticks.x = element_blank(),
                                 axis.line.y = element_line(linetype="blank"),
                                 axis.line.x = element_line(linetype="blank"),
                                 legend.title.align = 0,
                                 legend.text.align = 0)

  # guide_colorbar or guide_legend?
  #if (nonames) {
  #  p.m_index <- p.m_index + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  #} else {
  #  p.m_index <- p.m_index + theme(axis.text.y = element_text(color="black"))
  #}
  p.m_index <- p.m_index + theme(axis.text.y = element_text(color="black"))
  p.m_index <- p.m_index + theme(axis.text.x = element_text(color="black"))


  p.m_index
}

legend_extractor <- function (tmp) {
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

ylabl_extractor <- function(tmp) {
  yl <- which(grepl('axis.title.y.left', sapply(tmp$grobs, function(x) x$name)))
  ylabb <- tmp$grobs[[yl]]
  ylabb
}


xlabb_extractor <- function(tmp) {
  xl <- which(grepl('axis.title.x.bottom', sapply(tmp$grobs, function(x) x$name)))
  xlabb <- tmp$grobs[[xl]]
  xlabb
}


### Create ggplot2 objects 1: plots

## permanova

p.permar2c <- ggplot(adon.res.df.sub1, aes(x = reorder(Name, Variable.R2), y= Variable.R2, fill=Category))
p.permar2c <- p.permar2c + geom_bar(stat="identity", color='black') + coord_flip() + theme_classic(fig.textsize)
p.permar2c <- p.permar2c + ylab(fig.permar2.text)

p.permar2c <- p.permar2c + theme(axis.text.y = element_text(color="black"))
p.permar2c <- p.permar2c + theme(axis.text.x = element_text(color="black"))
p.permar2c <- p.permar2c + scale_fill_my_catcol() #scale_fill_manual(name="Phenotype category", values=get_catcol_palette())


p.permar2c <- p.permar2c + theme(legend.position="right",
                                 legend.direction="vertical",
                                 legend.box.margin=margin(0,0,0,0,"pt"),
                                 legend.margin=margin(0,0,0,0,"pt"),
                                 legend.justification="left",
                                 legend.title.align = 0,
                                 legend.text.align = 0)

print(p.permar2c)

g.permar2c <- ggplot_gtable(ggplot_build(p.permar2c))

g.permar2c.xlabb <- xlabb_extractor(g.permar2c)

g.permar2c.legend <- legend_extractor(g.permar2c)

p.permar2c.cl <- p.permar2c + theme(axis.text.x = element_blank(), axis.text.y=element_blank(),
                                    axis.title.x = element_blank(), axis.title.y=element_blank(),
                                    legend.position = "none")


g.permar2c.cl <- ggplot_gtable(ggplot_build(p.permar2c.cl))

## shannnon


p.shannon.sp0 <- single_col_hm_bar(shannon.df, md, 'Regression coefficient \nin linear model \nfor Shannon index',
                                   signf.text="Alpha diversity \nsignificance",
                                   amaxb=NULL, size=fig.textsize)

p.shannon.sp0 <- p.shannon.sp0 + scale_y_continuous(breaks=c(0.5), labels=c(expression(alpha)))
print(p.shannon.sp0)
g.shannon.sp0 <- ggplot_gtable(ggplot_build(p.shannon.sp0))
g.shannon.sp0.legend <- legend_extractor(g.shannon.sp0)
g.shannon.sp0.yaxis <- g.shannon.sp0$grobs[[3]]
p.shannon.sp0.cl <- p.shannon.sp0 + theme(legend.position = "none", axis.text.x = element_blank(),
				  axis.text.y=element_blank())
g.shannon.sp0.cl <- ggplot_gtable(ggplot_build(p.shannon.sp0.cl))


# PC axes


## arrange to single Grob

gs4 <- list(g.shannon.sp0.yaxis, g.shannon.sp0.cl$grobs[[6]], g.permar2c.cl$grobs[[6]], g.shannon.sp0.legend,
            g.shannon.sp0$grobs[[7]], g.permar2c$grobs[[7]], g.permar2c.xlabb, g.permar2c.cl$grobs[[3]],
            g.permar2c.legend)

g4 <- arrangeGrob(grobs = gs4, layout_matrix = rbind(c(1,2,8,3,NA),
                                                     c(1,2,8,3,4),
                                                     c(1,2,8,3,9),
                                                     c(1,2,8,3,NA),
                                                     c(NA,5,NA,6,NA),
                                                     c(NA,NA,NA,7,NA)),
                  widths=c(4,0.7,0.5,5,5), heights=c(3,12,8,21,1,2))


# TODO: automatize creation of gs4 grob list
# and creatiiong of layout table. should be simple enough.


ggsave(file=paste0(pltdir,'main.pdf'), plot=g4, height=fig.height, width=11)
ggsave(file=paste0(pltdir,'main.png'), plot=g4, height=fig.height, width=11)


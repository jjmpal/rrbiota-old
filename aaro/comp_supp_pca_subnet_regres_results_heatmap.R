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

source('code/phyloseq/generic_functions.R')
source("code/phyloseq/taxa_names_underscored.R")
source('code/phyloseq/netw_plot_funs.R')
source('code/phyloseq/lm_comp_functions.R')

source('code/phyloseq/netw_plot_funs.R')


source('code/phyloseq/edgeR_functions.R')
source('code/phyloseq/edgeR_plot_functions.R')

maindir <- "/csc/fr_metagenome/microbiome_scratch/scratch/"
datadir <- file.path(maindir, "data_private/")

pltdir <- 'docs/assets/PC_subnet_pheno_association/'
if (!dir.exists(pltdir)) {
  dir.create(pltdir)
}

pseqdate <- "2018-12-21"
phfg.pseq <- readRDS(paste0(datadir, "phfinrisk_genus_ok_drop50k_", pseqdate, ".RDs"))
pseq.rel.core <- core(microbiome::transform(phfg.pseq, "compositional"), detection=0.1/100, prevalence=1/100)

pseq.rel.core.rel.clr <- transform(transform(pseq.rel.core, "compositional"), "clr")


# latest sample data
phfg.pseq.all <- readRDS(paste0(datadir, "phfinrisk_genus_all_drop50k_", pseqdate, ".RDs"))

data("phfinrisk_metadatadesc", package="finriskmetagcommon")
names.df <- filter(phfinrisk_metadatadesc,
                   Ignored.Covariate.in.Cross.Sectional.Analysis.Aaro20181115 == 0)
phfg.pseq <- drop_unused_samplevars(phfg.pseq.all, names.df)

# read subnet info + rel abundances of subnets:

# from comp_fig2_plots_d50k.R
n.components <- readRDS(paste0(datadir, "spieceasi_0_", pseqdate,
                              "_components_fig2.RDs"))

n.relabunds.core.genera <- abundances(pseq.rel.core)
n.relabunds.core.genera.m <- melt(n.relabunds.core.genera) %>% transmute(OTU=Var1,
                                                                         SampleID=Var2,
                                                                         rel.abud=value)


n.relabunds.compos <- do.call("rbind",
                              lapply(unique(n.components$componentID), function (cid) {
                                       taxa.cid <- filter(n.components, componentID==cid)[,"OTU"]
                                       colSums(n.relabunds.core.genera[taxa.cid,])}))

n.relabunds.compos.clr <- transform(transform(n.relabunds.compos, "compositional"), "clr")


n.relabunds.compos.clr.m <- melt(n.relabunds.compos.clr) %>% transmute(ClusterNo=Var1,
                                                                     SampleID=Var2,
                                                                     CLR=value)

# PCA regression results

Z.prcomp.reg.res <- readRDS(paste0(datadir, "Z_prcomp_reg_res_drop50k_", pseqdate, ".RDs"))

# NOTE: if we state p-values, I believe they need be adjusted for npheno * naxis, not on per-axis base

# we choose to show the first 4 axis because first 4 are needed to explain more than 10% of the variation
# (together they explain 11% of the variation)

tmplist <- lapply(c(1,2,3,4,5), function (pcid) {
                tmp <- Z.prcomp.reg.res[[pcid]]
                tmp[['PC']] <- pcid
                tmp
              })

Z.df <-  do.call(rbind.data.frame, tmplist)
Z.df['P.adj.all'] <- p.adjust(Z.df[['P.vals']], method="BH")


# subnet regression results

res.subnets <- readRDS(paste0(datadir, "lm_stats_reslist_subnets_20181227.RDs"))
subnet_ids <- readRDS(paste0(datadir, "lm_stats_reslist_subnets_20181227_subnetids.RDs"))

tmplist.sn <- lapply(c(1,2,3,4), function (sid) {
                tmp <- res.subnets[[sid]]
                tmp[['SubnetId']] <- sid
                tmp
              })

sn.df <- do.call(rbind.data.frame, tmplist.sn)
sn.df['P.adj.all'] <- p.adjust(sn.df[['P.vals']], method="BH")

# PERMANOVA results (we do not show them but use them to select the variables that are shown)

ident <- "20181212_d50k_nperm50"

# TODO: also coefficient information should be made easily available
# (currently stored in temp results files under adonis_intermediate)
adon.res <- readRDS(paste0("/csc/fr_metagenome/microbiome_scratch/scratch/",
                           "data_private/adonis_result_", ident, ".RDs"))


adon.res.df <- adon.res[[1]] # only the variate of interest
adon.res.all.df <- readRDS(paste0("/csc/fr_metagenome/microbiome_scratch/scratch/",
                           "data_private/adonis.result.", ident, ".all.df.RDs"))# all output of aov.tab
rm(adon.res)


# PERMANOVA: merge covariate category into df

if (!exists('phfinrisk_metadatadesc')) {
  data("phfinrisk_metadatadesc", package="finriskmetagcommon")
}
names.df <- filter(phfinrisk_metadatadesc, Ignored.Covariate.in.Cross.Sectional.Analysis.Aaro20181115==0)

md <- names.df

adon.res.df <- merge(adon.res.df, md[,c('Covariate','Category','Name')],
                       by='Covariate')


source('code/phyloseq/adonis_add_r2.R')
adon.res.df <- adonis_add_r2(adon.res.df, adon.res.all.df)

adon.res.df.sig <- filter(adon.res.df, Qval < 0.05)
adon.res.df.sub1 <- adon.res.df.sig[order(adon.res.df.sig$Variable.R2, decreasing=TRUE)[1:40],]

my_df_helper <- function(mydf, betastr="Beta") {
  mydf.hmap <- mydf[mydf[['Covariate']] %in% adon.res.df.sub1[['Covariate']],]
  mydf.hmap <- merge(mydf.hmap, names.df[,c("Covariate", "Name", "Category")], by="Covariate")
  mydf.hmap[betastr] <- mydf.hmap['Estimate']
  mydf.hmap <- merge(mydf.hmap, adon.res.df.sub1[,c('Covariate', 'Variable.R2')], by="Covariate")
  mydf.hmap
}

sn.df.hmap <- my_df_helper(sn.df, betastr="Beta.Subnet")
Z.df.hmap <- my_df_helper(Z.df, betastr="Beta.PC")

p1 <- plot_maaslin_heatmap_final(sn.df.hmap, fname="SubnetId", pname="Name", vname="Beta.Subnet", xlab.text = "Subnet", signif="P.adj.all", cluster=FALSE,
                                 v.title.name = "Effect size (subnet)",
                                ystr="reorder(Name, Variable.R2)", sign.name="Significance.Subnet", no.sign.legend=TRUE)

p1; ggsave(paste0(pltdir, 'SubnetPlot.png'))

p2 <- plot_maaslin_heatmap_final(Z.df.hmap, fname="PC", pname="Name", vname="Beta.PC", xlab.text = "PC", signif="P.adj.all", cluster=FALSE,
                                 v.title.name = "Effect size (PC)",
                                 ystr="reorder(Name, Variable.R2)", sign.name="Significance.PC",
                                 no.sign.legend=TRUE)
p1 <- p1 + theme(legend.position="right",
                                 legend.direction="vertical",
                                 legend.box.margin=margin(0,0,0,0,"pt"),
                                 legend.margin=margin(0,0,0,0,"pt"),
                                 legend.justification=c(0.0,0.5),
                                 legend.title.align = 0,
                                 legend.text.align = 0)

p2 <- p2 + theme(legend.position="right",
                                 legend.direction="vertical",
                                 legend.box.margin=margin(0,0,0,0,"pt"),
                                 legend.margin=margin(0,0,0,0,"pt"),
                                 legend.justification=c(0.0,0.5),
                                 legend.title.align = 0,
                                 legend.text.align = 0)




g.p1 <- ggplot_gtable(ggplot_build(p1))
g.p2 <- ggplot_gtable(ggplot_build(p2))


legend_extractor <- function (tmp) {
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

p1.legend <- legend_extractor(g.p1)

p1cl <- p1 + theme(legend.position="none")

p2; ggsave(paste0(pltdir, 'PCPlot.png'))


p2.legend <- legend_extractor(g.p2)

p2cl <- p2 + theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.position = "none",
                   plot.background = element_blank())

tmp2 <- arrangeGrob(grobs = list(p1.legend, p2.legend),
                    layout_matrix= rbind(c(1),
                                         c(2)),
                    widths=c(1), heights=c(2,2))
pall <- plot_grid(p1cl, NULL, p2cl,NULL, tmp2, rel_widths=c(1, -0.4, 0.75, 0, 0.7), nrow=1)
pall
ggsave(paste0(pltdir, 'main.png'))

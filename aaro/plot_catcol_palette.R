library(ggsci)
library(dplyr)

# the main function is scale_fill_my_catcol
# which can be used as the usual ggplot scale_fill_*() functions

get_catcol_palette <- function() {
  if (!exists('phfinrisk_metadatadesc')) {
    data("phfinrisk_metadatadesc", package="finriskmetagcommon")
  }
  my_catcol_palette <- c("#00aaaa",
                         "#225522",
                         "#cccccc",
                         "#777777",
                         "#bb9977",
                         "#00aa00")
  allcats <- unique(phfinrisk_metadatadesc$Category)
  propercats <- allcats[!allcats %in% c(NA, "Unknown")] 
  names(my_catcol_palette) <- propercats
  my_catcol_palette
}


get_catcol_palette_grgr <- function() {
  if (!exists('phfinrisk_metadatadesc')) {
    data("phfinrisk_metadatadesc", package="finriskmetagcommon")
  }


  pal <- c("#00aa00", "#000000","#eeeeee")
  my_catcol_palette <-  colorRampPalette(pal)(n_distinct(phfinrisk_metadatadesc$Category))
  names(my_catcol_palette) <- unique(phfinrisk_metadatadesc$Category)
  my_catcol_palette
}


get_catcol_palette_gray <- function() {
  if (!exists('phfinrisk_metadatadesc')) {
    data("phfinrisk_metadatadesc", package="finriskmetagcommon")
  }
  my_catcol_palette <- gray(seq(0.0, 0.9, length=n_distinct(phfinrisk_metadatadesc$Category)))
  names(my_catcol_palette) <- unique(phfinrisk_metadatadesc$Category)
  my_catcol_palette
}


scale_fill_my_catcol <- function() {
  #scale_fill_manual(name="Metadata category", values=get_catcol_palette())
  #scale_fill_uchicago(name="Metadata category")
  scale_fill_startrek(name="Metadata category")
}

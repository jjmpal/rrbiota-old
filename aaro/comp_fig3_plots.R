library(scales)
library(gridExtra)

fig.height <- 4.0 + nrow(diversity) * 0.33
fig.permar2.text <- "R2 of variable \nin Bray-Curtis distance" # TODO latexify expression?
fig.textsize <- 15

diversity <- diversity %>%
    mutate(Qstar = ifelse(p.value < 0.5, '*', ' '))

amaxb <- diversity %>% pull(estimate) %>% abs %>% max %>% signif(., digits = 2)

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


plot.beta.diversity <- ggplot(diversity, aes(x = reorder(Name, R2), y= R2)) +
    geom_bar(stat="identity", color='black') +
    coord_flip() +
    theme_classic(fig.textsize) +
    ylab(fig.permar2.text) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black"),
          legend.position="right",
          legend.direction="vertical",
          legend.box.margin=margin(0,0,0,0,"pt"),
          legend.margin=margin(0,0,0,0,"pt"),
          legend.justification="left",
          legend.title.align = 0,
          legend.text.align = 0)

plot.gtable.beta.diversity <- ggplot_gtable(ggplot_build(plot.beta.diversity))
plot.xlab.beta.diversity <- xlabb_extractor(gtable.beta.diversity)

plot.beta.diversity <- plot.beta.diversity  +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x  =  element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
plot.gtable.beta.diversity.stripped <- ggplot_gtable(ggplot_build(plot.beta.diversity))



plot.alpha.diversity <- ggplot(diversity, aes(x=reorder(Name, R2),  y = 1, fill = estimate)) +
    geom_bar(stat="identity", color="black") +
    coord_flip() +
    scale_fill_gradientn(colours = c("darkblue", "blue", "white", "red", "darkred"),
                         values = rescale(c(-amaxb,-0.1, 0, 0.1, amaxb)),
                         name = 'Regression coefficient \nin linear model \nfor Shannon index',
                         limits = c(-1.05*amaxb, 1.05*amaxb),
                         na.value = "black") +
    theme_classic(20) +
    geom_point(aes(Name, y=0.5, shape=Qstar), show.legend=TRUE, color='black', size=20) +
    scale_shape_manual(name="",
                       values=c('*'='*', ' '=' '),
                       labels=c("*"="significant at FDR 0.05", ' '=' '),
                       breaks=c("*", ' ')) +
    xlab("") +
    ylab("") +
    theme(legend.position="right",
          legend.title=element_text(size = 20),
          legend.direction="vertical",
          legend.text=element_text(size = floor(0.75*20)),
          legend.box.margin=margin(0,0,0,0,"pt"),
          legend.margin=margin(0,0,0,0,"pt"),
          legend.justification="left",
          axis.line.y = element_line(linetype="blank"),
          axis.line.x = element_line(linetype="blank"),
          legend.title.align = 0,
          legend.text.align = 0,
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black"))+
    scale_y_continuous(breaks=c(0.5), labels=c(expression(alpha)))

plot.alpha.diversity.gtable <- ggplot_gtable(ggplot_build(plot.alpha.diversity))
plot.alpha.diversity.legend <- legend_extractor(plot.alpha.diversity.gtable)
plot.alpha.diversity.yaxis <- plot.alpha.diversity.gtable$grobs[[3]]

plot.alpha.diversity <- plot.alpha.diversity +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y=element_blank())

plot.gtable.alpha.diversity.stripped <- ggplot_gtable(ggplot_build(plot.alpha.diversity))

plot.gtable.alpha.diversity.stripped[[3]]

gs4 <- list(plot.alpha.diversity.yaxis,
            plot.gtable.alpha.diversity.stripped$grobs[[6]],
            plot.gtable.beta.diversity.stripped$grobs[[6]],
            plot.alpha.diversity.legend,
            plot.alpha.diversity.gtable$grobs[[7]],
            plot.gtable.beta.diversity$grobs[[7]],
           plot.xlab.beta.diversity)
#           plot.gtable.beta.diversity.stripped[[3]])

arrangeGrob(grobs = gs4,
            layout_matrix = rbind(
                c(1,2,8,3,NA),
                c(1,2,8,3,4),
                c(1,2,8,3,NA),
                c(1,2,8,3,NA),
                c(NA,5,NA,6,NA),
                c(NA,NA,NA,7,NA)),
            widths=c(4,0.7,0.5,5,5),
            heights=c(3,12,8,21,1,2))

ggsave(file = "alphabeta.pdf", plot=g4, height=fig.height, width=11)



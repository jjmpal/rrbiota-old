plot.alpha <- function(diversity,
                       alphabreaks = c(-Inf, -1, -0.1, -0.01, 0)) {
    diversity <- diversity %>%
        mutate(estimate_fac =  cut(estimate, breaks = alphabreaks))
    plot.alpha.diversity <- ggplot(diversity, aes(x=reorder(Name, R2),
                                                  y = 1,
                                                  fill = estimate_fac)) +
        geom_bar(stat="identity", color="black") +
        coord_flip() +
        scale_fill_brewer(palette="Blues",
                          name = 'Regression coefficient \nin linear model \nfor Shannon index',
                          drop=FALSE,
                          direction = -1) +
        #scale_color_discrete_sequential(palette = "Blues", nmax = 6, order = 2:6) +
#        scale_fill_gradientn(colours = c("blue", "white"),
#                             breaks = c(-colorlimit, 0),
#                             name = 'Regression coefficient \nin linear model \nfor Shannon index',
#                             limits = c(-colorlimit, 0),
#                             na.value = "black") +
        theme_classic(20) +
        geom_point(aes(Name, y=0.5, shape=Qstar), show.legend=FALSE, color='black', size=20) +
        scale_shape_manual(name="",
                           values=c('*'='*', ' '=' '),
                           labels=c("*"="significant at\nFDR 0.05", ' '=' '),
                           breaks=c("*", ' ')) +
        xlab("") +
        ylab("") +
        theme(legend.position="right",
              legend.title=element_text(size = 14),
              legend.direction="vertical",
              legend.text=element_text(size = floor(0.75*20)),
              legend.box.margin=margin(t = 0, unit = "cm"),
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
    return(list(plot = plot.gtable.alpha.diversity.stripped$grobs[[6]],
                yaxis = plot.alpha.diversity.yaxis,
                xaxis = plot.alpha.diversity.gtable$grobs[[7]],
                legend = plot.alpha.diversity.legend))
}

plot.beta <- function(diversity,
                      fig.permar2.text = "R2 of variable \nin Bray-Curtis distance",
                      fig.textsize = 15,
                      ymax = 0.0016) {
    diversity <- diversity %>% mutate(color = as.factor(ifelse(R2.p < 0.05, 1, 0)))
    plot.beta.diversity <- ggplot(diversity, aes(x = reorder(Name, R2), y= R2, color = color)) +
        geom_bar(stat="identity") +
        scale_color_manual(values = c("0" = "black", "1" = "red")) +
        coord_flip() +
        theme_classic(fig.textsize) +
        ylab(fig.permar2.text) +
        scale_y_continuous(limits = c(0, ymax)) +
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
    plot.xlab.beta.diversity <- xlabb_extractor(plot.gtable.beta.diversity)

    plot.beta.diversity <- plot.beta.diversity  +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x  =  element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    
    plot.gtable.beta.diversity.stripped <- ggplot_gtable(ggplot_build(plot.beta.diversity))

    return(list(plot = plot.gtable.beta.diversity.stripped$grobs[[6]],
                xaxis = plot.gtable.beta.diversity$grobs[[7]],
                xlab = plot.xlab.beta.diversity))
}

myforestplot <- function(data,
                         xlab = "",
                         ylab = "",
                         interceptone = FALSE,
                         subtitle = "",
                         my_y_scale = scale_y_continuous()) {
    ggplot(data = data,
           aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
        geom_point() +
        { if (interceptone) geom_hline(aes(fill=fake_salt), yintercept=1, linetype=2) } +
        geom_errorbar(width=0.1)+ 
        ggplot2::theme_minimal()  +
        theme(plot.title=element_text(size=16),
              axis.title=element_text(size=10),
              strip.background = element_blank()) +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.x = element_line(),
              axis.ticks.y = element_line(),
              axis.line = element_line(colour = "black")) +
        scale_x_discrete(labels=seq(3)) +
        my_y_scale + 
        xlab(xlab) +
        ylab(ylab) +
        labs(subtitle = subtitle)
}


gvlma.table.plot <- function(model) {
    { sink("/dev/null");
        gvlma.model <- gvlma::gvlma(model)
        gvlma.modek.summary <- summary(gvlma.model)
        ; sink(); }
    ret <- list("table" = gvlma.modek.summary,
                 "plot" = gvlma.model)
    return(ret)
}

plot.diversities <- function(diversities) {
    alpha.diversity.min <- plot.alpha(diversities$min, alphabreaks = seq(-1,0,0.2))
    alpha.diversity.max <- plot.alpha(diversities$max, alphabreaks = seq(-1,0,0.2))

    beta.diversity.min <- plot.beta(diversities$min)
    beta.diversity.max <- plot.beta(diversities$max)
    
    gs4 <- list(alpha.diversity.min$yaxis,
                alpha.diversity.min$plot,
                beta.diversity.min$plot,
                alpha.diversity.min$legend,
                alpha.diversity.min$xaxis,
                beta.diversity.min$xaxis,
                beta.diversity.min$xlab,
                text_grob("Age- and sex-adjusted model", size = 18),
                alpha.diversity.max$yaxis,
                alpha.diversity.max$plot,
                beta.diversity.max$plot,  
                alpha.diversity.max$xaxis,
                beta.diversity.max$xaxis,
                beta.diversity.max$xlab,
                text_grob("Multivariable-adjusted model", size = 18))
 
    g4 <- arrangeGrob(grobs = gs4,
                      layout_matrix = rbind(
                          c(8,    8,    8, NA,   NA, NA, NA),
                          c(NA, 1,    2,    NA, 3, NA, NA),
                          c(NA, 1,    2,    NA, 3, NA, 4),
                          c(NA, NA, 5,    NA,6, NA, NA),
                          c(NA, NA, NA, NA,7, NA, NA),
                          c(15,    15,    15, NA,   NA, NA, NA),
                          c(NA, 9,    10,  NA,  11, NA, NA),
                          c(NA, 9,    10,  NA,  11, NA, NA),
                          c(NA, NA, 12,   NA, 13, NA, NA),
                          c(NA, NA, NA, NA,14, NA, NA)),
                      widths=c(0.7, 4, 0.7, 0.1, 5, 0.7, 3),
                      heights=rep(c(1.5, 1, 6, 1, 2), 2))
    return(g4)
}

deseqhaetmap <- function(df, legend.name = "Log Base Two\nFold Change") {
    ggplot(df,
           aes(x = Feature,
               y = reorder(Name, deseqplotorder(Name)),
               fill = log2FoldChange)) +
        geom_tile(aes(fill = 1)) +
        geom_tile(colour="white", size=0.25) +
        scale_fill_distiller(palette = "RdBu",
                             name = legend.name,
                             limits=c(-1, 1),
                             breaks = c(-1, -0.5, 0.5, 0, 1),
                             na.value="grey30") +
        coord_fixed() +
        xlab("") +
        ylab("") +
        theme_classic() +  
        theme(axis.text.x = element_text(angle=90, size=10, hjust=1.0,vjust=0.5),
              axis.text.y = element_text(size=10),
              legend.title=element_text(size=10),
              legend.position = c(-0.09, -1))
}


deseqplotorder <- function(x) {
    y <- seq(length(x))
    y[x == "Systolic BP"] = 4
    y[x == "Diastolic BP"] = 3
    y[x == "Pulse pressure"] = 2
    y[x == "Mean arterial pressure"] = 1
    y[x == "Hypertension"] = 0
    y
}

saltplot <- function(dset) {
    ggplot(data = dset, aes(x = NA., y = g_Lactobacillus)) +
        geom_point(aes(color = SEX), size = 0.05, position = "jitter") +
        geom_line(data = lm_ribbon(dset)) +
        geom_ribbon(data = lm_ribbon(dset), aes(ymin=conf.low, ymax=conf.high), alpha=0.15) +
        scale_x_continuous(limits=c(20, 205), expand = c(0, 0)) +
        scale_y_continuous(limits=c(4.5, 8), expand = c(0, 0)) +
        scale_colour_manual(name = "SEX",
                            labels = c("Male", "Female"),
                            breaks=c("0", "1"),
                            values = c("blue", "red")) +
        xlab("24-hour urinary sodium") +
        ylab("CLR-transformed Lactobacillus abundance") +
        guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.4))) +
        theme_classic() +
        theme(plot.title = element_blank(),
              legend.title = element_blank(),
              legend.position = c(0.9,0.95))
}

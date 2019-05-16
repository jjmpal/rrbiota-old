plot.alpha <- function(diversity,  colorlimit = 0.25) {
    plot.alpha.diversity <- ggplot(diversity, aes(x=reorder(Name, R2),  y = 1, fill = estimate)) +
        geom_bar(stat="identity", color="black") +
        coord_flip() +
        scale_fill_gradientn(colours = c("blue", "white"),
                             breaks = c(-colorlimit, 0),
                             name = 'Regression coefficient \nin linear model \nfor Shannon index',
                             limits = c(-colorlimit, 0),
                             na.value = "black") +
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
                      fig.textsize = 15) {
    plot.beta.diversity <- ggplot(diversity, aes(x = reorder(Name, R2), y= R2)) +
        geom_bar(stat="identity", color="black") +
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
    plot.xlab.beta.diversity <- xlabb_extractor(plot.gtable.beta.diversity)

    plot.beta.diversity <- plot.beta.diversity  +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x  =  element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    plot.gtable.beta.diversity.stripped <- ggplot_gtable(ggplot_build(plot.beta.diversity))

      return(list(plot =plot.gtable.beta.diversity.stripped$grobs[[6]],
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
              axis.title=element_text(size=12),
              strip.background = element_blank()) +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.x = element_line(),
              axis.ticks.y = element_line(),
              axis.line = element_line(colour = "black")) +
        scale_x_discrete(labels=c("Q1", "Q2", "Q3", "Q4")) +
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

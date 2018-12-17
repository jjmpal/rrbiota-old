local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.uib.no/"
    options(repos = r)
})

if (FALSE) {
    install.packages("packrat")
    packrat::init("~/microbiome/articletwo")
    packrat::canUseGitHubDownloader()
}

getOption("repos")
install.packages("BH")
install.packages("dplyr")

install.packages("ggplot2")
install.packages("tidyr")

install.packages("Rcpp")
install.packages("rmarkdown")

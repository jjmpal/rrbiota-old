#local({r <- getOption("repos")
#       r["CRAN"] <- "http://cran.r-project.org" 
#       options(repos=r)
#})

myinstall.packages <- function(...) {
    if("gtools" %in% rownames(installed.packages()) == FALSE) {
        install.packages("gtools")
    }
    list.of.packages <- c(...)
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if (length(new.packages) == 0) { return(TRUE) }
    for (package in new.packages) {
        message(sprintf("Installing: %s", package))
        myinstall.packages(gtools::getDependencies(package))
        install.packages(package)
    }
}


#### -- Packrat Autoloader (version 0.5.0) -- ####
source("packrat/init.R")
#### -- End Packrat Autoloader -- ####

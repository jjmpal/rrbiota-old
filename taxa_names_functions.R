# Helper functions for handling FR02 taxa names

#' Remove some not-very-nice characters from taxa_names(pseq)
#'
#' spaces, hyphens, commas to underscores; parentheses to empty string
#'
#' @param pseq a phyloseq object
#'
#' @return processed taxa_names
#'
#' @export
#'
#' @importFrom phyloseq taxa_names
#'
#' @examples
#' \dontrun{
#' taxa_names(pseq) <- taxa_names_underscored(pseq)
#' }
taxa_names_underscored <- function (pseq) {
  taxanames <- taxa_names(pseq)
  taxanames <- characters_underscored(taxanames)
  taxanames
}


#' Remove some not-very-nice characters from character vector
#'
#' Helper function, intended main interface for pseqs is
#' taxa_names_underscored.
#' spaces, hyphens, commas to underscores; parentheses to empty string
#'
#' @param taxanames eg output of taxa_names(pseq)
#'
#' @return processed taxa_names
#'
#' @export
#'
#' @examples
#' \dontrun{
#' taxa_names(pseq) <- taxa_names_underscored(pseq)
#' }
characters_underscored <- function (taxanames) {
  taxanames <- unlist(lapply(taxanames, function(x) gsub(" ", "_", x)))
  taxanames <- unlist(lapply(taxanames, function(x) gsub("-", "_", x)))
  taxanames <- unlist(lapply(taxanames, function(x) gsub(",", "_", x)))
  taxanames <- unlist(lapply(taxanames, function(x) gsub("\\(", "", x)))
  taxanames <- unlist(lapply(taxanames, function(x) gsub("\\)", "", x)))
}


#' Return underscored version of a single taxa
#'
#' @param x character string
#'
#' @return underscored character string
#'
#' @export
taxa2underscore <- function(x) {
  gsub("\\)", "", gsub("\\(", "", gsub(",", "_", gsub("-", "_",  gsub(" ", "_", x)))))
}


#' Return 'pretty' taxa names for printing and plotting
#'
#' @param pseq a phyloseq object
#'
#' @return a vector of 'pretty' taxa_names
#'
#' @export
#'
#' @importFrom microbiome abundances
my_pretty_network_vertex_names <- function(pseq) {
  # TODO: sync this function with other underscoring functions

  tmpnames <- gsub("\\(BacteriaPlasmid\\)", "\nPlsmd.", colnames(t(abundances(pseq))))
  tmpnames <- gsub("\\(Bacteria\\)", "", tmpnames)
  tmpnames <- gsub("\\(Archaea\\)", "\nArch.", tmpnames)
  tmpnames <- gsub("\\(Viruses\\)", "\nVirus", tmpnames)
  tmpnames <- trimws(tmpnames)
  tmpnames <- characters_underscored(tmpnames)
  names(tmpnames) <-colnames(t(abundances(pseq)))
  tmpnames
}

#' @rdname my_pretty_network_vertex_names
#' @export
my_pretty_print_taxa_names <- my_pretty_network_vertex_names

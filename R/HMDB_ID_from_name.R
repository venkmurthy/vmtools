#' Query Human Metabolome Database by metabolite name
#'
#' This function searches the HMDB to find the canonical HMDB ID for a given metabolite name. It does look through
#' the synonym list as well. Importantly, it takes the *first* match among the names returned by the HMDB search engine
#' and if none, the first match among the synonyms.
#'
#' @param met.name string containing metabolite name
#' @return string with canonical HMDB ID
#' @keywords HMDB
#' @examples
#' HMDB_ID_from_one_name("valine")
#' HMDB_ID_from_one_name("ectoine")
#'
#' test.list <- c("valine-d8","phenylalanine-d8","glycine","alanine","serine","threonine","methionine","glutamate",
#' "asparagine","glutamine","histidine","arginine","lysine","valine","leucine","isoleucine","phenylalanine",
#' "tyrosine","tryptophan","proline","hydroxyproline","ornithine","citrulline","taurine","GABA",
#' "dimethylglycine","ADMA","SDMA","NMMA","2-aminoisobutyric acid","kynurenic acid","1-methylhistamine",
#' "N-carbamoyl-beta-alanine","thiamine","niacinamide","betaine")
#'
#' t0 <- Sys.time()
#' HMDB_ID_from_name(test.list[1:20])
#' Sys.time() -t0
#' @export
library(dplyr)
library(rvest)
library(xml2)
library(purrr)


HMDB_ID_from_name <- function(met.names,max.depth=25) {
  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"
  xml.url <- "https://hmdb.ca/metabolites/%s.xml"

  # simplify list
  met.names <- unlist(met.names)

  # Safe version of read_xml
  safe_read_xml <- function(u) {
    con <- file(u)
    on.exit(purrr::safely(close)(con))
    purrr::safely(read_xml)(con)
  }
  safe_read_html <- function(u) {
    con <- file(u)
    on.exit(purrr::safely(close)(con))
    purrr::safely(read_html)(con)
  }

  # Initialize output list
  out.ids <- rep("",length(met.names))
  names(out.ids) <- met.names

  # Loop over all elements
  for (i in seq_along(met.names)) {
    # Pull out one name
    x <- met.names[i]

    # if the element is NA then out.id is NA and move to the next
    if (is.na(x) | nchar(x)==0) {
      out.ids[i] <- NA
      next
    }

    found.id <- FALSE
    search.page <- 1
    ids.checked <- 0

    while (!found.id & ids.checked<=max.depth) {
      u <- sprintf(search.url,search.page,URLencode(tolower(x),reserved=TRUE))

      pause.length <- 0.2
      repeat {
        h <- safe_read_html(u)

        if (is.null(h$error)) {
          h <- h$result
          break
        } else {
          Sys.sleep(pause.length)
          pause.length <- pause.length * 1.5
        }
      }
      hmdb.ids <- h %>% html_nodes("div.result-link") %>% html_nodes("a") %>% html_text()

      for (j in seq_along(hmdb.ids)) {
        pause.length <- 0.2

        repeat {
          xml.entry <- safe_read_xml(sprintf(xml.url,hmdb.ids[j]))

          if (is.null(xml.entry$error)) {
            xml.entry <- xml.entry$result
            break
          } else {
            Sys.sleep(pause.length)
            pause.length <- pause.length * 1.5
          }
        }

        xml.names <- xml.entry %>% xml_find_first("//metabolite/name") %>% xml_text()
        xml.syns <- xml.entry %>% xml_find_all("//metabolite/synonyms/synonym") %>% xml_text()

        if (tolower(gsub(" ","",gsub("-","",x))) %in% tolower(c(xml.names,xml.syns))) {
          out.ids[i] <- hmdb.ids[j]
          found.id <- TRUE
          break
        }
      }

      ids.checked <- ids.checked + length(hmdb.ids)
    }

    if(!found.id) { out.ids[i] <- NA }
  }

  return(out.ids)
}



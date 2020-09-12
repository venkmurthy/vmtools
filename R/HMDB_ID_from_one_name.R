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
#' @export
library(dplyr)
library(rvest)
library(xml2)
library(purrr)

HMDB_ID_from_one_name <- function(met.name) {

  if(is.na(met.name)) return(NA)

  # Convert input met name to lower case
  q <- tolower(met.name)

  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"
  xml.url <- "https://hmdb.ca/metabolites/%s.xml"

  # Safe version of read_xml
  safe_read_xml <- purrr::safely(read_xml)

  # Initialize variables
  i <- 1
  hmdb.ids <- c()

  # Loop, downloading search results and extracting HMDB IDs
  repeat {
    p <- read_html(sprintf(search.url,i,URLencode(q,reserved=TRUE))) %>% html_nodes("div.result-link") %>% html_nodes("a") %>% html_text()

    hmdb.ids <- c(hmdb.ids, p)

    if (length(p)==0 | i >= 25) break
    else i <- i + 1
  }

  # Download all XMLs
  hmdb.xmls <- vector(mode="list",length=length(hmdb.ids))
  for (i in 1:length(hmdb.ids)) {

    # Set up a pause before reloading
    pause.length <- 0.05

    repeat {
      x <- safe_read_xml(sprintf(xml.url,hmdb.ids[i]))

      if (is.null(x$error)) {
        hmdb.xmls[[i]] <- x$result
        break
      } else {
        Sys.sleep(pause.length)
        pause.length <- pause.length * 2
      }
    }
  }

  # Pull out all names and synonyms
  hmdb.names <- lapply(hmdb.xmls, FUN=function(x) { x %>% xml_find_first("//metabolite/name") %>% xml_text( )})
  hmdb.syns <- lapply(hmdb.xmls, FUN=function(x) { x %>% xml_find_all("//metabolite/synonyms/synonym") %>% xml_text( )})

  # Initialize variable
  match.id <- NA

  # See if query metabolite matches any of the names
  match.names <- which(q == tolower(hmdb.names))

  if (length(match.names)>=1) {
    # If it matches a name, use that ID
    match.id <- hmdb.ids[match.names[1]]
  } else {
    # Loop over all synonym lists and see if any of those match
    i <- 1

    while(is.na(match.id) & i <= length(hmdb.syns)) {
      match.names <- which(q == tolower(hmdb.syns[[i]]))
      if (length(match.names)>=1) {
        match.id <- hmdb.ids[i]
      } else {
        i <- i + 1
      }
    }
  }

  return(match.id)
}



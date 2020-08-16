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
#' HMDB_ID_from_name("valine")
#' HMDB_ID_from_name("ectoine")
#' @export

library(dplyr)
library(rvest)
library(xml2)

HMDB_ID_from_name <- function(met.name) {
  # Convert input met name to lower case
  q <- tolower(met.name)

  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"
  xml.url <- "https://hmdb.ca/metabolites/%s.xml"

  # Initialize variables
  i <- 1
  hmdb.ids <- c()

  # Loop, downloading search results and extracting HMDB IDs
  repeat {
    p <- read_html(sprintf(search.url,i,q)) %>% html_nodes("div.result-link") %>% html_nodes("a") %>% html_text()

    hmdb.ids <- c(hmdb.ids, p)

    if (length(p)==0 | i >= 100) break
    else i <- i + 1
  }

  # Download all XMLs
  hmdb.xmls <- lapply(hmdb.ids, FUN=function(x) { Sys.sleep(0.05); read_xml(sprintf(xml.url,x)) } )

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


#' Query Human Metabolome Database by metabolite name
#'
#' This function searches the HMDB to find the canonical HMDB ID for a given ID.
#'
#' @param id string with ID starting with "HMDB" or number with just the numeric part
#' @return string with canonical HMDB ID
#' @keywords HMDB
#' @examples
#' HMDB_ID_from_ID(22)
#' HMDB_ID_from_ID("HMDB0006022")
#' @export
HMDB_ID_from_ID <- function(id) {
  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"

  # If id is numeric, convert to string, if string make sure it starts with HMDB
  if (is.numeric(id)) {
    id <- sprintf("HMDB%07i",id)
  } else if (is.character(id) & substr(id,1,4)=="HMDB") {

  } else {
    return(NA)
  }

  Sys.sleep(0.4)

  url <- httr::HEAD(sprintf(search.url,1,id))$all_headers[[1]]$headers$location
  url.parts <- unlist(strsplit(url, "/"))

  return(url.parts[length(url.parts)])
}


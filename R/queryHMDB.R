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
    p <- read_html(sprintf(search.url,i,URLencode(q,reserved=TRUE))) %>% html_nodes("div.result-link") %>% html_nodes("a") %>% html_text()

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
#' @param id vector of strings with ID starting with "HMDB" or numbers with just the numeric part
#' @return vector of strings with canonical HMDB IDs
#' @keywords HMDB
#' @examples
#' HMDB_ID_from_ID(22)
#' HMDB_ID_from_ID("HMDB0006022")
#' @export
HMDB_ID_from_ID <- function(ids) {
  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"

  # Simplify list
  ids <-unlist(ids)

  # Initialize output list
  # out.ids <- vector(mode="character",length=length(ids))
  out.ids <- rep("",length(ids))
  names(out.ids) <- ids

  # Loop over all elements
  for (i in seq_along(ids)) {
    # Pull out one ID
    x <- ids[i]

    # It is a character string that starts with "HMDB" then fine, likewise if it is a number. Otherwise NA
    if(is.character(x) & substr(x,1,4)=="HMDB") { out.ids[i] <- x
    } else if (is.numeric(x)) { out.ids[i] <- sprintf("HMDB%07i",x)
    } else out.ids[i] <- NA

    # If the ID is good above,
    if (!is.na(out.ids[i])) {

      # Initialize good.q
      good.q <- FALSE

      # Set up a pause before reloading
      pause.length <- 0.05

      # If the pulled headers are not good, and we haven't waited too long already
      while (!good.q & pause.length <10) {
        # Pause
        Sys.sleep(pause.length)

        # Increase next pause
        pause.length <- pause.length * 2

        # Retry pulling headers
        h <- httr::HEAD(sprintf(search.url,1,out.ids[i]))$all_headers

        # Did we get two sets of headers?
        good.q <- (length(h) >= 2)

        # If we got two sets of headers
        if (good.q) {
          # check that there are actually headers
          good.q <- ("headers" %in% names(h[[2]]))

          # If there are headers, is there a location header?
          if (good.q) (good.q <- ("location" %in% names(h[[2]][["headers"]])))
        }
      }

      # pull out the location header
      out.ids[[i]] <- h[[2]][["headers"]][["location"]]
      }
    }

  # Split out the HMDB ID part of the URL
  out.ids <- sapply(strsplit(out.ids,"\\/"),FUN=function(x) x[length(x)])

  return(out.ids)
}



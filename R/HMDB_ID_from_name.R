#' Query Human Metabolome Database by metabolite name
#'
#' This function searches the HMDB to find the canonical HMDB ID for a given metabolite name. It does look through
#' the synonym list as well. Importantly, it takes the *first* match among the names returned by the HMDB search engine
#' and if none, the first match among the synonyms.
#'
#' @param met.name string containing metabolite name
#' @param max.depth number indicating how many entries to be checked (may be slightly more than this,
#' depending on entries per page)
#' @param max.tries number of times to retry loading a page
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
#' "N-carbamoyl-beta-alanine","thiamine","niacinamide","betaine","PE(P-36:4)/PE(O-36:5)","PE(P-36:2)/PE(O-36:3)",
#' "PE(P-38:6)/PE(O-38:7)","PE(P-38:4)/PE(O-38:5)","PE(P-40:6)/PE(O-40:7)")
#'
#' t0 <- Sys.time()
#' HMDB_ID_from_name(test.list)
#' Sys.time() -t0
#'
#' HMDB_ID_from_name("C5 carnitine")
#' @export
library(purrr)

HMDB_ID_from_name <- function(met.names,max.depth=25,max.tries=5) {
  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"
  xml.url <- "https://hmdb.ca/metabolites/%s.xml"

  # simplify list
  met.names <- unlist(met.names)

  # Safe version of read_xml
  safe_read_xml <- function(u) {
    con <- file(u)
    on.exit(purrr::safely(close)(con))
    purrr::quietly(purrr::safely(xml2::read_xml,quiet=TRUE))(con)$result
  }
  safe_read_html <- function(u) {
    con <- file(u)
    on.exit(purrr::safely(close)(con))
    purrr::quietly(purrr::safely(xml2::read_html,quiet=TRUE))(con)$result
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

    # Initialize variables
    found.id <- FALSE
    search.page <- 1
    ids.checked <- 0

    # Keep looking up entries while not found and max depth of parsing not exceeded
    while (!found.id & ids.checked<=max.depth) {
      u <- sprintf(search.url,search.page,URLencode(tolower(x),reserved=TRUE))

      # Repeatedly try to read a specific page until success
      pause.length <- 5
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

      # Extract HMDB IDs
      hmdb.ids <- h %>% rvest::html_nodes("div.result-link") %>% rvest::html_nodes("a") %>% rvest::html_text()

      # For each HMDB ID
      for (j in seq_along(hmdb.ids)) {
        pause.length <- 5

        # Repeatedly try to pull XML entry until success
        num.tries <- 0
        found.xml <- NA
        repeat {
          xml.entry <- safe_read_xml(sprintf(xml.url,hmdb.ids[j]))

          if (is.null(xml.entry$error)) {
            xml.entry <- xml.entry$result
            found.xml <- TRUE
            break
          } else {
            num.tries <- num.tries + 1
            if (num.tries > max.tries) {
              found.xml <- FALSE
              break
            } else {
              Sys.sleep(pause.length)
              pause.length <- pause.length * 1.5
            }
          }
        }

        if (found.xml==TRUE) {
          # Pull out primary name and synonyms
          xml.names <- xml.entry %>% xml2::xml_find_first("//metabolite/name") %>% xml2::xml_text()
          xml.syns <- xml.entry %>% xml2::xml_find_all("//metabolite/synonyms/synonym") %>% xml2::xml_text()

          # Combine name and synonyms, create variations with hyphen vs. dash and make all lower case
          names.and.syns <- c(xml.names,xml.syns)
          names.and.syns <- tolower(unique(c(names.and.syns,gsub(" ","-",names.and.syns),gsub("-"," ",names.and.syns))))

          # Condense multiple spaces into one in the search key
          while(grepl("  ",x)) {
            gsub("  "," ",x)
          }

          # Convert search key to lower case
          x <- tolower(x)

          # Check if search key is in the names/synonyms list, if found, add it to output list and move on
          if (x %in% names.and.syns) {
            out.ids[i] <- hmdb.ids[j]
            found.id <- TRUE
            break
          }
        }
      }

      # Increment counter of number of IDs checked
      ids.checked <- ids.checked + length(hmdb.ids)
    }

    # If ID is not found, set to NA
    if(!found.id) { out.ids[i] <- NA }
  }

  return(out.ids)
}



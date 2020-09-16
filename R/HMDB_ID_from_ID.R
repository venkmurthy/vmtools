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
library(dplyr)
library(rvest)
library(xml2)
library(purrr)

HMDB_ID_from_ID <- function(ids) {
  # Set up some search constants
  search.url <- "https://hmdb.ca/unearth/q?button=&page=%i&query=%s&searcher=metabolites"

  # Simplify list
  ids <-unlist(ids)

  # Initialize output list
  out.ids <- rep("",length(ids))
  names(out.ids) <- ids

  # Loop over all elements
  for (i in seq_along(ids)) {
    # Pull out one ID
    x <- ids[i]

    # It is a character string that starts with "HMDB" then fine, likewise if it is a number. Otherwise NA
    if (is.na(x)) {
      out.ids[i] <- NA
      next
    } else if (is.character(x) & substr(x,1,4)=="HMDB") { out.ids[i] <- x
    } else if (is.numeric(x)) { out.ids[i] <- sprintf("HMDB%07i",x)
    } else {
      out.ids[i] <- NA
      next
    }

    # Initialize good.q
    good.q <- FALSE

    # Set up a pause before reloading
    pause.length <- 0.05

    # Repeatedly load, with longer and longer pauses if we fail
    repeat {
      # Retry pulling headers
      h <- httr::HEAD(sprintf(search.url,1,out.ids[i]))$all_headers

      # Find which is last element of chain which is status 302 (Found)
      list.302 <- which(sapply(h,FUN=function(x) x$status)==302)
      h.final <- ifelse(length(list.302)>0,max(list.302),NA)

      # If we found a record
      if (!is.na(h.final) & is.finite(h.final)) {

        # pull out the location header
        out.ids[[i]] <- h[[h.final]][["headers"]][["location"]]

        break
      }

      # Pause
      Sys.sleep(pause.length)

      # Increase next pause
      pause.length <- pause.length * 1.5
    }
  }

  # Split out the HMDB ID part of the URL
  out.ids <- sapply(strsplit(out.ids,"\\/"),FUN=function(x) x[length(x)])

  return(out.ids)
}



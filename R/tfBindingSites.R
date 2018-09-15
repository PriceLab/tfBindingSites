#' @importFrom methods new is
#' @import BiocGenerics
#' @import FimoClient
#' @import RSQLite
#' @import GenomicRanges
#'
#' @name tfBindingSites-class
#' @rdname tfBindingSites-class
#' @aliases tfBindingSites
#' @exportClass tfBindingSites

.tfBindingSites <- setClass("tfBindingSites",
                       slots = c(
                          tf="character",
                          fimoClient="FimoClientClass",
                          quiet="logical",
                          state="environment"
                          ),
                       )


#------------------------------------------------------------------------------------------------------------------------
#' Create a tfBindingSitesCtor object
#'
#' @description
#' once configured, this object reports motif match to supplied (ChIP-seq) sequences - recommended
#'
#' @rdname tfBindingSitesCtor
#' @aliases tfBindingSitesCtor
#'
#' @param tf character, a transcription factor gene symbol name
#' @param quiet do or do not print progress information
#'
#' @return An object of the tfBindingSites class
#'
#'
#' @export

tfBindingSitesCtor <- function(tf, quiet=TRUE)
{
   FIMO_HOST <- "localhost"
   FIMO_PORT <- 5558
   fimoClient <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=quiet)

   state <- new.env(parent=emptyenv())
   obj <- .tfBindingSites(tf=tf, fimoClient=fimoClient, quiet=quiet, state=state)

   obj

} # tfBindingSites ctor
#------------------------------------------------------------------------------------------------------------------------

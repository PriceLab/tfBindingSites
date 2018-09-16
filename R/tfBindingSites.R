#' @importFrom methods new is
#' @import BiocGenerics
#' @import FimoClient
#' @import RSQLite
#' @import GenomicRanges
#' @import trena
#'
#' @name tfBindingSites-class
#' @rdname tfBindingSites-class
#' @aliases tfBindingSites
#' @exportClass tfBindingSites

.tfBindingSites <- setClass("tfBindingSites",
                       slots = c(
                          tf="character",
                          genome="BSgenome",
                          fimoClient="FimoClientClass",
                          quiet="logical",
                          state="environment",
                          db="SQLiteConnection"
                          ),
                       )

#------------------------------------------------------------------------------------------------------------------------
setGeneric('getChipSeqHits', signature='obj', function (obj, chromosome, min, max) standardGeneric ('getChipSeqHits'))
setGeneric('getFimoHits',   signature='obj',  function (obj, tbl.regions, pvalThreshold) standardGeneric ('getFimoHits'))
setGeneric('getMotifMatcherHits', signature='obj', function(obj, tbl.regions, pfms, matchThreshold)
              standardGeneric('getMotifMatcherHits'))
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

tfBindingSitesCtor <- function(tf, genome, chipseq.db.file, fimoHost, fimoPort, quiet=TRUE)
{
   fimoClient <- FimoClient(fimoHost, fimoPort, quiet=quiet)

   if(!quiet)
      printf("connecting to SQLite ChIP-seq db.file '%s'", chipseq.db.file)

   db <- dbConnect(dbDriver("SQLite"), chipseq.db.file)
   state <- new.env(parent=emptyenv())

   obj <- .tfBindingSites(tf=tf, genome=genome, fimoClient=fimoClient, quiet=quiet, state=state, db=db)

   obj

} # tfBindingSites ctor
#------------------------------------------------------------------------------------------------------------------------
#' get a table of ChIP-seq hits of the specified tf, in the specified region
#'
#' @rdname getChipSeqHits
#' @aliases getChipSeqHits
#'
#' @param obj An object of class tfBindingSites
#' @param chromosome character
#' @param min numeric
#' @param max numeric
#'
#' @return A data.frame
#'
#' @export

setMethod('getChipSeqHits', 'tfBindingSites',

     function (obj, chromosome, min, max){
        query <- sprintf("select * from chipseq where chr='%s' and start >= %d and end <= %d and tf = '%s'",
                          chromosome, min, max, obj@tf)
        tbl <- dbGetQuery(obj@db, query)
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
#' get a table of FIMO hits in the regions specified in the incoming data.frame
#'
#' @rdname getFimoHits
#' @aliases getFimoHits
#'
#' @param obj An object of class tfBindingSites
#' @param tbl a data.frame with (at least) chr, start, end columns
#'
#' @return a data.frame
#'
#' @export

setMethod('getFimoHits', 'tfBindingSites',

     function (obj, tbl.regions, pvalThreshold){
        stopifnot(all(c("chr", "start", "end") %in% colnames(tbl.regions)))
        seqs <- with(tbl.regions, as.list(as.character(getSeq(obj@genome, chr, start, end))))
        names(seqs) <- sprintf("region.%02d", seq_len(nrow(tbl.regions)))
        tbl <- requestMatch(obj@fimoClient, seqs, pvalThreshold)
        if(nrow(tbl) > 0){
           tbl$chr <- tbl.regions$chr[1]
           tbl$pvalScore <- -log10(tbl$p.value)
           for(r in 1:nrow(tbl.regions)){
              chrom.start <- tbl.regions[r, "start"]
              rows.oi <- grep(sprintf("region.%02d", r), tbl$sequence.name)
              if(!obj@quiet) printf(" %d rows.oi for row %d", rows.oi, r)
              if(length(rows.oi) > 0){
                 tbl$start[rows.oi] <- tbl$start[rows.oi] + chrom.start;
                 tbl$stop[rows.oi] <- tbl$stop[rows.oi] + chrom.start;
                 } # ir rows.oi
              } # for r
           tbl <- tbl[, c("chr", "start", "stop", "pvalScore", "strand", "motif", "sequence.name", "matched.sequence")]
           }
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
#' get a table of bioconductor (MotifMatcher) hits in the regions specified in the incoming data.frame
#'
#' @rdname getMotifMatcherHits
#' @aliases getMotifMatcherHits
#'
#' @param obj An object of class tfBindingSites
#' @param tbl.regions a data.frame with (at least) chr, start, end columns
#' @parm pfms a list of position frequency matrices, from e.g. MotifDb
#'
#' @return a data.frame
#'
#' @export

setMethod('getMotifMatcherHits', 'tfBindingSites',

   function(obj, tbl.regions, pfms, matchThreshold=0.7){
      if(colnames(tbl.regions)[1] == "chr")
         colnames(tbl.regions)[1] <- "chrom"
      mm <- MotifMatcher("hg38", as.list(pfms))
      tbl <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=as.integer(100 * matchThreshold))
      tbl
      })

#------------------------------------------------------------------------------------------------------------------------

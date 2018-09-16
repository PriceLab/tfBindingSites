library(tfBindingSites)
library(RUnit)
library(BSgenome.Hsapiens.UCSC.hg38)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_getChipSeqHits()
   test_runFIMO()
   test_motifMatcher()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   tfbs <- tfBindingSitesCtor("CREM",
                              genome=BSgenome.Hsapiens.UCSC.hg38,
                              chipseq.db.file <- "~/s/data/public/human/remap-2018/remap-all.sqlite",
                              fimoHost="localHost",
                              fimoPort=5558,
                              quiet=TRUE)

   checkEquals(is(tfbs), "tfBindingSites")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_getChipSeqHits <- function()
{
   printf("--- test_getChipSeqHits")

   tf <- "CREM"
   tfbs <- tfBindingSitesCtor(tf,
                              genome=BSgenome.Hsapiens.UCSC.hg38,
                              chipseq.db.file= "~/s/data/public/human/remap-2018/remap-all.sqlite",
                              fimoHost="localHost",
                              fimoPort=5558,
                              quiet=TRUE)

   tbl.cs.crem <- getChipSeqHits(tfbs, "chr1", 167517115, 167525862)
   checkEquals(dim(tbl.cs.crem), c(2, 7))
   checkEquals(tbl.cs.crem$tf,  c(tf, tf))

    # save(tbl.cs.crem, file="../extdata/tbl.ts.crem.2.chr1.167517115.167525862.RData")

} # test_getChipSeqHits
#------------------------------------------------------------------------------------------------------------------------
test_runFIMO <- function()
{
   printf("--- test_runFIMO")

   tf <- "CREM"
   tfbs <- tfBindingSitesCtor(tf,
                              genome=BSgenome.Hsapiens.UCSC.hg38,
                              chipseq.db.file="~/s/data/public/human/remap-2018/remap-all.sqlite",
                              fimoHost="localHost",
                              fimoPort=5558,
                              quiet=TRUE)

   load(system.file(package="tfBindingSites", "extdata", "tbl.ts.crem.2.chr1.167517115.167525862.RData"))
   cutoff <- 0.005
   tbl.fimo <- getFimoHits(tfbs, tbl.cs.crem, pvalThreshold=cutoff)
   dim(tbl.fimo)
   checkEquals(dim(tbl.fimo), c(6, 8))
   checkTrue(all(tbl.fimo$pvalScore >= -log10(cutoff)))

} # test_runFIMO
#------------------------------------------------------------------------------------------------------------------------
test_motifMatcher <- function()
{
   printf("--- test_runMotifMatcher")

   tf <- "CREM"
   tfbs <- tfBindingSitesCtor(tf,
                              genome=BSgenome.Hsapiens.UCSC.hg38,
                              chipseq.db.file="~/s/data/public/human/remap-2018/remap-all.sqlite",
                              fimoHost="localHost",
                              fimoPort=5558,
                              quiet=TRUE)

   load(system.file(package="tfBindingSites", "extdata", "tbl.ts.crem.2.chr1.167517115.167525862.RData"))

   pfms <- MotifDb["Hsapiens-HOCOMOCOv10-CREM_HUMAN.H10MO.C"]

     # downstream constraints on pfms
   # MotifMatcher, .matchPwmForwardAndReverse(sequence, pfms[[i]], names(pfms)[i], min.match.percentage, quiet)
   # pfms[[i]] must be strictly a numeric matrix
   # name(pfms)[i] must be just a string

   cutoff <- 0.7
   tbl.bioc <- getMotifMatcherHits(tfbs, tbl.cs.crem, pfms, matchThreshold=cutoff)
   checkEquals(dim(tbl.bioc), c(5,13))
   checkTrue(all(tbl.bioc$motifRelativeScore >= cutoff))

} # test_motifMatcher
#------------------------------------------------------------------------------------------------------------------------
# explore the crem binding sites in the 9MB chr1 region explored here:
#   ~/github/trenaShinyApps/uw/survey/UM0463F-chr1/display.R
#
view_hits <- function()
{

   library(igvR)
   igv <- igvR()
   setGenome(igv, "hg38")
   showGenomicRegion(igv, "chr1:167517115-167525862")
   showGenomicRegion(igv, "chr1:162,671,057-172,522,660")

   tf <- "CREM"
   tfbs <- tfBindingSitesCtor(tf,
                              genome=BSgenome.Hsapiens.UCSC.hg38,
                              chipseq.db.file="~/s/data/public/human/remap-2018/remap-all.sqlite",
                              fimoHost="localHost",
                              fimoPort=5558,
                              quiet=TRUE)
   tbl.cs.crem <- getChipSeqHits(tfbs, "chr1", 167517115, 167525862)
   dim(tbl.cs.crem)
   track <- DataFrameAnnotationTrack("cs", tbl.cs.crem, "blue", trackHeight=25)
   displayTrack(igv, track)
   track <- DataFrameAnnotationTrack("peaks", tbl.cs.crem[, c("chr", "peakStart", "peakEnd", "name")])
   displayTrack(igv, track)

   tbl.fimo <- getFimoHits(tfbs, tbl.cs.crem, 10e-3)
   dim(tbl.fimo)
   tbl.wig <- tbl.fimo[, c("chr", "start", "stop", "pvalScore")]
   colnames(tbl.wig) <- c("chr", "start", "end", "value")
   track <- DataFrameQuantitativeTrack("fimo", tbl.wig, "darkRed", trackHeight=25, autoscale=TRUE)
   displayTrack(igv, track)

   pfms <- MotifDb["Hsapiens-HOCOMOCOv10-CREM_HUMAN.H10MO.C"]
   tbl.bioc <- getMotifMatcherHits(tfbs, tbl.cs.crem, pfms, matchThreshold=0.7)
   dim(tbl.bioc)
   tbl.wig <- tbl.bioc[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
   colnames(tbl.wig) <- c("chr", "start", "end", "value")
   tbl.wig$value <- as.integer(100 * tbl.wig$value)
   track <- DataFrameQuantitativeTrack("bioc", tbl.wig, "purple", trackHeight=25, autoscale=TRUE)
   displayTrack(igv, track)

      # look at the left binding site
   showGenomicRegion(igv, "chr1:167,517,518-167,518,017")
      # at the right
   showGenomicRegion(igv, "chr1:167,524,554-167,525,051")

   # conclusions
   # fimo 167,517,808-167,517,817, score 3.12
   # bioc 167,517,807-167,517,816, score 70     - forgive off-by-one error
   #
   #   actual sequence:                 CCATGAGGTGG
   #   consensusString(pfms[[1]])       C??TGACGTCA
   #   agreement                        |||||| ||
   #
   #   Hsapiens-HOCOMOCOv10-CREM_HUMAN.H10MO.C
   #          1    2    3 4    5    6    7    8    9   10   11
   #     A 0.03 0.47 0.27 0 0.03 0.97 0.03 0.23 0.10 0.03 0.83
   #     C 0.57 0.20 0.40 0 0.00 0.03 0.90 0.00 0.03 0.70 0.10
   #     G 0.27 0.33 0.27 0 0.97 0.00 0.07 0.77 0.00 0.23 0.03
   #     T 0.13 0.00 0.07 1 0.00 0.00 0.00 0.00 0.87 0.03 0.03


} # view_hits
#------------------------------------------------------------------------------------------------------------------------

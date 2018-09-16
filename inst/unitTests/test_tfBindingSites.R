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
                              quiet=FALSE)

   load(system.file(package="tfBindingSites", "extdata", "tbl.ts.crem.2.chr1.167517115.167525862.RData"))
   tbl.cs <- getFimoHits(tfbs, tbl.cs.crem)
   #checkEquals(dim(tbl.cs), c(2, 7))
   #checkEquals(tbl.ts$tf <- tf)

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
                              quiet=FALSE)

   load(system.file(package="tfBindingSites", "extdata", "tbl.ts.crem.2.chr1.167517115.167525862.RData"))

   pfms <- MotifDb["Hsapiens-HOCOMOCOv10-CREM_HUMAN.H10MO.C"]

     # downstream constraints on pfms
   # MotifMatcher, .matchPwmForwardAndReverse(sequence, pfms[[i]], names(pfms)[i], min.match.percentage, quiet)
   # pfms[[i]] must be strictly a numeric matrix
   # name(pfms)[i] must be just a string


   tbl.cs <- getMotifMatcherHits(tfbs, tbl.cs.crem, pfms, matchThreshold=0.7)
   checkTrue(nrow(tbl.cs) >= 4)

} # test_motifMatcher
#------------------------------------------------------------------------------------------------------------------------

library(tfBindingSites)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   tfbs <- tfBindingSitesCtor("CREM")
   checkEquals(is(tfbs), "tfBindingSites")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------

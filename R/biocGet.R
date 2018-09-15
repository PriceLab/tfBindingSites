biocGet <-
function (pkgName) 
{
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkgName)
}

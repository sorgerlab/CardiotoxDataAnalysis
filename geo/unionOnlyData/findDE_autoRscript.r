library(Biobase)
library(GEOquery)
library(limma)


# load series and platform data from GEO
options('download.file.method.GEOquery' = 'wget')
#options('download.file.extra.GEOquery' = '-L')

gset <- getGEO("GSE12260", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
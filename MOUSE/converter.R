library(Matrix)

path = "/home/nicola/Scrivania/DATASET_MOUSE"
setwd(path)


# generate single-cell RNA seq data
count <- read.table(file = "GSM3395913_L34690.tsv", sep = '\t', header=T, row.names = 1)
gbm <- t(count)

# save sparse matrix
sparse.gbm <- Matrix(gbm , sparse = T )
head(sparse.gbm)
## Market Exchange Format (MEX) format
writeMM(obj = sparse.gbm, file="matrix.mtx")

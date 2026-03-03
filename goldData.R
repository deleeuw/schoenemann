gold <-
structure(list(V1 = c(13.5, 1, 13.5, 16.5, 12, 9.5, 2, 9.5, 5.5, 
5.5, 15, 3, 5.5, 11, 16.5, 8, 5.5), V2 = c(13, 17, 6.5, 3, 4, 
13, 15.5, 1, 5, 11, 13, 8.5, 6.5, 10, 15.5, 8.5, 2), V3 = c(17, 
13, 12, 14, 11, 15, 3, 4, 1, 7.5, 5, 9, 10, 7.5, 16, 6, 2), V4 = c(15L, 
6L, 15L, 17L, 15L, 13L, 3L, 11L, 4L, 1L, 2L, 10L, 5L, 12L, 7L, 
9L, 8L), V5 = c(16L, 6L, 13L, 17L, 14L, 12L, 3L, 9L, 4L, 7L, 
2L, 8L, 1L, 15L, 10L, 11L, 5L), V6 = c(17, 10, 12.5, 11, 15, 
16, 5, 6, 2, 1, 8, 9, 7, 14, 12.5, 3, 4), V7 = c(17L, 9L, 15L, 
13L, 16L, 14L, 2L, 6L, 1L, 4L, 5L, 11L, 8L, 10L, 12L, 7L, 3L), 
    V8 = c(16L, 4L, 14L, 13L, 12L, 15L, 2L, 9L, 5L, 7L, 3L, 6L, 
    1L, 17L, 10L, 11L, 8L)), class = "data.frame", row.names = c(NA, 
-17L))

goldMatrix <- as.matrix(gold)
goldComplete <- smacofCompleteRectangular(goldMatrix)
goldDist <- as.dist(goldComplete)
weightMatrix <- 1 - diag(25)
weightMatrix[1:17, 1:17] <- 0
weightMatrix[17 + 1:8, 17 + 1:8] <- 0
weightDist <- as.dist(weightMatrix)
goldData <- makeMDSData(goldDist, weightDist)

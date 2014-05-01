totalData <- read.csv(file = "general_statistics_1000_1.csv", header = TRUE, sep = ',')
numCrosses <- totalData[1]
numInd <- totalData[2]
portionChrom <- totalData[5]
portionGenome <- totalData[6]
modData = list(numCrosses = numCrosses, numInd = numInd, portionChrom = portionChrom, portionGenome = portionGenome)


data <- read.csv(file = "~/GitHub/CrossSim/general_statistics_1000_1.csv", header = TRUE, sep = ',')
data <- data[order(data$Number.of.Back.Crosses), ]
qplot(data$Number.of.Back.Crosses, data$Percent.Selected.Chromosome, data = data, geom = "jitter")

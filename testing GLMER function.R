library(flowReMix)
cummean <- function(x) cumsum(x) / 1:length(x)
data(rv144)
#set.seed(502)
set.seed(504)
par(mfrow = c(1, 1), mar = rep(4, 4))
data <- rv144
leaves <- unique(data$population)
selected_populations = c(1:3, 7)
data <- subset(data, population %in% leaves[selected_populations])
data$population=factor(data$population)
data <- subset(data, stim != "sebctrl")
data$treatment <- as.numeric(data$stim == "env")
data$ptid <- as.numeric(data$ptid)
data$ptid[data$vaccine == "VACCINE"] <- data$ptid[data$vaccine == "VACCINE"] * 10^4
data$prop <- data$count / data$parentcount
data$population <- as.factor(data$population)
data <- data[order(data$population, data$ptid, data$stim, decreasing = FALSE), ]










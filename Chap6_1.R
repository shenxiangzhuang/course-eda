raw_matrix <- matrix(data = 
                       c(25.3, 32.1, 38.8, 25.4, 
                         25.3, 29.0, 31.0, 21.1, 
                         18.2, 18.8, 19.3, 20.3,
                         18.3, 24.3, 15.7, 24.0,
                         16.3, 19.0, 16.8, 17.5), ncol = 5)

raw_data <- data.frame(raw_matrix,
                       row.names = c("东北","中北", "南部", "西部"))
colnames(raw_data) <- c("<=8", "9-11", "12", "13-15", ">=16")

nrow <- dim(raw_matrix)[1]
ncol <- dim(raw_matrix)[2]
# add new col and new row
data0 <- raw_data
data0$newMed <- apply(raw_data, 1, median)
data0$lastTime = rep(0, nrow)
data0[nrow+1, ] = rep(0, ncol+2)
data0[nrow+2, ] = rep(0, ncol+2)
rownames(data0) <- c(rownames(data0)[1:nrow], "newMed", "lastTime")

# 第一次行平滑
data1 <- data0
data1$lastTime <- data1$newMed

subRowMed <- function(xs){
  xsMed <- median(xs)
  x <- xs - xsMed
}

resid_matrix <- t(apply(data1[c(1:nrow), c(1:ncol)], 1, subRowMed))
data1[c(1:nrow), c(1:ncol)] <- resid_matrix
data1["newMed", c(1:(ncol+2))] <- c(apply(resid_matrix, 2, median), 0,
                                    median(data1$lastTime[c(1:nrow)]))
# 第一次列平滑
data1_1 <- data1
data1_1["lastTime", ] <- data1["lastTime", ] + data1["newMed", ]

subColMed <- function(ys){
  ysMed <- median(ys)
  y <- ys - ysMed
}

resid_matrix <- apply(data1_1[c(1:nrow), c(1:ncol)], 2, subColMed)
data1_1[c(1:nrow), c(1:ncol)] <- resid_matrix

data1_1$newMed[c(1:nrow)] <- apply(resid_matrix, 1, median)
data1_1$lastTime[c(1:nrow)] <- data1_1$lastTime[c(1:nrow)] - data1_1["lastTime", "lastTime"]

#data1_1$newMed["lastTime"] <- median()

# 第二次行平滑
data2 <- data1_1
data2$lastTime <- data1_1$newMed + data1_1$lastTime

subRowMed <- function(xs){
  xsMed <- median(xs)
  x <- xs - xsMed
}

resid_matrix <- t(apply(data2[c(1:nrow), c(1:ncol)], 1, subRowMed))
data2[c(1:nrow), c(1:ncol)] <- resid_matrix

data2["newMed", c(1:(ncol+2))] <- c(apply(resid_matrix, 2, median), 0,
                                    median(data2$lastTime[c(1:nrow)]))

# 第二次列平滑
data2_2 <- data2

data2_2["lastTime", ] <- data2["lastTime", ] + data2["newMed", ]

subColMed <- function(ys){
  ysMed <- median(ys)
  y <- ys - ysMed
}

resid_matrix <- apply(data2_2[c(1:nrow), c(1:ncol)], 2, subColMed)
data2_2[c(1:nrow), c(1:ncol)] <- resid_matrix

data2_2$newMed[c(1:nrow)] <- apply(resid_matrix, 1, median)

data2_2["lastTime", ] <- data2["newMed", ] + data2["lastTime", ]

data2_2["NewEffect", ] <- data2["lastTime", ] - median(as.matrix(data2_2["lastTime", c(1:ncol)]))
data2_2["NewEffect", "lastTime"] <- data2_2["lastTime", "lastTime"] + median(as.matrix(data2_2["lastTime", c(1:ncol)]))
data2_2

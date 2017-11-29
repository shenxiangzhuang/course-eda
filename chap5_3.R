# 以儿童身高数据为例
age_data <- c(109, 113, 115, 116, 119, 120,
              121, 124, 126, 129, 130, 133,
              134, 135, 137, 139, 141, 142)

height_data <- c(137.6, 147.8, 136.8, 140.7, 132.7, 145.4,
                 135.0, 133.0, 148.5, 148.3, 147.5, 148.8,
                 133.2, 148.7, 152.0, 150.6, 165.3, 149.9)

raw_data <- data.frame(x=age_data, y=height_data)

# Wald
Wald.get.group <- function(raw_data){
  # sort by x
  data = raw_data[order(raw_data$x), ]
  
  n <- length(data$x)
  # if x_m == x_m+1 then remove both of them
  if(data$x[n/2] == data$x[n/2+1] &&
     n%%2 == 0){
    data = data[c(1:(n/2-1), c((n/2+2):n)), ]
  }
  # whatever n is even number or odd number
  group1 <- data[c(1:floor(n/2)), ]
  group2 <- data[c((floor(n/2)+1):n),]
  return(list(group1=group1, group2=group2))
}

Wald.get.para <- function(groups){
  b_W <- (sum(groups$group2$y) - sum(groups$group1$y))/(sum(groups$group2$x) - sum(groups$group1$x))
  x_bar <- mean(c(mean(groups$group1$x), mean(groups$group2$x)))
  y_bar <- mean(c(mean(groups$group1$y), mean(groups$group2$y)))
  a_W <- y_bar - b_W*x_bar
  return(list(a=a_W, b=b_W))
}

# get together all the steps
Wald.run <- function(raw_data){
  groups <- Wald.get.group(raw_data)
  paras <- Wald.get.para(groups)
  print(paras)
  
  plot(raw_data, col='blue')
  abline(a = paras$a, b = paras$b, col='green')
}

Wald.run(raw_data)

# Nair and Shrivastava 方法

NS.get.group <- function(raw_data, group_perc = 1/3){
  # sort by x
  data = raw_data[order(raw_data$x), ]
  n <- length(data$x)
  
  group_size = floor(n*group_perc)
  
  groupL <- data[c(1:group_size), ]
  groupU <- data[c(n-group_size:n),]
  return(list(groupL=groupL, groupU=groupU))
}

NS.get.para <- function(groups){
  xL_bar <- mean(groups$groupL$x)
  xU_bar <- mean(groups$groupU$x)
  x_bar <- mean(c(xL_bar, xU_bar))
  
  yL_bar <- mean(groups$groupL$y)
  yU_bar <- mean(groups$groupU$y)
  y_bar <- mean(c(yL_bar, yU_bar))
  
  b_NS <- (yU_bar - yL_bar)/(xU_bar - xL_bar)
  a_NS <- yL_bar - b_NS*xL_bar
  
  return(list(a=a_NS, b=b_NS))
}

# get together all the steps
NS.run <- function(raw_data){
  groups <- NS.get.group(raw_data)
  paras <- NS.get.para(groups)
  print(paras)
  
  plot(raw_data, col='blue')
  abline(a = paras$a, b = paras$b, col='green')
}

NS.run(raw_data)


# Barlett
Barlett.get.group <- function(raw_data, group_perc = 1/3){
  # sort by x
  data = raw_data[order(raw_data$x), ]
  n <- length(data$x)
  
  group_size = floor(n*group_perc)
  
  groupL <- data[c(1:group_size), ]
  groupU <- data[c(n-group_size:n),]
  return(list(groupL=groupL, groupU=groupU))
}

Barlett.get.para <- function(groups){
  xL_bar <- mean(groups$groupL$x)
  xU_bar <- mean(groups$groupU$x)
  x_bar <- mean(c(xL_bar, xU_bar))
  
  yL_bar <- mean(groups$groupL$y)
  yU_bar <- mean(groups$groupU$y)
  y_bar <- mean(c(yL_bar, yU_bar))
  
  b_Barl <- (yU_bar - yL_bar)/(xU_bar - xL_bar)
  a_Barl <- y_bar - b_NS*x_bar
  
  return(list(a=a_Barl, b=b_Barl))
}

# get together all the steps
Barlett.run <- function(raw_data){
  groups <- NS.get.group(raw_data)
  paras <- NS.get.para(groups)
  print(paras)
  
  plot(raw_data, col='blue')
  abline(a = paras$a, b = paras$b, col='green')
}

Barlett.run(raw_data)



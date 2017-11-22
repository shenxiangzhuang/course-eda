# # 三组耐抗线的进阶实现[克服振荡]

# 分成左，中，右三个组，使得组的大小尽可能相等
get.groups <- function(raw_data){
  # in: raw_data
  # out: group data[list]
  
  # sort
  sort_data <- raw_data[with(raw_data, order(x)), ]
  length_data <- length(raw_data$x)
  # group size
  k <- floor(length_data/3)
  rem_size <- length_data %% 3
  if(rem_size == 0){
    group_data <- list(left=sort_data[c(1:k), ],
                       mid=sort_data[c((k+1):(2*k)), ],
                       right=sort_data[c((2*k+1):(3*k)), ])
    return(group_data)
  }
  else if (rem_size == 1){
    group_data <- list(left=sort_data[c(1:k), ],
                       mid=sort_data[c((k+1):(2*k+1)), ],
                       right=sort_data[c((2*k+2):(3*k+1)), ])
    return(group_data)
  }
  else if(rem_size == 2){
    group_data <- list(left=sort_data[c(1:k+1), ],
                       mid=sort_data[c((k+2):(2*k+1)), ],
                       right=sort_data[c((2*k+2):(3*k+2)), ])
    return(group_data)
  }
}


# 分别获取各组的总括点
get.median.points <- function(group_data){
  median.points <- list(
    left = c(median(group_data$left$x), median(group_data$left$y)),
    mid = c(median(group_data$mid$x), median(group_data$mid$y)),
    right = c(median(group_data$right$x), median(group_data$right$y))
  )
  return(median.points)
}


# 根据总括点获取a, b
get.a.b <- function(median.points){
  xL <- median.points$left[1]
  yL <- median.points$left[2]
  
  xM <- median.points$mid[1]
  yM <- median.points$mid[2]
  
  xR <- median.points$right[1]
  yR <- median.points$right[2]
  
  b0 <- (yR-yL)/(xR-xL)
  a0 <- (1/3)*((yL-b0*(xL-xM)) + yM + (yR-b0*(xR-xM)))
  
  return(c(a0, b0, xM))
}



# data

x_bad <- c(-4, -3, -2, -1, 0, 1, 2, 3, 12)
y_bad <- c(0, 0, 0, 0, 0, 0, -5, 5, 1)
raw_data <- data.frame(x=x_bad, y=y_bad)

group_data <- get.groups(raw_data)
median.points <- get.median.points(group_data)

# 原始数据图&初始总括点
plot(raw_data$x, raw_data$y, pch=4)
# 总括点
medians_x <- unlist(median.points)[c(1,3,5)]
medians_y <- unlist(median.points)[c(2,4,6)]
points(medians_x, medians_y, pch=16, col="green")
# 分组线
abline(v=sort(group_data$left$x,decreasing = TRUE)[1], lty=2, col="green")
abline(v=sort(group_data$mid$x,decreasing = TRUE)[1], lty=2, col="green")

# 传入第一次迭代的模型参数
a0_b0_xM <- get.a.b(median.points)
a0 <- a0_b0_xM[1]
b0 <- a0_b0_xM[2]
xM <- a0_b0_xM[3]
# 初始化迭代参数
iter <- 0
b_abs_delta <- 1

# 迭代限制与收敛条件
MAX_ITER = 100
B_ABS_DELTA = 0.01
b_s <- c(b0)
a_s <- c(a0)

while(iter<MAX_ITER && b_abs_delta > B_ABS_DELTA){
  print(b0)
  # 残差代替y,得到新的分组数据
  left_r <- group_data$left$y-(a0+b0*(group_data$left$x-xM))
  mid_r <- group_data$mid$y-(a0+b0*(group_data$mid$x-xM))
  right_r <- group_data$right$y-(a0+b0*(group_data$right$x-xM))
  
  group_data <- list(
    left = data.frame(x=group_data$left$x, y=left_r),
    mid = data.frame(x=group_data$mid$x, y=mid_r),
    right = data.frame(x=group_data$right$x, y=right_r)
  )
  # 计算新的总括点
  median.points <- get.median.points(group_data)
  # 计算新的参数
  a0_b0_xM <- get.a.b(median.points)
  a0 <- a0_b0_xM[1]
  b0 <- a0_b0_xM[2]
  xM <- a0_b0_xM[3]
  
  a_s <- c(a_s, a0)
  b_s <- c(b_s, b0)
  # 迭代控制
  b_abs_delta <- abs(b0)  # 这里注意是变化量
  iter <- iter + 1
  
}

# 迭代结果整理
aFinal <- sum(a_s)
bFinal <- sum(b_s)
fit_model <- function(x){aFinal + bFinal*(x-xM)}


# 耐抗线拟合效果
y_result <- fit_model(raw_data$x)
residual <- raw_data$y - y_result
sum_square_residual <- sum(residual^2)
opar <- par(no.readonly = T)
par(mfrow=c(1,2))
plot(raw_data$x, residual, 
     main = paste("Residual=",as.character(sum_square_residual)),
     pch=4, ylim = c(-20,20))


# 和最小二乘法的比较
ols_fit_model <- lm(y~x, data=raw_data)
ols_residual <- ols_fit_model$residuals
sum_square_residual_ols <- sum(ols_residual^2)
plot(raw_data$x, ols_residual,
     main = paste("Residual=",as.character(sum_square_residual_ols)),
     pch=4, ylim = c(-20,20))


# 拟合直线的比较

x_loc <- c(-5:15)
y_robust <- aFinal+bFinal*(x_loc-xM)
par(opar)
plot(raw_data, ylim=range(y_robust), col='blue')
points(medians_x, medians_y, pch=16, col="green")

lines(x_loc, y_robust, col='red')
abline(ols_fit_model, col='blue')
legend(-4,14, 
       legend = c('Robust', 'Ols'), 
       col = c('red', 'blue'),
       lty =1)

# a_s, b_s
par(mfrow=c(1,2))
plot(a_s[1:20], main = "a_s")
plot(b_s[1:20], main = "b_s")


# 克服振荡现象，安全的迭代算法
# 一般振荡的b是异号的，我们最近选取相邻的异号的两个b
get.two.b <- function(b_s){
  for(b_i in b_s){
    for(b_j in b_s[c(2:length(b_s))]){
      if (b_i * b_j < 0)
        return(c(b_i, b_j))
    }
  }
}

b01 <- get.two.b(cumsum(b_s))
b0 <- b01[1]
b1 <- b01[2]

# delta_r_b0, delta_r_b1
# 对原始数据分组的操作
group_data <- get.groups(raw_data)
median.points <- get.median.points(group_data)

# function: get_delta_r_bx 
get.delta.r.bx <- function(bx){
  return(median(group_data$right$y - (group_data$right$x-xM)*bx) - 
           median(group_data$left$y - (group_data$left$x-xM)*bx))
}
delta_data_b0 <- get.delta.r.bx(b0)
delta_data_b1 <- get.delta.r.bx(b1)


# 线性迭代内插求出b2
# 插值函数
get.b2 <- function(b0, b1, delta_data_b0, delta_data_b1){
  b2 <- b1-delta_data_b1*(b1-b0)/(delta_data_b1-delta_data_b0)
  return(b2)
}

# 迭代前的一些数据
DELTA_R_B2_MIN = 0.0001
b2 = get.b2(b0, b1, delta_data_b0, delta_data_b1)
delta_r_b2 <- get.delta.r.bx(b2)
# 开始迭代
while(delta_r_b2 > DELTA_R_B2_MIN){
  if(b0*delta_r_b2>0) b0 <- b2
  else b1 <- b2
  b2 <- get.b2(b0, b1, delta_data_b0, delta_data_b1)
  delta_r_b2 <- get.delta.r.bx(b2)
}

# 根据b2计算截距项a
a <- median(raw_data$y - b2*(raw_data$x - xM))

# 最终拟合图查看
par(opar)
plot(raw_data)
points(medians_x, medians_y, pch=16, col="green")

abline(a, b2, col='red')
abline(ols_fit_model, col='blue')
legend(-4,5, 
       legend = c('Robust', 'Ols'), 
       col = c('red', 'blue'),
       lty =1)


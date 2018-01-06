data <- read.table("/home/shensir/Documents/MyPrograming/R/EDA-Course-WIth-R/paper/brain.txt")
data <- cbind.data.frame(data$V3, data$V2)
colnames(data) <- c("BodyW", "BrainW")

X <- data$BodyW
y <- data$BrainW

X.log <- log(data$BodyW)
y.log <- log(data$BrainW)

opar <- par(no.readonly = TRUE)
# for raw data
# fivenum
fivenum(X)
fivenum(y)
# log-data
fivenum(X.log)
fivenum(y.log)

# stem plot
library(aplpack)
stem.leaf(X)
stem.leaf(y)
stem.leaf.backback(X, y)

# log-data
stem.leaf(X.log)
stem.leaf(y.log)
stem.leaf.backback(X.log, y.log)

# boxplot
par(mfrow=c(1,2))
boxplot(X, main="X")
boxplot(y, main="y")

# log-data
par(mfrow=c(1,2))
boxplot(X.log, main="X.log")
boxplot(y.log, main="y.log")

# data plot
plot(X,y)
plot(X.log, y.log)


# plot
plot(data)
par(opar)
plot(data, log = "xy")


# 正态性检验
# qqnorm

par(mfrow=c(1,2))
qqnorm(X, main = "Normal Q-Q Plot of X")
qqline(X)
qqnorm(y, main = "Normal Q-Q Plot of y")
qqline(y)
# log-data
par(mfrow=c(1,2))
qqnorm(X.log, main = "Normal Q-Q Plot of X.log")
qqline(X.log)
qqnorm(y.log, main = "Normal Q-Q Plot of y.log")
qqline(y.log)

# transformation
ratio.X = max(data$BodyW)/min(data$BodyW)
ratio.y = max(data$BrainW)/min(data$BrainW)


# letter value  想要寻求变为对称的幂指数
#======================================
library(lvplot)
library(dplyr)

get.index.bylv <- function(data){
  letter.values <- lvplot::lvtable(data,k=7)
  letter.values <- as.data.frame(letter.values)
  M <- letter.values$LV[7]  # median
  
  letter.values.L <- rev(letter.values$LV[c(1:6)])
  letter.values.U <- letter.values$LV[c(8:13)]
  letter.values.New <- cbind.data.frame(letter.values.L,
                                        letter.values.U)
  rownames(letter.values.New) <- c("F", "E", "D", "C", "B", "A")
  colnames(letter.values.New) <- c("Low","Up")
  # columns calculation
  
  letter.values.New <- mutate(letter.values.New, subM = (Low+Up)/2-M,
                              numDiv = ((Up-M)^2+(M-Low)^2)/(4*M),
                              pvalue = 1-subM/numDiv,
                              index = 1-pvalue)
  # get index
  index.proper <- median(letter.values.New$index)
  return(index.proper)
}

x.index.bylv <- get.index.bylv(X)  # 0.03 -> log
y.index.bylv <- get.index.bylv(y)  # 0.01 -> log
boxplot(X.log, y.log)

#======================================

# letter value plot
letter.values <- lvplot::lvtable(X ,k=7)

lvplot::LVboxplot(X,col=blues9, k=7)
lvplot::LVboxplot(y,col=blues9, k=7)

lvplot::LVboxplot(X.log, col=blues9, k=7)
lvplot::LVboxplot(y.log, col=blues9, k=7)



# 为直线性寻求变换的幂指数
#======================================
library(MASS)
get.e.index <- function(X, y){
  Xm <- median(X)
  ym <- median(y)
  y1 <- y-ym
  X1 <- X-Xm
  f1 <- lqs(y1 ~ X1, method = "lms")
  C <- as.numeric(f1$coefficients[2])
  
  X2 <- C^2*(X-Xm)^2/(2*ym)
  y2 <- y-ym-C*(X-Xm)
  f2 <- lqs(y2~X2, method = "lms")
  plot(X2, y2)
  abline(a=as.numeric(f2$coefficients[1]), b=as.numeric(f2$coefficients[2]))
  e.index <- as.numeric(f2$coefficients[2])  
}

# 调换X，y的位置，寻求使得X变换后与y成直线关系
x.e.index <- get.e.index(y.log, X)  # 0 -> log
#======================================



# 相关分析
#cor(X, y, method = "pearson")
cor(X.log, y.log, method = "pearson")

# regression

fit <- lm(y.log~X.log)
summary(fit)
plot(fit)

# 增强的Q-Q Plot
library(car)
qqPlot(fit, labels=row.names(data), id.method="identify",
       simulate=TRUE, main="Q-Q Plot")

# 增强的标准残差图
residplot <- function(fit, nbreaks=10) {
  z <- rstudent(fit)
  hist(z, breaks=nbreaks, freq=FALSE,
       xlab="Studentized Residual",
       main="Distribution of Errors")
  rug(jitter(z), col="brown")
  curve(dnorm(x, mean=mean(z), sd=sd(z)),
        add=TRUE, col="blue", lwd=2)
  lines(density(z)$x, density(z)$y,
        col="red", lwd=2, lty=2)
  legend("topright",
         legend = c( "Normal Curve", "Kernel Density Curve"),
         lty=1:2, col=c("blue","red"), cex=.7)
}
residplot(fit)
# D-W test
durbinWatsonTest(fit)

# Acf
acf(fit$residuals)

# HOMOSCEDASTICITY[方差齐性]
# Breusch-Pagan test
ncvTest(fit)



# Unusual observations
# outlier
outlierTest(fit)

# High-leverage points
hat.plot <- function(fit) {
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  plot(hatvalues(fit), main="Index Plot of Hat Values")
  abline(h=c(2,3)*p/n, col="red", lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}
# 注意这里在交互式的环境中是可以手动选择特定点的，这里在rmd可能无法正常显示
hat.plot(fit) 

# Influential observations
cutoff <- 4/(nrow(data)-length(fit$coefficients)-2)
plot(fit, which=4, cook.levels=cutoff)
abline(h=cutoff, lty=2, col="red")

# combine the information from outlier, leverage, and influence plots
library(car)

influencePlot(fit, id.method="identify", 
              main="Influence Plot",
              sub="Circle size is proportional to Cook's distance")


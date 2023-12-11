## -----------------------------------------------------------------------------
## 生成时间序列数据
X <- matrix(0,1000,1)
Y <- matrix(0,1000,1)
for(i in 1:1000){   

X[i+1] <- X[i]+rnorm(1)
}

for(i in 1:1000){   

Y[i+1] <- Y[i]+rnorm(1)+0.3
}

plot(X,col=2,type="l")
lines(Y,col=3,type="l")
legend("topright",pch=c(15,15),legend=c("δ=0","δ=0.3"),col=c(2,3),bty="n")


## -----------------------------------------------------------------------------
my.sample <- function(x,size,replace =TRUE ,prob= NULL){
  ## 计算累计概率分布函数,记作cp
   x.size<-length(x)
   if(is.null(prob) )  cp <- c(seq(1/x.size,1,length =x.size))
    else cp<- cumsum(prob)
   ## 使用逆变换法 
    u<-runif(size)
    i<- findInterval(u,cp)
    return( x[i+1])
}
## 测试集
    set.seed(2023)
    my.sample(letters,10)
    my.sample(c(1:10),10)


## ----pressure-----------------------------------------------------------------
#定义抽样数
N<- 1000
#抽样
U<- runif(N)
library(VGAM)
x<-qlaplace(U)


# 结果包括积分上限的估计值
hist(x,probability= TRUE )
plot(density(x),col= "blue",main="density f(x)")
f <- function(x) 0.5*exp(-(abs(x))) 
lap.density<- f(seq(-6,6,length=N))
lines(seq(-6,6,length=N),lap.density,col= "red")
legend("topright",legend=c("simulated","True "),col=c("blue","red"),lty=c(1,1))

## -----------------------------------------------------------------------------
n <- 1e4;j<-k<-0;y <- numeric(n)
while (k < n) {
u <- runif(1)
j <- j + 1
x <- runif(1) #random variate from g(.)
if (x^2 * (1-x) > u) {
#we accept x
k <- k + 1
y[k] <- x
}
}
hist(y,probability= TRUE )
g <- function(x)
   12*x^2*(1-x) 
lap.density<- g(seq(0,1,length=n))
lines(seq(0,1,length=n),lap.density,col= "red",main="density")
legend("topright",legend=c("density"),col=c("red"),lty=c(1))


## -----------------------------------------------------------------------------
N<- 30000
u<-matrix(runif(N,-1,1),nrow = N/3,ncol = 3)
x<-numeric(length = N/3)
for (i in 1:(N/3)) {
  if(abs(u[i,3])>=max(abs(u[i,]))) x[i]<- u[i,2]
  else x[i]<- u[i,3]
}
hist(x,prob= TRUE)
lines(density(x),col= "red",main="density")
legend("topright",legend=c("density"),col=c("red"),lty=c(1))

## -----------------------------------------------------------------------------
## \rho =1 
set.seed(1006)
pihat1<- numeric(100)
pihat2<- numeric(100)
pihat3<- numeric(100)
rho <- 1
d <- 1
N <- 1e6
for (i in 1:100) {

X <- runif(N,0,d/2)
Y <- runif(N,0,pi/2)
pihat1[i] <- 2*rho/mean(1/2*sin(Y)>X)
}
rho <- 0.5
for (i in 1:100) {
X <- runif(N,0,d/2)
Y <- runif(N,0,pi/2)
pihat2[i] <- 2*rho/mean(1/2*sin(Y)>X)
}
rho <- 0.25
for (i in 1:100) {
X <- runif(N,0,d/2)
Y <- runif(N,0,pi/2)
pihat3[i] <- 2*rho/mean(1/2*sin(Y)>X)
}
## output the variance respectively，rho 1,0.5,0.25.
var(pihat1)
var(pihat2)
var(pihat3)

## -----------------------------------------------------------------------------
N <- 1e6

thetamc<- numeric(100)
thetaant <- numeric(100)
  set.seed(1006)
for ( i  in 1:100) {
  u<- runif(N)
  v<- u[1:(N/2)]
  thetamc[i]<-1/N*mean(exp(u))
  thetaant[i]<- 2/N * mean(exp(v)+exp(1-v))
}
 1-var(thetaant)/var(thetamc)

## -----------------------------------------------------------------------------
set.seed(20231015)
g <- function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>=1)
}
x = seq(1,8,0.01)
gs <- c(expression(g(x)),expression(f1(x)),expression(f2(x)))
par(mfrow=c(1,2))
# figure of g, f1, f2
plot(x, g(x), type="l", ylab="", ylim=c(0,0.7), lwd = 1, col=1)
lines(x, dnorm(x), lwd=1, col=2)
lines(x, dgamma(x,3,2), lwd=1, col=3)
legend("topright", legend = gs, lty=1, lwd=1, inset = 0.02,col=1:3)

# figure of g/f1, g/f2
plot(x, g(x)/dnorm(x), type="l", ylab="", ylim=c(0,5), lwd = 1, col=2)
lines(x, g(x)/dgamma(x,3,2), lwd=1, col=3)
legend("topright", legend = gs[-1], lty=1, lwd=1, inset = 0.02,col=2:3)

## -----------------------------------------------------------------------------
set.seed(20231015)
m = 1e4
theta <- se <- numeric(2)
# using f1
x <- rnorm(m) 
fg <- g(x) / dnorm(x)
theta[1] <- mean(fg)
se[1] <- sd(fg)
# using f2
x <- rgamma(m,3,2) 
fg <- g(x) / dgamma(x,3,2)
theta[2] <- mean(fg)
se[2] <- sd(fg)
rbind(theta, se)

## -----------------------------------------------------------------------------
se^2

## -----------------------------------------------------------------------------
alpha <- 3  
beta <- 1   

# Number of samples
N <- 100000

sample <- rgamma(N, shape = alpha, rate = beta)
est <- (mean((sample^2 * gamma(alpha)) / beta^alpha))

# Print the estimate
cat("Monte Carlo Estimate:", est, "\n")



## -----------------------------------------------------------------------------
set.seed(20231015)
m <- 1e6
estimation <- sd <- 0
g <- function(x){
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
estimation <- mean(fg)
sd <- sd(fg)

## -----------------------------------------------------------------------------
M <- 10000
k <- 5 
m <- M/k # replicates per stratum
T1 <- T2 <- numeric(k)
estimation_2 <- numeric(2)
fj <- function(x)5*exp(-x)/(1-exp(-1)) # fj(x)
g <- function(x)exp(-x)/(1+x^2)

## -----------------------------------------------------------------------------
set.seed(20231015)
for(j in 1:k){
  u = runif(m)
  x = -log(1-(1-exp(-1))*(u+j-1)/5) # inverse transform method
  T1[j] = mean(g(x)/fj(x))
  T2[j] = var(g(x)/fj(x))
}
estimation_2[1] = sum(T1) 
estimation_2[2] = sum(T2)
round(c(estimation,estimation_2[1]),4)
round(c(sd,sqrt(estimation_2[2])),5)

## -----------------------------------------------------------------------------
cl <- 0.95
s <- 20
n <- 10000
true_mean <- 2  

suc_intervals <- 0

for (i in 1:n) {
  x <- rchisq(s, df = 2)
  
  sample_mean <- mean(x)
  sample_sd <- sd(x)
  
  moe <- qt((1 + cl) / 2, df = s - 1) * (sample_sd / sqrt(s))
  interval <- c(sample_mean - moe, sample_mean + moe)
  
  if (true_mean >= interval[1] && true_mean <= interval[2]) {
    suc_intervals <- suc_intervals + 1
  }
}

c_p <- suc_intervals / n

cat("Estimated Coverage Probability:", c_p, "\n")




## -----------------------------------------------------------------------------
n <- 10000
alpha <- 0.05

true_mean_chi_sq <- 1
true_mean_uniform <- 1
true_mean_exponential <- 1

rej_chi_sq <- 0
rej_uniform <- 0
rej_exponential <- 0

# Perform Monte Carlo simulations for each case
for (i in 1:n) {
  
  sample_chi_sq <- rchisq(30, df = 1)  
  sample_uniform <- runif(30, min = 0, max = 2) 
  sample_exponential <- rexp(30, rate = 1)  
  
  t_test_chi_sq <- t.test(sample_chi_sq, mu = true_mean_chi_sq, alternative = "two.sided")
  t_test_uniform <- t.test(sample_uniform, mu = true_mean_uniform, alternative = "two.sided")
  t_test_exponential <- t.test(sample_exponential, mu = true_mean_exponential, alternative = "two.sided")
  
  if (t_test_chi_sq$p.value < alpha) {
    rej_chi_sq <- rej_chi_sq + 1
  }
  
  if (t_test_uniform$p.value < alpha) {
    rej_uniform <- rej_uniform + 1
  }
  
  if (t_test_exponential$p.value < alpha) {
    rej_exponential <- rej_exponential + 1
  }
}

t1e_chi_sq <- rej_chi_sq / n
t1e_uniform <- rej_uniform / n
t1e_exponential <- rej_exponential / n

# Print the results
cat("Type I Error Rate for χ²(1) Data:", t1e_chi_sq, "\n")
cat("Type I Error Rate for Uniform(0,2) Data:", t1e_uniform, "\n")
cat("Type I Error Rate for Exponential(1) Data:", t1e_exponential, "\n")


## -----------------------------------------------------------------------------
   set.seed(20231022)
    library(boot); library(MASS)
    lambda <- function(x,i) 1/mean(x[i])
    n <- 20
    x <- rexp(n,2)
    obj <- boot(data=x,statistic=lambda,R=1e4)
    round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
            se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
    n <- 10
    x <- rexp(n,2)
    obj <- boot(data=x,statistic=lambda,R=1e4)
    round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
            se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
    n <- 5
    x <- rexp(n,2)
    obj <- boot(data=x,statistic=lambda,R=1e4)
    round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
            se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot)
law1<-as.matrix(law)

boot.t.ci <-
function(x, B = 500, R = 100, level = .95, statistic){
#compute the bootstrap t CI
x <- as.matrix(x); n <- nrow(x)
stat <- numeric(B); se <- numeric(B)
boot.se <- function(x, R, f) {
#local function to compute the bootstrap
#estimate of standard error for statistic f(x)
x <- as.matrix(x); m <- nrow(x)
th <- replicate(R, expr = {
i <- sample(1:m, size = m, replace = TRUE)
f(x[i, ])
})
return(sd(th))
}
for (b in 1:B) {
j <- sample(1:n, size = n, replace = TRUE)
y <- x[j, ]
stat[b] <- statistic(y)
se[b] <- boot.se(y, R = R, f = statistic)
}
stat0 <- statistic(x)
t.stats <- (stat - stat0) / se
se0 <- sd(stat)
alpha <- 1 - level
Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
names(Qt) <- rev(names(Qt))
CI <- rev(stat0 - Qt * se0)
}

b.cor <- function(x,i) cor(x[i,1],x[i,2])

ci <- boot.t.ci(law1, statistic = b.cor, B=2000, R=200)
print(ci)

## -----------------------------------------------------------------------------
m<- 1000
set.seed(20231028)
library(boot)
b.mean<- function(x,i) mean(x[i])
data<-as.matrix(aircondit) 
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)

  bootout<- boot(data= data,statistic = b.mean ,R =1000)
  ci <- boot.ci(bootout,type=c("norm","basic","perc","bca"))
  print(ci)


## -----------------------------------------------------------------------------
library(bootstrap)
m<-nrow(scor)
Cov<-cov(scor)
theta=max(eigen(Cov)$values)/sum(diag(Cov))
thetajack<- numeric(m)
for (i in 1:m){
  cov<- cov(scor[-i,])
  sum<-sum(diag(cov))
  thetajack[i]<- max(eigen(cov)$values)/sum
}
jackmean<- mean(thetajack)
bias <- (m-1)*(jackmean-theta)
variance<- m/(m-1)* sum((thetajack-jackmean)^2)
print(c(bias,variance))

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
#in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- matrix(0,nrow = n*(n-1),ncol=2)
# for n-fold cross validation
# fit models on leave-two-out samples
for (k in 1:(n-1)) {
  for (l in 1:(n-1)) {
    
  
y1 <- magnetic[-k]
x1 <- chemical[-k]
y <- y1[-l]
x <- x1[-l]
J1 <- lm(y ~ x)
yhat11 <- J1$coef[1] + J1$coef[2] * chemical[k]
yhat12 <- J1$coef[1] + J1$coef[2] * x1[l]
e1[(k-1)*n+l,1] <- magnetic[k] - yhat11
e1[(k-1)*n+l,2] <- y1[l] - yhat12
J2 <- lm(y ~ x + I(x^2))
yhat21 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
yhat22 <- J2$coef[1] + J2$coef[2] * x1[l] +
J2$coef[3] * x1[l]^2
e2[(k-1)*n+l,1] <- magnetic[k] - yhat21
e2[(k-1)*n+l,2] <- y1[l] - yhat22

J3 <- lm(log(y) ~ x)
logyhat31 <- J3$coef[1] + J3$coef[2] * chemical[k]
logyhat32 <- J3$coef[1] + J3$coef[2] * x1[l]
yhat31 <- exp(logyhat31)
yhat32 <- exp(logyhat32)
e3[(k-1)*n+l,1] <- magnetic[k] - yhat31
e3[(k-1)*n+l,2] <- y1[l] - yhat32

J4 <- lm(log(y) ~ log(x))
logyhat41 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
logyhat42 <- J4$coef[1] + J4$coef[2] * log(x1[l])
yhat41 <- exp(logyhat41)
yhat42 <- exp(logyhat42)
e4[(k-1)*n+l,1] <- magnetic[k] - yhat41
e4[(k-1)*n+l,2] <- y1[l] - yhat42
}
}
c(mean(colMeans(e1^2)), mean(colMeans(e2^2)), mean(colMeans(e3^2)), mean(colMeans(e4^2)))


## -----------------------------------------------------------------------------
summary( lm(magnetic~ chemical+I(chemical^2)))

## -----------------------------------------------------------------------------
library(CDFt)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:26
D <- numeric(R) #storage for replicates
options(warn = -1)
D0 <- CramerVonMisesTwoSamples(x, y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = 14, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
D[i] <- CramerVonMisesTwoSamples(x1, y1)
}
p <- mean(c(D0, D) >= D0)
options(warn = 0)
p

## -----------------------------------------------------------------------------
set.seed(20231105)
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
R <- 999 #number of replicates
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
X <- x - mean(x) #centered by sample mean
Y <- y - mean(y)
z<- c(X,Y)
K <- 1:50
D <- numeric(R) #storage for replicates
options(warn = -1)
D0 <- count5test(X, Y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = 20, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
D[i] <- count5test(x1, y1)
}
p <- ifelse(D0 == 0, mean( c(D0, D) > D0 ), mean( c(D0, D) >= D0 ))
options(warn = 0)
p


## -----------------------------------------------------------------------------
set.seed(20231112)
my.root<- function (N,b1,b2,b3,f0)    {
x1 <- rpois(N,1);x2 <- rexp(N,1); x3 <- sample(0:1,N,replace=TRUE)
g <- function(alpha){
tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3); p <- 1/(1+tmp)
mean(p) - f0
}
solution <- uniroot(g,c(-1000,1000))
round(unlist(solution),5)[1]   }
my.root(N=1e6,b1=0,b2=1,b3=-1,f0=0.1)
my.root(N=1e6,b1=0,b2=1,b3=-1,f0=0.01)
my.root(N=1e6,b1=0,b2=1,b3=-1,f0=0.001)
my.root(N=1e6,b1=0,b2=1,b3=-1,f0=0.0001)
f<-seq(0,1,length=100)
vecroot<- Vectorize(my.root)
a<-vecroot(N=1e6,b1=0,b2=1,b3=-1,f0=f)
plot(-log(f),a,type="l")

## -----------------------------------------------------------------------------
set.seed(20231112)
 lapla <- function(x) 1/2*exp(-abs(x))
 
 rw.Metropolis <- function(sigma, x0, N) {
   x <- numeric(N)
   x[1] <- x0
   u <- runif(N)
   k <- 0
   for (i in 2:N) {
     y <- rnorm(1, x[i-1], sigma)
     if (u[i] <= (lapla(y) / lapla(x[i-1])))
       x[i] <- y else {
         x[i] <- x[i-1]
         k <- k + 1
       }
   }
   return(list(x=x, k=k/N))
 }
 
 N <- 2000
 sigma <- c(.05, .5,1, 2, 5,15)
 x0 <- 0
 rw1 <- rw.Metropolis(sigma[1], x0, N)
 rw2 <- rw.Metropolis(sigma[2], x0, N)
 rw3 <- rw.Metropolis(sigma[3], x0, N)
 rw4 <- rw.Metropolis(sigma[4], x0, N)
 rw5 <- rw.Metropolis(sigma[5], x0, N)
 rw6 <- rw.Metropolis(sigma[6], x0, N)
 plot(rw1$x,type="l")
 plot(rw2$x,type="l")
 plot(rw3$x,type="l")
 plot(rw4$x,type="l")
 plot(rw5$x,type="l")
 plot(rw6$x,type="l")
print(c(rw1$k, rw2$k, rw3$k, rw4$k,rw5$k,rw6$k))

## -----------------------------------------------------------------------------
 N<-2000
 
 
   x<-y <- numeric(N)
   x[1]<-y[1] <-0 
   u <- runif(N)
   k <- 0
   for (i in 2:N) {
     x0<- rnorm(1,0.9*y[i-1],sqrt(0.19))
     y0 <- rnorm(1,0.9*x[i-1], sqrt(0.19))
     
     x[i] <- x0 
     y[i] <- y0
   
   }
 
 li<- lm(y~x)
 plot(li)

## -----------------------------------------------------------------------------
set.seed(20231112)
## 函数
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}
f <- function(x, sigma) {
if (any(x < 0)) return (0)
stopifnot(sigma > 0)
return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}
Rai.chain<-function(sigma, m, x1) {
sigma <- 4
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rchisq(1, df = xt)
num <- f(y, sigma) * dchisq(xt, df = y)
den <- f(xt, sigma) * dchisq(y, df = xt)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1 #y is rejected
}
}
x
}
sigma <- .2 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 20000#length of chains
b <- 10000 #burn-in length
#choose overdispersed initial values
x0 <- c(-10, -5,5, 10)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- Rai.chain(sigma, n, x0[i])
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",
xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R",ylim=c(1,3))
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
library(boot)
#enter the payoff matrix
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
B<- A+2
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}
sB <- solve.game(B)
round(cbind(sB$x, sB$y), 7)
C<- A+B
sC <- solve.game(C)
round(cbind(sC$x, sC$y), 7)


## -----------------------------------------------------------------------------
library(microbenchmark)
library(Rcpp)
 gibbs_r <- function(N, thin , n,a,b) {
    mat <- matrix(nrow = N, ncol = 2)
    x <- y <- 0
    for (i in 1:N) {
      for (j in 1:thin) {
        x <- rbinom(1, n, y)
        y <- rbeta(1,x+a,n-x+b)
      }
      mat[i, ] <- c(x, y)
    }
    mat
 }
 dir_cpp <- 'C:/Users/11146/Desktop/studying/Rlab/statcomp/hw9/'
 sourceCpp(paste0(dir_cpp,"gibbsC.cpp"))
 ts <- microbenchmark(
 gibbs_r(100, 10 , 10 , 5,5),
 gibbs_cpp(100, 10 , 10,5,5)
 )
 summary(ts)


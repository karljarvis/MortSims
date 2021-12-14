par(mfrow=c(1,2))
# Option 2: inverse square
x <- seq(1,5,0.1)
y <- 1/(x^2)
plot(y, main = "Inverse Square", type = "l")

# negative exponential movement function (option 5 in CDPOP)
a <- 1
b <- 1
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, main = "Negative Exponential", type = "l")




# negative exponent comparing effects of parameters a and b
par(mfrow=c(3,4))
a <- 1
b <- 1
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=1 b=1')

a <- 2
b <- 1
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=2 b=1')

a <- 4
b <- 1
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=4 b=1')

a <- 8
b <- 1
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=8 b=1')

a <- 1
b <- 2
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=1 b=2')

a <- 1
b <- 4
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=1 b=4')

a <- 1
b <- 8
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=1 b=8')

a <- 2
b <- 2
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=2 b=2')

a <- 4
b <- 2
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=4 b=2')

a <- 8
b <- 2
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=8 b=2')

a <- 2
b <- 4
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=2 b=4')

a <- 2
b <- 8
x <- seq(1,5, 0.1)
y <- a*10^(-b*x)
plot(y, ylim=c(0,0.5), xlim=c(0,20), type = "l", main = 'a=2 b=8')


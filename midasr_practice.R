
library(midasr)
x = 1:12
fmls(x, k = 2, m = 3)
mls(x, 0:2, m = 3)
mls(x, 2:4, m = 3)

set.seed(1001)
n = 250   ## number of low frequency observations
trend = 1:n
x = rnorm(4*n)
z = rnorm(12*n)

## nealmon: exponential weight
fn_x <- nealmon(p = c(1, -0.5), d = 8)     ## coefficients for x
fn_z <- nealmon(p = c(2, 0.5, -0.1), d = 17)  ## coefficients for z
y <- 2 + 0.1 * trend + mls(x, 0:7, 4) %*% fn_x + mls(z, 0:16, 12) %*% fn_z + rnorm(n)

## unrestricted model, the following two models are equivalent
eq_u = lm(y~trend + mls(x, k = 0:7, m = 4) + mls(z, k = 0:16, m = 12))
eq_u = midas_r(y~trend + mls(x, k = 0:7, m = 4) + mls(z, k = 0:16, m = 12), start = NULL)

## restricted model
eq_r = midas_r(y~trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16,12,nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
# try a different start value
# eq_r = midas_r(y~trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16,12,nealmon), start = list(x = c(0.1, 0.3), z = c(0.1, 0.2, -0.3)))
summary(eq_r)
coef(eq_r)
coef(eq_r,midas=TRUE)
plot_midas_coef(eq_r)

eq_r = midas_r(y~trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16,12,nealmon), start = list(x = c(1, -0.5, 0.6), z = c(2, 0.5, -0.1)))
summary(eq_r)
coef(eq_r)

## use optim 
eq_r2 = midas_r(y~trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16,12,nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
                Ofunction = "optim", method = "Nelder-Mead")
summary(eq_r2)
coef(eq_r2)

eq_r2 = midas_r(y~trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16,12,nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
                Ofunction = "nls")
#summary(eq_r2)
#coef(eq_r2)

eq_r2 = update(eq_r2, Ofunction = "nls")

eq_r2 = midas_r(y~trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16,12,nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
                Ofunction = "optim", method = "Nelder-Mead")
eq_r2$opt
eq_r2$convergence   ## 1 is not convergent, 0 is convergent
eq_r2$coefficients

deriv_tests(eq_r, tol = 1e-06)
coef(eq_r,midas=TRUE)

## adequacy testing of restrictions
## H0: constraint is adequate
## H1: constraint is not adequate
hAh_test(eq_r)
hAhr_test(eq_r)
hAh_test(eq_r2)
hAhr_test(eq_r2)

## forecasting
newx = rnorm(4)
newz = rnorm(12)
forecast(eq_r2, newdata = list(x = newx, z = newz, trend = 251))

eq_f = midas_r(y ~ trend + mls(x, 4+0:7, 4, nealmon) + mls(z, 12 + 0:16, 12, nealmon),
               start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)) )
summary(eq_f)
hAh_test(eq_f)
hAhr_test(eq_f)
forecast(eq_f, newdata = list(x = rep(NA, 4), z = rep(NA, 12), trend = 251))

cbfc = select_and_forecast(y~trend+mls(x, 0, 4) + mls(z, 0, 12),
    from = list(x = c(4,8,12), z = c(12,24,36)),
    to = list(x = rbind(c(14,19), c(18,23), c(22,27)),
              z = rbind(c(22,27), c(34,39), c(46,51))),
    insample = 1:200, outsample = 201:250,
    weights = list(x = c("nealmon", "almonp"), z = c("nealmon", "almonp")),
    wstart = list(nealmon = rep(1,3), almonp = rep(1,3)),
    IC = "AIC", seltype = "restricted", ftype = "fixed",
    measures = c("MSE", "MAPE", "MASE"),
    fweights = c("EW", "BICW", "MSFE", "DMSFE")
    )
cbfc$accuracy$individual   ## best forecasting equations
cbfc$accuracy$average  ## the out-of-sample precision of forecasting combinations for each forecasting horizon

mod1 = midas_r(y ~ trend + mls(x, 4:14, 4, nealmon) + mls(z, 12:22, 12, nealmon),
               start = list(x = c(10,1,-0.1), z = c(2, -0.1)))
hAh_test(mod1)
hAhr_test(mod1)
avgf = average_forecast(list(mod1), data = list(y = y, x = x, z = z, trend = trend),
                        insample = 1:200, outsample = 201:250, type = "fixed", 
                        measures = c("MSE", "MAPE", "MASE"),
                        fweights = c("EW", "BICW", "MSFE", "DMSFE"))
## avgf$x - y[1:200,1]


data("USqgdp", package = "midasr")
data("USpayems", package = "midasr")
y = window(USqgdp, end = c(2011,2))
x = window(USpayems, end = c(2011,7))

yg = diff(log(y))*100
xg = diff(log(x))*100
nx = ts(c(NA,xg,NA,NA), start = start(x), frequency = 12)
ny = ts(c(rep(NA,33), yg, NA), start = start(x),frequency = 4)

xx = window(nx, start = c(1985,1), end = c(2009,3))
yy = window(ny, start = c(1985,1), end = c(2009,1))
beta0 = midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbeta), start = list(xx = c(1.7,1,5)))
hAh_test(beta0)
hAhr_test(beta0)
betan = midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbetaMT), start = list(xx = c(2, 1, 5, 0)))
hAh_test(betan)
hAhr_test(betan)

fulldata = list(xx = window(nx, start = c(1985,1), end = c(2011,6)),
                yy = window(ny, start = c(1985,1), end = c(2011,2)))
insample = 1:length(yy)
outsample = (1:length(fulldata$yy))[-insample]
vgf <- average_forecast(list(beta0, betan), data = fulldata, insample = insample, outsample = outsample)

data("rvsp500", package = "midasr")

harstep <- function(p, d, m) {
    if (d != 20)
        stop("HAR(3)-RV process requires 20 lags")
    out <- rep(0, 20)
    out[1] <- p[1] + p[2] / 5 + p[3] / 20
    out[2:5] <- p[2] / 5 + p[3] / 20
    out[6:20] <- p[3] / 20
    out 
}
spx2_rvol <- 100 * sqrt(252 * as.numeric(rvsp500[, "SPX2.rv"]))
mh = midas_r(rv ~ mls(rv, 1:20, 1, harstep), data = list(rv = spx2_rvol), start = list(rv = c(1, 1, 1)))
summary(mh)
hAh_test(mh)
hAhr_test(mh)

mr <- midas_r(rv ~ mls(rv, 1:20, 1, nealmon), data = list(rv = spx2_rvol),
              start = list(rv = c(0, 0, 0)), weight_gradients = list())
summary(mr)
hAhr_test(mr)





coef(eq_r)
coef(eq_r,midas=TRUE)
gamma0 = 1.3533616
gamma1 = -0.5074774

s = 0:7
gamma0*exp(gamma1* s)/sum(exp(gamma1*c(0:7)))
nealmon(p = c(gamma0, gamma1), d = 8)


eq_r2 = midas_r(y~trend + mls(x, 0:7, 4, almonp) + mls(z, 0:16,12,nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
coef(eq_r2)
coef(eq_r2, midas = TRUE)

gamma0 = 0.42776283
gamma1 = -0.05490551
s = 1:8
gamma0+gamma1*s

almonp(p = c(gamma0, gamma1), d = 8)
i = 1:8
poly(i, degree = 1, raw=TRUE)%*%gamma1+gamma0



AIC(eq_r)
BIC(eq_r2)






length(y)
length(x)
head(x,20)

eq_r = midas_r(y~ mls(x, 0:3, 4, nealmon), start = list(x=c(1,1)) )
coef(eq_r)
coef(eq_r, midas = TRUE)
head(eq_r$fitted.values)
a = mls(x, 0:3, 4, nealmon)[1,]
sum(a* (coef(eq_r, midas = TRUE)[2:5]) )+coef(eq_r, midas = TRUE)[1]
head(eq_r$fitted.values)

sum(rev(a)*(coef(eq_r, midas = TRUE)[2:5]))+coef(eq_r, midas = TRUE)[1]


a = mls(x, 0:3, 4, nealmon)[50,]
sum(a* (coef(eq_r, midas = TRUE)[2:5]) )+coef(eq_r, midas = TRUE)[1]
head(eq_r$fitted.values)

-0.6229437*0.6589577+(1.0915017)*0.5717220+(-0.1435595)*0.4960349+(-0.5573113)*0.4303677+14.5100010





head(mls(x, 0:2, 4, nealmon))
eq_r = midas_r(y~ mls(x, 0:2, 4, nealmon), start = list(x=c(1,1)) )
coef(eq_r)
coef(eq_r, midas = TRUE)
head(eq_r$fitted.values)
head(y)
head(x,20)

-0.6229437*0.5060015+1.0915017*0.5195755+(-0.1435595)*0.5335137+14.5427610








library(midasr)
x = 1:12

set.seed(1001)
n = 250   ## number of low frequency observations
trend = 1:n
x = rnorm(4*n)
z = rnorm(4*n)
w = rnorm(4*n)

## nealmon: exponential weight
fn_x <- nealmon(p = c(1, -0.5), d = 4)     ## coefficients for x
fn_z <- almonp(p = c(2, 0.5), d = 4)  ## coefficients for z
fn_w = nbeta(p = c(0.1, 0.5, 0.4), d = 4)  
y <- 2 + mls(x, 0:3, 4) %*% fn_x + mls(z, 0:3, 4) %*% fn_z + mls(w, 0:3, 4) %*% fn_w + rnorm(n)

## restricted model
eq_r = midas_r(y~trend + mls(x, 0:3, 4, nealmon) + mls(z, 0:3, 4, almonp) + mls(w, 0:2, 4, nbeta), 
               start = list(x = c(1, 1), z = c(1,1), w = c(1,1,1) ), Ofunction = "optim")
summary(eq_r)


for(i in 1:1000){
    set.seed(i+6000)
    message("seed: ", i)
    xst = runif(2)
    zst = runif(2)
    wst = runif(3)
    eq_r = midas_r(y~trend + mls(x, 0:3, 4, nealmon) + mls(z, 0:3, 4, almonp) + mls(w, 0:3, 4, nbeta), 
                   start = list(x = xst, z = zst, w = wst ), Ofunction = "optim")
    #summary(eq_r)
    hAh_test(eq_r)
    hAhr_test(eq_r)
    message("----------------------------------------------------------------------------------------")
}



set.seed(5)
xst = runif(2)
zst = runif(2)
wst = runif(3)
eq_r = midas_r(y~trend + mls(x, 0:3, 4, nealmon) + mls(z, 0:3, 4, almonp) + mls(w, 0:2, 4, nbeta), 
               start = list(x = xst, z = zst, w = wst ), Ofunction = "optim")
summary(eq_r)
hAh_test(eq_r)
hAhr_test(eq_r)





for(i in 1:1000){
    set.seed(i)
    message("seed: ", i)
    xst = runif(2)
    zst = runif(2)
    wst = runif(3)
    eq_r = midas_r(y~trend + mls(x, 0:3, 4, nealmon) + mls(z, 0:3, 4, almonp), 
                   start = list(x = xst, z = zst), Ofunction = "optim")
    summary(eq_r)
    hAh_test(eq_r)
    hAhr_test(eq_r)
    message("----------------------------------------------------------------------------------------")
}






R2 <- function(z) {
    r <- z$residuals
    f <- z$fitted.values
    mss <- if (attr(z$terms, "intercept")) 
        sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
    ans <- list(r.squared=NULL,adj.r.squared=NULL)
    
    n <- length(r)
    p <- length(coef(z))
    rdf <- n-p
    df.int <- if (attr(z$terms, "intercept")) 1L
    else 0L
    
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
    
    if(!is.null(z$unrestricted)) { ansu <- R2unres(z)
    }
    else {
        ansu <- NULL
    }
    
    out <- list(restricted=ans,unrestricted=ansu)
    class(out) <- "R2_midas_r"
    out
}

print.R2_midas_r <- function(x,digits = max(3L, getOption("digits") - 3L)) {
    cat("Multiple R-squared: ", formatC(x$restricted$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$restricted$adj.r.squared, digits = digits),"\n")
    if(!is.null(x$unrestricted)) {
        cat("MIDAS-U:\n")
        cat("Multiple R-squared: ", formatC(x$unrestricted$r.squared, digits = digits))
        cat(",\tAdjusted R-squared: ", formatC(x$unrestricted$adj.r.squared, digits = digits),"\n")
    }
}


R2unres <- function(x) {
    z <- x$unrestricted
    r <- z$residuals
    f <- z$fitted.values
    mss <- if ("`(Intercept)`" %in% names(coef(z))) 
        sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
    ans <- list(r.squared=NULL,adj.r.squared=NULL)
    
    n <- length(r)   
    rdf <- z$df.residual
    df.int <- if ("`(Intercept)`" %in% names(coef(z))) 1L
    else 0L
    
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
    ans
}


R2(eq_r)



fit = lm(y~x[1:250])
plot(fit)

fit$residuals/sd(fit$residuals)
summary(abs(rstandard(fit)-fit$residuals/sd(fit$residuals)))

library(ggplot2)
df = data.frame(fitted = eq_r$fitted.values, resid = eq_r$residuals, 
                stdresid = eq_r$residuals/sd(eq_r$residuals))

p1<-ggplot(df, aes(fitted, resid))+geom_point()
p1<-p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
p1<-p1+xlab("Fitted values")+ylab("Residuals")
p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_bw()+theme(plot.title = element_text(hjust=0.5))
p1

p2<-ggplot(df, aes(qqnorm(stdresid)[[1]], stdresid))+geom_point(na.rm = TRUE)
p2<-p2+geom_abline(slope=1,intercept=0)+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
p2<-p2+ggtitle("Normal Q-Q")+theme_bw()+theme(plot.title = element_text(hjust=0.5))
p2


p3<-ggplot(df, aes(fitted, sqrt(abs(stdresid))))+geom_point(na.rm=TRUE)
p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
p3<-p3+ggtitle("Scale-Location")+theme_bw()+theme(plot.title = element_text(hjust=0.5))
p3


library(ggplot2)
df = data.frame(LGD = runif(60), label = rep(c("historical", "predicted"), each=30), year = rep(1987:2016,2))
p = ggplot(data=df, aes(x=year, y=LGD, group=label)) + geom_line(aes(color=label))+
    xlim(c(1987, 2022))   #+scale_color_manual(values=c("maroon2", "cyan3"))
dfb = data.frame(year = 2017:2021, LGD = runif(5), label = rep("base", 5))
dfm = data.frame(year = 2017:2021, LGD = runif(5), label = rep("adverse", 5))
dfs = data.frame(year = 2017:2021, LGD = runif(5), label = rep("severe", 5))
p = p + geom_line(data = dfb, aes(x = year, y = LGD), color = "red")
p = p + geom_line(data = dfm, aes(x = year, y = LGD), color = "blue")
p = p + geom_line(data = dfs, aes(x = year, y = LGD), color = "yellow")
p = p + theme_bw()
p = p + scale_color_manual(values = c("maroon2", "cyan3", "red", "blue", "yellow"),
                           labels = c("true DD", "predicted DD", "base", "adverse", "severe"),
                           guide = guide_legend(override.aes = list(alpha = 1, size = 3)))
p




df = data.frame(LGD = runif(60), label = rep(c("historical", "predicted"), each=30), year = rep(1987:2016,2))
df1 = data.frame(LGD = runif(30), year = 1987:2016)
p = ggplot()+xlim(c(1987, 2022))
p = p + geom_line(data = df1, aes(x = year, y = LGD), color = "maroon2")
df2 = data.frame(LGD = runif(30), year = 1987:2016)
p = p + geom_line(data = df2, aes(x = year, y = LGD), color = "cyan3")
dfb = data.frame(year = 2017:2021, LGD = runif(5), label = rep("base", 5))
dfm = data.frame(year = 2017:2021, LGD = runif(5), label = rep("adverse", 5))
dfs = data.frame(year = 2017:2021, LGD = runif(5), label = rep("severe", 5))
p = p + geom_line(data = dfb, aes(x = year, y = LGD), color = "red")
p = p + geom_line(data = dfm, aes(x = year, y = LGD), color = "blue")
p = p + geom_line(data = dfs, aes(x = year, y = LGD), color = "yellow")
p = p + theme_bw()
p = p + scale_color_manual(values = c("maroon2", "cyan3", "red", "blue", "yellow"),
                           labels = c("true DD", "predicted DD", "base", "adverse", "severe"))
p




historical = runif(30)
pred_hist = runif(30)
pred_base = runif(5)
pred_adverse = runif(5)
pred_severe = runif(5)

df = data.frame(LGD = c(historical, pred_hist, pred_base, pred_adverse, pred_severe), 
                label = c(rep("historical", 30), rep("predicted", 30), rep("base", 5), rep("adverse", 5), rep("severe", 5)),
                year = c(1987:2016,1987:2016,2017:2021,2017:2021,2017:2021),
                type = c(rep("dashed",30), rep("solid",45)))
p = ggplot(data=df, aes(x=year, y=LGD, group=label)) + geom_line(aes(color=label, linetype = type))+
    xlim(c(1987, 2022))+scale_colour_discrete("guide") +
    scale_color_manual(name = "guide", values=c("maroon2", "cyan3", "red", "blue", "yellow"))+
    scale_linetype_manual(name="guide",values= c('dashed', 'solid', 'solid', "solid","solid"))+
    theme_bw()+theme(legend.title=element_blank())
p


#########################################################################################################
historical = runif(30)
pred_hist = runif(30)
pred_base = runif(5)
pred_adverse = runif(5)
pred_severe = runif(5)

df = data.frame(LGD = c(historical, pred_hist, pred_base, pred_adverse, pred_severe), 
                label = c(rep("historical", 30), rep("predicted", 30), rep("base", 5), rep("adverse", 5), rep("severe", 5)),
                year = c(1987:2016,1987:2016,2017:2021,2017:2021,2017:2021))
df$label = factor(df$label, levels=c("historical","predicted", "base", "adverse","severe"))
p = ggplot(data=df, aes(x=year, y=LGD, group=label)) + geom_line(aes(linetype = label, col = label))+
    xlim(c(1987, 2022))+
    scale_color_manual(name = "", values=c("maroon2", "cyan3", "red", "blue", "yellow"))+
    scale_linetype_manual(name="",values= c("solid", "dashed", "solid", "solid","solid"))+
    theme_bw()#+theme(legend.title=element_blank())
p
##############################################################################################################



scale_linetype_manual(values=c("twodash", "dotted"))+


scale_color_manual(values = c("darkgrey", "red", "blue"),
                   labels = c("not sign", "neg", "pos"),
                   guide = guide_legend(override.aes = list(alpha = 1, size = 3)))




plot(eq_r$fitted.values, eq_r$residuals)
str(eq_r)










library(tikzDevice)
library(ggplot2)
library(latex2exp)

n <- 250
m <- 139

beamer.parms = list(paperwidth   = 364.19536/72,  # converts `pt` to `in`
                    paperheight  = 273.14662/72,
                    textwidth    = 307.28987/72,
                    textheight   = 269.14662/72)

##########################################################################
tikz(file='tikZfigs/binomial.sampling.distribution.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.8, 
     height = beamer.parms$paperheight * 0.65)


theta.null <- 0.5

x <- seq(0, n)
y = dbinom(x, size = n, prob = theta.null)
expected.value <- n*theta.null
z <- abs(expected.value-x) >= abs(expected.value-m)

Df <- data.frame(x = x,
                 y = y,
                 z = z)

ggplot(Df,
       mapping=aes(x = x, y = y, fill=z)) + 
  geom_col(width = 0.65) +
  xlim(sapply(c(1e-5, 1-1e-5), function(x) qbinom(x, size=n, prob=theta.null))) + 
  theme_classic() + 
  guides(fill=FALSE) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  xlab('Observed number of Heads') +
  ylab('Probability')

dev.off()

############################################################################
tikz(file='tikZfigs/binomial.likelihood.2.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.8, 
     height = beamer.parms$paperheight * 0.6)

standardized.likelihood <- function(x, N, n){
  dbeta(x, n+1, N-n+1)/dbeta(n/N, n+1, N-n-1)
}

find.roots <- function(n, N, z=1/4) {
  
  f <- function(x, N, n, z) {
    standardized.likelihood(x, N, n) - z
  }
  
  left.root <- uniroot(f, 
                       c(0, n/N), 
                       tol = 0.0001, 
                       N = N, 
                       n=n,
                       z=z)
  
  right.root <- uniroot(f, 
                        c(n/N, 1), 
                        tol = 0.0001, 
                        N = N, 
                        n=n,
                        z=z)
  
  list(x = left.root$root,
       y = standardized.likelihood(left.root$root, N, n),
       xend = right.root$root,
       yend = standardized.likelihood(right.root$root, N, n))
}

theta <- seq(0, 1, by = 0.001)
ll <- standardized.likelihood(theta, n, m)

Df <- data.frame(theta,
                 ll = standardized.likelihood(theta, n, m))

likelihood.interval <- 1/4

roots <- find.roots(m, n, likelihood.interval)

ggplot(Df,
       mapping=aes(x = theta, y = ll)) + 
  geom_line() +
  theme_classic() + 
  ylab('$P(\\theta)$') +
  xlab('$\\theta$') +
  geom_segment(aes(x = roots$x,
                   y = roots$y,
                   xend = roots$xend,
                   yend = roots$yend),
               col='red') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=rel(0.5)),
        axis.title=element_text(size=rel(0.7)))

dev.off()
############################################################################


########################################################################

tikz(file='tikZfigs/poisson.likelihood.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.65)

par(cex.lab=1.0,
    cex.axis=1.0,
    mar = c(3.0,2.5,0,0) )

n <- 10
S <- 107
curve((exp(-n*x)*x^S)/(exp(-n*(S/n))*(S/n)^S),
      from = 5.0,
      to = 20,
      xlab='',
      ylab='',
      col='salmon4',
      axes=F,
      n=1001)

x.ticks <- seq(5, 20, 5)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=2.0)

y.ticks <- seq(0, 1.0, by = 0.2)
axis(2, at=y.ticks, labels=rep('', length(y.ticks)), col.axis="black", las=1)
mtext('$\\textrm{P}(D\\vert\\lambda)$', side=2, line=1.5)

dev.off()


##################################################################
tikz(file='tikZfigs/gamma.posterior.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(3,4,1,0) )

curve(dgamma(x, shape=108, scale=1/(10+1/100)),
      from = 0.0,
      to = 20,
      col='salmon4',
      #ylim=c(0, 0.04),
      xlab='',
      ylab='',
      axes=F,
      n=1001)

x.ticks <- seq(0, 100, 10)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=2)

y.ticks <- seq(0, 0.4, by = 0.1)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(\\lambda\\vert D)$', side=2, line=3.0)

dev.off()

source('../code/gamma.R')

tikz(file='tikZfigs/gamma.posterior.hpd.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(3,4,1,0) )

kappa <- 108
theta <- 1/(10+1/100)

HPD <- gamma.hpd.interval(kappa, theta)
hpd.interval <- HPD$hpd.interval
p.star <- HPD$p.star

curve(dgamma(x, shape=kappa, scale=theta),
      from = 0.0,
      to = 20,
      col='salmon4',
      #ylim=c(0, 0.04),
      xlab='',
      ylab='',
      axes=F,
      n=1001)

segments(hpd.interval[1], p.star, hpd.interval[2], p.star, lwd=1, col='salmon4')

x.ticks <- seq(0, 100, 10)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=2)

y.ticks <- seq(0, 0.4, by = 0.1)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(\\lambda\\vert D)$', side=2, line=3.0)

dev.off()

# Also shows the hpd interval.
gamma.plot.hpd.2 <- function(alpha, beta){
  
  kappa <- 108
  theta <- 1/(10+1/100)
  
  HPD <- gamma.hpd.interval(kappa, theta)
  hpd.interval <- HPD$hpd.interval
  p.star <- HPD$p.star
  
  x <- seq(hpd.interval[1],
           hpd.interval[2],
           length.out = 1001)
  
  y <- dgamma(x, kappa, scale=theta)
  x <- c(hpd.interval[1], x, hpd.interval[2])
  y <- c(0, y, 0)
  
  plot.new()
  plot(NULL, NULL,
       ylab=expression("P" * (lambda)),
       xlab=expression(lambda),
       xlim=c(0,20),
       ylim=c(0, 1.1 * dgamma((kappa-1)*theta, kappa, scale=theta))
  )
  
  polygon(x, y, col="lavender", lty=0)
  
  curve(dgamma(x, kappa, scale=theta), 
        ylab='',
        xlab='',
        from=0, 
        to=20, 
        n=1001,
        add=T)
  
  #segments(hpd.interval[1], p.star, hpd.interval[2], p.star, lwd=1)
  
}

tikz(file='tikZfigs/poisson.posterior.predictive.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

V <- posterior.predictive(x)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,1,0) )

x.new <- seq(0, 25)
p <- dnbinom(x.new, size=V$r, prob = 1-V$q)
plot(x.new, 
     p, 
     ylim=c(0, 0.11),
     ylab='',
     xlab='',
     type='h', 
     col='chocolate',
     axes=F)

axis(1, at=x.new, labels=x.new, col.axis="black", las=1)
mtext('$k$', side=1, line=3)

y.ticks <- seq(0, 0.1, by=0.02)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(k\\vert D)$', side=2, line=4)

dev.off()



############################################################

tikz(file='tikZfigs/poisson.distribution.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.6)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,1,0) )

x.new <- seq(0, 15)
p <- dpois(x.new, lambda=5)
plot(x.new, 
     p, 
     ylim=c(0, 0.18),
     ylab='',
     xlab='',
     type='h', 
     col='chocolate',
     axes=F)

axis(1, at=x.new, labels=x.new, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=3)

y.ticks <- seq(0, 0.15, by=0.05)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(x=k\\vert \\lambda)$', side=2, line=4)

dev.off()

tikz(file='tikZfigs/gamma.distributions.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.6)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,0,0) )

colors <- c('chocolate',
            'coral4',
            'cadetblue4')

curve(dgamma(x, shape=2.0, scale=10),
      from = 0.0,
      to = 100,
      col=colors[1],
      ylim=c(0, 0.04),
      xlab='',
      ylab='',
      axes=F,
      n=1001)

curve(dgamma(x, shape=3.0, scale=12),
      from = 0.0,
      to = 100,
      col=colors[2],
      n=1001, add=T)

curve(dgamma(x, shape=1.0, scale=100),
      from = 0.0,
      to = 100,
      col=colors[3],
      n=1001, add=T)

legend(50, 0.035, 
       legend=c('$\\kappa=2$, $\\theta=10$', '$\\kappa=3$, $\\theta=12$', '$\\kappa=1.0$, $\\theta=100$'),
       lty=1, 
       col=colors, 
       bty='n')

x.ticks <- seq(0, 100, 10)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=3)

y.ticks <- seq(0, 0.04, by = 0.01)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(\\lambda)$', side=2, line=4.0)

dev.off()


# Contrôle croisé avec le package de référence dlnm (Gasparrini, 2011).
# install.packages("dlnm")
library(dlnm); library(splines)
d <- read.csv("serie journaliere simulee.csv")
d$date <- as.Date(d$date); d$temps <- seq_len(nrow(d))
n_annees <- length(unique(format(d$date, "%Y")))
cb <- crossbasis(d$tmax, lag = 7, argvar = list(fun = "ns", df = 3),
                 arglag = list(fun = "ns", knots = logknots(7, 2)))
m <- glm(admissions ~ cb + ns(temps, df = 7 * n_annees) + factor(jour_semaine),
         family = quasipoisson, data = d)
p99 <- quantile(d$tmax, 0.99)
pred <- crosspred(cb, m, cen = median(d$tmax), at = p99)
print(c(rr = pred$allRRfit, bas = pred$allRRlow, haut = pred$allRRhigh))

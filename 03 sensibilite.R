# Analyses de sensibilité et diagnostics
source("00 fonctions.R")
d <- read.csv("serie journaliere simulee.csv")
d$date <- as.Date(d$date); d$temps <- seq_len(nrow(d))
n_annees <- length(unique(format(d$date, "%Y")))
ref <- median(d$tmax); p99 <- quantile(d$tmax, 0.99)

ajuster <- function(lag_max = 7, df_an = 7, nd_var = 3) {
  cb <- base_croisee(d$tmax, lag_max = lag_max, nd_var = nd_var)
  m <- glm(admissions ~ cb + ns(temps, df = df_an * n_annees) + factor(jour_semaine),
           family = quasipoisson, data = d)
  r <- rr_cumule(m, cb, p99, ref)
  list(m = m, r = r)
}

scenarios <- expand.grid(lag_max = c(7, 14, 21), df_an = c(4, 7, 10))
scenarios <- cbind(scenarios, t(sapply(seq_len(nrow(scenarios)), function(i) {
  r <- ajuster(scenarios$lag_max[i], scenarios$df_an[i])$r
  round(c(rr = r$rr, bas = r$bas, haut = r$haut), 3)
})))
write.csv(scenarios, "analyses de sensibilite.csv", row.names = FALSE)
print(scenarios)

# Diagnostic : autocorrélation des résidus du modèle principal
principal <- ajuster()
res <- residuals(principal$m, type = "deviance")
lb <- Box.test(res, lag = 14, type = "Ljung-Box")
cat("Ljung-Box (14 lags) : p =", signif(lb$p.value, 3), "\n")

png("3 sensibilite et diagnostics.png", width = 1800, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))
lab <- paste0("lag ", scenarios$lag_max, " / ", scenarios$df_an, " ddl")
o <- seq_len(nrow(scenarios))
plot(scenarios$rr, o, xlim = range(c(scenarios$bas, scenarios$haut)), yaxt = "n",
     pch = 19, col = "#2b5d63", xlab = "RR cumulé au P99", ylab = "",
     main = "Sensibilité du résultat", cex.main = 1)
segments(scenarios$bas, o, scenarios$haut, o, col = "#2b5d63", lwd = 2)
axis(2, at = o, labels = lab, las = 1, cex.axis = 0.6); abline(v = 1, col = "grey50")
vrai_p99 <- exp(0.035 * max(p99 - 29, 0)) / exp(0.035 * max(ref - 29, 0))
abline(v = vrai_p99, col = "#c0583e", lty = 2, lwd = 2)
legend("bottomright", bty = "n", lty = 2, lwd = 2, col = "#c0583e", legend = "Vrai effet", cex = 0.7)
acf(res, lag.max = 30, main = "ACF des résidus de déviance", cex.main = 1)
dev.off()

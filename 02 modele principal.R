# Modèle principal : quasi-Poisson + base croisée (DLNM), lags 0 à 7
# (choix a priori : les effets de la chaleur sont surtout immédiats ;
# 14 et 21 jours sont testés en sensibilité)
source("00 fonctions.R")
d <- read.csv("serie journaliere simulee.csv")
d$date <- as.Date(d$date)
n_annees <- length(unique(format(d$date, "%Y")))
d$temps <- seq_len(nrow(d))

cb <- base_croisee(d$tmax, lag_max = 7, nd_var = 3)
m <- glm(admissions ~ cb + ns(temps, df = 7 * n_annees) + factor(jour_semaine),
         family = quasipoisson, data = d)

phi <- summary(m)$dispersion
reference <- median(d$tmax)          # référence : température médiane
p99 <- quantile(d$tmax, 0.99)
grille <- seq(floor(min(d$tmax)), ceiling(max(d$tmax)), by = 0.2)
courbe <- rr_cumule(m, cb, grille, reference)
lagc   <- rr_lag(m, cb, p99, reference)
vrai   <- read.csv("vrai effet cumule.csv")
vrai   <- vrai[vrai$x >= min(grille) & vrai$x <= max(grille), ]
# La vraie courbe est recentrée sur la même référence
vrai$rr <- vrai$rr / exp(0.035 * max(reference - 29, 0))

rr_p99 <- rr_cumule(m, cb, p99, reference)
vrai_p99 <- exp(0.035 * max(p99 - 29, 0)) / exp(0.035 * max(reference - 29, 0))

res <- data.frame(indicateur = c("Dispersion (phi)", "Température de référence (°C)",
                                 "P99 de Tmax (°C)", "RR cumulé au P99 (estimé)",
                                 "IC 95 % bas", "IC 95 % haut", "RR cumulé au P99 (vrai)"),
                  valeur = round(c(phi, reference, p99, rr_p99$rr, rr_p99$bas,
                                   rr_p99$haut, vrai_p99), 3))
write.csv(res, "resultats principaux.csv", row.names = FALSE)
print(res)

# Figure 1 : courbe exposition-réponse cumulée
png("1 courbe exposition reponse.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1), family = "sans")
plot(courbe$x, courbe$rr, type = "n", log = "y", ylim = range(c(courbe$bas, courbe$haut)),
     xlab = "Température maximale (°C)", ylab = "Risque relatif cumulé (lags 0-7)",
     main = "Chaleur et admissions : effet cumulé estimé (DLNM)", cex.main = 1)
polygon(c(courbe$x, rev(courbe$x)), c(courbe$bas, rev(courbe$haut)),
        col = adjustcolor("#3f7f86", 0.25), border = NA)
lines(courbe$x, courbe$rr, col = "#2b5d63", lwd = 2.5)
lines(vrai$x, vrai$rr, col = "#c0583e", lwd = 2, lty = 2)
abline(h = 1, col = "grey50"); abline(v = c(reference, p99), col = "grey70", lty = 3)
legend("topleft", bty = "n", lwd = c(2.5, 2, 8), lty = c(1, 2, 1),
       col = c("#2b5d63", "#c0583e", adjustcolor("#3f7f86", 0.25)),
       legend = c("Estimation", "Vrai effet simulé", "IC 95 %"))
dev.off()

# Figure 2 : structure des décalages au P99
png("2 structure des decalages.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1))
plot(lagc$lag, lagc$rr, type = "n", ylim = range(c(lagc$bas, lagc$haut)),
     xlab = "Décalage (jours)", ylab = "Risque relatif",
     main = sprintf("Effet d'un jour à %.1f °C (P99), selon le décalage", p99), cex.main = 1)
polygon(c(lagc$lag, rev(lagc$lag)), c(lagc$bas, rev(lagc$haut)),
        col = adjustcolor("#3f7f86", 0.25), border = NA)
lines(lagc$lag, lagc$rr, col = "#2b5d63", lwd = 2.5, type = "b", pch = 19)
abline(h = 1, col = "grey50")
dev.off()



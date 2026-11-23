# Étude de simulation : le modèle retrouve-t-il le vrai effet de façon fiable ?
# Une seule série simulée ne suffit pas à juger une méthode. On répète donc
# 200 fois : nouvelle série (même mécanisme que « 01 simuler donnees.R »,
# graine différente) -> même modèle que « 02 modele principal.R » -> RR cumulé
# au P99 de la série, comparé au vrai RR. On mesure le biais, l'erreur
# quadratique moyenne (RMSE) et la couverture de l'IC à 95 %.
source("R/00 fonctions.R")

simuler_serie <- function(graine) {
  set.seed(graine)
  dates <- seq(as.Date("2015-01-01"), as.Date("2023-12-31"), by = "day")
  n <- length(dates); jour_an <- as.numeric(format(dates, "%j"))
  saison <- 27 + 4 * sin(2 * pi * (jour_an - 300) / 365.25)
  bruit <- as.numeric(arima.sim(list(ar = 0.7), n = n, sd = 1.1))
  tmax <- round(saison + bruit + 0.04 * (seq_len(n) / 365.25), 1)
  log_rr <- stats::filter(0.035 * pmax(tmax - 29, 0), c(0.4, 0.3, 0.2, 0.1), sides = 1)
  log_rr[is.na(log_rr)] <- 0
  jour_sem <- factor(weekdays(dates))
  effet_sem <- c(0.06, 0.04, 0.02, 0, -0.02, -0.10, -0.12)[as.integer(jour_sem)]
  log_base <- log(25) + 0.08 * cos(2 * pi * jour_an / 365.25) + effet_sem +
              0.015 * (seq_len(n) / 365.25)
  data.frame(date = dates, tmax = tmax, jour_semaine = jour_sem,
             admissions = rnbinom(n, mu = exp(log_base + log_rr), size = 40))
}

ajuster <- function(d) {
  d$temps <- seq_len(nrow(d)); n_annees <- length(unique(format(d$date, "%Y")))
  cb <- base_croisee(d$tmax, lag_max = 7, nd_var = 3)
  m <- glm(admissions ~ cb + ns(temps, df = 7 * n_annees) + factor(jour_semaine),
           family = quasipoisson, data = d)
  ref <- median(d$tmax); p99 <- unname(quantile(d$tmax, 0.99))
  r <- rr_cumule(m, cb, p99, ref)
  vrai <- exp(0.035 * (max(p99 - 29, 0) - max(ref - 29, 0)))
  c(estime = r$rr, bas = r$bas, haut = r$haut, vrai = vrai)
}

n_rep <- 200
res <- t(sapply(seq_len(n_rep), function(i) ajuster(simuler_serie(10000 + i))))
res <- as.data.frame(res)
res$couvert <- res$bas <= res$vrai & res$vrai <= res$haut
write.csv(round(res, 4), "resultats/etude de simulation (200 repetitions).csv", row.names = FALSE)

err <- log(res$estime) - log(res$vrai)
bilan <- data.frame(
  indicateur = c("Répétitions", "Vrai RR au P99 (moyenne)", "RR estimé (moyenne géométrique)",
                 "Biais relatif du log-RR (%)", "RMSE du log-RR",
                 "Couverture de l'IC à 95 % (%)", "IC 95 % de la couverture : bas (%)",
                 "IC 95 % de la couverture : haut (%)", "Largeur moyenne de l'IC (rapport haut/bas)",
                 "Répétitions avec un effet significatif (puissance, %)"),
  valeur = c(n_rep, mean(res$vrai), exp(mean(log(res$estime))),
             100 * mean(err) / mean(log(res$vrai)), sqrt(mean(err^2)),
             100 * mean(res$couvert),
             100 * binom.test(sum(res$couvert), n_rep)$conf.int,
             mean(res$haut / res$bas), 100 * mean(res$bas > 1)))
bilan$valeur <- round(bilan$valeur, 3)
write.csv(bilan, "resultats/etude de simulation (bilan).csv", row.names = FALSE)
print(bilan)

# Figure 4 : les 200 estimations et leurs IC, triées, face au vrai effet
o <- order(res$estime)
png("figures/4 etude de simulation.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1))
plot(seq_len(n_rep), res$estime[o], type = "n", log = "y",
     ylim = range(c(res$bas, res$haut)), xlab = "Répétition (triée par estimation)",
     ylab = "RR cumulé au P99", cex.main = 1,
     main = sprintf("200 séries simulées : couverture de l'IC à 95 %% = %.1f %%", 100 * mean(res$couvert)))
segments(seq_len(n_rep), res$bas[o], seq_len(n_rep), res$haut[o],
         col = ifelse(res$couvert[o], adjustcolor("#3f7f86", 0.5), "#c0583e"))
points(seq_len(n_rep), res$estime[o], pch = 19, cex = 0.35, col = "#2b5d63")
abline(h = mean(res$vrai), col = "#c0583e", lwd = 2, lty = 2); abline(h = 1, col = "grey50")
legend("topleft", bty = "n", lty = c(1, 1, 2), lwd = c(2, 2, 2),
       col = c(adjustcolor("#3f7f86", 0.5), "#c0583e", "#c0583e"),
       legend = c("IC 95 % contenant le vrai effet", "IC 95 % ne le contenant pas", "Vrai effet (moyenne)"))
dev.off()

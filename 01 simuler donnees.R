# Simulation d'une série journalière réaliste (2015-2023) :
# température maximale et admissions hospitalières, avec un effet
# de la chaleur CONNU, pour vérifier que le modèle le retrouve.

set.seed(2026)
dates <- seq(as.Date("2015-01-01"), as.Date("2023-12-31"), by = "day")
n <- length(dates)
jour_an <- as.numeric(format(dates, "%j"))

# Température : saison + persistance (AR1) + réchauffement de 0,04 °C par an
saison <- 27 + 4 * sin(2 * pi * (jour_an - 300) / 365.25)
bruit <- as.numeric(arima.sim(list(ar = 0.7), n = n, sd = 1.1))
tendance <- 0.04 * (seq_len(n) / 365.25)
tmax <- round(saison + bruit + tendance, 1)

# Vrai effet : nul sous 29 °C, puis log-RR qui croît de 0,035 par °C,
# réparti sur les lags 0 à 3 (40 %, 30 %, 20 %, 10 %)
poids_lag <- c(0.4, 0.3, 0.2, 0.1)
effet_jour <- 0.035 * pmax(tmax - 29, 0)
log_rr <- stats::filter(effet_jour, poids_lag, sides = 1)
log_rr[is.na(log_rr)] <- 0

# Niveau de base : saisonnalité propre, effet du jour de la semaine, tendance
jour_sem <- factor(weekdays(dates))
effet_sem <- c(0.06, 0.04, 0.02, 0, -0.02, -0.10, -0.12)[as.integer(jour_sem)]
log_base <- log(25) + 0.08 * cos(2 * pi * jour_an / 365.25) + effet_sem +
            0.015 * (seq_len(n) / 365.25)

# Surdispersion : loi binomiale négative (taille 40)
admissions <- rnbinom(n, mu = exp(log_base + log_rr), size = 40)

d <- data.frame(date = dates, tmax = tmax, admissions = admissions,
                jour_semaine = jour_sem)
write.csv(d, "serie journaliere simulee.csv", row.names = FALSE)

# Vraie courbe exposition-réponse cumulée, pour comparaison
vrai <- data.frame(x = seq(20, 37, 0.1))
vrai$rr <- exp(0.035 * pmax(vrai$x - 29, 0))
write.csv(vrai, "vrai effet cumule.csv", row.names = FALSE)
cat("Série simulée :", n, "jours ; moyenne", round(mean(admissions), 1),
    "admissions par jour\n")

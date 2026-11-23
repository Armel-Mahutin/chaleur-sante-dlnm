# Fonctions utilitaires : base croisée (cross-basis) d'un modèle DLNM,
# écrite avec les seuls packages de base de R (stats, splines).
# Le package dlnm (Gasparrini, 2011) fait la même chose ; « 99 controle dlnm.R »
# compare les deux implémentations.

library(splines)

# Matrice des valeurs décalées : colonne l = exposition au jour t - l
matrice_decalee <- function(x, lag_max) {
  n <- length(x)
  sapply(0:lag_max, function(l) c(rep(NA, l), x[seq_len(n - l)]))
}

# Base croisée : produit tensoriel d'une base sur l'exposition (var)
# et d'une base sur les décalages (lag)
base_croisee <- function(x, lag_max, nd_var = 3, noeuds_lag = NULL) {
  bvar <- ns(x, df = nd_var)                      # base de l'exposition
  if (is.null(noeuds_lag))
    noeuds_lag <- exp(seq(log(1), log(lag_max), length.out = 4))[2:3]
  blag <- ns(0:lag_max, knots = noeuds_lag, intercept = TRUE)
  X <- matrice_decalee(x, lag_max)
  # Pour chaque fonction de base de l'exposition, on applique la base des lags
  cb <- do.call(cbind, lapply(seq_len(ncol(bvar)), function(v) {
    Bv <- apply(X, 2, function(col) predict(bvar, col)[, v])
    Bv %*% blag
  }))
  colnames(cb) <- paste0("cb", seq_len(ncol(cb)))
  attr(cb, "bvar") <- bvar
  attr(cb, "blag") <- blag
  attr(cb, "lag_max") <- lag_max
  cb
}

# Risque relatif cumulé (somme sur les lags) pour une grille d'expositions,
# centré sur une valeur de référence, avec IC 95 % par la méthode delta
rr_cumule <- function(modele, cb, grille, reference) {
  bvar <- attr(cb, "bvar"); blag <- attr(cb, "blag")
  idx  <- grep("^cb", names(coef(modele)))
  beta <- coef(modele)[idx]
  V    <- vcov(modele)[idx, idx]
  somme_lag <- colSums(blag)                       # cumul sur tous les lags
  Bx  <- predict(bvar, grille)
  Bref <- predict(bvar, reference)
  # Ligne de contraste pour chaque valeur de la grille
  L <- t(apply(Bx, 1, function(b) kronecker(b - Bref, somme_lag)))
  est <- as.vector(L %*% beta)
  se  <- sqrt(rowSums((L %*% V) * L))
  data.frame(x = grille, rr = exp(est),
             bas = exp(est - 1.96 * se), haut = exp(est + 1.96 * se))
}

# Courbe lag-réponse pour une exposition donnée (contre la référence)
rr_lag <- function(modele, cb, valeur, reference) {
  bvar <- attr(cb, "bvar"); blag <- attr(cb, "blag")
  idx  <- grep("^cb", names(coef(modele)))
  beta <- coef(modele)[idx]; V <- vcov(modele)[idx, idx]
  d <- predict(bvar, valeur) - predict(bvar, reference)
  L <- t(sapply(seq_len(nrow(blag)), function(l) kronecker(d, blag[l, ])))
  est <- as.vector(L %*% beta); se <- sqrt(rowSums((L %*% V) * L))
  data.frame(lag = 0:(nrow(blag) - 1), rr = exp(est),
             bas = exp(est - 1.96 * se), haut = exp(est + 1.96 * se))
}

# PARTIE 1----
# Script des cas limites----
# Script de trajectoires de population sans dimension - figure 2
#Installer et charger les bibliothèques nécessaires
install.packages("deSolve")
install.packages("ggplot2")
library(deSolve)
library(ggplot2)

#Modèle d'équations différentielles
model <- function(t, state, parameters) {
  P_auto <- state[1]  #Population autochtone
  P_allo <- state[2]  #Population allochthone
  K <- parameters$K   #Capacité de soutien
  I <- parameters$I   #Immigration
  Gamma <- parameters$Gamma   #Taux de regénération intrinsèque (facteur de croissance)
  
  #Normalisation des paramètres par la capacité de soutien
  P_auto_norm <- P_auto / K
  P_allo_norm <- P_allo / K
  
  #Équations différentielles du modèle, basées sur les équations 3 de l'article
  #On additionne les deux populations pour le facteur de densité-dépendance pour que notre échelle des y soit sur 1, et non 2
  dP_auto <- (P_auto_norm + P_allo_norm + I/Gamma) * (1 - (P_auto_norm + P_allo_norm))  
  dP_allo <- I/Gamma * (1 - (P_auto_norm + P_allo_norm))  
  
  list(c(dP_auto, dP_allo))
}

#Paramètres du modèle
parameters <- list(K = 1, I = 1)  

#Conditions initiales (ces conditions sont les mêmes que dans le graphique de l'article)
state <- c(P_auto = 0.01, P_allo = 0.01)

#Temps de simulation
times <- seq(0, 5, by = 0.1)

#Résoudre l'ODE pour différents Gamma
parameters$Gamma <- 10   #Cas intermédiaire
out_intermediate <- ode(y = state, times = times, func = model, parms = parameters)

parameters$Gamma <- 0.1   #Immigration dominante (Γ → 0)
out_immigration <- ode(y = state, times = times, func = model, parms = parameters)

parameters$Gamma <- 100   #Croissance intrinsèque dominante (Γ → ∞)
out_intrinsic <- ode(y = state, times = times, func = model, parms = parameters)

#Convertir les résultats en dataframes
df_intermediate <- as.data.frame(out_intermediate)
df_immigration <- as.data.frame(out_immigration)
df_intrinsic <- as.data.frame(out_intrinsic)

#Calculer la population totale
df_intermediate$Total <- df_intermediate$P_auto + df_intermediate$P_allo
df_intrinsic$Total <- df_intrinsic$P_auto + df_intrinsic$P_allo
df_immigration$Total <- df_immigration$P_auto + df_immigration$P_allo

#Tracer la figure avec inversion des couleurs des lignes pleines
ggplot() +
  #Courbes en pointillés pour les cas limites (pas modifiées)
  geom_line(data = df_intrinsic, aes(x = time, y = Total, color = "Γ → ∞"), linetype = "dashed", size = 1) +
  geom_line(data = df_immigration, aes(x = time, y = Total, color = "Γ → 0"), linetype = "dashed", size = 1) +
  
  #Courbes solides pour Γ = 10 
  geom_line(data = df_intermediate, aes(x = time, y = Total, color = "Total"), size = 1) +
  geom_line(data = df_intermediate, aes(x = time, y = P_allo, color = "Autochthonous"), size = 1) +  
  geom_line(data = df_intermediate, aes(x = time, y = P_auto, color = "Allochthonous"), size = 1) +  
  
  #Personnalisation avec couleurs corrigées
  labs(x = "Time (τ)", y = expression(P^"*"), color = "Population Type") +
  theme_minimal() +
  ggtitle("Nondimensional population recovery trajectories") +
  scale_color_manual(values = c("Γ → ∞" = "steelblue", "Γ → 0" = "goldenrod", 
                                "Total" = "black", "Allochthonous" = "goldenrod", "Autochthonous" = "steelblue"))


#Script graphiques avant-récif 3a ----
#Aller chercher dans les bibliothèques
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)


#Paramètres du modèle
params <- list(
  rb = 0.05,   #Taux de croissance autochtone de la population backreef
  rf = 0.15,  #Taux de croissance autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  
  #Coefficients d’interactions (ζ(γβ))
  zeta_bb = 0.01,   #Influence de backreef sur backreef
  zeta_fb = 0.002,    #Influence de forereef sur backreef
  zeta_ff = 0.02,    #Influence de forereef sur forereef AUTO
  zeta_bf = 0.1,   #Influence de backreef sur forereef ALLO (négligeable)
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion?
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

# MODELE 2 - Équations différentielles selon les équations 4
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dFallo  <- sigma_f * zeta_bf * pBtotal + If * (1 - pFtotal / Kf)
    dFauto  <- pFtotal * (rf + sigma_f * zeta_ff) * (1 - pFtotal / Kf)
    dFtotal <- dFallo + dFauto
    dBallo  <- sigma_f * zeta_fb * pFtotal + Ib * (1 - pBtotal / Kb)
    dBauto  <- pBtotal * (rb + sigma_f * zeta_bb) * (1 - pBtotal / Kb)
    dBtotal <- dBallo + dBauto
    list(c(dFallo, dFauto, dFtotal, dBallo, dBauto, dBtotal))
  })
}
#Conditions initiales  MODELE 2 de nos populations
state <- c(pFallo = 0.01, pFauto = 0.01, pFtotal = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années
times <- seq(0, 7, by = 0.1)

#Résolution numérique du système
out <- ode(y = state, times = times, func = model, parms = params)
out <- as.data.frame(out)

#Mise en forme pour ggplot
plot_data <- out %>%
  select(time, pFallo,pFauto,pFtotal) %>%
  rename(
    Total = pFtotal,
    Allochthonous = pFallo,
    Autochthonous = pFauto
  ) %>%
  pivot_longer(-time, names_to = "Source", values_to = "Value")

#Graphique
ggplot(plot_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.4) +
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochthonous" = "goldenrod",
      "Autochthonous" = "steelblue"  
    )
  ) +
  labs(
    title = "Croissance de la population corallienne (Pf)",
    x = "Temps (années)",
    y = "Couverture corallienne proportion"
  ) +
  theme_minimal(base_size = 14)





# Script graphiques avant-récif 3c----
#Aller chercher dans les bibliothèques
library(ggplot2)
library(dplyr)
library(tidyr)

#Données modélisées
years <- 2011:2018   #Années de l'étude
n <- length(years)

Pf <- seq(0.05, 0.6, length.out = n)   #Population forereef
Pb <- seq(0.2, 0.5, length.out = n)   #Population backreef
Pf[1] <- 0   #Forcer à commencer à 0
Pb[1] <- 0   #Force à commencer à 0

Kf <- 0.8   #Capacité de charge de 80%
rf <- 0.12   #Taux de croissance asexuée 
sigma_f <- 0.05   #Taux de reproduction larvaire (estimé)

zeta_ff <- seq(0.2, 0.8, length.out = n)   #Autorecrutement (forereef vers forereef)
zeta_bf <- seq(0.01, 0.01, length.out = n)   #Recrutement d'échange (forereef vers backreef)
If <- c(0.0075, 0.0095, 0.0085, 0.008, 0.007, 0.006, 0.005, 0.002)   #Immigration externe (données modélisées)

#Calcul des composantes dPf/dt 
asexual <- Pf * rf   #Croissance asexuée
self_recruit <- Pf * sigma_f * zeta_ff   #Autorecrutement (autochtone)
local_exchange <- Pb * sigma_f * zeta_bf   #Échange local
immigration <- If   #Immigration (allochtone)

density_factor <- (1 - Pf / Kf)   #Facteur de densité-dépendance

#Ajustement de chaque composante par le facteur de densité-dépendance
asexual_adj <- asexual * density_factor
self_recruit_adj <- self_recruit * density_factor
local_exchange_adj <- local_exchange * density_factor
immigration_adj <- immigration * density_factor

#Création du dataframe pour le graphique
df <- data.frame(
  Year = rep(years, 4),
  Component = rep(c("Asexual Growth", "Self-Recruitment", "Local Exchange", "Immigration"), each = n),
  Value = c(asexual_adj, self_recruit_adj, local_exchange_adj, immigration_adj)
)

#Mise en forme des données en pourcentage de la capacité de charge par an
df$PercentK <- 100 * df$Value / Kf
scaling_factor <- 2.5   #Amplification des courbes pour un meilleur effet visuel
df$PercentK <- df$PercentK * scaling_factor

#Graphique lissé
ggplot(df, aes(x = Year, y = PercentK, linetype = Component)) +
  geom_smooth(se = FALSE, method = "loess", size = 1.2, color = "black") +
  scale_linetype_manual(
    values = c("solid", "dotted", "dotdash", "dashed"),
    breaks = c("Asexual Growth", "Immigration", "Self-Recruitment", "Local Exchange")
  ) +
  labs(
    title = "Components contributing to changes in population growth (LTER 1 Forereef)",
    y = expression("%K / yr"),
    x = NULL,
    linetype = NULL
  ) +
  coord_cartesian(ylim = c(0, 11)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    plot.title = element_text(hjust = 0.4, face = "bold", size = 8)
  )





# Script graphiques après-récifs 4a----
#Différences: ajout de mortalité et changements états initiaux. 
params <- list(
  mortal <- 0.3,  #taux de mortalité du backreef. 
  
  rb = 0.05,   #Taux de croissance autochtone de la population backreef
  rf = 0.15,  #Taux de croissance autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  
  #Coefficients d’interactions (ζ(γβ))
  zeta_bb = 0.01,   #Influence de backreef sur backreef
  zeta_fb = 0.002,    #Influence de forereef sur backreef #faible
  zeta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  zeta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef  tiré du graphique 4c
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion?
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

# MODELE 2 - Équations différentielles selon les équations 4
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dFallo  <- sigma_f * zeta_bf * pBtotal + If * (1 - pFtotal / Kf)
    dFauto  <- pFtotal * (rf + sigma_f * zeta_ff) * (1 - pFtotal / Kf)
    dFtotal <- dFallo + dFauto
    dBallo  <- sigma_f * zeta_fb * pFtotal + Ib * (1 - pBtotal / Kb)
    dBauto  <- pBtotal * (rb + sigma_f * zeta_bb) * (1 - pBtotal / Kb)
    dBtotal <- dBallo + dBauto - mortal * pBtotal           ####ici, ajout du paramètre de mortalité
    list(c(dFallo, dFauto, dFtotal, dBallo, dBauto, dBtotal))
  })
}
#Conditions initiales  MODELE 2 de nos populations
state <- c(pFallo = 0.01, pFauto = 0.01, pFtotal = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années
times <- seq(0, 6, by = 0.1)

#Résolution numérique du système
out <- ode(y = state, times = times, func = model, parms = params)
out <- as.data.frame(out)

#graphique qui reproduit la figure 4a
ggplot(out, aes(x = time, y = pBtotal)) +
  geom_line(color = "red", size = 1.4) +
  labs(
    title = "Couverture corallienne totale du backreef",
    x = "Temps (années)",
    y = "Proportion de couverture (pBtotal)"
  ) +
  theme_minimal(base_size = 14)


#Si on veux afficher vraiment toutes les courbes.(garder en tête que la mortalité est appliquée après la somme des courbes)
backreef_data <- out %>%
  select(time, pBallo, pBauto, pBtotal) %>%
  rename(
    Allochthonous = pBallo,
    Autochthonous = pBauto,
    Total = pBtotal
  ) %>%
  pivot_longer(-time, names_to = "Source", values_to = "Value")

# Graphique
ggplot(backreef_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.4) +
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochthonous" = "goldenrod",
      "Autochthonous" = "steelblue"
    )
  ) +
  labs(
    title = "Dynamique de la population corallienne — Backreef",
    subtitle = "Contributions autochtone, allochone et totale",
    x = "Temps (années)",
    y = "Couverture corallienne (proportion)",
    color = "Source"
  ) +
  theme_minimal(base_size = 14)





# Script graphiques après-récifs 4c----
#Aller chercher les bibliothèques nécessaires
library(ggplot2)
library(dplyr)
library(tidyr)

#Paramètres
years <- 2011:2018   #Années de l'étude
n <- length(years)

# Populations
Pf <- seq(0.05, 0.6, length.out = n)  #Population forereef
Pb <- seq(0.2, 0.5, length.out = n)   #Population backreef
Pf[1] <- 0   #Forcer à commencer à 0
Pb[1] <- 0

#Paramètres biologiques
Kb <- 0.8   #Capacité de soutien (backreef)
rb <- 0.1   #Croissance asexuée (backreef)
sigma_f <- 0.05   #Taux de reproduction larvaire

zeta_bb <- seq(0.2, 0.7, length.out = n)  #Autorecrutement backreef vers backreef
zeta_fb <- rep(0.01, n)  #Échange forereef vers backreef

#Calcul des composantes dPb/dt
intrinsic <- Pb * (rb + sigma_f * zeta_bb)  #Croissance et autorecrutement
exchange <- Pf * sigma_f * zeta_fb          #Recrutement d’échange local
immigration <- rep(0.0007, n)   #Immigration

#Facteur de densité-dépendance
density_factor <- (1 - Pb / Kb)

#Ajustement de chaque composante par le facteur de densité-dépendance
intrinsic_adj <- intrinsic * density_factor
exchange_adj <- exchange * density_factor
immigration_adj <- immigration   #L'immigragtion reste constante peu importe la densité de la population

#Variation totale
total_change <- intrinsic_adj + exchange_adj + immigration_adj

#Estimation de la courbe de mortalité
mortalite <- -0.003 + 0.0022 * (years - min(years)) / ((years - min(years)) + 4)

#Création du dataframe pour le graphique
df_backreef <- data.frame(
  Year = rep(years, 2),
  Component = rep(c("Immigration", "Mortality"), each = n),
  Value = c(immigration_adj, mortalite)
)

#Mise en forme des données en pourcentage de la capacité de charge par an
df_backreef$PercentKb <- 100 * df_backreef$Value / Kb

#Graphique lissé
ggplot(df_backreef, aes(x = Year, y = PercentKb, linetype = Component)) +
  geom_smooth(se = FALSE, method = "loess", size = 1.2, color = "black") +
  scale_linetype_manual(values = c("Immigration" = "dashed", "Mortality" = "solid")) +
  labs(
    title = "Components contributing to population change (LTER 1 Backreef)",
    y = expression("%K"[b] * " / yr"),
    x = NULL,
    linetype = NULL
  ) +
  coord_cartesian(ylim = c(min(df_backreef$PercentKb), max(df_backreef$PercentKb))) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10)
  )





# PARTIE 2----
# Script graphique peak=1 period=7----
#fonction qui représente la variation de mortalité aux 7 ans
times <- seq(0,30, by=1)

#voir photo chatgpt pour comprendre que pour une vrai chute de 99%, il faut peak=5.
taux_perte <- function(t, peak = 1, period = 7, pulse_width = 1) {   #peak=amplitude de la perturbation. C'est une perte continue de 0,50/ans appliqué sur [width] ans.
  if ((t %% period) < pulse_width) { 
    return(peak)
  } else {
    
                  #donc si le temps est pas dans la période de perturbation, on retourne 0
    return(0) 
  }
}
#visualisation de cette fonction
plot(times, sapply(times, taux_perte), type = "l",
     ylab = "Taux de perte", xlab = "Temps", col = "red", lwd = 2)


#Paramètres du modèle
params <- list(
  mortal <- 0.6,  #taux de mortalité constante du backreef pck c difficile la vie.
  
  rb = 0.05,   #Taux de croissance autochtone de la population backreef
  rf = 0.15,  #Taux de croissance autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  
  #Coefficients d’interactions (ζ(γβ))
  zeta_bb = 0.01,   #Influence de backreef sur backreef
  zeta_fb = 0.002,    #Influence de forereef sur backreef
  zeta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  zeta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion (basé sur étendue max observée, tiré de l'article)
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

# MODELE 2 - Équations différentielles selon les équations 4 ici, j'ai enlevé pftotal, et je l'ai remplacé par la somme. 
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    pertes_t <- taux_perte(time)
    
    dFallo  <- (sigma_f * zeta_bf * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * zeta_ff) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * zeta_fb * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * zeta_bb) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
    dBtotal <- (dBallo + dBauto - mortal * pBtotal) 
    
    list(c(dFallo, dFauto, dBallo, dBauto, dBtotal))
  })
}
#Conditions initiales  MODELE 2 de nos populations
state <- c(pFallo = 0.01, pFauto = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années
times <- seq(0, 30, by = 0.1)

#Résolution numérique du système
out <- ode(y = state, times = times, func = model, parms = params)
out <- as.data.frame(out)
out$pFtotal <- out$pFauto + out$pFallo  #calcul de la somme

#Mise en forme pour ggplot
plot_data <- out %>%
  select(time, pFallo,pFauto,pFtotal) %>%
  rename(
    Total = pFtotal,
    Allochthonous = pFallo,
    Autochthonous = pFauto
  ) %>%
  pivot_longer(-time, names_to = "Source", values_to = "Value")

#Graphique
ggplot(plot_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.4) +
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochthonous" = "goldenrod",
      "Autochthonous" = "steelblue"  
    )
  ) +
  labs(
    title = "Suivi de la population corallienne (Forereef)",
    subtitle = "Perturbation réelle de ±63% aux 7 ans, immigration standard",
    x = "Temps (années)",
    y = "Couverture corallienne"
  ) +
  theme_minimal(base_size = 12)




# Script graphique peak=1 period=2----
#fonction qui représente la variation de mortalité aux 7 ans
times <- seq(0,30, by=1)

#voir photo chatgpt pour comprendre que pour une vrai chute de 99%, il faut peak=5.
taux_perte <- function(t, peak = 1, period = 2, pulse_width = 1) {   #peak=amplitude de la perturbation. C'est une perte continue de 0,50/ans appliqué sur [width] ans.
  if ((t %% period) < pulse_width) { 
    return(peak)
  } else {
    
    #donc si le temps est pas dans la période de perturbation, on retourne 0
    return(0) 
  }
}
#visualisation de cette fonction
plot(times, sapply(times, taux_perte), type = "l",
     ylab = "Taux de perte", xlab = "Temps", col = "red", lwd = 2)


#Paramètres du modèle
params <- list(
  mortal <- 0.6,  #taux de mortalité constante du backreef pck c difficile la vie.
  
  rb = 0.05,   #Taux de croissance autochtone de la population backreef
  rf = 0.15,  #Taux de croissance autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  
  #Coefficients d’interactions (ζ(γβ))
  zeta_bb = 0.01,   #Influence de backreef sur backreef
  zeta_fb = 0.002,    #Influence de forereef sur backreef
  zeta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  zeta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion (basé sur étendue max observée, tiré de l'article)
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

# MODELE 2 - Équations différentielles selon les équations 4 ici, j'ai enlevé pftotal, et je l'ai remplacé par la somme. 
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    pertes_t <- taux_perte(time)
    
    dFallo  <- (sigma_f * zeta_bf * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * zeta_ff) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * zeta_fb * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * zeta_bb) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
    dBtotal <- (dBallo + dBauto - mortal * pBtotal) 
    
    list(c(dFallo, dFauto, dBallo, dBauto, dBtotal))
  })
}
#Conditions initiales  MODELE 2 de nos populations
state <- c(pFallo = 0.01, pFauto = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années
times <- seq(0, 30, by = 0.1)

#Résolution numérique du système
out <- ode(y = state, times = times, func = model, parms = params)
out <- as.data.frame(out)
out$pFtotal <- out$pFauto + out$pFallo  #calcul de la somme

#Mise en forme pour ggplot
plot_data <- out %>%
  select(time, pFallo,pFauto,pFtotal) %>%
  rename(
    Total = pFtotal,
    Allochthonous = pFallo,
    Autochthonous = pFauto
  ) %>%
  pivot_longer(-time, names_to = "Source", values_to = "Value")

#Graphique
ggplot(plot_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.4) +
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochthonous" = "goldenrod",
      "Autochthonous" = "steelblue"  
    )
  ) +
  labs(
    title = "Suivi de la population corallienne (Forereef)",
    subtitle = "Perturbation réelle de ±63% aux 2 ans, immigration standard",
    x = "Temps (années)",
    y = "Couverture corallienne"
  ) +
  theme_minimal(base_size = 12)

# Script graphique peak=5 period=7----
#fonction qui représente la variation de mortalité aux 7 ans
#fonction qui représente la variation de mortalité aux 7 ans
times <- seq(0,30, by=1)

#voir photo chatgpt pour comprendre que pour une vrai chute de 99%, il faut peak=5.
taux_perte <- function(t, peak = 5, period = 7, pulse_width = 1) {   #peak=amplitude de la perturbation. C'est une perte continue de 0,50/ans appliqué sur [width] ans.
  if ((t %% period) < pulse_width) { 
    return(peak)
  } else {
    
    #donc si le temps est pas dans la période de perturbation, on retourne 0
    return(0) 
  }
}
#visualisation de cette fonction
plot(times, sapply(times, taux_perte), type = "l",
     ylab = "Taux de perte", xlab = "Temps", col = "red", lwd = 2)


#Paramètres du modèle
params <- list(
  mortal <- 0.6,  #taux de mortalité constante du backreef pck c difficile la vie.
  
  rb = 0.05,   #Taux de croissance autochtone de la population backreef
  rf = 0.15,  #Taux de croissance autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  
  #Coefficients d’interactions (ζ(γβ))
  zeta_bb = 0.01,   #Influence de backreef sur backreef
  zeta_fb = 0.002,    #Influence de forereef sur backreef
  zeta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  zeta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion (basé sur étendue max observée, tiré de l'article)
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

# MODELE 2 - Équations différentielles selon les équations 4 ici, j'ai enlevé pftotal, et je l'ai remplacé par la somme. 
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    pertes_t <- taux_perte(time)
    
    dFallo  <- (sigma_f * zeta_bf * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * zeta_ff) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * zeta_fb * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * zeta_bb) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
    dBtotal <- (dBallo + dBauto - mortal * pBtotal) 
    
    list(c(dFallo, dFauto, dBallo, dBauto, dBtotal))
  })
}
#Conditions initiales  MODELE 2 de nos populations
state <- c(pFallo = 0.01, pFauto = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années
times <- seq(0, 30, by = 0.1)

#Résolution numérique du système
out <- ode(y = state, times = times, func = model, parms = params)
out <- as.data.frame(out)
out$pFtotal <- out$pFauto + out$pFallo  #calcul de la somme

#Mise en forme pour ggplot
plot_data <- out %>%
  select(time, pFallo,pFauto,pFtotal) %>%
  rename(
    Total = pFtotal,
    Allochthonous = pFallo,
    Autochthonous = pFauto
  ) %>%
  pivot_longer(-time, names_to = "Source", values_to = "Value")

#Graphique
ggplot(plot_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.4) +
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochthonous" = "goldenrod",
      "Autochthonous" = "steelblue"  
    )
  ) +
  labs(
    title = "Suivi de la population corallienne (Forereef)",
    subtitle = "Perturbation réelle de ±99% aux 7 ans, immigration standard",
    x = "Temps (années)",
    y = "Couverture corallienne"
  ) +
  theme_minimal(base_size = 12)

# Script graphique peak=5 period=2----
#fonction qui représente la variation de mortalité aux 7 ans
times <- seq(0,30, by=1)

#voir photo chatgpt pour comprendre que pour une vrai chute de 99%, il faut peak=5.
taux_perte <- function(t, peak = 5, period = 2, pulse_width = 1) {   #peak=amplitude de la perturbation. C'est une perte continue de 0,50/ans appliqué sur [width] ans.
  if ((t %% period) < pulse_width) { 
    return(peak)
  } else {
    
    #donc si le temps est pas dans la période de perturbation, on retourne 0
    return(0) 
  }
}
#visualisation de cette fonction
plot(times, sapply(times, taux_perte), type = "l",
     ylab = "Taux de perte", xlab = "Temps", col = "red", lwd = 2)


#Paramètres du modèle
params <- list(
  mortal <- 0.6,  #taux de mortalité constante du backreef pck c difficile la vie.
  
  rb = 0.05,   #Taux de croissance autochtone de la population backreef
  rf = 0.15,  #Taux de croissance autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  
  #Coefficients d’interactions (ζ(γβ))
  zeta_bb = 0.01,   #Influence de backreef sur backreef
  zeta_fb = 0.002,    #Influence de forereef sur backreef
  zeta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  zeta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion (basé sur étendue max observée, tiré de l'article)
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

# MODELE 2 - Équations différentielles selon les équations 4 ici, j'ai enlevé pftotal, et je l'ai remplacé par la somme. 
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    pertes_t <- taux_perte(time)
    
    dFallo  <- (sigma_f * zeta_bf * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * zeta_ff) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * zeta_fb * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * zeta_bb) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
    dBtotal <- (dBallo + dBauto - mortal * pBtotal) 
    
    list(c(dFallo, dFauto, dBallo, dBauto, dBtotal))
  })
}
#Conditions initiales  MODELE 2 de nos populations
state <- c(pFallo = 0.01, pFauto = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années
times <- seq(0, 30, by = 0.1)

#Résolution numérique du système
out <- ode(y = state, times = times, func = model, parms = params)
out <- as.data.frame(out)
out$pFtotal <- out$pFauto + out$pFallo  #calcul de la somme

#Mise en forme pour ggplot
plot_data <- out %>%
  select(time, pFallo,pFauto,pFtotal) %>%
  rename(
    Total = pFtotal,
    Allochthonous = pFallo,
    Autochthonous = pFauto
  ) %>%
  pivot_longer(-time, names_to = "Source", values_to = "Value")

#Graphique
ggplot(plot_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.4) +
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochthonous" = "goldenrod",
      "Autochthonous" = "steelblue"  
    )
  ) +
  labs(
    title = "Suivi de la population corallienne (Forereef)",
    subtitle = "Perturbation réelle de ±99% aux 2 ans, immigration standard",
    x = "Temps (années)",
    y = "Couverture corallienne"
  ) +
  theme_minimal(base_size = 12)

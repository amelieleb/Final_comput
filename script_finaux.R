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
  dP_allo <- I/Gamma * (1 - (P_auto_norm + P_allo_norm))  #On prend seulement I/Gamma, car cette population est seulement immigratrice
  
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
# ALLER CHERCHER LES PACKAGES NÉCESSAIRES
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)


# DÉFINITION DES PARAMÈTRES DU MODÈLE
params <- list(
  
  #Paramètre de reproduction asexuée
  rb = 0.05,        #Reproduction asexuée de la population backreef
  rf = 0.15,        #Reproduction asexuée de la population forereef
  
  
  #Paramètre de reproduction sexuée
  sigma = 0.5,    #Intensité de la croissance larvaire
  fecondity = 1,    #Fécondité
  zeta = 1,         #Probabilité de maturation larvaire
  
  #Coefficients d’interactions (γβ)
  gamma_beta_bb = 0.01,   #Influence de backreef sur backreef
  gamma_beta_fb = 0.002,  #Influence de forereef sur backreef
  gamma_beta_ff = 0.02,   #Influence de forereef sur forereef
  gamma_beta_bf = 0.1,    #Influence de backreef sur forereef
  
  #Paramètre d'immigration nette (apport externe constant)
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,   #Immigration du forereef
  
  #Capacités de soutien en proportion
  Kb = 0.64,   #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

#Conditions initiale du modèle
state <- c(pFallo = 0.01, pFauto = 0.01, pFtotal = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années (2011 à 2018)
times <- seq(0, 7, by = 0.1)




# CRÉATION DU MODÈLE DE CROISSANCE SELON LES ÉQUATIONS 4 DE L'ARTICLE
model_proportion_forereef <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 1. Séparation de l'équation de la croissance de sous-population avant-récif en 3 équations:
    # Croissance du au recrutement autochtone et à la reproduction asexuée (dFauto/dt)
    dFauto  <- (pFtotal * (rf + (sigma * fecondity * zeta * gamma_beta_ff))) * (1 - pFtotal / Kf)
    
    # Croissance du au recrutement allochtone et l'immigration nette (dFallo/dt)
    dFallo  <- ((pBtotal * (sigma * fecondity * zeta * gamma_beta_bf)) + If) * (1 - pFtotal / Kf)
    
    # Croissance totale de l'avant-récif (dB/dt)
    dFtotal <- dFallo + dFauto
    
    
    #2. Sépration de l'équation de la croissance de sous-population après-récif en 3 équations:
    # Croissance du au recrutement autochtone et à la reproduction asexuée (dBauto/dt)
    dBauto  <- (pBtotal * (rb + (sigma * fecondity * zeta * gamma_beta_bb))) * (1 - pBtotal / Kb)
    
    # Croissance du au recrutement allochtone et l'immigration nette (dBallo/dt)
    dBallo  <- ((pFtotal * (sigma * fecondity * zeta * gamma_beta_fb)) + Ib) * (1 - pBtotal / Kb)
    
    # Croissance totale de l'après-récif (dB/dt)
    dBtotal <- dBallo + dBauto
    
    
    list(c(dFallo, dFauto, dFtotal, dBallo, dBauto, dBtotal))
  })
}



# RÉSOLUTION NUMÉRIQUE DU MODÈLE
out <- ode(y = state, times = times, func = model_proportion_forereef, parms = params)
out <- as.data.frame(out)



# GRAPHIQUE
#Mise en forme pour ggplot
plot_data <- out %>%
  select(time, pFallo,pFauto,pFtotal) %>%
  
  #Renommer les paramètres utilisées
  rename(
    Total = pFtotal,
    Allochtone = pFallo,
    Autochtone = pFauto
  ) %>%
  
  #Mettre en pourcentage les valeurs modélisées
  mutate(
    Total = Total*100,
    Allochtone = Allochtone*100,
    Autochtone = Autochtone*100
  ) %>%
  
  pivot_longer(-time, names_to = "Source", values_to = "Value")


#Graphique
ggplot(plot_data, aes(x = time, y = Value, color = Source)) +
  geom_line(size = 1.7) +
  
  #Définition de l'axe des x
  scale_x_continuous(
    breaks = c(1,3,5,7),
    labels= c(2012, 2014, 2016, 2018)
  ) +
  
  #Définition de l'axe des y
  scale_y_continuous(
    limits=c(0,75),
    breaks = c(0, 25, 50, 75)
  ) +
  
  #Définition des couleurs
  scale_color_manual(
    values = c(
      "Total" = "black",
      "Allochtone" = "goldenrod",
      "Autochtone" = "steelblue"  
    )
  ) +
  
  #Définition des titres
  labs(
    title = "Croissance de la sous-population avant-récif corallienne (Pf)",
    x = NULL,
    y = "Couverture corallienne (%)",
    color = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  
  #Ajustement de certains paramètres visuels
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left","top"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(10, 30, 10, 10)
  ) +
  
  #Aligner l'axe des x avec le 0 de l'axe des y
  coord_cartesian(expand = FALSE)




# Script graphiques avant-récif 3c----

# ALLER CHERCHER LES PACKAGES NÉCESSAIRES
library(ggplot2)
library(dplyr)
library(tidyr)


# DÉFINITION DES PARAMÈTRES UTILISÉS

#Paramètres de temps
years <- 2011:2018   #Années de l'étude
n <- length(years)

#Paramètre de taille de population
Pf <- seq(0.05, 0.6, length.out = n)   #Population forereef
Pb <- seq(0.2, 0.5, length.out = n)    #Population backreef
Pf[1] <- 0                             #Forcer à commencer à 0
Pb[1] <- 0                             #Force à commencer à 0

#Paramètre de capacité de soutien
Kf <- 0.8            #Capacité de charge de 80%

#Paramètres de reproduction
rf <- 0.12           #Reproduction asexuée 
sigma_f <- 0.05      #Intensité de la croissance larvaire

#Coefficients d'interaction (γβ)
gamma_beta_ff <- seq(0.2, 0.8, length.out = n)   #forereef vers forereef
gamma_beta_bf <- seq(0.01, 0.01, length.out = n) #forereef vers backreef

#Immigration externe (données modélisées)
If <- c(0.0075, 0.0095, 0.0085, 0.008, 0.007, 0.006, 0.005, 0.002)



#Calcul des 4 composantes  de dPf/dt et ajustement selon le facteur de densité-dépendance
#Facteur de densité-dépendance
density_factor <- (1 - Pf / Kf)

#Croissance asexuée
asexual <- Pf * rf
asexual_adj <- asexual * density_factor

#Auto-recrutement (autochtone)
self_recruit <- Pf * sigma_f * gamma_beta_ff
self_recruit_adj <- self_recruit * density_factor

#Échange local
local_exchange <- Pb * sigma_f * gamma_beta_bf
local_exchange_adj <- local_exchange * density_factor

#Immigration nette (allochtone)
immigration <- If 
immigration_adj <- immigration * density_factor



# GRAPHIQUE
#Création du dataframe pour le graphique
df <- data.frame(
  Year = rep(years, 4),
  Component = rep(c("Croissance asexuée", "Auto-recrutement", "Échange local", "Immigration"), each = n),
  Value = c(asexual_adj, self_recruit_adj, local_exchange_adj, immigration_adj)
)

#Mise en forme des données en pourcentage de la capacité de charge par an
df$PercentK <- 100 * df$Value / Kf

#Amplification des courbes pour un meilleur effet visuel
scaling_factor <- 2.5

df$PercentK <- df$PercentK * scaling_factor


#Graphique
ggplot(df, aes(x = Year, y = PercentK, linetype = Component)) +
  geom_smooth(se = FALSE, method = "loess", size = 1.2, color = "black") +
  
  #Définition du type de lignes utilisées
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash"),
    breaks = c("Croissance asexuée", "Immigration", "Auto-recrutement", "Échange local")
  ) +
  
  #Définition de l'axe des y
  scale_y_continuous(
    limits = c(0, 11),
    breaks = c(0, 3, 7, 11)
  )+
  
  #Définition des titres
  labs(
    title = "Composants contribuant aux changements de la croissance de la population",
    y = expression("%K / yr"),
    x = NULL,
    linetype = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  
  #Ajustement de certains paramètres visuels
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.4, face = "bold", size = 20),
    plot.margin = margin(10, 30, 10, 10)
  )+
  
  #Aligner l'axe
  coord_cartesian(expan = FALSE)
  





# Script graphiques après-récifs 4a----
# ALLER CHERCHER LES PACKAGES NÉCESSAIRES
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)


# DÉFINITION DES PARAMÈTRES DU MODÈLE
#Différences: ajout de mortalité et changements états initiaux.

params <- list(
  
  #Taux de mortalité du backreef
  mortal <- 0.3,   
  
  #Paramètre de reproduction asexuée
  rb = 0.05,        #Reproduction asexuée de la population backreef
  rf = 0.15,        #Reproduction asecuée de la population forereef
  
  #Paramètre de reproduction sexuée
  sigma = 0.5,    #Intensité de la croissance larvaire 
  fecondity = 1,    #Fécondité
  zeta = 1,         #Probabilité de maturation larvaire
  
  #Coefficients d’interactions (γβ)
  gamma_beta_bb = 0.01,   #Influence de backreef sur backreef
  gamma_beta_fb = 0.002,  #Influence de forereef sur backreef
  gamma_beta_ff = 0.2,    #Influence de forereef sur forereef 
  gamma_beta_bf = 0.1,    #Influence de backreef sur forereef 
  
  #Paramètre d'immigration nette (apport externe constant)
  Ib = 0.08,     #Immigration du backreef (tiré du graphique 4c)
  If = 0.09,     #Immigration du forereef
  
  #Capacités de soutien en proportion
  Kb = 0.64,   #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)

#Conditions initiales du modèle
state <- c(pFallo = 0.01, pFauto = 0.01, pFtotal = 0.01, pBallo = 0.3, pBauto = 0.1, pBtotal = 0.4)

#Intervalle de temps en années (2011 à 2018)
times <- seq(0, 7, by = 0.1)





# CRÉATION DU MODÈLE DE CROISSANCE SELON LES ÉQUATIONS 4 DE L'ARTICLE + MORTALITÉ
model_proportion_backreef <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 1. Séparation de l'équation de la croissance de sous-population avant-récif en 3 équations:
    # Croissance du au recrutement autochtone et à la reproduction asexuée (dFauto/dt)
    dFauto  <- (pFtotal * (rf + (sigma * fecondity * zeta * gamma_beta_ff))) * (1 - pFtotal / Kf)
    
    # Croissance du au recrutement allochtone et l'immigration nette (dFallo/dt)
    dFallo  <- ((pBtotal * (sigma * fecondity * zeta * gamma_beta_bf)) + If) * (1 - pFtotal / Kf)
    
    # Croissance totale de l'avant-récif (dB/dt)
    dFtotal <- dFallo + dFauto
    
    
    #2. Sépration de l'équation de la croissance de sous-population après-récif en 3 équations:
    # Croissance du au recrutement autochtone et à la reproduction asexuée (dBauto/dt)
    dBauto  <- (pBtotal * (rb + (sigma * fecondity * zeta * gamma_beta_bb))) * (1 - pBtotal / Kb)
    
    # Croissance du au recrutement allochtone et l'immigration nette (dBallo/dt)
    dBallo  <- ((pFtotal * (sigma * fecondity * zeta * gamma_beta_fb)) + Ib) * (1 - pBtotal / Kb)
    
    # Croissance totale de l'après-récif (dB/dt)
    # Ajout de la mortalité ici
    dBtotal <- dBallo + dBauto - (mortal * pBtotal)
    
    
    list(c(dFallo, dFauto, dFtotal, dBallo, dBauto, dBtotal))
  })
}



# RÉSOLUTION NUMÉRIQUE DU MODÈLE
out <- ode(y = state, times = times, func = model_proportion_backreef, parms = params)
out <- as.data.frame(out)



# GRAPHIQUE
#Ajustement des données en pourcentage
out$pBtotal_pourcentage <- out$pBtotal*100

#graphique qui reproduit la figure 4a
ggplot(out, aes(x = time, y = pBtotal_pourcentage)) +
  geom_line(size = 1.7) +
  
  #Définition de l'axe des x
  scale_x_continuous(
    breaks = c(1,3,5,7),
    labels= c(2012, 2014, 2016, 2018)
  ) +
  
  #Définition de l'axe des y
  scale_y_continuous(
    limits=c(0,50),
    breaks = c(0, 20, 50)
  ) +
  
  #Définition des titres
  labs(
    title = "Couverture corallienne de la sous-population après-récif (Pb)",
    x = NULL,
    y = "Couverture corallienne (%)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  #Ajustement de certains paramètres visuels
  theme(
    legend.position = c(0.80, 0.95),
    legend.justification = c("left","top"),    
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(10, 30, 10, 10)
  ) +
  
  #Aligner l'axe des x avec le 0 de l'axe des y
  coord_cartesian(expand = FALSE)


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
# ALLER CHERCHER LES PACKAGES NÉCESSAIRES
library(ggplot2)
library(dplyr)
library(tidyr)

# DÉFINITION DES PARAMÈTRES UTILISÉS

#Paramètres de temps
years <- 2011:2018   #Années de l'étude
n <- length(years)

# Paramètres de taille de population
Pf <- seq(0.05, 0.6, length.out = n)  #Population forereef
Pb <- seq(0.2, 0.5, length.out = n)   #Population backreef
Pf[1] <- 0                            #Forcer à commencer à 0
Pb[1] <- 0                            #Forcer à commencer à 0

#Paramètre de capacité de soutien
Kb <- 0.8            #Capacité de charge de 80%

#Paramètres de reproduction
rb <- 0.1            #Reproduction asexuée
sigma_f <- 0.05      #Intensité de la croissance larvaire

#Coefficients d'interaction (γβ)
gamma_beta_bb <- seq(0.2, 0.7, length.out = n)  #backreef vers backreef
gamma_beta_fb <- rep(0.01, lenght.out = n)      #forereef vers backreef



# CALCUL DES COMPOSANTES dPb/dt ET AJUSTEMENT SELON LE FACTEUR DE DENSITÉ-DÉPENDANCE
#Facteur de densité-dépendance
density_factor <- (1 - Pb / Kb)

#Croissance asexuée et auto-recrutement (autochtone)
intrinsic <- Pb * (rb + sigma_f * gamma_beta_bb)
intrinsic_adj <- intrinsic * density_factor

#Échange local
exchange <- Pf * sigma_f * gamma_beta_fb
exchange_adj <- exchange * density_factor

#Immigration nette (allochtone)
immigration <- rep(0.0007, n)
immigration_adj <- immigration 

#Variation totale
total_change <- intrinsic_adj + exchange_adj + immigration_adj

#Estimation de la courbe de mortalité
mortalite <- -0.003 + 0.0022 * (years - min(years)) / ((years - min(years)) + 4)



# GRAPHIQUE
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
  
  #Définition du type de lignes utilisés
  scale_linetype_manual(values = c("Mortality" = "solid", "Immigration" = "dashed")) +
  
  #Définition de l'axe des y
  scale_y_continuous(
    limits = c(-0.4, 0.1),
    breaks = c(-0.4, -0.1, 0.1)
  )+
  
  #Définition des titres
  labs(
    title = "Composants contribuant aux changements de la croissance de la population dans le backreef",
    y = expression("%K / yr"),
    x = NULL,
    linetype = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  
  #Ajustement de certains paramètres visuels
  theme(
    legend.position = c(0.80, 0.20),
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.margin = margin(10, 30, 10, 10)
  )+
  
  #Aligner l'axe des x avec le 0 de l'axe des y
  coord_cartesian(expan=FALSE)





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
  mortal = 0.6,  #taux de mortalité constante du backreef
  
  rb = 0.05,   #Reproduction autochtone de la population backreef
  rf = 0.15,  #Reproduction autochtone de la population forereef
  sigma_f = 0.5,   #Intensité de la croissance larvaire (829 polypes/colonie)
  fecondity = 1,   #Fécondité
  zeta = 1,   #Probabilité de maturation larvaire
  
  #Coefficients d’interactions (γβ)
  gamma_beta_bb = 0.01,   #Influence de backreef sur backreef
  gamma_beta_fb = 0.002,    #Influence de forereef sur backreef
  gamma_beta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  gamma_beta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
  #Apports externes constants
  Ib = 0.08,   #Immigration du backreef
  If = 0.09,       #Immigration du forereef
  
  #Capacités de soutien en proportion (basé sur étendue max observée, tiré de l'article)
  Kb = 0.64,     #Capacité de soutien du backreef 0.64 dans article (observé)
  Kf = 0.8     #Capacité de soutien du forereef  0.8 dans article
)



model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    pertes_t <- taux_perte(time)
    
    dFallo  <- (sigma_f * gamma_beta_bf * fecondity * zeta * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * gamma_beta_ff * fecondity * zeta) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * gamma_beta_fb * fecondity * zeta * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * gamma_beta_bb * fecondity * zeta) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
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
  fecondity = 1,   #Fécondité
  zeta = 1,   #Probabilité de maturation larvaire
  
  #Coefficients d’interactions (γβ)
  gamma_beta_bb = 0.01,   #Influence de backreef sur backreef
  gamma_beta_fb = 0.002,    #Influence de forereef sur backreef
  gamma_beta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  gamma_beta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
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
    
    dFallo  <- (sigma_f * gamma_beta_bf * fecondity * zeta * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * gamma_beta_ff * fecondity * zeta) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * gamma_beta_fb * fecondity * zeta * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * gamma_beta_bb * fecondity * zeta) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
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
#Fonction qui représente la variation de mortalité aux 7 ans
#Fonction qui représente la variation de mortalité aux 7 ans
times <- seq(0,30, by=1)

#Pour une vrai chute de 99%, il faut peak=5.
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
  fecondity = 1,   #Fécondité
  zeta = 1,   #Probabilité de maturation larvaire
  
  #Coefficients d’interactions (γβ)
  gamma_beta_bb = 0.01,   #Influence de backreef sur backreef
  gamma_beta_fb = 0.002,    #Influence de forereef sur backreef
  gamma_beta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  gamma_beta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
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
    
    dFallo  <- (sigma_f * gamma_beta_bf * fecondity * zeta * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * gamma_beta_ff * fecondity * zeta) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * gamma_beta_fb * fecondity * zeta * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * gamma_beta_bb * fecondity * zeta) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
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
  fecondity = 1,   #Fécondité
  zeta = 1,   #Probabilité de croissance larvaire
  
  #Coefficients d’interactions (ζ(γβ))
  gamma_beta_bb = 0.01,   #Influence de backreef sur backreef
  gamma_beta_fb = 0.002,    #Influence de forereef sur backreef
  gamma_beta_ff = 0.2,    #Influence de forereef sur forereef AUTO
  gamma_beta_bf = 0.1,   #Influence de backreef sur forereef ALLO
  
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
    
    dFallo  <- (sigma_f * gamma_beta_bf * fecondity * zeta * pBtotal + If * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFallo))
    dFauto  <- ((pFauto + pFallo) * (rf + sigma_f * gamma_beta_ff * fecondity * zeta) * (1 - (pFauto + pFallo) / Kf)) - (pertes_t * (pFauto))
    
    dBallo  <- (sigma_f * gamma_beta_fb * fecondity * zeta * (pFauto + pFallo) + Ib * (1 - pBtotal / Kb)) - (pertes_t * pBallo)
    dBauto  <- (pBtotal * (rb + sigma_f * gamma_beta_bb * fecondity * zeta) * (1 - pBtotal / Kb)) - (pertes_t * pBauto)
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



# modele en r

model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dPb <-   Pb (Rb + omeg * Fe * Mp *(YBbb))   +    omeg * Fe * Mp * (YBfb) * Pf + Ib   *   (1- (Pb/Kb))
    dPf <-   Pb (Rf + omeg * Fe * Mp *(YBff))   +    omeg * Fe * Mp * (YBbf) * Pb + If   *   (1- (Pf/Kf))                                                        
                                                                               
                                                                               
                                                                               
  }

parameters <- c(
  Rb   =
  Rf   =
  Mp   = 0.5
  Fe   = fixe
  omeg = fixe
  Ib   =
  If   =
  Kb   = fixe
  Kf   = fixe
)

state <- c(Pb = 0.01, Pf = 0.01)
times <- seq(0, 50, by = 0.1) 


out <- ode(y = state, times = times, func = model, parms = parameters)
out <- as.data.frame(out)




#FIGURE 2
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
  scale_color_manual(values = c("Γ → ∞" = "blue", "Γ → 0" = "goldenrod", 
                                "Total" = "black", "Allochthonous" = "goldenrod", "Autochthonous" = "blue"))



#FIGURES 3c & 3d

#---------- Figure 3c ----------
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
rf <- 0.12   #Taux de croissance asexuée (estimée)
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

#Cération du dataframe pour le graphique
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


#---------- Figure 3d ----------
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
rf <- 0.085   #Taux de croissance asexuée (estimée)
sigma_f <- 0.05   #Taux de reproduction larvaire (estimé)

zeta_ff <- seq(0.2, 0.8, length.out = n)   #Autorecrutement (forereef vers forereef)
zeta_bf <- seq(0.01, 0.01, length.out = n)   #Recrutement d'échange (forereef vers backreef)
If <- c(0.005, 0.0065, 0.006, 0.006, 0.0058, 0.0053, 0.005, 0.002)   #Immigration externe (données modélisées)

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

#Cération du dataframe pour le graphique
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
    title = "Components contributing to changes in population growth (LTER 2 Forereef)",
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
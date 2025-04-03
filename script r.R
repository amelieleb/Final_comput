# FIGURE 2 DANS L'INTRODUCTION DE L'ARTICLE (FIGURE DE BASE)
# Installer et charger les bibliothèques nécessaires
install.packages("deSolve")
install.packages("ggplot2")
library(deSolve)
library(ggplot2)

# Modèle d'équations différentielles
model <- function(t, state, parameters) {
  P_auto <- state[1]  # Population autochtone
  P_allo <- state[2]  # Population allochthone
  K <- parameters$K
  I <- parameters$I
  Gamma <- parameters$Gamma
  
  # Normalisation
  P_auto_norm <- P_auto / K
  P_allo_norm <- P_allo / K
  
  # Équations différentielles
  dP_auto <- (P_auto_norm + P_allo_norm + I/Gamma) * (1 - (P_auto_norm + P_allo_norm))  
  dP_allo <- I * (1 - (P_auto_norm + P_allo_norm)) / Gamma  
  
  list(c(dP_auto, dP_allo))
}

# Paramètres du modèle
parameters <- list(K = 1, I = 1)  

# Conditions initiales
state <- c(P_auto = 0.01, P_allo = 0.01)

# Temps
times <- seq(0, 5, by = 0.1)

# Résoudre l'ODE pour différents Gamma
parameters$Gamma <- 10   # Cas intermédiaire
out_intermediate <- ode(y = state, times = times, func = model, parms = parameters)

parameters$Gamma <- 0.1   # Immigration dominante (Γ → 0)
out_immigration <- ode(y = state, times = times, func = model, parms = parameters)

parameters$Gamma <- 100   # Croissance intrinsèque dominante (Γ → ∞)
out_intrinsic <- ode(y = state, times = times, func = model, parms = parameters)

# Convertir les résultats en dataframes
df_intermediate <- as.data.frame(out_intermediate)
df_immigration <- as.data.frame(out_immigration)
df_intrinsic <- as.data.frame(out_intrinsic)

# Calculer la population totale
df_intermediate$Total <- df_intermediate$P_auto + df_intermediate$P_allo
df_intrinsic$Total <- df_intrinsic$P_auto + df_intrinsic$P_allo
df_immigration$Total <- df_immigration$P_auto + df_immigration$P_allo

# Tracer la figure avec inversion des couleurs des lignes pleines
ggplot() +
  # Courbes en pointillés pour les cas limites (pas modifiées)
  geom_line(data = df_intrinsic, aes(x = time, y = Total, color = "Γ → ∞"), linetype = "dashed", size = 1) +
  geom_line(data = df_immigration, aes(x = time, y = Total, color = "Γ → 0"), linetype = "dashed", size = 1) +
  
  # Courbes solides pour Γ = 10 
  geom_line(data = df_intermediate, aes(x = time, y = Total, color = "Total"), size = 1) +
  geom_line(data = df_intermediate, aes(x = time, y = P_allo, color = "Autochthonous"), size = 1) +  
  geom_line(data = df_intermediate, aes(x = time, y = P_auto, color = "Allochthonous"), size = 1) +  
  
  # Personnalisation avec couleurs corrigées
  labs(x = "Time (τ)", y = expression(P^"*"), color = "Population Type") +
  theme_minimal() +
  ggtitle("Nondimensional population recovery trajectories") +
  scale_color_manual(values = c("Γ → ∞" = "blue", "Γ → 0" = "goldenrod", 
                                "Total" = "black", "Allochthonous" = "goldenrod", "Autochthonous" = "blue"))
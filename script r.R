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
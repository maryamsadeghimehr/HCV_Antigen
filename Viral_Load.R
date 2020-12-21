Viral_Load <- function(bl, t)
{
#require(deSolve)


##############################################################
c = 6.3 # viral clearance rate
s = 10 ^ 5 # Target cell production rate
dT = 0.01 # Target cell death rate
beta = 10 ^ (-8) # infectivity
lambda = 1 # effector cell recruitment rate
mu = 2 # death rate of effector cells
dE = 4 # Antigen-induced exhaustion rate of effector cells
kD = 2.7 * 10 ^ 4 # hill function scalling for effector cell exhaustion
bE = 1 # Antigen-induced proliferation rate of effector cell
kB = 10 ^ 3 # hill function scaling for effector cell proliferation
##############################################################

func_VL <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    f= 0
    dT = 10^5 - 0.01 * T - 10 ^ (-8) * V * I
    dV = 4.4 * (1 - f) * I - 6.3 * V
    dI = 10 ^ (-8) * V * T - 0.01 * I - 0.244 * E * I
    dE = 1 + I * E / (I + 10 ^ 3) - 4 * I * E /((2.7 * 10 ^ 4) + I) - 2 * E
    list(c(dT, dV, dI, dE))
  })
}
##########################################################
#times      <- seq(0, 100, by = 0.3344482)
times = t
parameters <- c()

###############################################################

##### Initial points #########################################
#V0 = as.numeric(bl["Viral0"])
V0 = 34
T0 = s / (dT + beta * V0)
I0 = beta * V0 * T0 / s
E0 = lambda * (mu + dE * I0 /(kD + I0) - bE * I0 /(kB + I0))
##############################################################

state <- c(T = T0, V = V0, I = I0, E = E0)

out <- ode(y = state, times = times, func = func_VL, parms = parameters)
viral_load <- c * out[,3]/out[,4]
return(viral_load)
}


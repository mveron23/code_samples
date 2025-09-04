#' ---
#' title: "Final Project: Part II"
#' author: "Melanie Veron"
#' date: "February 9, 2018"
#' output: pdf_document
#' ---
#'
#' In this part of the project, the task is to build and simulate a version of
#' the Ross-MacDonald malaria model and examine the threshold for persistence of
#' malaria.
#' 
#' Formulae:
#' 
#' $$\frac {dS_{H}} {dt} = \mu N_{H} - \frac {\beta_{V} I_{V} S_{H}} {N_{H}} + 
#' \sigma I_{H} - \mu S_{H}$$
#' 
#' $$\frac {dI_{H}} {dt} = \frac {\beta_{V} I_{V} S_{H}} {N_{H}} - \sigma 
#' I_{H} - \mu I_{H}$$
#' 
#' $$\frac {dS_{V}} {dt} = c N_{V} - \frac {\beta_{H} S_{V} I_{H}} {N_{H}} - c 
#' S_{V}$$
#' 
#' $$\frac {dI_{V}} {dt} = \frac {\beta_{H} S_{V} I_{H}} {N_{H}} - c I_{V},$$
#' 
#' where $\mu$ is the host mortality rate, $\sigma$ is the host recovery rate, 
#' $\beta_{V} = b p_{HV}$ (where $b$ is the vector bite rate and $p_{HV}$ is the
#' probability of transmission from vector to host), $c$ is the vector mortality 
#' rate, and $\beta_{H} = b \frac {N_{V}} {N_{H}} p_{VH}$ (where $N_{V}$ and 
#' $N_{H}$ are the sizes of the vector and host populations and $p_{VH}$ is the 
#' probability of transmission from host to vector).

# load in R library to do integration and other sources files
source("epi-helper.R")
source("epi-ode-SIS.R")

# set up initial population parameters
pop.size.host <- 1000
pop.size.vector <- 10000
initial.infecteds.host <- 1
initial.infecteds.vector <- 1
initial.susceptibles.host <- pop.size.host - initial.infecteds.host
initial.susceptibles.vector <-
  pop.size.vector - initial.infecteds.vector

# set host and vector lifespans (L)
lifespan.host <- 50 * 52
lifespan.vector <- 1

# set host infectious period (D)
infectious.period <- 4

# set recovery rate (sigma)
sigma <- 1 / infectious.period

# set natural mortality rates for host (mu) and vector (c)
mu <- 1 / lifespan.host
c <- 1 / lifespan.vector

# set vector bite rate (b) and transmission probabilities from vector to
# host (pHV) and from host to vector (pVH)
bite.rate <- 0.5 * 7
pHV <- 0.5
pVH <- 0.8

# set vector and host transmission rates (beta)
beta.vector <- bite.rate * pHV
beta.host <- ((bite.rate * pop.size.vector) / pop.size.host) * pVH

# set simulation times
start.time <- 0
end.time <- 1
step <- 0.01

# set initial populations data frame
initial.populations <- data.frame(
  time = start.time,
  susceptibles.host = initial.susceptibles.host,
  infecteds.host = initial.infecteds.host,
  susceptibles.vector = initial.susceptibles.vector,
  infecteds.vector = initial.infecteds.vector
)

# integrate the differential equations
final.populations <- run.integration(
  ode.deterministic.SIS,
  initial.populations,
  end.time,
  timestep = step,
  beta.host = beta.host,
  beta.vector = beta.vector,
  mu = mu,
  c = c,
  sigma = sigma
)

# plot the results
plot.populations(final.populations)

#'
#' Formula to calculate $R_{0}$:
#'
#' $$R_{0} = \sqrt {\frac {b^{2} p_{HV} p_{VH}} {c (\mu + \sigma)} \frac {N_{V}}
#' {N_{H}}}$$

R0 <-
  sqrt((bite.rate ^ 2 * pHV * pVH) / (c * (mu + sigma)) * (pop.size.vector /
                                                             pop.size.host))
R0

R0 <- 0.99

pop.size.host <- 1000
pop.size.vector <-
  (((R0 ^ 2) * c * (mu + sigma)) / ((bite.rate ^ 2) * pHV * pVH)) * pop.size.host
pop.size.vector

pop.size.vector <- 10 ^ 4
pop.size.host <-
  (((bite.rate ^ 2) * pHV * pVH) / ((R0 ^ 2) * c * (mu + sigma))) * pop.size.vector
pop.size.host

#'
#' The $R_{0}$, assuming the original parameters, is 14. To eliminate malaria,
#' we would need to reduce mosquito density to around 50 mosquitos per square
#' kilometer. At a mosquito density of 10,000 per square kilometer, the human
#' population would need to increase to at least approximately 199,672
#' individuals per square kilometer to achieve an $R_{0}$ of less than 1.

#'
#' Scenario 1: $R_{0} < 1$ by Reducing the Mosquito Population Size
#' ----------------------------------------------------------------

# set up initial population parameters
pop.size.host <- 1000
pop.size.vector <-
  (((R0 ^ 2) * c * (mu + sigma)) / ((bite.rate ^ 2) * pHV * pVH)) *
  pop.size.host
initial.infecteds.host <- 1
initial.infecteds.vector <- 1
initial.susceptibles.host <- pop.size.host - initial.infecteds.host
initial.susceptibles.vector <-
  pop.size.vector - initial.infecteds.vector

# set host and vector lifespans (L)
lifespan.host <- 50 * 52
lifespan.vector <- 1

# set host infectious period (D)
infectious.period <- 4

# set recovery rate (sigma)
sigma <- 1 / infectious.period

# set natural mortality rates for host (mu) and vector (c)
mu <- 1 / lifespan.host
c <- 1 / lifespan.vector

# set vector bite rate (b) and transmission probabilities from vector to
# host (pHV) and from host to vector (pVH)
bite.rate <- 0.5 * 7
pHV <- 0.5
pVH <- 0.8

# set vector and host transmission rates (beta)
beta.vector <- bite.rate * pHV
beta.host <- ((bite.rate * pop.size.vector) / pop.size.host) * pVH

# set simulation times
start.time <- 0
end.time <- 1
step <- 0.01

# set initial populations data frame
initial.populations <- data.frame(
  time = start.time,
  susceptibles.host = initial.susceptibles.host,
  infecteds.host = initial.infecteds.host,
  susceptibles.vector = initial.susceptibles.vector,
  infecteds.vector = initial.infecteds.vector
)

# integrate the differential equations
final.populations <- run.integration(
  ode.deterministic.SIS,
  initial.populations,
  end.time,
  timestep = step,
  beta.host = beta.host,
  beta.vector = beta.vector,
  mu = mu,
  c = c,
  sigma = sigma
)

# plot the results
plot.populations(final.populations)

#'
#' When the mosquito population is brought down to below critical density (i.e.,
#' when $R_{0} < 1$), there are zero infected hosts or vectors by the end time.

#'
#' Scenario 2: $R_{0} < 1$ by Reducing the Human Population Size
#' -------------------------------------------------------------

# set up initial population parameters
pop.size.vector <- 10000
pop.size.host <-
  (((bite.rate ^ 2) * pHV * pVH) / ((R0 ^ 2) * c * (mu + sigma))) *
  pop.size.vector
initial.infecteds.host <- 1
initial.infecteds.vector <- 1
initial.susceptibles.host <- pop.size.host - initial.infecteds.host
initial.susceptibles.vector <-
  pop.size.vector - initial.infecteds.vector

# set host and vector lifespans (L)
lifespan.host <- 50 * 52
lifespan.vector <- 1

# set host infectious period (D)
infectious.period <- 4

# set recovery rate (sigma)
sigma <- 1 / infectious.period

# set natural mortality rates for host (mu) and vector (c)
mu <- 1 / lifespan.host
c <- 1 / lifespan.vector

# set vector bite rate (b) and transmission probabilities from vector to
# host (pHV) and from host to vector (pVH)
bite.rate <- 0.5 * 7
pHV <- 0.5
pVH <- 0.8

# set vector and host transmission rates (beta)
beta.vector <- bite.rate * pHV
beta.host <- ((bite.rate * pop.size.vector) / pop.size.host) * pVH

# set simulation times
start.time <- 0
end.time <- 1
step <- 0.01

# set initial populations data frame
initial.populations <- data.frame(
  time = start.time,
  susceptibles.host = initial.susceptibles.host,
  infecteds.host = initial.infecteds.host,
  susceptibles.vector = initial.susceptibles.vector,
  infecteds.vector = initial.infecteds.vector
)

# integrate the differential equations
final.populations <- run.integration(
  ode.deterministic.SIS,
  initial.populations,
  end.time,
  timestep = step,
  beta.host = beta.host,
  beta.vector = beta.vector,
  mu = mu,
  c = c,
  sigma = sigma
)

# plot the results
plot.populations(final.populations)

#' 
#' Again, endemic equilibrium fails when $R_{0} < 1$, this time due to bringing
#' the host population size above the critical density threshold. Keeping the 
#' host population high and the vector population low is key to controlling
#' vector-borne disease transmission.

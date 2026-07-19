library(dplyr)
library(deSolve)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales) # Required for label_comma() and label_number()

seir_full_model <- function(time, compartments, parameters) {
  with(as.list(c(compartments, parameters)), {
    N  <- S + E + I + R + V
    
    dS <- (par_mu * N) - (par_beta * S * I / N) - (par_mu * S) - (par_psi * S)
    dE <- (par_beta * S * I / N) - (par_sigma * E) - (par_mu * E)
    dI <- (par_sigma * E) - (par_gamma * I)        - (par_mu * I)
    dR <- (par_gamma * I) - (par_mu * R)
    dV <- (par_psi * S)   - (par_mu * V)
    dC <- par_rho * par_sigma * E  # cumulative reported cases
    
    list(c(dS, dE, dI, dR, dV, dC))
  })
}

# --- Shared ---
sigma_val <- 1 / 8     # latent rate σ  (8-day latent period)
gamma_val <- 1 / 5      # recovery rate γ (5-day infectious period)
rho_val   <- 0.15

# --- UK ---
N_uk      <- 66805830
R0_uk     <- 15
mu_uk     <- 1 / (82 * 365)     # vital dynamics: 82-yr life expectancy
# β derived from R0 for SEIR with demography:
#   R0 = β·σ / [(σ + μ)(γ + μ)]  →  β = R0·(σ+μ)·(γ+μ) / σ
beta_uk   <- R0_uk * (sigma_val + mu_uk) * (gamma_val + mu_uk) / sigma_val
psi_uk    <- 0.00001 # 0.001%

# --- Saudi Arabia ---
N_sa      <- 30906342
R0_sa     <- 13
mu_sa     <- 1 / (76 * 365)     # vital dynamics: 76-yr life expectancy
beta_sa   <- R0_sa * (sigma_val + mu_sa) * (gamma_val + mu_sa) / sigma_val
psi_sa    <- 0.00001 # 0.001%

# Package parameters for each country
params_uk <- c(par_beta  = beta_uk,  par_sigma = sigma_val, par_gamma = gamma_val,
               par_mu    = mu_uk,    par_psi   = psi_uk,    par_rho   = rho_val)

params_sa <- c(par_beta  = beta_sa,  par_sigma = sigma_val, par_gamma = gamma_val,
               par_mu    = mu_sa,    par_psi   = psi_sa,    par_rho   = rho_val)

run_seir <- function(N_pop, params, vac) {
  unvac <- 1 - vac
  compartments  <- c(S = unvac * N_pop - 1, E = 0, I = 1, R = vac * N_pop, V = 0, C = 0)
  times <- seq(0, 50 * 365, by = 1)   # ← 50 years
  
  ode(y      = compartments,
      times  = times,
      func   = seir_full_model,
      parms  = params,
      method = "lsoda") |>
    as.data.frame()
}

# vac = fraction with pre-existing immunity at model start
out_uk <- run_seir(N_uk, params_uk, vac = 0.89)
out_sa <- run_seir(N_sa, params_sa, vac = 0.98)

# =====================================================================
# INITIAL PLOT: Long-term Dynamics
# =====================================================================
# 1. Assign the plot to the variable 'fig_initial'
fig_initial <- bind_rows(
  out_uk |> mutate(country = "UK"),
  out_sa |> mutate(country = "Saudi Arabia")
) |>
  ggplot(aes(x = time / 365, y = I, colour = country)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~country, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = label_comma()) + 
  labs(x = "Years", y = "Infectious (I)", title = "Simulated Measles Dynamics") +
  theme_classic() +
  theme(
    # Explicitly enlarge text elements
    legend.text       = element_text(size = 13),
    axis.title        = element_text(size = 15, face = "bold"),
    axis.text =       element_text(size = 12),
    legend.position = "none",
    # Add grid lines back to classic theme
    panel.grid.major = element_line(colour = "#e5e5e5", linewidth = 0.5),
    panel.grid.minor = element_line(colour = "#f0f0f0", linewidth = 0.25)
  )
# 2. Print it to your screen
print(fig_initial)
# =====================================================================
# Save Initial Plot: Long-term Dynamics
ggsave(
  filename = "figure0_initial_dynamics.png", 
  plot = fig_initial, 
  width = 6, 
  height = 4,       # Slightly taller since it's a faceted plot (two stacked charts)
  dpi = 300,
  bg = "white"
)

# =====================================================================
# Post-process: derive year and daily new_cases from ODE output
# new_cases = dC/dt = rho * sigma * E
# =====================================================================
out_uk$year      <- out_uk$time / 365
out_uk$new_cases <- rho_val * sigma_val * out_uk$E

out_sa$year      <- out_sa$time / 365
out_sa$new_cases <- rho_val * sigma_val * out_sa$E

# =====================================================================
# Plot aesthetics
# =====================================================================
col_uk     <- "#1f77b4"   # blue  — United Kingdom
col_sa     <- "#2ca02c"   # green — Saudi Arabia
k_lab      <- label_number(scale = 1e-3, suffix = "k")
base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor  = element_blank(),
        legend.background = element_rect(fill = "white", colour = NA))

# =====================================================================
# FIGURE 1  — flip: UK 89% vs Saudi 98%  (imp = 0, 10 years)
# =====================================================================
out_uk$country <- "United Kingdom — 89% immune"
out_sa$country <- "Saudi Arabia — 98% immune"
flip <- rbind(out_uk, out_sa)

pk <- out_uk[which.max(out_uk$new_cases), ]

fig1 <- ggplot(flip, aes(year, new_cases, colour = country)) +
  geom_area(data = out_uk, aes(year, new_cases), fill = col_uk, alpha = 0.08, inherit.aes = FALSE) +
  geom_line(linewidth = 1.1) +
  
  # UK Peak Point
  annotate("point", x = pk$year, y = pk$new_cases, colour = col_uk, size = 2.5) +
  
  # FIXED UK LABEL: Shifted right (x + 1.5) and adjusted vertical alignment
  annotate("text", x = pk$year + 1.5, y = pk$new_cases,
           label = sprintf("~%s,000 cases/day\nat the peak", round(pk$new_cases / 1000)),
           colour = col_uk, hjust = 0, vjust = 0.7, fontface = "bold", size = 3.6) +
  
  # FIXED SA LABEL: Moved to empty upper space (x = 18, y = 7500)
  annotate("text", x = 18, y = 7500,
           label = "Initial outbreak prevented\n(measles dies out)",
           colour = col_sa, hjust = 0, fontface = "bold", size = 3.6) +
  
  # NEW: Arrow pointing from SA label down to the flat green line
  annotate("curve", x = 20, y = 6500, xend = 8, yend = 200,
           colour = col_sa, curvature = 0.2, arrow = arrow(length = unit(0.2, "cm")), alpha = 0.6) +
  
  scale_colour_manual(values = c("United Kingdom — 89% immune" = col_uk,
                                 "Saudi Arabia — 98% immune"   = col_sa)) +
  scale_y_continuous(labels = k_lab) +
  labs(x      = "Years after measles is introduced",
       y      = "New measles cases per day",
       colour = NULL) +
  base_theme + 
  theme(
    legend.position = c(0.78, 0.9),
    # Explicitly enlarge text elements
    legend.text       = element_text(size = 13),
    axis.title        = element_text(size = 15, face = "bold"),
    axis.text         = element_text(size = 12),
    # Force grid lines to render
    panel.grid.major = element_line(colour = "#e5e5e5", linewidth = 0.5),
    panel.grid.minor = element_line(colour = "#f0f0f0", linewidth = 0.25)
  )

print(fig1)
# Save Figure 1: Making the physical image file horizontally wider
ggsave(
  filename = "figure1_uk_vs_sa_measles.png", 
  plot = fig1, 
  width = 8,       # Increased from 10 to 14 to stretch it out horizontally
  height = 5,       # Slightly reduced height to emphasize the length
  dpi = 300,        
  bg = "white"      
)


# =====================================================================
# FIGURE 2: SA - Lower Initial Immunity, Higher psi, with HIT (First 2 Years)
# =====================================================================
# 1. Calculate Saudi Arabia's Herd Immunity Threshold (HIT)
hit_sa <- 1 - (1 / R0_sa) # 1 - 1/13 ≈ 0.923

# 2. Modify parameters: Higher continuous vaccination (psi)
# Increasing psi from 0.00001 to 0.002 (vaccinating 0.2% of susceptibles daily)
params_sa_high <- params_sa
params_sa_high["par_psi"] <- 0.002

# 3. Run model: Lower initial immunity (85% instead of 98%)
out_sa_new <- run_seir(N_sa, params_sa_high, vac = 0.85)

# 4. Post-process to get new cases and total immunity fraction
out_sa_new$year      <- out_sa_new$time / 365
out_sa_new$new_cases <- rho_val * sigma_val * out_sa_new$E
# Immunity fraction = (Recovered + Vaccinated) / Total Population
out_sa_new$immune_frac <- (out_sa_new$R + out_sa_new$V) / N_sa

# 5. Determine scaling coefficient for the dual y-axis
max_cases <- max(out_sa_new$new_cases)
coeff <- max_cases / 1.05 

# 6. Build the Plot
col_imm <- "#9467bd" # Purple for immunity

fig2 <- ggplot(out_sa_new, aes(x = year)) +
  # Case wave (Primary Y-axis)
  geom_area(aes(y = new_cases), fill = col_sa, alpha = 0.15) +
  geom_line(aes(y = new_cases, colour = "New Cases"), linewidth = 1) +
  
  # Immunity trajectory (Secondary Y-axis, scaled up by coeff)
  geom_line(aes(y = immune_frac * coeff, colour = "Immunity Level"), linewidth = 1.2) +
  
  # Herd Immunity Threshold line
  geom_hline(yintercept = hit_sa * coeff, linetype = "dashed", colour = "darkred", linewidth = 0.8) +
  
  # HIT Annotation (Moved to the right side, safely below the lines and right-aligned)
  annotate("text", x = 1.95, y = hit_sa * coeff - (0.02 * coeff), 
           label = sprintf("Herd Immunity Threshold (~%.1f%%)", hit_sa * 100), 
           colour = "darkred", hjust = 1, vjust = 1, fontface = "bold", size = 3.6) +
  
  # Scales and Aesthetics
  scale_colour_manual(values = c("New Cases" = col_sa, "Immunity Level" = col_imm)) +
  scale_x_continuous(breaks = seq(0, 2, by = 0.5)) + # Cleaner breaks for a 2-year window
  scale_y_continuous(
    name = "New measles cases per day",
    labels = k_lab,
    # Secondary axis scales the values back down by 'coeff' and formats as %
    sec.axis = sec_axis(~ . / coeff, name = "Population Immunity", labels = label_percent())
  ) +
  # Zoom in on the first 2 years without dropping data points that calculate the lines
  coord_cartesian(xlim = c(0, 2)) + 
  labs(
    x = "Years after measles is introduced",
    title = "Saudi Arabia: 85% Initial Immunity + Aggressive Vaccination Campaign",
    subtitle = "First 2 Years",
    colour = NULL
  ) +
  base_theme +
  theme(
    legend.position = "top",
    # Explicitly enlarge text elements
    legend.text       = element_text(size = 13),
    axis.title        = element_text(size = 15, face = "bold"),
    axis.text         = element_text(size = 12),
    # Color the right axis text to match the immunity line for clarity
    axis.title.y.right = element_text(colour = col_imm, face = "bold"),
    axis.text.y.right = element_text(colour = col_imm),
    # Force grid lines to render
    panel.grid.major = element_line(colour = "#e5e5e5", linewidth = 0.5), 
    panel.grid.minor = element_line(colour = "#f0f0f0", linewidth = 0.25)
  )

print(fig2)

ggsave(
  filename = "figure2_sa_hit_2years.png", 
  plot = fig2, 
  width = 7, 
  height = 4, 
  dpi = 300,
  bg = "white"
)

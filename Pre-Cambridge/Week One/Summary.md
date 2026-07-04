# Data Visualization and Statistical Inference in R

This document summarizes week one, which contains R programming, data visualization, and statistical modeling.

## Day 1: Getting Started, Data Types, and Plotting

### DA1: Getting started & Data types & structures
Programming allows you to move away from point-and-click activities, forcing you to think about the explicit steps in your analysis and greatly improving reproducibility. R is built on a variety of distinct data structures, with the simplest collection of data being a 1D vector consisting of elements of the exact same data type. Values are typically assigned to objects using the `<-` operator, and creating sequences or grouped collections of variables is handled by functions like `c()` and `seq()`. In real-world data, observations may be absent, which R represents as `NA`. When dealing with missing data, most numerical functions (like `mean()`) will return `NA` unless specifically instructed to ignore them using the `na.rm = TRUE` argument.

### DA2: Data & plotting
Tabular data frames are the most prevalent way to store information, with rows representing unique observations and columns containing variables. To import and handle this data robustly, R uses functions like `read_csv()` from the `tidyverse` suite. Once the data is loaded, R’s Grammar of Graphics (`ggplot2`) builds visualizations layer by layer. To create a plot, you must specify the dataset, map the variables inside the `aes()` helper function, and define the geometric shape using a `geom_*` function such as `geom_point()` or `geom_boxplot()`.

```R
# Day 1: Extensive R Example - Basics, I/O, and Plotting

# Load required libraries
library(tidyverse)

# --- DA1: Data types & structures ---
# Creating basic vectors and handling NA values
patient_ages <- c(23, 45, 67, 34, NA, 82)
mean_age <- mean(patient_ages, na.rm = TRUE)
print(mean_age)
> [1] 50.2

# Generating number sequences using seq()
even_numbers <- seq(2, 21, by = 2)
even_numbers
> [1]  2  4  6  8 10 12 14 16 18 20

# --- DA2: Working with tabular data ---
# Loading tabular data from a CSV file
surveys <- read_csv("~/Course_Materials/Week1_1_Intro_to_R/data/surveys.csv")

# Quick diagnostics of the tabular data structure
str(surveys)
> str(surveys)
spc_tbl_ [35,549 × 9] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
 $ record_id      : num [1:35549] 1 2 3 4 5 6 7 8 9 10 ...
 $ month          : num [1:35549] 7 7 7 7 7 7 7 7 7 7 ...
 $ day            : num [1:35549] 16 16 16 16 16 16 16 16 16 16 ...
 $ year           : num [1:35549] 1977 1977 1977 1977 1977 ...
 $ plot_id        : num [1:35549] 2 3 2 7 3 1 2 1 1 6 ...
 $ species_id     : chr [1:35549] "NL" "NL" "DM" "DM" ...
 $ sex            : chr [1:35549] "M" "M" "F" "M" ...
 $ hindfoot_length: num [1:35549] 32 33 37 36 35 14 NA 37 34 20 ...
 $ weight         : num [1:35549] NA NA NA NA NA NA NA NA NA NA ...

head(surveys)
> head(surveys)
# A tibble: 6 × 9
  record_id month   day  year plot_id species_id sex   hindfoot_length weight
      <dbl> <dbl> <dbl> <dbl>   <dbl> <chr>      <chr>           <dbl>  <dbl>
1         1     7    16  1977       2 NL         M                  32     NA
2         2     7    16  1977       3 NL         M                  33     NA
3         3     7    16  1977       2 DM         F                  37     NA
4         4     7    16  1977       7 DM         M                  36     NA
5         5     7    16  1977       3 DM         M                  35     NA
6         6     7    16  1977       1 PF         M                  14     NA

summary(surveys)
> summary(surveys)
   record_id         month             day             year         plot_id    
 Min.   :    1   Min.   : 1.000   Min.   : 1.00   Min.   :1977   Min.   : 1.0  
 1st Qu.: 8888   1st Qu.: 4.000   1st Qu.: 9.00   1st Qu.:1984   1st Qu.: 5.0  
 Median :17775   Median : 6.000   Median :16.00   Median :1990   Median :11.0  
 Mean   :17775   Mean   : 6.478   Mean   :15.99   Mean   :1990   Mean   :11.4  
 3rd Qu.:26662   3rd Qu.:10.000   3rd Qu.:23.00   3rd Qu.:1997   3rd Qu.:17.0  
 Max.   :35549   Max.   :12.000   Max.   :31.00   Max.   :2002   Max.   :24.0  
                                                                               
     species_id           sex        hindfoot_length     weight      
 Length   :35549   Length   :35549   Min.   : 2.00   Min.   :  4.00  
 N.unique :   48   N.unique :    2   1st Qu.:21.00   1st Qu.: 20.00  
 N.blank  :    0   N.blank  :    0   Median :32.00   Median : 37.00  
 Min.nchar:    2   Min.nchar:    1   Mean   :29.29   Mean   : 42.67  
 Max.nchar:    2   Max.nchar:    1   3rd Qu.:36.00   3rd Qu.: 48.00  
 NAs      :  763   NAs      : 2511   Max.   :70.00   Max.   :280.00  
                                     NAs    :4111    NAs    :3266  

# Counting unique values & handling missing data
num_species <- length(unique(surveys$species_id))
> num_species
[1] 49

missing_weights <- sum(is.na(surveys$weight))
> missing_weights
[1] 3266

# --- DA2: Plotting data ---
# Building a multi-layered ggplot
my_plot <- ggplot(data = surveys, mapping = aes(x = weight, y = hindfoot_length, colour = sex)) +
  geom_point(alpha = 0.5, size = 1.5) +  # Scatter points with transparency to handle overlap
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) + # Add linear regression line
  facet_wrap(vars(species_id)) + # Sub-panels for each species
  scale_colour_brewer(palette = "Dark2") + # Apply color-blind friendly palette
  labs(
    title = "Weight vs Hindfoot Length across Species",
    x = "Weight (grams)",
    y = "Hindfoot Length (mm)"
  ) +
  theme_minimal()

# Save the plot
ggsave("images/day1_survey_plot.png", plot = my_plot, width = 10, height = 8, dpi = 300)
```
![](Pre-Cambridge/Images/day1_survey_plot.png)
Species variation: Species like NL and DS are significantly heavier (stretching further right on the x-axis) compared to species like PF or RM.

Sexual dimorphism: By looking at the separation or overlap of the teal and orange lines/points within a single panel, you can tell if males and females of a specific species differ in size.

## Day 2: Manipulating, Reshaping, and Combining Data

### DA3: Manipulating columns & rows with chaining
As analyses grow, dealing with large datasets requires efficient filtering and transformation. Instead of creating numerous intermediate objects, R utilizes the pipe operator (`|>` or `%>%`), which chains functions together sequentially (read conceptually as "and then"). You can select specific variables with `select()`, filter observations fulfilling logical conditions (like `>`, `<=`, `==`) with `filter()`, and generate new columns utilizing `mutate()`. When investigating properties of specific subgroups within your dataset, you apply the "split-apply-combine" method by using `group_by()` followed by `summarise()` or `count()`.

### DA4: Organise, combine, and clean
Data arrives in varied shapes: "long" format (where each measured variable has its own column) is typically preferred for plotting, while "wide" format can be useful for cross-variable mathematical operations. You can shift between these structures using `pivot_longer()` and `pivot_wider()`. Joining disparate tables requires a shared identifier key, allowing operations like `left_join()`, `inner_join()`, or `full_join()` to safely aggregate data without manual copy-pasting. Real-world data often suffers from inconsistent naming and encoding, which can be remedied using functions like `clean_names()` from the `janitor` package and `case_when()` for programmatic text replacement.

```R
# Day 2: Extensive R Example - Data Wrangling Pipeline

library(tidyverse)
library(janitor)

# Load disparate datasets
infections <- read_csv("data/infections.csv")
hospital_info <- read_csv("data/hospital_info.csv")

# --- DA4: Cleaning, Styling & Arranging ---
# Clean column names to lowercase and underscores
infections <- clean_names(infections)

# Fix encoding issues in the data using case_when() and ensure booleans are logical
infections_cleaned <- infections |>
  mutate(
    # Recode inconsistent country names
    country = case_when(
      country %in% c("United kingdom", "U.K.", "UK") ~ "United Kingdom",
      TRUE ~ country
    ),
    # Force generic text flags to strict logical TRUE/FALSE
    icu_admission = as.logical(tolower(icu_admission) %in% c("yes", "y", "true"))
  )

# --- DA4: Combining data ---
# Left join hospital information to the infections table using the 'hospital' key
infections_joined <- left_join(infections_cleaned, hospital_info, by = "hospital")

# --- DA3: Chaining Operations & Grouped Operations ---
# Create an analytical pipeline: filter, group, summarize, and count
infections_summary <- infections_joined |>
  filter(!is.na(crp_level), age_group != "Unknown") |>  # Remove missing CRP and Unknown ages
  mutate(crp_log = log(crp_level)) |>                  # DA3: Create new column
  group_by(hospital_name, age_group, icu_admission) |> # DA3: Grouping
  summarise(
    mean_crp = mean(crp_level, na.rm = TRUE),
    max_crp  = max(crp_level, na.rm = TRUE),
    n_patients = n(),                                  # Count patients per group
    .groups = "drop"                                   # Always drop grouping when done
  ) |>
  arrange(desc(mean_crp))                              # DA3: Order rows by severity

# --- DA4: Reshaping data ---
# Pivot wider to compare ICU vs Non-ICU directly side-by-side
infections_wide <- infections_summary |>
  select(hospital_name, age_group, icu_admission, mean_crp) |>
  pivot_wider(
    names_from = icu_admission, 
    values_from = mean_crp,
    names_prefix = "icu_"
  )

print(head(infections_wide))
```

## Day 3: Statistical Inference & Categorical Predictors

### CS1: Statistical inference (One & Two Samples)
To determine if a sample mean significantly deviates from a known population value, we employ a one-sample t-test, assessing assumptions such as independence and a normal distribution using a `shapiro.test()` or interpreting Q-Q diagnostic plots. If the data is non-parametric (not normally distributed), the Wilcoxon signed-rank test is used to evaluate the median instead of the mean. When comparing two independent groups, we apply Student's t-test (or Welch’s t-test if variances are unequal, checked via Bartlett's or Levene's tests), formulating null hypotheses that the true group means are equal. If the two groups contain related data points (e.g., measurements taken from the same individual before and after an event), a paired t-test is required.

### CS2: Categorical predictors (ANOVA & Kruskal-Wallis)
When examining differences in continuous response variables spanning three or more groups, we utilize a one-way Analysis of Variance (ANOVA), which operates identically to a linear model under the hood. After constructing a linear model with `lm()`, diagnostic tools like `ggResidpanel` evaluate linearity, normality, and influential points via Cook's distance. Following a significant `anova()` result, post-hoc analysis via Tukey's Honest Significant Difference test (`tukey_hsd()`) isolates exactly which groups diverge. If the data fundamentally violates normality assumptions, the non-parametric Kruskal-Wallis test is leveraged instead, paired with Dunn's test for multiple pairwise comparisons.

```R
# Day 3: Extensive R Example - Statistical Inference

library(tidyverse)
library(rstatix)
library(ggResidpanel)

# --- CS1: One-sample & Two-sample tests ---
rivers <- read_csv("data/CS1-twosample.csv")

# Assess normality visually with QQ plot and mathematically with Shapiro-Wilk
rivers |>
  group_by(river) |>
  summarise(shapiro_p = shapiro.test(length)$p.value)

# Test for homogeneity of variance (Levene's test for non-normal tolerance)
rivers |> levene_test(length ~ river)

# Implement Student's Two-Sample t-test (assuming equal variance based on tests)
t_test_result <- rivers |>
  t_test(length ~ river, alternative = "two.sided", var.equal = TRUE)
print(t_test_result)

# If non-parametric was required, we'd use Mann-Whitney U:
# wilcox_test(length ~ river, alternative = "two.sided", data = rivers)

# --- CS2: ANOVA (Categorical Predictors > 2 groups) ---
oystercatcher <- read_csv("data/CS2-oystercatcher-feeding.csv")

# 1. Define the linear model
lm_oystercatcher <- lm(feeding ~ site, data = oystercatcher)

# 2. Check assumptions graphically
resid_panel(lm_oystercatcher, plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE)

# 3. Implement ANOVA
anova_results <- anova(lm_oystercatcher)
print(anova_results)

# 4. Implement Tukey's post-hoc test to find pair-wise significance
tukey_results <- tukey_hsd(lm_oystercatcher)
print(tukey_results)

# --- CS2: Kruskal-Wallis (Non-parametric ANOVA) ---
spidermonkey <- read_csv("data/CS2-spidermonkey.csv")

# Kruskal-Wallis test for non-normal continuous distributions across multiple groups
kw_result <- kruskal_test(aggression ~ familiarity, data = spidermonkey)
print(kw_result)

# Follow up with Dunn's test for non-parametric post-hoc comparisons
dunn_res <- dunn_test(aggression ~ familiarity, data = spidermonkey, p.adjust.method = "holm")
print(dunn_res)
```

## Day 4: Continuous Predictors & Two Predictors

### CS3: Continuous predictors
Correlation assesses the directional relationship between variables but does not imply causation, typically measured by Pearson's r (for linear relations) or Spearman's rank (robust to outliers). To test for an actual predictive capability, we use simple linear regression to calculate the intercept ($\beta_0$) and the slope ($\beta_1$) of the line of best fit, testing if the slope significantly differs from zero. Models must be vetted via diagnostic residual plots to ensure they are properly specified and relationships are truly linear.

### CS4: Two predictors
When systems are influenced by multiple variables (e.g., one categorical and one continuous, or two categorical factors), modeling requires two predictors. Here, it is vital to check for statistical *interactions* using the `*` operator in R (e.g., `y ~ categorical * continuous`), which investigates if the effect of one predictor intrinsically depends on the level of the other. For instance, visual lines of best fit that are not parallel heavily indicate a significant interaction effect.

```R
# Day 4: Extensive R Example - Continuous and Multiple Predictors

library(tidyverse)
library(rstatix)
library(ggResidpanel)
library(broom)

# --- CS3: Correlations & Simple Linear Regression ---
us_arrests <- read_csv("data/CS3-usarrests.csv")

# Compute Pearson's correlation matrix
cormat <- us_arrests |> 
  select_if(is.numeric) |> 
  cor_test(method = "pearson")
print(cormat)

# Build a simple linear regression: Does assault predict murder?
lm_arrests <- lm(murder ~ assault, data = us_arrests)
resid_panel(lm_arrests, plots = "all") # Check assumptions

# Determine significance using ANOVA table
anova(lm_arrests)

# Visualize the regression line
ggplot(us_arrests, aes(x = assault, y = murder)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_minimal()

# --- CS4: Two Predictors (Linear Regression with grouped data & Interaction) ---
tree_light <- read_csv("data/CS4-treelight.csv")

# Define a model featuring interaction: depth * species 
# Evaluates if light reduction by depth depends on the woodland species
lm_light_interaction <- lm(light ~ depth * species, data = tree_light)

# Check diagnostic plots before interpreting
resid_panel(lm_light_interaction, plots = c("resid", "qq", "ls", "cookd"))

# Run the ANOVA to check for significant interaction (depth:species term)
interaction_test <- anova(lm_light_interaction)
print(interaction_test)

# If interaction is NOT significant, fit an additive model (+) to calculate 
# lines with different intercepts but the same slope
lm_light_additive <- lm(light ~ depth + species, data = tree_light)
summary(lm_light_additive)

# Use broom::augment to easily append fitted predicted values for visualization
light_fitted <- augment(lm_light_additive, data = tree_light)

ggplot(light_fitted, aes(x = depth, y = light, color = species)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = .fitted), linewidth = 1.2) +
  labs(title = "Additive Linear Model of Light by Canopy Depth and Species")
```

## Day 5: Multiple Predictors, Model Comparisons, and Power Analysis

### CS5: Multiple predictors & Model comparisons
Expanding to three or more predictors often yields highly complex models involving multiple main effects and potentially multi-way interactions. To identify the "best" and most parsimonious model, we employ Backwards Stepwise Elimination utilizing the Akaike Information Criterion (AIC). By comparing the AIC scores of reduced models to a full model (dropping one term at a time), we eliminate non-significant variables—whenever a model drops in AIC or increases by less than 2, the simpler model is highly supported. R automates this rigorous simplification technique via the `step()` function.

### CS6: Statistical power
Because hypotheses testing risks returning false positives (Type I error) or false negatives (Type II error), rigorous experimental design utilizes power analysis. This requires a desired power threshold (often 80%), a significance level (often 0.05), and a pre-determined *effect size* measured biologically, often defined via Cohen's benchmarks (like Cohen's d for t-tests, or $f^2$ for linear models calculated from the $R^2$ variance). Incorporating these variables via the `pwr` package computes the precise minimal sample size needed, ensuring experiments are ethically and structurally sound.

```R
# Day 5: Extensive R Example - Model Selection and Power Analysis

library(tidyverse)
library(broom)
library(pwr)

# --- CS5: Multiple Linear Regression & Model Comparisons ---
pollution <- read_csv("data/CS5-pm2_5.csv")

# Define the full comprehensive model with an interaction
lm_full <- lm(pm2_5 ~ avg_temp + rain_mm + wind_m_s * location, data = pollution)
summary(lm_full)

# Perform Backwards Stepwise Elimination to find the minimal most predictive model
# The step() function iteratively removes factors checking for optimal AIC
lm_minimal <- step(lm_full, direction = "backward")

# Inspect the final model recommended by the AIC criterion
summary(lm_minimal)
anova(lm_minimal)

# --- CS6: Statistical Power Analysis ---
# Example 1: Planning a two-sample t-test
# We want to detect a "medium" effect size (Cohen's d = 0.5) at 80% power and 0.05 alpha
t_test_power <- pwr.t.test(
  d = 0.5, 
  sig.level = 0.05, 
  power = 0.8, 
  type = "two.sample", 
  alternative = "two.sided"
)
print(t_test_power) # Output 'n' denotes required sample size PER group

# Example 2: Retrospective power analysis on a Linear Model
lobsters <- read_csv("data/CS2-lobsters.csv")
lm_lobster <- lm(weight ~ diet, data = lobsters)

# Extract R-squared using broom's glance() function
r_sq <- glance(lm_lobster) |> pull(r.squared)

# Calculate Cohen's f^2 effect size for linear models (f2 = R2 / (1 - R2))
f2_effect <- r_sq / (1 - r_sq)

# Numerator Degrees of Freedom (u) = number of parameters - 1 
# Here, 3 diet groups - 1 = 2
u_df <- 2 

# Calculate the required denominator degrees of freedom (v) to achieve 80% power
lm_power <- pwr.f2.test(
  u = u_df, 
  f2 = f2_effect, 
  sig.level = 0.05, 
  power = 0.8
)
print(lm_power)

# Total required observations = u + v + 1
total_obs_needed <- ceiling(u_df + lm_power$v + 1)
message("Total observations needed for 80% power: ", total_obs_needed)
```

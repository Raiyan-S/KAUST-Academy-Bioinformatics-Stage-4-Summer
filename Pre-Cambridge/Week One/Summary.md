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
![Survey_Plot](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day1_survey_plot.png)
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

# Load the datasets
infections <- read_csv("~/Course_Materials/Week1_1_Intro_to_R/data/infections.csv")
hospital_info <- read_csv("~/Course_Materials/Week1_1_Intro_to_R/data/hospital_info.csv")
messy_data <- read_csv("~/Course_Materials/Week1_1_Intro_to_R/data/messy_data.csv")

# --- DA4: Cleaning, Styling & Arranging (messy_data) ---
# Clean column names to lowercase and underscores
messy_data
> messy_data
# A tibble: 100 × 8
   ID    Age         Gender Score Income.in.GBP country        employed.or.not notes
   <chr> <chr>       <chr>  <chr>         <dbl> <chr>          <chr>           <chr>
 1 id001 50          NA     3             36495 United kingdom FALSE           none 
 2 id002 34          f      1             31794 NA             y               NA   
 3 id003 33          f      4             37418 United Kingdom TRUE            none 
 4 id004 22          NA     2             33702 U.K.           no              ok   
 5 id005 twenty-five M      high          24622 United kingdom yes             NA   
 6 id006 unknown     NA     1             26982 NA             n               ok   
 7 id007 56          NA     high          55768 U.K.           n               NA   
 8 id008 33          NA     2             47478 United Kingdom FALSE           error
 9 id009 44          Female NA            23406 United Kingdom n               good 
10 id010 45          NA     4             62616 UK             y               ok 

messy_data <- clean_names(messy_data)
   id    age         gender score income_in_gbp country        employed_or_not notes
   <chr> <chr>       <chr>  <chr>         <dbl> <chr>          <chr>           <chr>
 1 id001 50          NA     3             36495 United kingdom FALSE           none 
 2 id002 34          f      1             31794 NA             y               NA   
 3 id003 33          f      4             37418 United Kingdom TRUE            none 
 4 id004 22          NA     2             33702 U.K.           no              ok   
 5 id005 twenty-five M      high          24622 United kingdom yes             NA   
 6 id006 unknown     NA     1             26982 NA             n               ok   
 7 id007 56          NA     high          55768 U.K.           n               NA   
 8 id008 33          NA     2             47478 United Kingdom FALSE           error
 9 id009 44          Female NA            23406 United Kingdom n               good 
10 id010 45          NA     4             62616 UK             y               ok   

# Fix encoding issues in the 'country' column using case_when()
messy_data_cleaned <- messy_data |>
  mutate(
    country = case_when(
      country %in% c("United kingdom", "U.K.", "UK") ~ "United Kingdom",
      TRUE ~ country
    )
  )
> messy_data_cleaned
# A tibble: 100 × 8
   id    age         gender score income_in_gbp country        employed_or_not notes
   <chr> <chr>       <chr>  <chr>         <dbl> <chr>          <chr>           <chr>
 1 id001 50          NA     3             36495 United Kingdom FALSE           none 
 2 id002 34          f      1             31794 NA             y               NA   
 3 id003 33          f      4             37418 United Kingdom TRUE            none 
 4 id004 22          NA     2             33702 United Kingdom no              ok   
 5 id005 twenty-five M      high          24622 United Kingdom yes             NA   
 6 id006 unknown     NA     1             26982 NA             n               ok   
 7 id007 56          NA     high          55768 United Kingdom n               NA   
 8 id008 33          NA     2             47478 United Kingdom FALSE           error
 9 id009 44          Female NA            23406 United Kingdom n               good 
10 id010 45          NA     4             62616 United Kingdom y               ok 

# --- DA4: Combining data (infections) ---
# Left join hospital information to the infections table using the 'hospital' key
infections
> infections
# A tibble: 1,400 × 11
   patient_id hospital   quarter infection_type vaccination_status age_group
   <chr>      <chr>      <chr>   <chr>          <chr>              <chr>    
 1 ID_0001    hospital_3 Q2      none           NA                 65+      
 2 ID_0002    hospital_3 Q2      viral          NA                 18 - 64  
 3 ID_0003    hospital_2 Q2      none           unknown            65+      
 4 ID_0004    hospital_2 Q3      fungal         unvaccinated       < 18     
 5 ID_0005    hospital_3 Q2      fungal         vaccinated         65+      
 6 ID_0006    hospital_5 Q3      none           vaccinated         65+      
 7 ID_0007    hospital_4 Q1      fungal         unvaccinated       18 - 64  
 8 ID_0008    hospital_1 Q1      NA             unvaccinated       18 - 64  
 9 ID_0009    hospital_2 Q1      viral          NA                 65+      
10 ID_0010    hospital_3 Q3      none           unvaccinated       NA

hospital_info
> hospital_info
# A tibble: 6 × 5
  hospital   hospital_name               location   bed_capacity teaching_hospital
  <chr>      <chr>                       <chr>             <dbl> <lgl>            
1 hospital_1 Royal London Hospital       London              784 TRUE             
2 hospital_2 Manchester General          Manchester          849 FALSE            
3 hospital_3 Bristol Royal Infirmary     Bristol             768 TRUE             
4 hospital_4 Edinburgh Medical Centre    Edinburgh           551 FALSE            
5 hospital_5 Cardiff University Hospital Cardiff             582 TRUE             
6 hospital_6 Leeds General Infirmary     Leeds               750 TRUE

infections_joined <- left_join(infections, hospital_info, by = "hospital")
> infections_joined
# A tibble: 1,400 × 15
   patient_id hospital   quarter infection_type vaccination_status age_group
   <chr>      <chr>      <chr>   <chr>          <chr>              <chr>    
 1 ID_0001    hospital_3 Q2      none           NA                 65+      
 2 ID_0002    hospital_3 Q2      viral          NA                 18 - 64  
 3 ID_0003    hospital_2 Q2      none           unknown            65+      
 4 ID_0004    hospital_2 Q3      fungal         unvaccinated       < 18     
 5 ID_0005    hospital_3 Q2      fungal         vaccinated         65+      
 6 ID_0006    hospital_5 Q3      none           vaccinated         65+      
 7 ID_0007    hospital_4 Q1      fungal         unvaccinated       18 - 64  
 8 ID_0008    hospital_1 Q1      NA             unvaccinated       18 - 64  
 9 ID_0009    hospital_2 Q1      viral          NA                 65+      
10 ID_0010    hospital_3 Q3      none           unvaccinated       NA 

# The joined columns are now 15 (11+5-1)
# --- DA3: Chaining Operations & Grouped Operations ---
# Create an analytical pipeline: filter, group, and summarize the infections data
infections_summary <- infections_joined |>
  filter(!is.na(crp_level)) |>                         # DA3: Filter out missing CRP values
  group_by(hospital_name, age_group, icu_admission) |> # DA3: Grouping by three variables
  summarise(
    mean_crp = mean(crp_level, na.rm = TRUE),          # Calculate average CRP
    max_crp  = max(crp_level, na.rm = TRUE),           # Calculate max CRP
    n_patients = n(),                                  # Count patients per group
    .groups = "drop"                                   # Always drop grouping when done
  ) |>
  arrange(desc(mean_crp))                              # DA3: Order rows by severity
infections_summary
> infections_summary
# A tibble: 67 × 6
   hospital_name               age_group icu_admission mean_crp max_crp n_patients
   <chr>                       <chr>     <lgl>            <dbl>   <dbl>      <int>
 1 Royal London Hospital       NA        FALSE             27.1    58.2          7
 2 Cardiff University Hospital NA        TRUE              27.0    46.1          3
 3 Edinburgh Medical Centre    18 - 64   NA                25.1    58.3          5
 4 Royal London Hospital       < 18      NA                24.9    26.6          3
 5 Cardiff University Hospital NA        NA                23.8    23.8          1
 6 Royal London Hospital       65+       TRUE              23.7    58.9         43
 7 Edinburgh Medical Centre    18 - 64   TRUE              23.5    53.1         17
 8 Manchester General          NA        TRUE              23.0    38.9          6
 9 NA                          < 18      TRUE              22.9    55.2          5
10 Cardiff University Hospital < 18      TRUE              22.8    55.8         21

# --- DA4: Reshaping data ---
# Pivot wider to compare ICU vs Non-ICU directly side-by-side
infections_wide <- infections_summary |>
  select(hospital_name, age_group, icu_admission, mean_crp) |>
  pivot_wider(
    names_from = icu_admission, 
    values_from = mean_crp,
    names_prefix = "icu_"
  )

# Preview the reshaped tabular data
print(head(infections_wide))
> print(head(infections_wide))
# A tibble: 6 × 5
  hospital_name               age_group icu_FALSE icu_TRUE icu_NA
  <chr>                       <chr>         <dbl>    <dbl>  <dbl>
1 Royal London Hospital       NA             27.1     16.6  20.5 
2 Cardiff University Hospital NA             14.5     27.0  23.8 
3 Edinburgh Medical Centre    18 - 64        21.2     23.5  25.1 
4 Royal London Hospital       < 18           20.1     18.8  24.9 
5 Royal London Hospital       65+            21.3     23.7  NA   
6 Manchester General          NA             18.4     23.0   8.28
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
rivers <- read_csv("~/Course_Materials/Week1_2_Core_Stats/data/CS1-twosample.csv")
rivers
> rivers
# A tibble: 68 × 2
   river   length
   <chr>    <dbl>
 1 Guanapo   19.1
 2 Guanapo   23.3
 3 Guanapo   18.2
 4 Guanapo   16.4
 5 Guanapo   19.7
 6 Guanapo   16.6
 7 Guanapo   17.5
 8 Guanapo   19.9
 9 Guanapo   19.1
10 Guanapo   18.8

# Assess normality visually with QQ plot and mathematically with Shapiro-Wilk
rivers |>
  group_by(river) |>
  summarise(shapiro_p = shapiro.test(length)$p.value)
# A tibble: 2 × 2
  river   shapiro_p
  <chr>       <dbl>
1 Aripo      0.0280
2 Guanapo    0.176 

# Test for homogeneity of variance (Levene's test for non-normal tolerance)
rivers |> levene_test(length ~ river)
# A tibble: 1 × 4
    df1   df2 statistic     p
  <int> <int>     <dbl> <dbl>
1     1    66      1.77 0.188

# Implement Student's Two-Sample t-test (assuming equal variance based on tests)
t_test_result <- rivers |>
  t_test(length ~ river, alternative = "two.sided", var.equal = TRUE)
print(t_test_result)
> print(t_test_result)
# A tibble: 1 × 8
  .y.    group1 group2     n1    n2 statistic    df        p
* <chr>  <chr>  <chr>   <int> <int>     <dbl> <dbl>    <dbl>
1 length Aripo  Guanapo    39    29      3.84    66 0.000275

# If non-parametric (non-normality) was required, we'd use Mann-Whitney U:
# wilcox_test(length ~ river, alternative = "two.sided", data = rivers)


# --- CS2: ANOVA (Categorical Predictors > 2 groups) ---
oystercatcher <- read_csv("data/CS2-oystercatcher-feeding.csv")
oystercatcher
> oystercatcher
# A tibble: 120 × 2
   site    feeding
   <chr>     <dbl>
 1 exposed    12.2
 2 exposed    13.1
 3 exposed    17.9
 4 exposed    13.9
 5 exposed    14.1
 6 exposed    18.4
 7 exposed    15.0
 8 exposed    10.3
 9 exposed    11.8
10 exposed    12.5

# 1. Define the linear model
lm_oystercatcher <- lm(feeding ~ site, data = oystercatcher)

# 2. Check assumptions graphically
resid_panel(lm_oystercatcher, plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE)
```
![Resid_Panel_Plot](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day3_oyster_plot.png)
```
# 3. Implement ANOVA
anova_results <- anova(lm_oystercatcher)
print(anova_results)
> print(anova_results)
Analysis of Variance Table

Response: feeding
           Df  Sum Sq Mean Sq F value    Pr(>F)    
site        2 1878.02  939.01  150.78 < 2.2e-16 ***
Residuals 117  728.63    6.23                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# 4. Implement Tukey's post-hoc test to find pair-wise significance
tukey_results <- tukey_hsd(lm_oystercatcher)
print(tukey_results)
> print(tukey_results)
# A tibble: 3 × 9
  term  group1  group2  null.value estimate conf.low conf.high    p.adj p.adj.signif
* <chr> <chr>   <chr>        <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
1 site  exposed partial          0     3.26     1.93      4.58 1.43e- 7 ****        
2 site  exposed shelte…          0     9.53     8.21     10.9  1.31e-14 ****        
3 site  partial shelte…          0     6.27     4.95      7.60 2.58e-14 ****  

# --- CS2: Kruskal-Wallis (Non-parametric ANOVA) [The dataset shows normality (I think), but we will assume not I guess]---
spidermonkey <- read_csv("data/CS2-spidermonkey.csv")
ggplot(spidermonkey, aes(sample = aggression)) +
  stat_qq() +                     # Plots the data points
  stat_qq_line(col = "red") +     # Adds the reference line
  theme_minimal() +               # Cleans the background visual
  labs(title = "ggplot2 Q-Q Plot")
```
![Kruskal_(Non)_Normality](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day3_kruskal.png)
```
spidermonkey |>
  group_by(familiarity) |>
  summarise(shapiro_p = shapiro.test(aggression)$p.value)
# A tibble: 3 × 2
  familiarity shapiro_p
  <chr>           <dbl>
1 high            0.429
2 low             0.264
3 none            0.915

# Kruskal-Wallis test for non-normal continuous distributions across multiple groups
kw_result <- kruskal_test(aggression ~ familiarity, data = spidermonkey)
print(kw_result)
> print(kw_result)
# A tibble: 1 × 6
  .y.            n statistic    df       p method        
* <chr>      <int>     <dbl> <int>   <dbl> <chr>         
1 aggression    21      13.6     2 0.00112 Kruskal-Wallis

# Follow up with Dunn's test for non-parametric post-hoc comparisons
dunn_res <- dunn_test(aggression ~ familiarity, data = spidermonkey, p.adjust.method = "holm")
print(dunn_res)
> print(dunn_res)
# A tibble: 3 × 9
  .y.        group1 group2    n1    n2 statistic        p    p.adj p.adj.signif
* <chr>      <chr>  <chr>  <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
1 aggression high   low        7     7      1.41 0.160    0.160    ns          
2 aggression high   none       7     7      3.66 0.000257 0.000771 ***         
3 aggression low    none       7     7      2.25 0.0245   0.0490   *  
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
us_arrests <- read_csv("~/Course_Materials/Week1_2_Core_Stats/data/CS3-usarrests.csv")

# Compute Pearson's correlation matrix
cormat <- us_arrests |> 
  select_if(is.numeric) |> 
  cor_test(method = "pearson")
print(cormat)
> print(cormat)
# A tibble: 16 × 8
   var1      var2        cor     statistic        p conf.low conf.high method 
   <chr>     <chr>     <dbl>         <dbl>    <dbl>    <dbl>     <dbl> <chr>  
 1 murder    murder     1          Inf     0          1          1     Pearson
 2 murder    assault    0.8          9.30  2.6 e-12   0.674      0.883 Pearson
 3 murder    urban_pop  0.07         0.483 6.31e- 1  -0.213      0.341 Pearson
 4 murder    robbery    0.56         4.73  2.03e- 5   0.338      0.728 Pearson
 5 assault   murder     0.8          9.30  2.6 e-12   0.674      0.883 Pearson
 6 assault   assault    1          Inf     0          1          1     Pearson
 7 assault   urban_pop  0.26         1.86  6.95e- 2  -0.0210     0.501 Pearson
 8 assault   robbery    0.67         6.17  1.36e- 7   0.475      0.796 Pearson
 9 urban_pop murder     0.07         0.483 6.31e- 1  -0.213      0.341 Pearson
10 urban_pop assault    0.26         1.86  6.95e- 2  -0.0210     0.501 Pearson
11 urban_pop urban_pop  1    464943848.    0          1          1     Pearson
12 urban_pop robbery    0.41         3.13  3   e- 3   0.150      0.619 Pearson
13 robbery   murder     0.56         4.73  2.03e- 5   0.338      0.728 Pearson
14 robbery   assault    0.67         6.17  1.36e- 7   0.475      0.796 Pearson
15 robbery   urban_pop  0.41         3.13  3   e- 3   0.150      0.619 Pearson
16 robbery   robbery    1          Inf     0          1          1     Pearson

# Build a simple linear regression: Does assault predict murder?
lm_arrests <- lm(murder ~ assault, data = us_arrests)
resid_panel(lm_arrests,
            plots = c("resid", "qq", "ls", "cookd"),
            smoother = TRUE) # Check assumptions
```
![Us_State_Plots](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day4_us_states_plot.png)
```
# Determine significance using ANOVA table
anova(lm_arrests)
> anova(lm_arrests)
Analysis of Variance Table

Response: murder
          Df Sum Sq Mean Sq F value    Pr(>F)    
assault    1 597.70  597.70  86.454 2.596e-12 ***
Residuals 48 331.85    6.91                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Visualize the regression line
ggplot(us_arrests, aes(x = assault, y = murder)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_minimal()
```
![Us_State_linear_Plots](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day4_us_states_linear_plot.png)
```
# --- CS4: Two Predictors (Linear Regression with grouped data & Interaction) ---
tree_light <- read_csv("~/Course_Materials/Week1_2_Core_Stats/data/CS4-treelight.csv")

# Define a model featuring interaction: depth * species 
# Evaluates if light reduction by depth depends on the woodland species
lm_light_interaction <- lm(light ~ depth * species, data = tree_light)

# Check diagnostic plots before interpreting
resid_panel(lm_light_interaction, plots = c("resid", "qq", "ls", "cookd"))
```
![Us_State_linear_Plots](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day4_light_plots.png)
```
# Run the ANOVA to check for significant interaction (depth:species term)
interaction_test <- anova(lm_light_interaction)
print(interaction_test)
> print(interaction_test)
Analysis of Variance Table

Response: light
              Df   Sum Sq  Mean Sq  F value    Pr(>F)    
depth          1 30812910 30812910 107.8154 2.861e-09 ***
species        1 51029543 51029543 178.5541 4.128e-11 ***
depth:species  1   218138   218138   0.7633    0.3932    
Residuals     19  5430069   285793                       
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# If interaction is NOT significant, fit an additive model (+) to calculate 
# lines with different intercepts but the same slope
lm_light_additive <- lm(light ~ depth + species, data = tree_light)
summary(lm_light_additive)
> summary(lm_light_additive)

Call:
lm(formula = light ~ depth + species, data = tree_light)

Residuals:
   Min     1Q Median     3Q    Max 
-842.4 -351.3 -216.3  318.1 1091.2 

Coefficients:
               Estimate Std. Error t value Pr(>|t|)    
(Intercept)     7962.03     231.36  34.415  < 2e-16 ***
depth           -262.17      39.92  -6.567 2.13e-06 ***
speciesConifer -3113.03     231.59 -13.442 1.78e-11 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 531.4 on 20 degrees of freedom
Multiple R-squared:  0.9354,	Adjusted R-squared:  0.929 
F-statistic: 144.9 on 2 and 20 DF,  p-value: 1.257e-12

# Use broom::augment to easily append fitted predicted values for visualization
light_fitted <- augment(lm_light_additive, data = tree_light)

ggplot(light_fitted, aes(x = depth, y = light, color = species)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = .fitted), linewidth = 1.2) +
  labs(title = "Additive Linear Model of Light by Canopy Depth and Species")
```
![Us_State_linear_Plots](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/day4_light_linear_plot.png)

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
> summary(lm_full)

Call:
lm(formula = pm2_5 ~ avg_temp + rain_mm + wind_m_s * location, 
    data = pollution)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.1247 -0.6800  0.0212  0.6520  3.2281 

Coefficients:
                        Estimate Std. Error t value Pr(>|t|)    
(Intercept)            18.182858   0.171992 105.719   <2e-16 ***
avg_temp                0.010451   0.007971   1.311    0.190    
rain_mm                -0.027880   0.048229  -0.578    0.563    
wind_m_s               -0.285450   0.026668 -10.704   <2e-16 ***
locationouter          -2.070843   0.190643 -10.862   <2e-16 ***
wind_m_s:locationouter -0.429455   0.037886 -11.335   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.024 on 724 degrees of freedom
Multiple R-squared:  0.8389,	Adjusted R-squared:  0.8378 
F-statistic: 754.1 on 5 and 724 DF,  p-value: < 2.2e-16

# Perform Backwards Stepwise Elimination to find the minimal most predictive model
# The step() function iteratively removes factors checking for optimal AIC
lm_minimal <- step(lm_full, direction = "backward")

# Inspect the final model recommended by the AIC criterion
summary(lm_minimal)
> summary(lm_minimal)

Call:
lm(formula = pm2_5 ~ wind_m_s + location + wind_m_s:location, 
    data = pollution)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.1884 -0.6896  0.0084  0.6548  3.2592 

Coefficients:
                       Estimate Std. Error t value Pr(>|t|)    
(Intercept)            18.24221    0.13259  137.58   <2e-16 ***
wind_m_s               -0.28509    0.02661  -10.71   <2e-16 ***
locationouter          -2.05975    0.19027  -10.83   <2e-16 ***
wind_m_s:locationouter -0.43182    0.03780  -11.43   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.024 on 726 degrees of freedom
Multiple R-squared:  0.8385,	Adjusted R-squared:  0.8378 
F-statistic:  1256 on 3 and 726 DF,  p-value: < 2.2e-16

anova(lm_minimal)
> anova(lm_minimal)
Analysis of Variance Table

Response: pm2_5
                   Df  Sum Sq Mean Sq F value    Pr(>F)    
wind_m_s            1  821.35  821.35  782.79 < 2.2e-16 ***
location            1 2995.68 2995.68 2855.01 < 2.2e-16 ***
wind_m_s:location   1  136.95  136.95  130.52 < 2.2e-16 ***
Residuals         726  761.77    1.05                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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
> print(t_test_power) # Output 'n' denotes required sample size PER group

     Two-sample t test power calculation 

              n = 63.76561
              d = 0.5
      sig.level = 0.05
          power = 0.8
    alternative = two.sided

NOTE: n is number in *each* group

# Example 2: Retrospective power analysis on a Linear Model
lobsters <- read_csv("data/CS2-lobsters.csv")
lm_lobster <- lm(weight ~ diet, data = lobsters)

# Extract R-squared using broom's glance() function
r_sq <- glance(lm_lobster) |> pull(r.squared)
> r_sq
[1] 0.1797219

# Calculate Cohen's f^2 effect size for linear models (f2 = R2 / (1 - R2))
f2_effect <- r_sq / (1 - r_sq)
> f2_effect
[1] 0.2190987

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
> print(lm_power)

     Multiple regression power calculation 

              u = 2
              v = 44.10317
             f2 = 0.2190987
      sig.level = 0.05
          power = 0.8

# Total required observations = u + v + 1
total_obs_needed <- ceiling(u_df + lm_power$v + 1)
message("Total observations needed for 80% power: ", total_obs_needed)
# Total observations needed for 80% power: 48
```
# URLs
- https://cambiotraining.github.io/data-analysis-in-r-and-python/
- https://cambiotraining.github.io/corestats/

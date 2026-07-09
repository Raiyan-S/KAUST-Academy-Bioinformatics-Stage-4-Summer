## Statistical Test Selection Guide

> **Always rely on plots first when checking assumptions, then use formal statistical tests to support your conclusions.**

### Assumption Checks

- **Normality:** Use a **QQ plot** (visual) and the **Shapiro–Wilk test** (formal).
- **Equal variance:** Use **Bartlett's test** when the data are approximately normal. If normality is questionable, use **Levene's test**, since it is more robust.

---

### One-sample data

**Used to compare the mean (or median) of a single sample to a known or hypothesized value.**

- Check normality.
- If the data are normally distributed, use a **one-sample t-test** (`t.test`) with `mu = ...`.
- Otherwise, use a **one-sample Wilcoxon signed-rank test** (`wilcox.test`) with `mu = ...`.

---

### Two independent samples

**Used to compare the means (or distributions) of two independent groups.**

- Check normality for both groups.
- If both groups are approximately normal:
  - Check equal variance:
    - **Bartlett's test** if normality is satisfied.
    - **Levene's test** if normality is questionable.
  - If the variance test has **p > 0.05**, use **Student's t-test** with `var.equal = TRUE`.
  - If the variance test has **p ≤ 0.05**, use **Welch's t-test** (`t.test`, default `var.equal = FALSE`).
- If the data are not approximately normal, use the **Wilcoxon rank-sum test (Mann–Whitney U test)** with `alternative = "two.sided"`.

---

### Paired samples

**Used to compare two related measurements from the same subjects or matched pairs (e.g., before vs. after treatment).**

- Check the normality of the **differences** between paired observations (not each sample separately).
- No Bartlett or Levene test is needed.
- If the differences are approximately normal, use a **paired t-test** with `paired = TRUE`.
- Otherwise, use a **paired Wilcoxon signed-rank test** with `paired = TRUE` and `alternative = "two.sided"`.

---

### One-way ANOVA

**Used when comparing the means of three or more independent groups for one continuous response variable and one categorical predictor.**

- Check normality of the **model residuals**:
  - Fit an `lm` model.
  - Inspect diagnostic plots (e.g., `resid_panel()`).
  - Perform the **Shapiro–Wilk test** on the residuals.
- Check equal variance:
  - **Bartlett's test** if residuals are approximately normal.
  - **Levene's test** if normality is questionable.
- If residuals are approximately normal and variances are equal (**p > 0.05**), use **one-way ANOVA** (`anova(lm)` or `aov()`).
- If residuals are approximately normal but variances are unequal (**p ≤ 0.05**), use **Welch's ANOVA** (`oneway.test(..., var.equal = FALSE)`).
- If the residuals are not approximately normal, use the **Kruskal–Wallis test**.

**Post-hoc tests**

- ANOVA only tells you whether **at least one group mean differs**.
- If ANOVA is significant:
  - Use **Tukey's HSD** (`tukey_hsd()` or `TukeyHSD()`) after standard one-way ANOVA.
  - **Do not use Tukey after Welch's ANOVA; use Games–Howell instead.**

---

### Kruskal–Wallis Test

**Used when comparing three or more independent groups but the assumptions for one-way ANOVA (especially normality) are not met.**

- Use `kruskal.test()`.
- If the result is significant, perform **Dunn's test** for pairwise comparisons (the non-parametric equivalent of Tukey's HSD).

---

### Linear Regression with Grouped Data

**Used when the response variable is continuous and there is at least one continuous predictor, with optional categorical predictors (groups).**

- Determines how continuous predictors affect the response while accounting for group differences.
- Interaction terms (e.g., `x * group`) can be included to test whether the relationship (slope) differs between groups.
- Can also be used for prediction.

---

### Two-way ANOVA

**Used when the response variable is continuous and there are two categorical predictors (factors).**

Tests:

- The main effect of Factor A.
- The main effect of Factor B.
- The interaction between Factors A and B (whether the effect of one factor depends on the level of the other factor).

Example:

- Plant height (continuous response) measured across different fertilizer types and genotypes (two categorical predictors).

Assumptions:

- Normal residuals.
- Equal variances.

---

### Choosing Between Linear Regression and Two-way ANOVA

- Use **Linear Regression** when **at least one predictor is continuous**.
- Use **Two-way ANOVA** when **all predictors are categorical**.
- Both are special cases of the general linear model (`lm()` in R).

---

### AIC (Akaike Information Criterion)

**Used to compare competing statistical models fitted to the same dataset.**

- Balances model fit and model complexity (number of parameters).
- Lower AIC values indicate a better trade-off between goodness of fit and simplicity.
- AIC does **not** test statistical significance; it is only used for model comparison.

General interpretation:

- **ΔAIC < 2:** Models have similar support.
- **ΔAIC between 2 and 7:** The model with the lower AIC has noticeably stronger support.
- **ΔAIC > 10:** The model with the lower AIC is strongly preferred.

# Quick Summary
### Decision Flow
### Decision Flow

```text
Start
│
├── One sample?
│   ├── Check normality
│   │   ├── Normal → One-sample t-test
│   │   └── Not normal → One-sample Wilcoxon signed-rank test
│
├── Two independent samples?
│   ├── Check normality
│   │   ├── Normal
│   │   │   ├── Equal variance?
│   │   │   │   ├── Yes → Student's t-test
│   │   │   │   └── No  → Welch's t-test
│   │   └── Not normal → Wilcoxon rank-sum (Mann–Whitney U) test
│
├── Paired samples?
│   ├── Check normality of the differences
│   │   ├── Normal → Paired t-test
│   │   └── Not normal → Paired Wilcoxon signed-rank test
│
├── Three or more independent groups?
│   ├── Check residual normality
│   │   ├── Normal
│   │   │   ├── Equal variance?
│   │   │   │   ├── Yes → One-way ANOVA
│   │   │   │   │          └── Significant? → Tukey's HSD
│   │   │   │   └── No  → Welch's ANOVA
│   │   │   │              └── Significant? → Games–Howell
│   │   └── Not normal → Kruskal–Wallis
│   │                     └── Significant? → Dunn's test
│
├── Continuous response with predictors?
│   ├── At least one predictor is continuous
│   │   └── Linear Regression (lm)
│   │       ├── Optional categorical predictors
│   │       ├── Optional interaction terms
│   │       └── Can be used for prediction
│   │
│   └── All predictors are categorical
│       └── Two-way ANOVA
│           ├── Tests main effect of Factor A
│           ├── Tests main effect of Factor B
│           └── Tests interaction (A × B)
│
└── Comparing multiple candidate models fitted to the same data?
    └── Compare AIC values
        └── Choose the model with the lowest AIC
```

# Libraries
library(readr)
library(skimr)
library(ggplot2)
library(tidyverse)
library(fastDummies)
library(corrplot)
library(survival)
library(ggfortify)
library(caret)
library(MASS)
library(riskRegression)
library(survminer)
library(gridExtra)
library(randomForestSRC)
library(pracma)
library(mice)
library(matrixStats)
library(splines)
library(glmnet)

# I. Understanding the context of the problem and dataset
df <- read_tsv("brca_metabric_clinical_data.tsv")
head(df)
summary(df)
skim(df)
glimpse(df)

# II. Exploratory data analysis

# 2.
# Numerical summaries
# Continuous variables
df %>%
  summarise(
    n = n (),
    mean_age = mean(`Age at Diagnosis`, na.rm = TRUE),
    medidan_tumor_size = median(`Tumor Size`, na.rm = TRUE),
    mean_overall_survival_months = mean(`Overall Survival (Months)`, na.rm = TRUE),
    mean_rfs_months = mean(`Relapse Free Status (Months)`, na.rm = TRUE)
  )

df %>%
  summarise(
    `Number of patients` = n(),
    `Mean age` = mean(`Age at Diagnosis`, na.rm = TRUE),
    `Median tumor size` = median(`Tumor Size`, na.rm = TRUE),
    `Mean overall survival (months)` = mean(`Overall Survival (Months)`, na.rm = TRUE),
    `Mean relapse-free survival (months)` = mean(`Relapse Free Status (Months)`, na.rm = TRUE)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# Frequency tables for categorical variables
table(df$`ER Status`)
table(df$`Tumor Stage`)

# Graphical Summaries
# Age Distribution
ggplot(df, aes(x = `Age at Diagnosis`)) +
  geom_histogram(binwidth = 5, fill = "grey70", color = "black") +
  labs(x = "Age at Diagnosis (years)", y = "Count") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Tumor Size Distribution
ggplot(df, aes(x = `Tumor Size`)) +
  geom_histogram(binwidth = 5, fill = "grey70", color = "black") +
  labs(x = 'Tumor Size (mm)', y = "Count") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Tumor Stage
ggplot(df, aes(x = `Tumor Stage`, fill = factor(`Overall Survival Status`))) +
  # geom_bar(fill = "darkorange") +
  geom_bar() +
  labs(x = "Stage", y = "Count")


df_long <- df %>%
  pivot_longer(
    cols = c(13, 16, 28),
    names_to = "Marker",
    values_to = "Status"
  )

ggplot(df_long, aes(x = Status)) +
  geom_bar(fill = "grey70", color = "black", na.rm = TRUE) +
  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +
  labs(
    x = NULL,
    y = "Number of Patients"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# 3. Explore dependence and correlations between variables
# Transform character variables to factor variables (categorical)
exclude <- c("Study ID", "Patient ID", "Sample ID")

# Transform character variables to factor variables
df[] <- lapply(names(df), function(col) {
  if (is.character(df[[col]]) && !(col %in% exclude)) {
    as.factor(df[[col]])
  } else {
    df[[col]]
  }
})

# drop columns with only 1 factor level
one_level_col <- names(df)[sapply(df, function(x) length(unique(na.omit(x))) == 1)]
df <- df[, !(names(df) %in% one_level_col)]

binary_var <- names(df)[sapply(df, function(x) length(unique(na.omit(x))) == 2)]
binary_var <- setdiff(binary_var, c("Overall Survival Status", "Relapse Free Status"))
ordinal_var <- c("Cellularity", "Neoplasm Histologic Grade", "HER2 status measured by SNP6", "Tumor Stage")

# Preserving order for binary variables
df$`Type of Breast Surgery` <- factor(df$`Type of Breast Surgery`, levels = c("BREAST CONSERVING", "MASTECTOMY"))
df$Chemotherapy  <- factor(df$Chemotherapy , levels = c("NO", "YES"))
df$`ER status measured by IHC` <- factor(df$`ER status measured by IHC`, levels = c("Negative", "Positve"))
df$`ER Status` <- factor(df$`ER Status` , levels = c("Negative", "Positive"))
df$`HER2 Status` <- factor(df$`HER2 Status` , levels = c("Negative", "Positive"))
df$`Hormone Therapy` <- factor(df$`Hormone Therapy` , levels = c("NO", "YES"))
df$`Inferred Menopausal State` <- factor(df$`Inferred Menopausal State` , levels = c("Pre", "Post"))
df$`Primary Tumor Laterality` <- factor(df$`Primary Tumor Laterality` , levels = c("Left", "Right"))
df$`PR Status` <- factor(df$`PR Status` , levels = c("Negative", "Positive"))
df$`Radio Therapy` <- factor(df$`Radio Therapy` , levels = c("NO", "YES"))

# Encode the event variables to 0/1
df$`Overall Survival Status` <- as.numeric(factor(df$`Overall Survival Status`, levels = c("0:LIVING", "1:DECEASED"))) - 1
df$`Relapse Free Status` <- as.numeric(factor(df$`Relapse Free Status` , levels = c("0:Not Recurred", "1:Recurred"))) - 1

# Preserving order for ordinal variables
df$Cellularity <- factor(df$Cellularity, levels = c("Low", "Moderate", "High"))
df$`Neoplasm Histologic Grade` <- factor(df$`Neoplasm Histologic Grade`, levels = c(1, 2, 3))
# For `HER2 status measured by SNP6` var, make UNDEFINED value to NA
df$`HER2 status measured by SNP6` <- ifelse(df$`HER2 status measured by SNP6` == "UNDEF", NA, as.character(df$`HER2 status measured by SNP6`))
df$`HER2 status measured by SNP6` <- factor(df$`HER2 status measured by SNP6`, levels = c("LOSS", "NEUTRAL", "GAIN"))
df$`Tumor Stage` <- factor(df$`Tumor Stage`, levels = c(0, 1, 2, 3, 4))

# Make a copy of datafram df
df_corr <- df

# Make binary and ordinal variables to numeric
for (col in binary_var) {
  df_corr[[col]] <- as.numeric(df_corr[[col]]) - 1
}

for (col in ordinal_var) {
  df_corr[[col]] <- as.numeric(df_corr[[col]])
}

# Transform norminal variables to dummies variables
col_to_exclude <- c("Patient ID", "Sample ID", "Relapse Free Status (Months)",
                    "Overall Survival (Months)", "Overall Survival Status",
                    "Relapse Free Status", "Patient's Vital Status", binary_var, ordinal_var)
df_dummy_norminal <-  dummy_cols(df_corr[, !(names(df_corr) %in% col_to_exclude)],
                remove_selected_columns = TRUE, ignore_na = TRUE)

df_binary <- df_corr[, binary_var]
df_ordinal <- df_corr[, ordinal_var]

df_dummy <- cbind(df_binary, df_ordinal, df_dummy_norminal)

# Correlation matrix using spearman method
cor_mat <- cor(df_dummy, use = 'pairwise.complete.obs', method = 'spearman')
corrplot(cor_mat, method = 'color', tl.cex = 0.38, tl.srt = 45)

# Find strong correlations
strong_corr_idx <- which(abs(cor_mat) > 0.6 & abs(cor_mat) < 0.99, arr.ind = TRUE)

# Create a data frame with variable names and correlation values
strong_corr <- data.frame(
  Var1 = rownames(cor_mat)[strong_corr_idx[,1]],
  Var2 = colnames(cor_mat)[strong_corr_idx[,2]],
  Correlation = cor_mat[strong_corr_idx]
)

# remove duplicates (Var1-Var2 and Var2-Var1)
strong_corr <- strong_corr[!duplicated(t(apply(strong_corr[,1:2], 1, sort))), ]


# 4. Graphically represent the survival functions in the subgroups defined by the categorical variables.
# For Tumor Stage
KM_fit_OS <- survfit(Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`) ~ df$`Tumor Stage`)
autoplot(KM_fit_OS, ylab = "Overall Survival (%)", xlab = "Time (Months)")

# For Type of Breast Surgery
KM_fit_OS <- survfit(Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`) ~ df$`Type of Breast Surgery`)
autoplot(KM_fit_OS, ylab = "Overall Survival (%)", xlab = "Time (Months)")

# For ER Status
KM_fit_OS <- survfit(Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`) ~ df$`ER Status`)
autoplot(KM_fit_OS, ylab = "Overall Survival (%)", xlab = "Time (Months)")

# For HER2 Status
KM_fit_OS <- survfit(Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`) ~ df$`HER2 Status`)
autoplot(KM_fit_OS, ylab = "Overall Survival (%)", xlab = "Time (Months)")

# 5. Log Rank test for each categorical variable
cat_var <- names(df)[sapply(df, is.factor)]
cat_var <- setdiff(cat_var, "Patient's Vital Status")
test_results <- data.frame(variable = character(), OS_p_value = numeric())

for (var in cat_var){
  # OS test
  os_test <- survdiff(Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`) ~ df[[var]])
  os_p <- os_test$pval

  test_results <- rbind(test_results, data.frame(
    variable = var,
    OS_p_value = os_p
    # RFS_p_value = rfs_p
  ))
}

test_results[order(test_results$OS_p_value), ]


# III. Statistical modelling and analysis
# Handling missing values
# Remove rows of outcome variables that have missing values
df_clean <- df[!is.na(df$`Overall Survival (Months)`) & !is.na(df$`Overall Survival Status`)
                  & !is.na(df$`Relapse Free Status (Months)`) & !is.na(df$`Relapse Free Status`), ]

# Calculate missing percentage per column
missing_pct <- colSums(is.na(df_clean)) / nrow(df_clean) * 100

# Thresholds
median_mode_thresh <- 5    # <5% missing
mice_thresh <- 30          # 5–30% missing

# Median/mode imputation for <5% missing
for (col in names(df_clean)) {
  pct <- missing_pct[col]

  if (pct > 0 & pct <= median_mode_thresh) {
    if (is.numeric(df_clean[[col]])) {
      # Median for numeric
      df_clean[[col]][is.na(df_clean[[col]])] <- median(df_clean[[col]], na.rm = TRUE)
      cat(col, "numeric imputed (median)\n")
    } else {
      # Mode for categorical
      mode_val <- names(sort(table(df_clean[[col]]), decreasing = TRUE))[1]
      df_clean[[col]][is.na(df_clean[[col]])] <- mode_val
      cat(col, "categorical imputed (mode)\n")
    }
  }
}

# MICE imputation for 5–30% missing
cols_mice <- names(missing_pct[missing_pct > median_mode_thresh & missing_pct <= mice_thresh])
if (length(cols_mice) > 0) {
  # Subset for MICE
  mice_data <- df_clean[, cols_mice]

  # Run MICE
  imputed <- mice(mice_data, m = 1, method = 'pmm', seed = 123)
  df_clean[, cols_mice] <- complete(imputed, 1)

  cat("MICE imputation done for columns:\n")
  print(cols_mice)
}

# Check final missing values
colSums(is.na(df_clean))

# 6. Train-test spit
set.seed(44)
train_index <- createDataPartition(df_clean$`Overall Survival Status`, p = 0.75, list = FALSE, times = 1)
train_data <- df_clean[train_index, ]
test_data <- df_clean[-train_index, ]

# Check percentage of censorship
prop.table(table(train_data$`Overall Survival Status`))
prop.table(table(test_data$`Overall Survival Status`))

# 7.1 Linear Cox Model
# Perform stepwise selection using AIC score for variable selection
full_cox <- coxph(Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ . - `Patient ID` - `Sample ID`
                  - `Relapse Free Status (Months)` - `Relapse Free Status` - `Patient's Vital Status`, data = train_data, x = TRUE, y = TRUE)

cox_linear <- stepAIC(full_cox, direction = 'both')
cox_linear

# Model Diagnostic
# Check for linearity with martingales residuals for continous variables
conti_covariates <- c("Age at Diagnosis", "Tumor Size")
res <- residuals(cox_linear, type = "martingale")
plots <- list()

for (var in conti_covariates) {
  col <- train_data[[var]]

  p <- ggplot(data.frame(x = col, y = res), aes(x = x, y = y)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(x = var, y = "Martingale Residuals") +
    theme_minimal()

  plots[[var]] <- p
}

do.call(grid.arrange, c(plots, ncol = 1))

# 7.2. Non-linear cox model
cox_nonlinear <- coxph(Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ ns(`Age at Diagnosis`, df = 3) + log(`Tumor Size`)
                       + `Lymph nodes examined positive` + `Nottingham prognostic index` + `Type of Breast Surgery`+ Chemotherapy
                       + `Pam50 + Claudin-low subtype` + `ER Status` + `HER2 Status` + `Inferred Menopausal State`
                       + `Radio Therapy` + `Tumor Stage`, data = train_data, x = TRUE, y = TRUE)

cox_nonlinear

# Check for time invariance via Schoenfeld residuals
# ggcoxzph(cox.zph(cox_nonlinear))
# par(mfrow = c(4, 3), mar = c(6, 4, 2, 1))
par(mfrow = c(4, 3), mar = c(3, 3, 2, 1))
plot(cox.zph(cox_nonlinear), col = "red")
cox.zph(cox_nonlinear)

# 8.
# Time points at which to evaluate the Brier score
times <- seq(0, max(train_data$`Overall Survival (Months)`), by = 12)

# Compute Brier score
brier_res <- Score(
  list("Linear Cox Model" = cox_linear, "Non Linear Cox Model" = cox_nonlinear),
  formula = Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ 1,
  data = test_data,
  times = times,
  metrics = "brier",
  cens.model = "km"
)

brier_df <- brier_res$Brier$score

total_time <- max(times)
ibs_table <- brier_df %>%
  group_by(model) %>%
  summarise(IBS = trapz(times, Brier) / total_time)
ibs_table

ggplot(brier_df, aes(x = times, y = Brier, color = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.2, color = NA) +
  labs(
    x = "Time",
    y = "Brier Score"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

# 9
new_patients <- data.frame(
  # varied covariates
  `Age at Diagnosis` = c(40, 55, 70),
  `Pam50 + Claudin-low subtype` = factor(c("LumA", "LumB", "Basal"), levels = levels(train_data$`Pam50 + Claudin-low subtype`)),
  Chemotherapy = factor(c("NO", "YES", "YES"), levels = levels(train_data$Chemotherapy)),
  `Tumor Stage` = factor(c(1, 3, 2), levels = levels(train_data$`Tumor Stage`)),

  # fix covariates
  `ER Status` = factor(c("Positive", "Positive", "Positive"), levels = levels(train_data$`ER Status`)),
  `HER2 Status` = factor(c("Negative", "Negative", "Negative"), levels = levels(train_data$`HER2 Status`)),
  `Inferred Menopausal State` = factor(c("Pre", "Pre", "Pre"), levels = levels(train_data$`Inferred Menopausal State`)),
  `Lymph nodes examined positive` = c(2, 2, 2),
  `Nottingham prognostic index` = c(2, 2, 2),
  `Radio Therapy` = factor(c("NO", "NO", "NO"), levels = levels(train_data$`Radio Therapy`)),
  `Tumor Size` = c(25, 25, 25),
  `Type of Breast Surgery` = factor(c("BREAST CONSERVING", "BREAST CONSERVING", "BREAST CONSERVING"), levels = levels(train_data$`Type of Breast Surgery`)),
   check.names = FALSE
)

# Compute survival curves for these patients using non-linear cox models
surv_nonlinear <- survfit(cox_nonlinear, newdata = new_patients)

ggsurvplot(
  surv_nonlinear,
  data = new_patients,
  legend.labs = c("Patient 1", "Patient 2", "Patient 3"),
  legend.title = "",
  palette = c("darkgreen", "orange", "red"),
  xlab = "Months",
  ylab = "Survival Probability",
  size = 0.9,      # line thickness
  conf.int = TRUE,
  ggtheme = theme_bw(base_size = 16)
)


# Survival probability at 12 years using non linear cox model
summary(surv_nonlinear, times = 144)$surv

# 10
# times in months
t0_months <- 7 * 12    # 84
t1_months <- 12 * 12   # 144

# estimate probabilities at times t0 and t1 for each new patient
S_t0_nonlinear <- summary(surv_nonlinear, times = t0_months)$surv
S_t1_nonlinear <- summary(surv_nonlinear, times = t1_months)$surv

# conditional survival P(T > t1 | T > t0) = S(t1)/S(t0)
cond_nonlinear <- S_t1_nonlinear / S_t0_nonlinear

results_nonlinear <- data.frame(
  Patients = paste0("Patient", 1:nrow(new_patients)),
  S_t0_months = c(S_t0_nonlinear),
  S_t1_months = c(S_t1_nonlinear),
  Conditional_Prob = c(cond_nonlinear)
)
print(results_nonlinear)

# 11. Random Forest Model
# Default
rsf_default <- rfsrc(
  Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ .,
  data = train_data[, !names(train_data) %in% c("Patient ID", "Sample ID", "Relapse Free Status (Months)",
                                                "Relapse Free Status", "Patient's Vital Status")],
  ntree = 200
)

print(rsf_default)

# Calibrated RSF
rsf_calibrated <- rfsrc(
  Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ .,
  data = train_data[, !names(train_data) %in% c("Patient ID", "Sample ID", "Relapse Free Status (Months)",
                                                "Relapse Free Status", "Patient's Vital Status")],
  ntree = 200,
  mtry = 7,
  nodesize = 30,
  nsplit = 20,
  importance = "permute",
)

print(rsf_calibrated)

vimp_result <- vimp(rsf_calibrated)
plot(vimp_result)

times <- seq(0, max(train_data$`Overall Survival (Months)`), by = 12)

brier_res <- Score(
  list("Linear Cox Model" = cox_linear, "Non-linear Cox Model" = cox_nonlinear,
       "RSF Default" = rsf_default, "RSF Calibrated" = rsf_calibrated),
  formula = Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ 1,
  data = test_data,
  times = times,
  # metrics = "brier",
)

brier_df <- brier_res$Brier$score

total_time <- max(times)
ibs_table <- brier_df %>%
  group_by(model) %>%
  summarise(IBS = trapz(times, Brier) / total_time)
ibs_table

ggplot(brier_df, aes(x = times, y = Brier, color = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.2, color = NA) +
  labs(
    x = "Time",
    y = "Brier Score"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())


# 13
# We can use IBS (Integrated Brier Score) to compare models
# Interpret variable importance in Random Forest model
plot(vimp(rsf_calibrated))

# 14.
# Add geonomic data
gene_df <- read.table("brca_metabric/data_mrna_illumina_microarray.txt", header = TRUE)
gene_df <- gene_df[, !names(gene_df) %in% "Entrez_Gene_Id"]

# Tranpose genomic data to have gene expression as columns
gene_df_t <- as.data.frame(t(gene_df))

# Remove HugoSymbol rows and change Patient ID structure to fit with clinical dataset
colnames(gene_df_t) <- gene_df$Hugo_Symbol
gene_df_t <- gene_df_t[-1, ]
gene_df_t <- cbind(`Patient ID` = rownames(gene_df_t), gene_df_t)
rownames(gene_df_t) <- NULL
gene_df_t$`Patient ID` <- gsub("\\.", "-", gene_df_t$`Patient ID`)

# Combine Clinical and Genomic Data
merge_data <- merge(df_clean, gene_df_t, by = "Patient ID")

# Remove Patient ID, Sample ID, and Patient's Vital Status columns
merge_data <- merge_data[, -c(1, 2)]
merge_data <- merge_data[, !names(merge_data) %in% "Patient's Vital Status"]

#remove one row where survival month = 0 to avoid error
merge_data <- merge_data[merge_data$`Overall Survival (Months)` != 0, ]

# Remove 16 missing value rows
merge_data <- na.omit(merge_data)

# Convert Genomic data to numeric
merge_data[32:ncol(merge_data)] <- lapply(merge_data[32:ncol(merge_data)], as.numeric)

gene_cols <- names(merge_data[32:ncol(merge_data)])

# create partition stratified on status
train_idx_merge <- createDataPartition(merge_data$`Overall Survival Status`, p = 0.75, list = FALSE)
train_merge <- merge_data[train_idx_merge , ]
test_merge  <- merge_data[-train_idx_merge , ]

# Check percentage of censorship
prop.table(table(train_merge$`Overall Survival Status`))
prop.table(table(test_merge$`Overall Survival Status`))

top_n_genes <- 500

# Extract gene matrices (patients x genes)
expr_train <- as.matrix(train_merge[, gene_cols, drop = FALSE])
expr_test  <- as.matrix(test_merge[,  gene_cols, drop = FALSE])

# compute variance on training set
gene_vars <- colVars(expr_train, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n_genes]

# Subset matrices to top genes
expr_train_top <- expr_train[, top_genes, drop = FALSE]
expr_test_top  <- expr_test[,  top_genes, drop = FALSE]

# Standardize using training mean/sd (no imputation needed)
train_mean <- colMeans(expr_train_top)
train_sd   <- apply(expr_train_top, 2, sd)

# avoid division by zero
train_sd[train_sd == 0] <- 1

scale_with_train <- function(mat, mean_vec, sd_vec) {
  t( (t(mat) - mean_vec) / sd_vec )
}

hv_gene_train <- scale_with_train(expr_train_top, train_mean, train_sd)
hv_gene_test  <- scale_with_train(expr_test_top,  train_mean, train_sd)

# combine high variance genes with clinical data
clinical_data_train <- train_merge[, -c(32:ncol(train_merge))]
clinical_data_train <- clinical_data_train[, -c(22 ,23, 26, 27)]
clinical_data_train <- model.matrix(~ . -1, data = clinical_data_train)

clinical_data_test <- test_merge[, -c(32:ncol(test_merge))]
clinical_data_test <- clinical_data_test[, -c(22 ,23, 26, 27)]
clinical_data_test <- model.matrix(~ . -1, data = clinical_data_test)

x_train <- as.matrix(cbind(clinical_data_train, hv_gene_train))
x_test <- as.matrix(cbind(clinical_data_test , hv_gene_test))

# Fit elastic-net Cox
set.seed(123)
y <- Surv(train_merge$`Overall Survival (Months)`, train_merge$`Overall Survival Status`)
elastic_cox <- cv.glmnet(x_train, y, family = "cox", alpha = 0.5)

plot(elastic_cox)

# Predict risk score for test set
risk_score <- predict(elastic_cox, newx = x_test, s = "lambda.min", type = "link")

y_test <- Surv(test_merge$`Overall Survival (Months)`, test_merge$`Overall Survival Status`)

# Compute concordance
c_index_obj <- concordance(y_test ~ risk_score)
c_index <- c_index_obj$concordance
c_index


# Compare concordance index with previous clinical data models
risk_cox_linear <- predict(cox_linear, newdata = test_merge, type = "risk")
risk_cox_nonlinear <- predict(cox_nonlinear, newdata = test_merge, type = "risk")
risk_rsf_default <- predict(rsf_default, newdata = test_merge)$predicted
risk_rsf_calibrated <- predict(rsf_calibrated, newdata = test_merge)$predicted

y_test <- Surv(test_merge$`Overall Survival (Months)`, test_merge$`Overall Survival Status`)

#  C-index for each
c_cox_linear <- concordance(y_test ~ risk_cox_linear)$concordance
c_cox_nonlinear <- concordance(y_test ~ risk_cox_nonlinear)$concordance
c_rsf_default <- concordance(y_test ~ risk_rsf_default)$concordance
c_rsf_calibrated <- concordance(y_test ~ risk_rsf_calibrated)$concordance
c_elastic <- concordance(y_test ~ risk_score)$concordance

# Combine into a table
cindex_table <- data.frame(
  Model = c("Cox Linear", "Cox Nonlinear", "RSF Default", "RSF Calibrated", "Elastic-net Cox"),
  C_index = c(c_cox_linear, c_cox_nonlinear, c_rsf_default, c_rsf_calibrated, c_elastic)
)

print(cindex_table)




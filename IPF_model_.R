library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(caret)
library(pROC)
library(randomForest)
library(glmnet)
library(org.Hs.eg.db)
library(ROSE)
library(limma)

#Step 1: Load and Preprocess the Data

#loading the dataset
raw_data <- read.csv("/Users/shivaniravindran/Library/CloudStorage/Box-Box/Class_files_SEM2/BMIII/Raw_count_IPF.csv")
head(raw_data)

metadata_mod <-  read.csv("/Users/shivaniravindran/Library/CloudStorage/Box-Box/Class_files_SEM2/BMIII/metadata_modified.csv")
head(metadata_mod)
# I wanted to set the gene names as rows
rownames(raw_data) <- raw_data$X
colnames(raw_data)
# since the gene name column has become the row names there is another column X with gene names , so removing it
raw_data <- raw_data[,-1]
head(raw_data)
# now setting the sample column which is the "title" column as the row names
rownames(metadata_mod) <- metadata_mod$title
#removing the redundant sample column "title"
metadata_mod <- metadata_mod[,-1]

#I also dont want the Donor_id column, so removing it too
metadata_mod <- metadata_mod[,-1]
head(metadata_mod)



#Step 3: Check if the Data Matches (No Data Leakage!)

# Ensure that metadata and raw data are aligned
all(colnames(raw_data) %in% rownames(metadata_mod))
all(colnames(raw_data) == rownames(metadata_mod))

# If they don't match, reorder them
raw_data <- raw_data[, order(colnames(raw_data))]
metadata_mod <- metadata_mod[order(rownames(metadata_mod)), ]

all(colnames(raw_data)==rownames(metadata_mod))



#Step 4: Data Normalization
# Log2-CPM normalization (using DESeq2's VST)
str(metadata_mod)
metadata_mod$Disease_status <- as.factor(metadata_mod$Disease_status)
dds <- DESeqDataSetFromMatrix(countData = raw_data, colData = metadata_mod, design = ~ Disease_status)
dds$Disease_status <- relevel(dds$Disease_status, ref = "Normal")
dds <- DESeq(dds)
#using DESeq2’s vst() or variance-stabilizing transformation for more robust normalization
vst_data <- vst(dds, blind = TRUE)
logcpm <- assay(vst_data)
head(logcpm)



#Step 5: Initial Data split (training and hold out test set)
X_full <- t(logcpm)  # samples x genes
y_full <- as.factor(metadata_mod$Disease_status)  # target variable (Disease vs Normal)
#just visualizing the number of disease and normal samples in the data
barplot(table(y_full), main="Class Distribution in Full Dataset", col=c("lightblue", "salmon"))

set.seed(1999)
train_indices <- createDataPartition(y_full, p = 0.7, list = FALSE)

X_train <- X_full[train_indices, ]
y_train <- y_full[train_indices]
X_test <- X_full[-train_indices, ]
y_test <-y_full[-train_indices]

cat("Training set size:", nrow(X_train), "samples\n")
cat("Test set size:", nrow(X_test), "samples\n")
print(table(y_train))
print(table(y_test))

#Step:6----Defining control paameters(including SMOTE for imbalance)-----
  
ctrl <- trainControl(method = "repeatedcv",
                     number=10,
                     repeats=3,
                     classProbs = TRUE,
                     summaryFunction=twoClassSummary,
                     savePredictions = "final",
                     sampling = "smote",
                     allowParallel = TRUE)
#Step:7 ----- Feature selection (LASSO + stability selection) within the training data --

# This entire block must ONLY use X_train and y_train.
# We will use this to select genes for BOTH RF and LASSO models that follow.
# IMPORTANT: You can embed feature selection directly into `caret`'s `train` function
# using `preProc` or `filter` options, but for more complex stability selection,
# doing it once on the *entire* training set (X_train) before the `train` call
# is a common and acceptable practice. The key is that it's NOT done on X_full.

cat("\nPerforming Stability Selection on training data...\n")
N_STABILITY <- 100 # Number of bootstrap samples for stability selection
sel_counts <- integer(ncol(X_train))
names(sel_counts) <- colnames(X_train)

for (i in seq_len(N_STABILITY)) {
  # Perform bootstrap sampling of the training data
  boot_idx <- sample(seq_len(nrow(X_train)), replace = TRUE)
  X_boot <- X_train[boot_idx, , drop = FALSE]
  y_boot <- y_train[boot_idx]
  
  # Fit LASSO on bootstrap sample
  # Note: `nfolds` should be smaller for smaller boot samples, or use specific lambda.
  # Here, we'll let cv.glmnet decide, but it might warn about fold sizes.
  fit_b <- cv.glmnet(x = X_boot, y = y_boot, family = "binomial", alpha = 1, nfolds = 5)
  coef_b <- coef(fit_b, s = "lambda.1se") # Use lambda.1se for more parsimonious models
  
  # Count selected genes
  # The intercept is the first coefficient, so we exclude it.
  sel_genes_boot <- rownames(coef_b)[which(as.vector(coef_b)[-1] != 0) + 1] # +1 to align with features
  sel_counts[sel_genes_boot] <- sel_counts[sel_genes_boot] + 1
}

sel_freq <- sel_counts / N_STABILITY
sel_freq <- as.data.frame(sel_freq) %>% arrange(desc(sel_freq))
# Define a threshold for stable genes (e.g., selected in 60% or more bootstraps)
stable_genes_threshold <- 0.60
stable_genes <- rownames(sel_freq)[which(sel_freq$sel_freq >= stable_genes_threshold)]

cat(length(stable_genes), "stable genes identified with frequency >= ", stable_genes_threshold, ":\n")
print(head(stable_genes))


#Step 8: Train Models on the Training Set (with selected stable genes) ---
# Now, `train` will perform its own cross-validation on X_train (using only stable genes)
# and apply SMOTE and pre-processing within each fold.

X_train_stable <- X_train[,stable_genes]
X_test_stable <- X_test[,stable_genes]

set.seed(1999)

############# Random Forest Model Training ##############
rf_fit <- train(x= X_train_stable, y= y_train, method ="rf",
                metric = "ROC",
                trControl = ctrl,
                tuneLength=3,
                ntree= 500,
                preProc = c("center","scale"))
print(rf_fit)

############# LASSO Model Training ##############
cat("\nTraining LASSO (glmnet) model...\n")
set.seed(1999) # For reproducibility of model training
lasso_fit <- train(x = X_train_stable, y = y_train,
                   method = "glmnet",
                   metric = "ROC",
                   trControl = ctrl,
                   tuneLength = 8, # Number of alpha/lambda combinations to try
                   preProc = c("center", "scale")) # Scaling *within* CV folds

print(lasso_fit)

#step 9: Now let's evaluate model performance on the test set ----
POS_CLASS <- "Disease"

################ Random Forest Evaluation on test ##############
cat("\n Evaluating Random forest on the test set:\n")
rf_preds_test_prob <-predict(rf_fit, newdata= X_test_stable, type = "prob")
rf_roc_test <- roc(response = y_test, predictor = rf_preds_test_prob[[POS_CLASS]], 
                   levels = c("Normal","Disease"), direction = "<")

cat("RF AUC on test set:", auc (rf_roc_test), "\n")

################ LASSPO Evaluation on test ##############
cat("\nEvaluating LASSO on Test Set:\n")
lasso_preds_test_prob <- predict(lasso_fit, newdata = X_test_stable, type = "prob")
lasso_roc_test <- roc(response = y_test, predictor = lasso_preds_test_prob[[POS_CLASS]],
                      levels = c("Normal", "Disease"), direction = "<")
cat("LASSO AUC on Holdout Test Set:", auc(lasso_roc_test), "\n")


# --- Step 9: Compare Training vs. Test Performance (for potential overfitting diagnosis) ---
cat("\n--- Training vs. Test Set Performance Comparison ---\n")

# Random Forest
rf_preds_train_prob <- predict(rf_fit, newdata = X_train_stable, type = "prob")
rf_roc_train <- roc(response = y_train, predictor = rf_preds_train_prob[[POS_CLASS]],
                    levels = c("Normal", "Disease"), direction = "<")
cat("RF AUC on Training Set (Internal CV mean):", max(rf_fit$results$ROC), "\n") # This is mean CV AUC
cat("RF AUC on Full Training Set (for comparison):", auc(rf_roc_train), "\n")
cat("RF AUC on Holdout Test Set:", auc(rf_roc_test), "\n")
cat("Difference (Train - Test):", auc(rf_roc_train) - auc(rf_roc_test), "\n")

# Plot Random Forest ROC curves
plot(rf_roc_train, col = "darkblue", lwd = 2, main = "Random Forest ROC Curves",
     legacy.axes = TRUE, xlab = "False Positive Rate", ylab = "True Positive Rate")
lines(rf_roc_test, col = "skyblue", lwd = 2, lty = 2)
abline(a = 0, b = 1, lty = 3, col = "gray") # Reference line for random classifier
legend("bottomright", legend = c(paste0("Train AUC: ", round(auc(rf_roc_train), 3)),
                                 paste0("Test AUC: ", round(auc(rf_roc_test), 3))),
       col = c("darkblue", "skyblue"), lwd = 2, lty = c(1, 2), bty = "n")



# LASSO
lasso_preds_train_prob <- predict(lasso_fit, newdata = X_train_stable, type = "prob")
lasso_roc_train <- roc(response = y_train, predictor = lasso_preds_train_prob[[POS_CLASS]],
                       levels = c("Normal", "Disease"), direction = "<")
cat("LASSO AUC on Training Set (Internal CV mean):", max(lasso_fit$results$ROC), "\n")
cat("LASSO AUC on Full Training Set (for comparison):", auc(lasso_roc_train), "\n")
cat("LASSO AUC on Holdout Test Set:", auc(lasso_roc_test), "\n")
cat("Difference (Train - Test):", auc(lasso_roc_train) - auc(lasso_roc_test), "\n")



Step 10: Check for Overfitting

# Check Random Forest OOB error
rf_oob <- rf_fit$finalModel$err.rate
cat("RF OOB Error:", mean(rf_oob[, "OOB"]), "\n")

# Permutation test for Random Forest
Nperm <- 500
obs_roc_rf <- auc(rf_roc_test)
perm_rocs <- numeric(Nperm)
for (i in seq_len(Nperm)) {
  y_perm <- sample(y_test)
  set.seed(1999 + i)
  y_test_perm <- sample(y_test) # Permute test labels
  perm_preds_rf_prob <- predict(rf_fit, newdata = X_test_stable, type ="prob")# Predict using the *trained* RF model on the *original* X_test_stable, but with *permuted* labels
  perm_roc_rf <- roc(response = y_test_perm, predictor = perm_preds_rf_prob[[POS_CLASS]],
                     levels=c("Normal", "Disease"), direction = "<")
  perm_rocs[i] <- auc(perm_roc_rf)
}

obs_roc_rf<- as.numeric(obs_roc_rf)
perm_roc_rf <- as.numeric(perm_rocs)
p_val <- (sum(perm_roc_rf >= obs_roc_rf, na.rm = TRUE) + 1) / (Nperm + 1)
cat("Permutation p-value (RF CV ROC):", p_val, "\n")


# Permutation test for Random Forest
Nperm <- 300
obs_roc_rf <- auc(roc_rf)
perm_rocs <- numeric(Nperm)
for (i in seq_len(Nperm)) {
  y_perm <- sample(y)
  set.seed(1000 + i)
  perm_fit <- train(x = X[, stable_genes], y = y_perm, method = "rf", metric = "ROC", trControl = ctrl, tuneLength = 2, ntree = 200)
  perm_rocs[i] <- max(perm_fit$results$ROC, na.rm = TRUE)
}
p_val <- (sum(perm_rocs >= obs_roc_rf, na.rm = TRUE) + 1) / (Nperm + 1)
cat("Permutation p-value (RF CV ROC):", p_val, "\n")




## Let's take a look at the stable genes by visualizing their expression patterns
stable_genes_df <- as.data.frame(X_full[,stable_genes])
stable_genes_df$Disease_status <- y_full # adding the disease status column
stable_genes_df$Sample_id <- rownames(X_full)

#we will convert to long format which is suitable for ggplot2
stable_genes.l <- stable_genes_df %>%
  pivot_longer(cols = -c(Disease_status, Sample_id),
               names_to = "Gene",
               values_to = "Expression")
# Plotting
top_10 <- min(10, length(stable_genes)) 
if(top_10>0){
  genes_to_plot <- stable_genes[1:top_10]
  ggplot(stable_genes.l %>% filter(Gene %in% genes_to_plot),
        aes(x= Disease_status, y= Expression, fill=Disease_status)) +
    geom_boxplot()+
    geom_jitter(width= 0.2, alpha =0.6)+
    facet_wrap(~Gene, scales ="free_y", ncol= 5)+
    labs(title="Expression of Top Stable Genes by Disease Status",
         x= "Normalized Expression(VST)",
         y= "Disease status")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

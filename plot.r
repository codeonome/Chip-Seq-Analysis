# STEP 1: Setup
# Set working directory 
setwd("/gpfs/data/khodadadilab/home/temp/Di-Stefano-Lab-Assignment/Task-2/data")

# Load required libraries
library(ggplot2)

# STEP 2: Load and Clean the Data
# Read the processed signal data
df <- read.delim("h3k27ac_CORRECTED_data.tab", header = TRUE, check.names = FALSE)

# Remove single quotes from column names
colnames(df) <- gsub("'", "", colnames(df))

# Ensure signal columns are numeric
df$`CTRL_rep1.bw` <- as.numeric(df$`CTRL_rep1.bw`)
df$`CTRL_rep2.bw` <- as.numeric(df$`CTRL_rep2.bw`)
df$`DDX6_rep1.bw` <- as.numeric(df$`DDX6_rep1.bw`)
df$`DDX6_rep2.bw` <- as.numeric(df$`DDX6_rep2.bw`)


# STEP 3: Compute log2(RPKM + 1) Values
# Calculate region lengths in kilobases
df$length_kb <- (df$end - df$start) / 1000

# Calculate RPKM for each replicate
df$CTRL_1_rpkm <- df$`CTRL_rep1.bw` / df$length_kb
df$CTRL_2_rpkm <- df$`CTRL_rep2.bw` / df$length_kb
df$DDX6_1_rpkm <- df$`DDX6_rep1.bw` / df$length_kb
df$DDX6_2_rpkm <- df$`DDX6_rep2.bw` / df$length_kb

# Log2-transform the RPKM values
df$CTRL_mean_logRPKM <- rowMeans(log2(cbind(df$CTRL_1_rpkm, df$CTRL_2_rpkm) + 1))
df$DDX6_mean_logRPKM <- rowMeans(log2(cbind(df$DDX6_1_rpkm, df$DDX6_2_rpkm) + 1))

# STEP 4: Perform Statistical Test

# Two-sample t-test on mean log2(RPKM + 1) values
ttest_result <- t.test(df$CTRL_mean_logRPKM, df$DDX6_mean_logRPKM)

# Format p-value for plot
pval <- signif(ttest_result$p.value, 3) 
pval_label <- paste0("p = ", pval)

# STEP 5: Prepare Plotting Data
# Combine data for ggplot2
df_plot <- data.frame(
  Condition = rep(c("CTRL", "DDX6"), each = nrow(df)),
  logRPKM = c(df$CTRL_mean_logRPKM, df$DDX6_mean_logRPKM)
)

# STEP 6: Plot (Violin + Boxplot + P-value)

ggplot(df_plot, aes(x = Condition, y = logRPKM, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 3) +
  annotate("text", x = 1.5, y = max(df_plot$logRPKM) + 0.1, label = pval_label, size = 5) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("CTRL" = "#1f77b4", "DDX6" = "#d62728")) +
  labs(
    title = "H3K27ac signal at 684 Super-Enhancers",
    y = "log₂(RPKM + 1)", x = ""
  )

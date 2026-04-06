# ==========================================
# 顶刊标准: 三重交叉验证 + 脑营养竞争指数验证
# ==========================================
library(randomForest)
library(glmnet)
library(caret)
library(pls)
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(ggpubr) # 用于画带P值的箱线图 (如果没有请先 install.packages("ggpubr"))

# 1. 读入并整理数据
data <- read.csv("Table_S1_Cleaned_Metabolomics_Matrix.csv", stringsAsFactors = FALSE)
data <- data %>%
  mutate(Group = case_when(
    grepl("-C-M-", SampleID) ~ "Healthy_Control",
    grepl("-MS-M-", SampleID) ~ "Risk_Exposure",
    TRUE ~ "Unknown"
  )) %>%
  mutate(Group = factor(Group, levels = c("Healthy_Control", "Risk_Exposure"))) %>%
  filter(Group != "Unknown")

X <- as.matrix(data[, 3:ncol(data)])
Y <- data$Group

cat("\n🚀 正在启动三重机器学习交叉验证...\n")

# ==========================================
# 算法 1: Random Forest
# ==========================================
set.seed(2026)
rf_model <- randomForest(x = X, y = Y, ntree = 500, importance = TRUE)
rf_imp <- as.data.frame(importance(rf_model))
rf_markers <- rownames(rf_imp)[order(rf_imp$MeanDecreaseGini, decreasing = TRUE)][1:15]

# ==========================================
# 算法 2: LASSO
# ==========================================
set.seed(2026)
Y_num <- ifelse(Y == "Risk_Exposure", 1, 0)
lasso_cv <- cv.glmnet(X, Y_num, family = "binomial", alpha = 1, nfolds = nrow(X))
lasso_coef <- coef(lasso_cv, s = "lambda.min")
lasso_markers <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_markers <- setdiff(lasso_markers, "(Intercept)")

# ==========================================
# 算法 3: PLS-DA (精准修复底层提取逻辑)
# ==========================================
set.seed(2026)
pls_fit <- plsda(x = X, y = Y, ncomp = 2)
# 直接获取 varImp 对象
pls_imp_raw <- varImp(pls_fit)
# 智能判断其结构并转换为纯数据框
if("importance" %in% names(pls_imp_raw)) {
  pls_imp_df <- as.data.frame(pls_imp_raw$importance)
} else {
  pls_imp_df <- as.data.frame(pls_imp_raw)
}
# 安全提取第一列并降序
importance_scores <- as.numeric(pls_imp_df[, 1])
pls_markers <- rownames(pls_imp_df)[order(importance_scores, decreasing = TRUE)][1:15]

# ==========================================
# 韦恩图绘制 (寻找黄金交集)
# ==========================================
venn_list <- list(`Random Forest` = rf_markers, `LASSO` = lasso_markers, `PLS-DA` = pls_markers)
golden_markers <- Reduce(intersect, venn_list)
cat("\n🏆 三大算法黄金标志物:\n")
print(golden_markers)

p_venn <- ggVennDiagram(venn_list, color = "black", lwd = 0.8, set_color = "black", set_size = 5) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Multi-Algorithm Validation")
ggsave("Venn_Diagram_Validation.pdf", plot = p_venn, width = 6, height = 5, device = "pdf")

# ==========================================
# 🚀 进阶临床护理价值: 计算“血脑屏障竞争指数” (强制修复 margin 版)
# ==========================================
library(ggpubr)
library(dplyr)
library(ggplot2)

cat("\n⚙️ 正在重新加载数据并强制修正临床分组...\n")
# 1. 重新读取数据
data <- read.csv("Cleaned_Data_ST004597_Merged.csv", stringsAsFactors = FALSE)

# 2. 强制修复分组 (彻底消灭单一的 C57 组)
data <- data %>%
  mutate(Group = case_when(
    grepl("-C-M-", SampleID) ~ "Healthy_Control",
    grepl("-MS-M-", SampleID) ~ "Risk_Exposure",
    TRUE ~ "Unknown"
  )) %>%
  filter(Group != "Unknown") %>%
  mutate(Group = factor(Group, levels = c("Healthy_Control", "Risk_Exposure")))

cat("📊 当前真实分组校验:\n")
print(table(data$Group)) 

# 3. 寻找列名并计算指数
bcaa_col <- grep("Acetylleucine", colnames(data), ignore.case = TRUE, value = TRUE)[1]
trp_col <- grep("Tryptophan", colnames(data), ignore.case = TRUE, value = TRUE)[1]

if(is.na(bcaa_col) | is.na(trp_col)) {
  cat("⚠️ 警告: 找不到对应的代谢物列名，请检查拼写。\n")
} else {
  # 计算竞争指数
  data$BCAA_Trp_Ratio <- data[[bcaa_col]] - data[[trp_col]]
  
  # 4. 绘制临床价值箱线图 (使用 Wilcoxon 非参数检验)
  p_ratio <- ggboxplot(data, x = "Group", y = "BCAA_Trp_Ratio", 
                       color = "Group", palette = c("#4DBBD5FF", "#E64B35FF"),
                       add = "jitter", shape = "Group", size = 1) +
    # 核心修复点1：用 method.args 传递 exact 参数
    stat_compare_means(method = "wilcox.test", 
                       method.args = list(exact = FALSE), 
                       label = "p.format", 
                       label.x.npc = "center", 
                       size = 5) +
    theme_classic(base_size = 14) +
    labs(title = "Blood-Brain Barrier Transport Competition Index",
         subtitle = "(N-Acetylleucine vs Tryptophan)",
         x = "", y = "Log2 Ratio (BCAA / Trp)") +
    # 核心修复点2：强制使用 ggplot2::margin 彻底解决报错
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(b = 10)),
          plot.subtitle = element_text(hjust = 0.5, face = "italic", color = "grey40", margin = ggplot2::margin(b = 15)))
  
  ggsave("BCAA_Tryptophan_Competition_Index.pdf", plot = p_ratio, width = 5, height = 6, device = "pdf")
  cat("✅ 竞争指数箱线图已成功生成并保存为: BCAA_Tryptophan_Competition_Index.pdf\n")
}

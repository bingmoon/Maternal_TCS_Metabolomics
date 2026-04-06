# ==========================================
# 护理科研项目: 早期预警标志物筛选与机器学习建模
# ==========================================
# 加载必备R包 (如果没有请先 install.packages(c("ggplot2", "ggrepel", "randomForest")))
library(dplyr)
library(ggplot2)
library(ggrepel)
library(randomForest)

# 1. 读入您刚才生成的黄金矩阵
data <- read.csv("Table_S1_Cleaned_Metabolomics_Matrix.csv", stringsAsFactors = FALSE)

# 2. 智能修复分组标签 (基于样本名解析真实临床分组)
cat("\n🔧 [Step 1] 正在修复并重构真实临床分组...\n")
data <- data %>%
  mutate(Group = case_when(
    grepl("-C-M-", SampleID) ~ "Healthy_Control",
    grepl("-MS-M-", SampleID) ~ "Risk_Exposure",
    TRUE ~ "Unknown"
  )) %>%
  mutate(Group = factor(Group, levels = c("Healthy_Control", "Risk_Exposure")))

cat("📊 当前真实分组情况:\n")
print(table(data$Group))

# 3. 差异分析 (T检验与Log2FC计算)
cat("\n⚙️ [Step 2] 正在进行差异表达分析...\n")
features <- colnames(data)[3:ncol(data)]
diff_results <- data.frame(Feature = features, log2FC = NA, pvalue = NA)

for (i in 1:length(features)) {
  feat <- features[i]
  val_control <- data[[feat]][data$Group == "Healthy_Control"]
  val_risk <- data[[feat]][data$Group == "Risk_Exposure"]
  
  # Log2 矩阵的差异倍数计算 (暴露组 - 对照组)
  log2fc <- mean(val_risk) - mean(val_control)
  # T检验
  ptest <- tryCatch(t.test(val_risk, val_control)$p.value, error = function(e) NA)
  
  diff_results$log2FC[i] <- log2fc
  diff_results$pvalue[i] <- ptest
}

# 标记显著性 (针对小样本预警探索，设置P<0.05，|log2FC|>0.3为阈值)
diff_results <- diff_results %>%
  mutate(Significance = case_when(
    pvalue < 0.05 & log2FC > 0.3 ~ "Up (Risk Marker)",
    pvalue < 0.05 & log2FC < -0.3 ~ "Down (Protective Marker)",
    TRUE ~ "Not Sig"
  ))

# 4. 绘制高颜值火山图
library(ggplot2)
library(ggrepel)
library(dplyr)

diff_results$pvalue[diff_results$pvalue == 0] <- 1e-10 

# 定义顶刊配色 (如 Nature Medicine 常用的经典红蓝配)
color_palette <- c("Up (Risk Marker)" = "#E64B35FF", 
                   "Down (Protective Marker)" = "#4DBBD5FF", 
                   "Not Sig" = "#CCCCCC")

# 核心绘图代码
p_volcano_top <- ggplot(diff_results, aes(x = log2FC, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.85, size = 3, stroke = 0.5) + # 优化散点质感
  scale_color_manual(values = color_palette) +
  
  # 添加阈值虚线 (精细化线条)
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "#333333", size = 0.5, alpha = 0.7) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "longdash", color = "#333333", size = 0.5, alpha = 0.7) +
  
  # 智能排版标签 (仅标注显著的代谢物，去除重叠，增加牵引线美感)
  geom_text_repel(data = filter(diff_results, Significance != "Not Sig"),
                  aes(label = Feature), 
                  size = 3.5, 
                  fontface = "italic", # 代谢物名称通常用斜体显得更专业
                  box.padding = 0.5, 
                  point.padding = 0.3, 
                  segment.color = 'grey50',
                  max.overlaps = 20) +
  
  # 顶刊极简主题 (theme_classic 去除无用背景)
  theme_classic(base_size = 14) + 
  theme(
    text = element_text(color = "black"), 
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(color = "black", size = 0.8), # 加粗坐标轴线
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "top", # 图例放上方，节省作图空间
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    # 【已修复】：强制指定调用 ggplot2 的 margin 函数
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10) 
  ) +
  labs(x = expression(bold(Log[2]~"Fold Change (Risk / Control)")), 
       y = expression(bold(-Log[10]~italic(P)-"value")))

# 导出为 PDF
ggsave("Volcano_Plot_TopTier.pdf", plot = p_volcano_top, width = 7, height = 6, device = "pdf")
cat("✅ 顶刊级别火山图已成功保存为: Volcano_Plot_TopTier.pdf\n")

# 5. 随机森林 (机器学习) 筛选核心预警 Panel
cat("\n🤖 [Step 3] 正在启动 Random Forest 机器学习筛选核心特征...\n")
set.seed(2026) # 固定种子确保复现
rf_model <- randomForest(Group ~ ., data = data[, -1], ntree = 500, importance = TRUE)

# 提取特征重要性得分 (MeanDecreaseGini 越高，预测能力越强)
rf_imp <- as.data.frame(importance(rf_model))
rf_imp$Feature <- rownames(rf_imp)
rf_imp <- rf_imp %>% arrange(desc(MeanDecreaseGini))

cat("\n🏆 随机森林预测出的 Top 5 最核心预警标志物 (Nurse Warning Panel):\n")
print(head(rf_imp[, c("Feature", "MeanDecreaseGini")], 5))

# 导出结果供写文章使用
write.csv(diff_results, "Table_S2_Differential_Analysis.csv", row.names = FALSE)
write.csv(rf_imp, "Table_S3_Random_Forest_Importance.csv", row.names = FALSE)
cat("\n🎉 分析全部完成！详细表格已保存为CSV文件。\n")

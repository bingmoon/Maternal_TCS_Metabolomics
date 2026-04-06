# ==========================================
# 护理科研项目: 母婴环境暴露与营养干预预测预警
# Target Dataset: ST004597 (PR002895)
# 任务: 正负离子双模式数据融合与严格质控清洗
# ==========================================

# 加载必要R包
library(jsonlite)
library(dplyr)
library(stringr)
library(impute) # BiocManager::install("impute")

# ==========================================
# Step 1: 从云端拉取 ST004597 临床分组信息
# ==========================================
cat("\n📥 [Step 1] 正在安全拉取 ST004597 临床分组信息...\n")

url_4597 <- "https://www.metabolomicsworkbench.org/rest/study/study_id/ST004597/factors"
raw_meta <- bind_rows(fromJSON(url_4597))

# 智能提取并标准化临床信息
meta_4597 <- raw_meta %>%
  mutate(
    mb_sample_id = trimws(ifelse("sample_id" %in% colnames(raw_meta), sample_id, mb_sample_id)),
    local_sample_id = trimws(local_sample_id),
    # 动态抓取主要分组因子 (通常在 factors 列中)
    RawGroup = str_match(factors, "([^:\\|]+):\\s*([^|\\s]+)")[,3] 
  ) %>%
  # 这里为了护理学研究目标进行重新定义映射（需根据实际返回的RawGroup名称微调）
  # 结合原摘要，推测分组包含：对照组、TCS暴露组(患病风险组)、TCS+Trp组(营养干预组)
  mutate(
    Group = case_when(
      grepl("Control|Vehicle|WT", RawGroup, ignore.case = TRUE) ~ "Healthy_Control",
      grepl("TCS", RawGroup, ignore.case = TRUE) & !grepl("Trp|Tryptophan", RawGroup, ignore.case = TRUE) ~ "Risk_Exposure",
      grepl("Trp|Tryptophan", RawGroup, ignore.case = TRUE) ~ "Nutrition_Intervention",
      TRUE ~ RawGroup # 保留未命中规则的原始命名
    )
  ) %>%
  dplyr::select(mb_sample_id, local_sample_id, Group)

cat("📊 [Step 1] 临床分组拉取成功！样本分布如下:\n")
print(table(meta_4597$Group))

# ==========================================
# Step 2: 正负离子双矩阵读取、对齐与融合
# ==========================================
cat("\n🚀 [Step 2] 正在加载并融合双离子模式表达矩阵...\n")

file_pos <- "MSdata_ST004597_1.txt"
file_neg <- "MSdata_ST004597_2.txt"
if(!file.exists(file_pos) | !file.exists(file_neg)) stop("❌ 找不到原始数据文件，请检查路径！")

# 内部函数：读取并初步处理单模式矩阵
process_raw_matrix <- function(file_path, mode_suffix) {
  df <- read.delim(file_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  # 使用代谢物名称作为行名，加上正负离子后缀防止重名冲突
  rownames(df) <- make.unique(paste0(as.character(df[[1]]), "_", mode_suffix))
  expr <- df[, -1]
  if("RefMet_name" %in% colnames(expr)) expr <- expr[, colnames(expr) != "RefMet_name"]
  return(expr)
}

expr_pos <- process_raw_matrix(file_pos, "POS")
expr_neg <- process_raw_matrix(file_neg, "NEG")

# 寻找正负离子模式共有的样本列
common_samples <- intersect(colnames(expr_pos), colnames(expr_neg))
cat(sprintf("   ✅ 正负模式共有样本数: %d\n", length(common_samples)))

# ==========================================
# Step 3: 临床信息匹配与数据质控清洗
# ==========================================
cat("\n⚙️ [Step 3] 正在进行严格的数据质控 (过滤、归一化、插补、Log2)...\n")

# 1. 事实匹配：找出临床信息与表达矩阵中完全一致的样本
if (any(common_samples %in% meta_4597$local_sample_id)) {
  match_col <- "local_sample_id"
} else if (any(common_samples %in% meta_4597$mb_sample_id)) {
  match_col <- "mb_sample_id"
} else {
  stop("❌ 严重错误：表达矩阵列名与云端 ID 无法匹配！")
}

# 提取交集样本并按临床信息排序
meta_matched <- meta_4597[meta_4597[[match_col]] %in% common_samples, ]
final_samples <- meta_matched[[match_col]]

# 2. 矩阵融合: 上下拼接正负离子矩阵
expr_merged <- rbind(expr_pos[, final_samples], expr_neg[, final_samples])
cat(sprintf("   ✅ 样本严格对齐: 匹配 %d 个样本, 融合后总特征数: %d\n", length(final_samples), nrow(expr_merged)))

# 3. 转换为数值矩阵并处理0值
expr_mat <- as.matrix(expr_merged)
class(expr_mat) <- "numeric"
expr_mat[expr_mat == 0] <- NA

# 4. 核心质控 1: 剔除缺失率 > 30% 的噪点特征
missing_rate <- rowMeans(is.na(expr_mat))
expr_mat_filtered <- expr_mat[missing_rate <= 0.3, ]
cat(sprintf("   ✅ 特征过滤 (缺失率 <= 30%%): 保留了 %d 个高质量代谢特征\n", nrow(expr_mat_filtered)))

# 5. 核心质控 2: TIC 归一化 (矫正样本总量偏差)
tic <- colSums(expr_mat_filtered, na.rm = TRUE)
expr_mat_tic <- sweep(expr_mat_filtered, 2, tic, "/") * median(tic, na.rm = TRUE)

# 6. 核心质控 3: KNN 插补 (填补少量缺失值)
set.seed(2026) # 固定种子，保证结果100%可复现
cat("   ⏳ 正在进行 KNN 缺失值插补 (K=10)...\n")
imputed_res <- impute.knn(expr_mat_tic, k = 10)$data

# 7. 核心质控 4: Log2 转换 (使特征服从正态分布)
min_val <- min(imputed_res[imputed_res > 0], na.rm = TRUE)
imputed_log2 <- log2(imputed_res + (min_val / 2))

# ==========================================
# Step 4: 组装黄金分析矩阵
# ==========================================
data_clean_4597 <- as.data.frame(t(imputed_log2))
data_clean_4597$SampleID <- rownames(data_clean_4597)

# 与分组信息合并
data_clean_4597 <- data_clean_4597 %>%
  left_join(meta_matched %>% dplyr::select(all_of(match_col), Group), 
            by = c("SampleID" = match_col)) %>%
  dplyr::select(SampleID, Group, everything())

# 为了后续机器学习和差异分析，将Group转化为Factor
# 注: 如果上面的case_when没有命中，这里会用到真实的RawGroup，您可以打印 unique(data_clean_4597$Group) 确认
if("Healthy_Control" %in% data_clean_4597$Group){
  data_clean_4597$Group <- factor(data_clean_4597$Group, 
                                  levels = c("Healthy_Control", "Risk_Exposure", "Nutrition_Intervention"))
} else {
  data_clean_4597$Group <- as.factor(data_clean_4597$Group)
}

cat("\n🎉 [完美收官] ST004597 双模式数据融合与清洗全部完成！\n")
cat("   📊 最终黄金矩阵维度: ", nrow(data_clean_4597), "个样本 x ", ncol(data_clean_4597) - 2, "个代谢特征。\n")

# 保存清洗后的黄金矩阵
write.csv(data_clean_4597, "Table_S1_Cleaned_Metabolomics_Matrix.csv", row.names = FALSE)
cat("   💾 已保存至本地文件: Cleaned_Data_ST004597_Merged.csv\n")


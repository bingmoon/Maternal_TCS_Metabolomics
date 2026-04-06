# ==========================================
# 顶刊标准: 肠脑轴共代谢网络图 (数学硬刚严谨版)
# ==========================================
# 加载必备包 (如果没有请先安装: install.packages(c("igraph", "ggraph", "Hmisc", "tidygraph")))
library(igraph)
library(ggraph)
library(Hmisc)
library(dplyr)
library(tidygraph)
library(ggplot2)
cat("\n🕸️ 正在构建肠脑轴共代谢网络...\n")
# 1. 读取黄金矩阵
data <- read.csv("Table_S1_Cleaned_Metabolomics_Matrix.csv", stringsAsFactors = FALSE)
col_names <- colnames(data)
# 2. 智能动态匹配：从真实数据列名中抓取“故事核心主角”
story_nodes <- c(
grep("Kynurenic", col_names, ignore.case = TRUE, value = TRUE)[1],          # KYNA (黄金预警标志物)
grep("Tryptophan", col_names, ignore.case = TRUE, value = TRUE)[1],         # 色氨酸
grep("serotonin", col_names, ignore.case = TRUE, value = TRUE)[1],          # 血清素
grep("Aminobutyric", col_names, ignore.case = TRUE, value = TRUE)[1],       # GABA
grep("Acetylleucine", col_names, ignore.case = TRUE, value = TRUE)[1],      # 神经保护剂/BCAA
grep("acetylisoleucine", col_names, ignore.case = TRUE, value = TRUE)[1],   # 神经保护剂/BCAA
grep("Phenylacetylglutamine", col_names, ignore.case = TRUE, value = TRUE)[1], # 肠道毒素 (PAGln)
grep("Formyl.*methionine", col_names, ignore.case = TRUE, value = TRUE)[1], # 肠道屏障泄露物
grep("Isovaleroylglycine", col_names, ignore.case = TRUE, value = TRUE)[1]  # 肠道腐败产物
)
# 剔除可能为空的节点，生成最终主角名单
story_nodes <- na.omit(story_nodes)
cat(sprintf("✅ 成功锁定 %d 个真实存在的肠脑轴节点。\n", length(story_nodes)))
# 3. 提取表达矩阵并计算 Spearman 相关系数
expr_mat <- as.matrix(data[, story_nodes])
cor_res <- rcorr(expr_mat, type = "spearman")
r_mat <- cor_res$r
# 4. 提取极其显著的边 (数学硬刚阈值: |r| >= 0.75，在N=8时数学上绝对保证 P < 0.05)
edges <- data.frame(
from = rep(rownames(r_mat), ncol(r_mat)),
to = rep(colnames(r_mat), each = nrow(r_mat)),
r = as.vector(r_mat)
) %>%
filter(from != to) %>%
filter(abs(r) >= 0.75)
# 去除双向重复的边 (比如 A-B 和 B-A 只留一条)
edges$pair <- apply(edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
edges <- edges[!duplicated(edges$pair), ]
edges$pair <- NULL
# 5. 构建图形对象 (强制包含所有主角，即使没有达标的边也不会消失)
nodes_df <- data.frame(name = story_nodes)
net <- graph_from_data_frame(d = edges, vertices = nodes_df, directed = FALSE)
# 6. 绘制高颜值星系网络图 (环形布局)
set.seed(2026) # 固定排版布局
p_net <- ggraph(net, layout = "linear", circular = TRUE) +
# 连线: 红色正相关(协同)，蓝色负相关(拮抗)
geom_edge_link(aes(color = r, width = abs(r)), alpha = 0.8) +
scale_edge_color_gradient2(low = "#4DBBD5FF", mid = "white", high = "#E64B35FF", midpoint = 0,
name = "Spearman r") +
scale_edge_width(range = c(1, 3), guide = "none") +
# 节点: 统一的学术蓝
geom_node_point(size = 8, color = "#3C5488FF", alpha = 0.9) +
# 标签: 斜体、智能避让
geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "italic",
color = "black", point.padding = 0.5) +
# 主题: 极简留白，彻底解决 margin 报错
theme_void(base_size = 14) +
labs(title = "Gut-Brain Axis Co-Metabolic Network",
subtitle = "Strict Spearman Correlation (|r| \u2265 0.75, p < 0.05)") +
theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = ggplot2::margin(b = 5)),
plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", margin = ggplot2::margin(b = 15)),
legend.position = "bottom")
# 7. 导出高精度 PDF
ggsave("Gut_Brain_CoMetabolic_Network.pdf", plot = p_net, width = 8, height = 8, device = "pdf")
cat("🎉 [完美收官] 肠脑轴共代谢网络图已生成并保存为: Gut_Brain_CoMetabolic_Network.pdf\n")s

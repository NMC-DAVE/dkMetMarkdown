## 中国区域：降水-温度二维双变量地图（CHELSA + biscale）

# 0. 包安装与加载 -------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# 按需安装 GitHub 包
if (!requireNamespace("rchelsa", quietly = TRUE)) remotes::install_github("inSileco/rchelsa")
if (!requireNamespace("biscale", quietly = TRUE)) remotes::install_github("chris-prener/biscale")

pacman::p_load(
  geodata, tidyverse, sf, terra,
  rchelsa, biscale, cowplot
)

# 1. 数据准备 -----------------------------------------------------------------
main_dir <- getwd()

# 下载 CHELSA 生物气候变量：BIO1（年均温）和 BIO12（年降水）
ids <- c(1, 12)

download_chelsa_data <- function(id, path) {
  # 兼容不同 rchelsa 版本的函数名：优先使用 get_chelsa_data
  if ("get_chelsa_data" %in% getNamespaceExports("rchelsa")) {
    rchelsa::get_chelsa_data(categ = "clim", type = "bio", id = id, path = path)
  } else if ("get_chelsea_data" %in% getNamespaceExports("rchelsa")) {
    rchelsa::get_chelsea_data(categ = "clim", type = "bio", id = id, path = path)
  } else {
    stop("rchelsa: 未找到 get_chelsa_data/get_chelsea_data 函数，请更新 rchelsa 包。")
  }
}

# 若本地已存在 CHELSA tif，则跳过下载
chelsa_files <- file.path(main_dir, c("CHELSA_bio10_01.tif", "CHELSA_bio10_12.tif"))
if (!all(file.exists(chelsa_files))) {
  invisible(lapply(ids, download_chelsa_data, path = main_dir))
} else {
  message("CHELSA tif 已存在，跳过下载。")
}

# 载入栅格（CHELSA 1981–2010 基期：bio10 文件名前缀）
temp <- terra::rast(file.path(main_dir, "CHELSA_bio10_01.tif"))  # 年均温（单位依版本可能缩放）
prec <- terra::rast(file.path(main_dir, "CHELSA_bio10_12.tif"))  # 年降水（mm/年）

# 将年降水转为“月均降水”（mm/月）
prec_monthly <- prec / 12

# 组合并命名图层
temp_prec <- c(temp, prec_monthly)
names(temp_prec) <- c("temperature", "precipitation")

# 2. 中国边界、裁剪与投影 -----------------------------------------------------
country_sf <- geodata::gadm(country = "CHN", level = 0, path = main_dir) |>
  sf::st_as_sf()

# 采用 Web Mercator（通用网络地图投影）
target_crs <- "EPSG:3857"

# 裁剪 + 掩膜
temp_prec_country <- terra::crop(temp_prec, country_sf, mask = TRUE)

# 重投影（双线性插值）
temp_prec_proj <- terra::project(temp_prec_country, target_crs, method = "bilinear")

# 转数据框以便 ggplot 绘制（包含 x, y 坐标）
temp_prec_df <- as.data.frame(temp_prec_proj, xy = TRUE, na.rm = TRUE)

# 3. 双变量分类与配色 ---------------------------------------------------------
breaks <- biscale::bi_class(
  temp_prec_df,
  x = temperature, y = precipitation,
  style = "fisher", dim = 3
)

pal <- "DkBlue"

theme_for_the_win <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(color = "grey10", hjust = .5, face = "bold", vjust = -1),
      plot.subtitle = ggplot2::element_text(hjust = .5, vjust = -1),
      plot.caption = ggplot2::element_text(size = 9, color = "grey20", hjust = .5, vjust = 1),
      plot.margin = grid::unit(c(0, 0, 0, 0), "lines")
    )
}

# 4. 绘制 2D 双变量图 + 图例 --------------------------------------------------
map <- ggplot2::ggplot(breaks) +
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = bi_class), show.legend = TRUE) +
  biscale::bi_scale_fill(pal = pal, dim = 3, flip_axes = TRUE, rotate_pal = FALSE) +
  ggplot2::labs(
    title = "中国：温度与降水（双变量）",
    subtitle = "1981–2010 年：年均温与月均降水",
    caption = "数据：CHELSA | 制图：Your Name",
    x = "", y = ""
  ) +
  ggplot2::coord_sf(crs = target_crs) +
  theme_for_the_win()

legend <- biscale::bi_legend(
  pal = pal, flip_axes = TRUE, rotate_pal = FALSE,
  dim = 3, xlab = "温度 (°C)", ylab = "月均降水 (mm)", size = 8
)

full_map <- cowplot::ggdraw() +
  cowplot::draw_plot(plot = map, x = 0, y = 0, width = 1, height = 1) +
  cowplot::draw_plot(plot = legend, x = .05, y = .13, width = .25, height = .25)

print(full_map)

# 5. 导出 PNG -----------------------------------------------------------------
ggplot2::ggsave(
  filename = "china_bivariate_2d.png",
  plot = full_map,
  width = 9, height = 8, dpi = 600, bg = "white"
)

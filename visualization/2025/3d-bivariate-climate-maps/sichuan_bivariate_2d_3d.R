## 中国四川区域：降水-温度二维/三维双变量地图（CHELSA + biscale + rayshader）

# 0) 依赖安装与加载 -------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# 按需安装 GitHub 包（rchelsa 与 biscale）
if (!requireNamespace("rchelsa", quietly = TRUE)) remotes::install_github("inSileco/rchelsa")
if (!requireNamespace("biscale", quietly = TRUE)) remotes::install_github("chris-prener/biscale")

pacman::p_load(
  geodata, tidyverse, sf, terra,
  rchelsa, biscale, elevatr, cowplot,
  rayshader
)

# 1) 数据准备（CHELSA 温度/降水）------------------------------------------------
main_dir <- getwd()

# BIO1（年均温）与 BIO12（年降水）
ids <- c(1, 12)

download_chelsa_data <- function(id, path) {
  # 兼容不同 rchelsa 版本函数名
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

# 载入 CHELSA 生物气候栅格（1981–2010 基期，bio10 前缀）
temp <- terra::rast(file.path(main_dir, "CHELSA_bio10_01.tif"))  # 年均温（°C 或缩放单位）
prec <- terra::rast(file.path(main_dir, "CHELSA_bio10_12.tif"))  # 年降水（mm/年）

# 转为“月均降水”（mm/月）
prec_monthly <- prec / 12

# 合并为两层并命名
temp_prec <- c(temp, prec_monthly)
names(temp_prec) <- c("temperature", "precipitation")

# 2) 行政区边界（中国四川）-----------------------------------------------------
# 获取中国一级行政区（省级）边界
adm1 <- geodata::gadm(country = "CHN", level = 1, path = main_dir) |>
  sf::st_as_sf()

# 在 NAME_1 字段中匹配四川（大小写不敏感）
name_col <- if ("NAME_1" %in% names(adm1)) "NAME_1" else names(adm1)[grepl("NAME_1", names(adm1), ignore.case = TRUE)][1]
if (is.na(name_col) || is.null(name_col)) stop("未找到省级名称字段（NAME_1）。")

sichuan_sf <- adm1[grepl("sichuan", tolower(adm1[[name_col]])), , drop = FALSE]
if (nrow(sichuan_sf) == 0) stop("未在 GADM level 1 中匹配到四川省（Sichuan）。")

# 3) 裁剪、DEM、重投影与网格对齐 -----------------------------------------------
# 目标投影：Web Mercator（EPSG:3857）
target_crs <- "EPSG:3857"

# 按四川边界先裁剪气候（原始坐标系）
temp_prec_sichuan <- terra::crop(temp_prec, sichuan_sf, mask = TRUE)

# DEM 缓存到本地，避免重复下载（原始坐标系）
dem_path <- file.path(main_dir, "dem_sichuan_z8.tif")
if (file.exists(dem_path)) {
  dem <- terra::rast(dem_path)
} else {
  dem <- elevatr::get_elev_raster(locations = sichuan_sf, z = 8, clip = "locations") |>
    terra::rast() |>
    terra::crop(sichuan_sf, mask = TRUE)
  terra::writeRaster(dem, dem_path, overwrite = TRUE)
}

# 将四川边界投影至目标坐标系
sichuan_proj <- sf::st_transform(sichuan_sf, target_crs)

# 分别重投影至目标坐标系
temp_prec_proj <- terra::project(temp_prec_sichuan, target_crs, method = "bilinear")
dem_proj <- terra::project(dem, target_crs)

# 使用投影后的四川边界统一裁剪与掩膜，确保范围完全一致
temp_prec_proj_crop <- terra::crop(temp_prec_proj, sichuan_proj, mask = TRUE)
dem_proj_crop <- terra::crop(dem_proj, sichuan_proj, mask = TRUE)

# 在相同坐标系下将气候格网对齐到 DEM 格网（双线性），确保像元网格一致
temp_prec_aligned <- terra::resample(temp_prec_proj_crop, dem_proj_crop, method = "bilinear")

# 转为数据框用于 ggplot（包含坐标）
temp_prec_df <- as.data.frame(temp_prec_aligned, xy = TRUE, na.rm = TRUE)

# 几何一致性断言（像元大小/范围/对齐）
terra::compareGeom(temp_prec_aligned, dem_proj_crop, stopOnError = TRUE)

# 统一绘图范围（基于 DEM 裁剪结果）
dem_ext <- terra::ext(dem_proj_crop)
xlim <- c(dem_ext$xmin, dem_ext$xmax)
ylim <- c(dem_ext$ymin, dem_ext$ymax)

# 4) 双变量分类、配色与绘图 -----------------------------------------------------
breaks <- biscale::bi_class(
  temp_prec_df,
  x = temperature,
  y = precipitation,
  style = "fisher",
  dim = 3
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

# 极简面板主题（供 3D 贴图/高度渲染，去除所有装饰）
theme_panel <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))
}

map <- ggplot2::ggplot(breaks) +
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = bi_class), show.legend = TRUE) +
  biscale::bi_scale_fill(pal = pal, dim = 3, flip_axes = TRUE, rotate_pal = FALSE) +
  ggplot2::labs(
    title = "中国四川：温度与降水（双变量）",
    subtitle = "1981–2010 年：年均温与月均降水",
    caption = "数据：CHELSA | 制图：Your Name",
    x = "", y = ""
  ) +
  ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  ggplot2::coord_equal() +
  theme_for_the_win()

legend <- biscale::bi_legend(
  pal = pal,
  flip_axes = TRUE,
  rotate_pal = FALSE,
  dim = 3,
  xlab = "温度 (°C)",
  ylab = "月均降水 (mm)",
  size = 8
)

full_map <- cowplot::ggdraw() +
  cowplot::draw_plot(plot = map, x = 0, y = 0, width = 1, height = 1) +
  cowplot::draw_plot(plot = legend, x = .05, y = .13, width = .25, height = .25)

print(full_map)

# 保存二维图
ggplot2::ggsave(
  filename = "sichuan_bivariate_2d.png",
  plot = full_map,
  width = 9, height = 8, dpi = 600, bg = "white"
)

# 5) 生成贴图与地形“纯面板”------------------------------------------------------
# 3D 贴图使用无标题/无图例的面板，避免 grob 引入的边距导致错位
map_texture <- ggplot2::ggplot(breaks) +
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = bi_class)) +
  biscale::bi_scale_fill(pal = pal, dim = 3, flip_axes = TRUE, rotate_pal = FALSE) +
  ggplot2::guides(fill = "none") +
  ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  ggplot2::coord_equal() +
  theme_panel()

dem_df <- dem_proj_crop |>
  as.data.frame(xy = TRUE, na.rm = TRUE)
names(dem_df)[3] <- "dem"

dem_panel <- ggplot2::ggplot(dem_df, ggplot2::aes(x = x, y = y, fill = dem)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_gradientn(colors = "white") +
  ggplot2::guides(fill = "none") +
  ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  ggplot2::coord_equal() +
  theme_panel()

# 6) 渲染 3D 场景 --------------------------------------------------------------
rayshader::plot_gg(
  ggobj = map_texture,
  ggobj_height = dem_panel,
  width = 7,
  height = 7,
  windowsize = c(800, 800),
  scale = 120,
  shadow = TRUE,
  shadow_intensity = 1,
  phi = 87, theta = 0, zoom = .60,
  multicore = TRUE
)

# 调整相机（视角/缩放）
rayshader::render_camera(zoom = .65)

# 7) 环境光与高质量渲染 --------------------------------------------------------
url <- "https://dl.polyhaven.org/file/ph-assets/HDRIs/hdr/4k/brown_photostudio_02_4k.hdr"
hdri_file <- basename(url)
if (!file.exists(hdri_file)) {
  download.file(url = url, destfile = hdri_file, mode = "wb")
} else {
  message("HDRI 已存在，跳过下载：", hdri_file)
}

rayshader::render_highquality(
  filename = "sichuan-bivariate-3d.png",
  preview = TRUE,
  light = FALSE,
  environment_light = hdri_file,
  intensity = 1,
  rotate_env = 90,
  parallel = TRUE,
  width = 1200, height = 1200,
  interactive = FALSE
)

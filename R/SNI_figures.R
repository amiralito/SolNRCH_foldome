#!/usr/bin/env Rscript
################################################################################
# SNI figures — v3 (manuscript style)
#
#   Fig 1A / S6A / S6B   benchmark heatmaps (+ *_2way: parameters clustered too)
#   Fig 2A / S7A / S7B   clade heatmaps     (+ *_2way)
#   Fig 2B               contact counts with the NB-ARC tree drawn in the same
#                        coordinate system (tips aligned to bars, clade supports)
#   Fig S4               per-sequence SNI values + LOESS, tree aligned underneath
#   ParamClust_*         parameter correlation heatmap + bootstrapped dendrogram
#   v1_vs_v2_*           diagnostics
#
# Typography: Helvetica Neue. With the font file available (--font, or the macOS
# system font) every PDF embeds it and SVG/PNG use it; otherwise PDFs reference
# "HelveticaNeue" by name and Greek symbols fall back to the Symbol font.
# All panels are drawn in ggplot2 and assembled with patchwork on shared
# coordinates, so dendrogram leaves, tree tips, tiles, numbers and labels align.
#
# Usage:
#   Rscript SNI_figures.R --base /path/to/analysis_v2 \
#       [--main main/resistosome_analysis_summary.csv] \
#       [--benchmark benchmark/resistosome_analysis_summary.csv] \
#       [--trees inputs/trees]   folder with the rooted NB-ARC tree and clades/*.tree
#       [--v1_main ...] [--v1_benchmark ...] [--out figures] [--nboot 1000] \
#       [--font /System/Library/Fonts/HelveticaNeue.ttc]
################################################################################

invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(ape); library(patchwork); library(svglite)
})

# ------------------------------------------------------------------ arguments
argv <- commandArgs(trailingOnly = TRUE)
opt <- list(base = ".", main = "main/resistosome_analysis_summary.csv",
            benchmark = "benchmark/resistosome_analysis_summary.csv",
            v1_main = NA, v1_benchmark = NA, out = "figures", trees = "inputs/trees", nboot = "1000")
i <- 1
while (i <= length(argv)) { opt[[sub("^--", "", argv[i])]] <- argv[i + 1]; i <- i + 2 }
rp <- function(p) if (is.na(p) || grepl("^/", p)) p else file.path(opt$base, p)
OUT <- rp(opt$out); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NBOOT <- as.integer(opt$nboot)

# ------------------------------------------------------------------ typography
# If the Helvetica Neue font file is available (opt --font, or installed), every
# output is rendered with it through cairo (PDF embeds the real font; Greek
# symbols as bold-italic Helvetica Neue glyphs, exactly as in the manuscript).
# Otherwise PDFs reference "HelveticaNeue" by name and Greek falls back to Symbol.
FAM <- "Helvetica Neue"
if (!is.null(opt$font)) opt$font <- rp(opt$font)
if (!is.null(opt$font) && file.exists(opt$font) && Sys.info()[["sysname"]] != "Darwin") {
  fd <- file.path(Sys.getenv("HOME"), ".fonts"); dir.create(fd, showWarnings = FALSE)
  dest <- file.path(fd, basename(opt$font))
  if (normalizePath(opt$font) != normalizePath(dest, mustWork = FALSE)) file.copy(opt$font, dest, overwrite = TRUE)
  system("fc-cache -f >/dev/null 2>&1")
}
REAL_FONT <- if (Sys.info()[["sysname"]] == "Darwin") {         # macOS: cairo uses the system font directly
  file.exists("/System/Library/Fonts/HelveticaNeue.ttc") || (!is.null(opt$font) && file.exists(opt$font))
} else isTRUE(any(grepl("HelveticaNeue|Helvetica Neue",
                        tryCatch(suppressWarnings(system("fc-match -f '%{file}' 'Helvetica Neue:bold:italic'",
                                                         intern = TRUE, ignore.stderr = TRUE)),
                                 error = function(e) ""))))
cat("Helvetica Neue font file:", if (REAL_FONT) "found (cairo, embedded)" else "not found (Type1 reference)", "\n")
if (!REAL_FONT) local({
  src <- system.file("afm", package = "grDevices"); d <- file.path(tempdir(), "afm")
  dir.create(d, showWarnings = FALSE)
  mk <- function(f, ps) {
    l <- readLines(gzfile(file.path(src, f)))
    l <- sub("^FontName .*", paste("FontName", ps), l); l <- sub("^FullName .*", paste("FullName", ps), l)
    l <- sub("^FamilyName .*", "FamilyName Helvetica Neue", l)
    o <- file.path(d, paste0(ps, ".afm")); writeLines(l, o); o
  }
  afm <- c(mk("Helvetica.afm.gz", "HelveticaNeue"), mk("Helvetica-Bold.afm.gz", "HelveticaNeue-Bold"),
           mk("Helvetica-Oblique.afm.gz", "HelveticaNeue-Italic"),
           mk("Helvetica-BoldOblique.afm.gz", "HelveticaNeue-BoldItalic"), file.path(src, "Symbol.afm.gz"))
  pdfFonts("Helvetica Neue" = Type1Font("HelveticaNeue", afm))
})
dev_pdf <- function(filename, width, height, bg = "white", ...) {
  if (REAL_FONT) cairo_pdf(filename, width = width, height = height, family = FAM, bg = bg)
  else pdf(file = filename, width = width, height = height, family = FAM, bg = bg, useDingbats = FALSE)
}
dev_png <- function(filename, width, height, res = 300, bg = "white", ...)
  png(filename, width = width, height = height, units = "in", res = res, bg = bg, type = "cairo", family = FAM)

save_fig <- function(p, base, w, h) {
  suppressWarnings({
    ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = dev_pdf, limitsize = FALSE)
    ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE,
           device = dev_png)
    ggsave(paste0(base, ".svg"), p, width = w, height = h, device = svglite, limitsize = FALSE)
  })
  s <- readLines(paste0(base, ".svg"), warn = FALSE)                  # live text in Helvetica Neue
  s <- gsub('font-family: "[^"]*"', 'font-family: "Helvetica Neue"', s)
  s <- gsub(" textLength='[^']*' lengthAdjust='[^']*'", "", s)
  writeLines(s, paste0(base, ".svg"))
}

TXT <- "#1f1f1f"; MUTED <- "#6b6b6b"
theme_sni <- function(size = 10) theme_void(base_family = FAM, base_size = size) +
  theme(text = element_text(colour = TXT), plot.margin = margin(0, 0, 0, 0))

# ------------------------------------------------------------------ constants
SCIENTIST_ONLY   <- c("ipTM_LCB", "Sum_CONTACTS", "S_PROTO", "L_APEX_mean")
OVERLAPPING      <- c("BSA_INT_PROTO", "SD_THETA_ROT", "D_APEX")
COSCIENTIST_ONLY <- c("D_MHD_P_mean", "THETA_APEX_mean", "H_ABS_mean", "MU_H_mean")
PARAM_ORDER <- c(SCIENTIST_ONLY, OVERLAPPING, COSCIENTIST_ONLY)
SCI_SET  <- c(SCIENTIST_ONLY, OVERLAPPING); COSI_SET <- c(OVERLAPPING, COSCIENTIST_ONLY)
PARAM_SOURCE <- setNames(c(rep("Scientist", 4), rep("Both", 3), rep("Co-scientist", 4)), PARAM_ORDER)
SOURCE_COLOURS <- c("Scientist" = "#a6a8ab", "Co-scientist" = "#3f6293")      # manuscript Fig 1 brackets
BY_SCI <- c(SCIENTIST_ONLY, OVERLAPPING); BY_COSI <- c(OVERLAPPING, COSCIENTIST_ONLY)
CATEGORY <- c(ipTM_LCB = "Confidence & interface", Sum_CONTACTS = "Confidence & interface",
              BSA_INT_PROTO = "Confidence & interface", SD_THETA_ROT = "Ring symmetry",
              S_PROTO = "Ring symmetry", D_APEX = "Ring symmetry",
              THETA_APEX_mean = "CC-domain helix", L_APEX_mean = "CC-domain helix",
              H_ABS_mean = "CC-domain helix", MU_H_mean = "CC-domain helix",
              D_MHD_P_mean = "Motif distance")
CATEGORY_COLOURS <- c("Confidence & interface" = "#003049", "Ring symmetry" = "#d52a28",
                      "CC-domain helix" = "#f78e1e", "Motif distance" = "#ffd641")
# Mathematical symbols, manuscript style (Fig 1): bold-italic symbol and subscript,
# summation sign upright, units roman.
G <- if (REAL_FONT) list(sigma = "bolditalic('σ')", theta = "θ", phi = "bolditalic('φ')", mu = "bolditalic('μ')",
                         sum = "'∑'") else list(sigma = "bold(sigma)", theta = "theta", phi = "bold(phi)", mu = "bold(mu)", sum = "symbol('\\345')")
bi <- function(x) sprintf("bolditalic('%s')", x)
SYM <- c(ipTM_LCB      = paste0(bi("ipTM"), "[", bi("LCB"), "]"),
         Sum_CONTACTS  = paste0(G$sum, "[", bi("CONTACTS"), "]"),
         S_PROTO       = paste0(bi("S"), "[", bi("PROTO"), "]"),
         L_APEX_mean   = paste0(bi("L"), "[", bi("APEX"), "]"),
         BSA_INT_PROTO = paste0(bi("BSA"), "[", bi("INTER-PROTO"), "]"),
         SD_THETA_ROT  = paste0(G$sigma, "[", if (REAL_FONT) bi("θ ROT") else "bold(theta)~bolditalic('ROT')", "]"),
         D_APEX        = paste0(bi("D"), "[", bi("APEX"), "]"),
         D_MHD_P_mean  = paste0(bi("D"), "[", bi("MHD-P"), "]"),
         THETA_APEX_mean = paste0(G$phi, "[", bi("APEX"), "]"),
         H_ABS_mean    = paste0(bi("H"), "[", bi("ABS"), "]"),
         MU_H_mean     = paste0(G$mu, "[", bi("H"), "]"))
UNIT <- c(S_PROTO = "~'(Å)'", BSA_INT_PROTO = "~'(Å'^2*')'", SD_THETA_ROT = "~'(°)'",
          D_APEX = "~'(Å)'", D_MHD_P_mean = "~'(Å)'", THETA_APEX_mean = "~'(°)'")
lab_expr <- function(p, units = TRUE) {
  s <- unname(SYM[p]); if (units) s <- paste0(s, ifelse(is.na(UNIT[p]), "", UNIT[p]))
  parse(text = s)
}
plain_label <- c(ipTM_LCB = "ipTM_LCB", Sum_CONTACTS = "Σ_CONTACTS", S_PROTO = "S_PROTO",
                 L_APEX_mean = "L_APEX", BSA_INT_PROTO = "BSA_INTER-PROTO", SD_THETA_ROT = "σθ_ROT",
                 D_APEX = "D_APEX", D_MHD_P_mean = "D_MHD-P", THETA_APEX_mean = "φ_APEX",
                 H_ABS_mean = "H_ABS", MU_H_mean = "μ_H")

CLAMP_COLS <- c("Sum_CONTACTS", "BSA_INT_PROTO", "ipTM_LCB", "CONTACTS_HBOND", "CONTACTS_SALT_BRIDGE",
                "CONTACTS_HYDROPHOBIC", "CONTACTS_DISULFIDE", "CONTACTS_VDW")
HEAT <- c("#1C39BB", "white", "#ED2939")
NA_FILL <- "#DCDCDC"; REF_RED <- "#E8202A"

BENCH_LABELS <- c(
  nbnrc2a = "NbNRC2a", nrc2 = "SlNRC2", nrc3 = "SlNRC3", nbnrc4c = "NbNRC4c", nbnrc3 = "NbNRC3",
  gos_00000000001_capsicum_annuum_helixer_chr00_010241_1 = "CaNRC2",
  gca_000710875_1_capsicum_annuum_helixer_cm002822_1_001604_1 = "CaNRC3",
  gcf_019175385_1_lycium_barbarum_helixer_nc_083346_1_002027_1 = "LbNRC4",
  gca_000715095_1_nicotiana_tabacum_helixer_awok01s169090_1_000001_1_cnl = "NtNRC4",
  rx = "Rx", rx4 = "Rx4", rpi_amr1e = "Rpi-amr1e", rpi_amr3 = "Rpi-amr3",
  rpi_nig3_1084a = "Rpi-nig3", rpa1 = "Rpa1", gpa2 = "Gpa2", cnl_11990 = "CNL-11990", bs2 = "Bs2")
BENCH_SENSORS <- c("rx", "rx4", "rpi_amr1e", "rpi_amr3", "rpi_nig3_1084a", "rpa1", "gpa2", "cnl_11990", "bs2")
BENCH_REFS <- c("NbNRC2a", "SlNRC3", "NbNRC4c")         # experimentally resolved (red labels)
MAIN_REFS  <- c("NRC2", "NRC3", "NRC4other")

# ------------------------------------------------------------------ data helpers
# IDs are matched after lower-casing and mapping every non-alphanumeric character to "_"
# (AF3 job names turn "." and "-" into "_", e.g. tree tip "...E86-69_chr01_001162.1_CNL").
norm_id <- function(x) tolower(gsub("[^A-Za-z0-9]", "_", x))
read_summary <- function(path) read_csv(path, show_col_types = FALSE) %>% mutate(ID_norm = norm_id(ID))
clamp <- function(df) df %>% mutate(across(any_of(CLAMP_COLS), ~ pmax(.x, 0)))
read_clades <- function(dir) {
  map_dfr(list.files(dir, "\\.tree$", full.names = TRUE), function(f)
    tibble(ID_original = read.tree(f)$tip.label, clade = tools::file_path_sans_ext(basename(f)))) %>%
    mutate(ID_norm = norm_id(ID_original)) %>% distinct(ID_norm, .keep_all = TRUE)
}
group_matrix <- function(df, group_col, params) {
  params <- intersect(params, names(df))
  gm <- df %>% group_by(.data[[group_col]]) %>%
    summarise(n_members = n(), across(all_of(params), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  raw <- gm %>% select(-n_members) %>% column_to_rownames(group_col) %>% as.matrix()
  raw[is.nan(raw)] <- NA
  raw <- raw[, apply(raw, 2, function(v) isTRUE(var(v, na.rm = TRUE) > 0)), drop = FALSE]
  list(raw = raw, z = scale(raw), n = setNames(gm$n_members, gm[[group_col]]))
}
fmt_val <- function(p, v) {
  d <- ifelse(p %in% c("H_ABS_mean", "S_PROTO", "MU_H_mean", "ipTM_LCB"), 2,
              ifelse(p == "BSA_INT_PROTO", 0, 1))
  ifelse(is.na(v), "", sprintf(paste0("%.", d, "f"), v))
}

# ------------------------------------------------------------------ dendrogram geometry
# Leaves at integer positions 1..n in hc$order; returns segment table (pos, height)
dend_segments <- function(hc) {
  n <- length(hc$order); pos <- numeric(n); pos[hc$order] <- seq_len(n)
  m <- hc$merge; h <- hc$height; nx <- numeric(nrow(m)); seg <- vector("list", nrow(m))
  gx <- function(k) if (k < 0) pos[-k] else nx[k]
  gh <- function(k) if (k < 0) 0 else h[k]
  for (j in seq_len(nrow(m))) {
    a <- m[j, 1]; b <- m[j, 2]; xa <- gx(a); xb <- gx(b); nx[j] <- (xa + xb) / 2
    seg[[j]] <- data.frame(p0 = c(xa, xb, xa), h0 = c(gh(a), gh(b), h[j]),
                           p1 = c(xa, xb, xb), h1 = c(h[j], h[j], h[j]))
  }
  list(seg = bind_rows(seg), node_pos = nx, node_h = h)
}

# ------------------------------------------------------------------ parameter annotation
# Labels (symbols) + slim category bar + two source tracks: a grey dot if proposed by
# the scientists, a blue dot if proposed by the AI co-scientist (both dots = both).
SRC_TITLE <- "Proposed by"
param_marks <- function(params, y, ylim, units = TRUE, label_size = 10.5, order0 = 3) {
  y <- unname(y); a <- tibble(param = params, y = y, cat = CATEGORY[params])
  dots <- bind_rows(tibble(y = y[params %in% BY_SCI], x = 1.82, who = "Scientist"),
                    tibble(y = y[params %in% BY_COSI], x = 2.22, who = "Co-scientist"))
  ggplot() +
    geom_tile(data = a, aes(x = 1.2, y = y, fill = cat), width = 0.3, height = 0.9) +
    geom_point(data = tibble(x = rep(c(1.82, 2.22), each = length(y)), y = rep(unname(y), 2)), aes(x, y),
               shape = 21, size = 1.45, colour = "#e2e2e2", fill = NA, stroke = 0.35) +
    geom_point(data = dots, aes(x, y, colour = who), size = 1.45) +
    scale_fill_manual(values = CATEGORY_COLOURS, name = NULL, breaks = names(CATEGORY_COLOURS),
                      guide = guide_legend(order = order0, nrow = 1, keywidth = unit(3, "mm"),
                                           keyheight = unit(3, "mm"))) +
    scale_colour_manual(values = SOURCE_COLOURS, name = SRC_TITLE, breaks = names(SOURCE_COLOURS),
                        guide = guide_legend(order = order0 + 1, nrow = 1, override.aes = list(size = 2))) +
    scale_x_continuous(limits = c(0.95, 2.45), expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0), breaks = a$y, labels = lab_expr(a$param, units)) +
    theme_sni() + theme(axis.text.y.left = element_text(size = label_size, hjust = 1, margin = margin(r = 4)),
                        legend.text = element_text(size = 7.5, colour = MUTED),
                        legend.title = element_text(size = 7.5, colour = MUTED))
}

# ------------------------------------------------------------------ SNI heatmap
# z, raw: groups x params ; returns invisible(list(hc, order))
sni_heatmap <- function(z, raw, params, file_base, col_labels = NULL, col_strip = NULL,
                        strip_cols = NULL, strip_name = NULL, ref_labels = NULL,
                        cluster_params = FALSE, legend_title = "z-score") {
  params <- intersect(PARAM_ORDER, intersect(params, colnames(z)))
  zz <- t(z[, params, drop = FALSE]); rr <- t(raw[, params, drop = FALSE])
  hc <- hclust(dist(t(zz)), method = "ward.D2")
  ng <- ncol(zz); np <- length(params)
  gx <- setNames(match(colnames(zz), colnames(zz)[hc$order]), colnames(zz))   # column -> x
  if (cluster_params) {
    r <- suppressWarnings(cor(t(zz), use = "pairwise.complete.obs")); r[is.na(r)] <- 0
    hr <- hclust(as.dist(1 - r), method = "ward.D2")
    py <- setNames(np + 1 - match(params, params[hr$order]), params)
  } else {                                             # Scientist / Both / Co-scientist blocks
    src <- factor(PARAM_SOURCE[params], levels = c("Scientist", "Both", "Co-scientist"))
    gap <- cumsum(c(0, diff(as.integer(src)) != 0)) * 0.3
    py <- setNames(np - seq_len(np) + 1 - gap + max(gap), params)
  }
  ytop <- max(py) + 0.5
  if (!is.null(col_strip)) ytop <- ytop + 0.62                 # annotation strip inside the panel
  ylim <- c(min(py) - 0.5, ytop); xlim <- c(0.5, ng + 0.5)
  cl <- if (is.null(col_labels)) colnames(zz) else unname(col_labels[colnames(zz)])
  cl_ord <- cl[hc$order]
  df <- expand.grid(param = params, group = colnames(zz), stringsAsFactors = FALSE) %>%
    mutate(x = gx[group], y = py[param], z = zz[cbind(param, group)], v = rr[cbind(param, group)],
           lab = fmt_val(param, v), zc = pmax(pmin(z, 3), -3))

  p_main <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = zc), colour = "white", linewidth = 0.7, width = 1, height = 1) +
    geom_text(aes(label = lab), family = FAM, size = 2.75, colour = TXT) +
    scale_fill_gradientn(colours = HEAT, limits = c(-3, 3), na.value = NA_FILL, name = legend_title,
                         breaks = c(-3, -2, -1, 0, 1, 2, 3),
                         guide = guide_colourbar(direction = "horizontal", title.position = "top",
                                                 barwidth = unit(38, "mm"), barheight = unit(2.2, "mm"),
                                                 ticks = FALSE, order = 1)) +
    scale_x_continuous(limits = xlim, expand = c(0, 0), breaks = seq_len(ng), labels = cl_ord,
                       position = "top") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    theme_sni() + theme(legend.title = element_text(size = 8, colour = MUTED),
                        legend.text = element_text(size = 7, colour = MUTED),
                        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10,
                                                       colour = ifelse(cl_ord %in% ref_labels, REF_RED, TXT),
                                                       margin = margin(b = 3)))
  if (!is.null(col_strip)) {
    sv <- unname(col_strip[colnames(zz)[hc$order]])
    if (is.numeric(sv)) {
      ramp <- colorRampPalette(c("#f2f2f2", "#3a3a3a"))(101)
      sfill <- ramp[1 + round(100 * (sv - min(sv)) / max(1e-9, diff(range(sv))))]
      sdf <- tibble(x = seq_len(ng), fill = sfill, key = sv)
      p_main <- p_main +
        geom_tile(data = sdf, aes(x = x, y = ytop - 0.27), fill = sdf$fill, height = 0.42, width = 1,
                  colour = "white", linewidth = 0.7, inherit.aes = FALSE) +
        geom_text(data = sdf, aes(x = x, y = ytop - 0.27, label = key), family = FAM, size = 2.2,
                  colour = ifelse(sv > median(sv), "white", TXT), inherit.aes = FALSE)
    } else {
      sdf <- tibble(x = seq_len(ng), v = sv)
      p_main <- p_main +
        geom_tile(data = sdf, aes(x = x, y = ytop - 0.27), fill = strip_cols[sdf$v], height = 0.42,
                  width = 1, colour = "white", linewidth = 0.7, inherit.aes = FALSE) +
        geom_point(data = sdf, aes(x = x, y = ytop - 0.27, shape = v), alpha = 0, inherit.aes = FALSE) +
        scale_shape_manual(values = c(15, 15), name = NULL,
                           guide = guide_legend(order = 2, nrow = 1,
                                                override.aes = list(size = 3.6, alpha = 1, colour = unname(strip_cols))))
    }
  }

  # left: parameter labels + category bar + source tracks
  p_left <- param_marks(params, py[params], ylim)

  # top: column dendrogram (+ optional annotation strip with column labels)
  ds <- dend_segments(hc)
  p_dend <- ggplot(ds$seg) + geom_segment(aes(x = p0, xend = p1, y = h0, yend = h1),
                                          linewidth = 0.35, colour = TXT, lineend = "square") +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) + theme_sni()
  top <- list(p_dend)
  heights <- c(2.2)
  # right: parameter dendrogram (two-way)
  right_main <- plot_spacer(); wr <- 0.05
  if (cluster_params) {
    dr <- dend_segments(hr)
    seg <- dr$seg %>% mutate(y0 = np + 1 - p0, y1 = np + 1 - p1)
    right_main <- ggplot(seg) + geom_segment(aes(x = h0, xend = h1, y = y0, yend = y1),
                                             linewidth = 0.35, colour = TXT, lineend = "square") +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.03))) + theme_sni()
    wr <- 2
  }
  nt <- length(top)
  rows <- list()
  for (k in seq_len(nt)) rows <- c(rows, list(plot_spacer(), top[[k]], plot_spacer()))
  rows <- c(rows, list(p_left, p_main, right_main))
  pw <- wrap_plots(rows, ncol = 3, widths = c(1.1, ng, wr), heights = c(heights, diff(ylim))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.justification = "left", legend.spacing.x = unit(1, "mm"),
          legend.spacing.y = unit(1.5, "mm"),
          legend.margin = margin(t = 2, b = 2), legend.box.margin = margin(t = 6))
  cell <- 0.36                                             # inches per tile
  W <- max(2.2 + (0.6 + ng + wr) * cell, 7.2)
  H <- 1.4 + (sum(heights) + diff(ylim)) * cell + 1.3
  save_fig(pw, file_base, W, H)
  write_csv(as.data.frame(t(rr)) %>% rownames_to_column("group"), paste0(file_base, "_raw_values.csv"))
  write_csv(as.data.frame(t(zz)) %>% rownames_to_column("group"), paste0(file_base, "_zscores.csv"))
  invisible(list(hc = hc, order = colnames(zz)[hc$order]))
}

# ------------------------------------------------------------------ parameter clustering
param_clustering <- function(df, file_base, label, nboot = NBOOT) {
  X <- df %>% select(any_of(PARAM_ORDER)) %>% as.data.frame()
  X <- X[, sapply(X, function(v) sum(!is.na(v)) > 2 && isTRUE(var(v, na.rm = TRUE) > 0)), drop = FALSE]
  ps <- colnames(X); np <- length(ps)
  r <- cor(X, use = "pairwise.complete.obs")
  set.seed(20260503)                                  # reproducible AU/BP
  pv <- pvclust::pvclust(as.matrix(X), method.dist = "correlation", method.hclust = "ward.D2",
                         nboot = nboot, quiet = TRUE, use.cor = "pairwise.complete.obs")
  hc <- pv$hclust                                       # same tree as hclust(as.dist(1 - r), "ward.D2")
  ord <- ps[hc$order]; pos <- setNames(seq_len(np), ord)
  ds <- dend_segments(hc)
  lim <- c(0.5, np + 0.5)

  # (1) bootstrapped dendrogram, horizontal: leaves (rows) carry symbols, category bar and
  #     source tracks; AU | BP printed just above each internal node, inside the panel.
  au <- round(100 * pv$edges$au); bp <- round(100 * pv$edges$bp)
  #     Leaves sit at the left, next to their labels; root at the right. Each AU | BP label
  #     is written inside its own fork (left of the node, between its two child branches),
  #     where no other branch or label can run. Leaf stems start at -w so labels of shallow
  #     nodes stay inside the panel. The root (always 100 | 100 by construction) is unlabelled.
  hmax <- max(ds$node_h); w <- 0.13 * hmax
  seg <- ds$seg %>% mutate(y0 = np + 1 - p0, y1 = np + 1 - p1, h0 = ifelse(h0 == 0 & p0 == p1, -w, h0))
  nodes <- tibble(y = np + 1 - ds$node_pos, h = ds$node_h, au = au, bp = bp) %>%
    mutate(lab = ifelse(is.na(bp), paste0(au), paste0(au, " | ", bp)))
  inner <- head(nodes, -1)                                   # all merges except the root
  ylab <- setNames(np + 1 - seq_len(np), ord)
  p_tree <- ggplot() +
    geom_segment(data = seg, aes(x = h0, xend = h1, y = y0, yend = y1), linewidth = 0.45,
                 colour = TXT, lineend = "square") +
    geom_point(data = nodes, aes(h, y), size = 0.9, colour = TXT) +
    geom_text(data = inner, aes(h, y, label = lab), family = FAM, size = 2.35, colour = MUTED,
              hjust = 1, nudge_x = -0.012 * hmax) +
    scale_y_continuous(limits = lim, expand = c(0, 0)) +
    scale_x_continuous(limits = c(-w, hmax * 1.02), expand = c(0, 0)) +
    theme_sni()
  p_lab <- param_marks(ord, ylab[ord], lim, units = FALSE, label_size = 11, order0 = 1)
  pd <- wrap_plots(list(p_lab, p_tree), ncol = 2, widths = c(0.9, 4.2)) +
    plot_layout(guides = "collect") +
    plot_annotation(caption = sprintf("%s  ·  distance 1 - Pearson r, Ward.D2  ·  node labels: AU | BP (%%), %d bootstraps",
                                      label, nboot),
                    theme = theme(plot.caption = element_text(family = FAM, size = 7.5, colour = MUTED, hjust = 0),
                                  plot.margin = margin(10, 14, 6, 10))) &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.justification = "left", legend.spacing.y = unit(1.5, "mm"), legend.margin = margin(t = 2, b = 2))
  save_fig(pd, paste0(file_base, "_dendrogram"), 7.4, 1.6 + np * 0.36 + 0.9)

  # (2) correlation heatmap ordered by the same tree, dendrogram on top
  cd <- expand.grid(a = ord, b = ord, stringsAsFactors = FALSE) %>%
    mutate(x = pos[a], y = np + 1 - pos[b], v = r[cbind(a, b)], lab = formatC(v, format = "f", digits = 2))
  p_h <- ggplot(cd, aes(x, y)) +
    geom_tile(aes(fill = v), colour = "white", linewidth = 0.7) +
    geom_text(aes(label = lab), family = FAM, size = 2.5, colour = TXT) +
    scale_fill_gradientn(colours = HEAT, limits = c(-1, 1), name = "Pearson r",
                         guide = guide_colourbar(direction = "horizontal", title.position = "top",
                                                 barwidth = unit(32, "mm"), barheight = unit(2.2, "mm"),
                                                 ticks = FALSE, order = 1)) +
    scale_x_continuous(limits = lim, expand = c(0, 0), breaks = seq_len(np), labels = lab_expr(ord, FALSE),
                       position = "bottom") +
    scale_y_continuous(limits = lim, expand = c(0, 0)) +
    theme_sni() + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
                        legend.title = element_text(size = 8, colour = MUTED),
                        legend.text = element_text(size = 7, colour = MUTED))
  p_l <- param_marks(ord, setNames(np + 1 - seq_len(np), ord), lim, units = FALSE, order0 = 2)
  p_t <- ggplot(ds$seg) + geom_segment(aes(x = p0, xend = p1, y = h0, yend = h1), linewidth = 0.35,
                                       colour = TXT, lineend = "square") +
    scale_x_continuous(limits = lim, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) + theme_sni()
  pw <- wrap_plots(list(plot_spacer(), p_t, p_l, p_h), ncol = 2, widths = c(0.75, np), heights = c(2, np)) +
    plot_layout(guides = "collect") +
    plot_annotation(caption = sprintf("%s  ·  n = %d NLRs  ·  ordered by 1 - r, Ward.D2", label, nrow(X)),
                    theme = theme(plot.caption = element_text(family = FAM, size = 7.5, colour = MUTED, hjust = 0.5))) &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.justification = "left", legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1.5, "mm"),
          legend.margin = margin(t = 2, b = 2), legend.box.margin = margin(t = 6))
  save_fig(pw, paste0(file_base, "_correlation"), max(7.2, 2.2 + (0.6 + np) * 0.42), 2.2 + (2 + np) * 0.42 + 1.4)
  write_csv(as.data.frame(r) %>% rownames_to_column("parameter"), paste0(file_base, "_correlation.csv"))
  memb <- character(nrow(hc$merge))
  for (j in seq_len(nrow(hc$merge))) memb[j] <- paste(unlist(lapply(hc$merge[j, ], function(k)
    if (k < 0) hc$labels[-k] else strsplit(memb[k], ";")[[1]])), collapse = ";")
  write_csv(tibble(node = seq_along(memb), members = memb, height = hc$height,
                   AU = round(100 * pv$edges$au), BP = round(100 * pv$edges$bp)),
            paste0(file_base, "_bootstrap.csv"))
  invisible(hc)
}

# ------------------------------------------------------------------ phylogeny geometry
# Rectangular tree with tips at x = 1..N (ladderized, same order as ggtree right = TRUE,
# top-to-bottom), depth upward from the root (y = 0) — tips point up, toward the data.
tree_layout <- function(tr) {
  tr <- ladderize(tr, right = TRUE); n <- Ntip(tr)
  tip_seq <- rev(tr$edge[tr$edge[, 2] <= n, 2])        # tip node ids, left -> right
  x <- numeric(n + tr$Nnode); x[tip_seq] <- seq_len(n)
  dep <- node.depth.edgelength(tr)
  internal <- (n + 1):(n + tr$Nnode)
  # assign internal x as midpoint of children, processing deepest first
  ord_int <- internal[order(-dep[internal])]
  for (v in ord_int) { ch <- tr$edge[tr$edge[, 1] == v, 2]; x[v] <- (min(x[ch]) + max(x[ch])) / 2 }
  edges <- tibble(parent = tr$edge[, 1], child = tr$edge[, 2])
  vseg <- edges %>% transmute(x0 = x[child], y0 = dep[parent], x1 = x[child], y1 = dep[child])
  hseg <- edges %>% group_by(parent) %>% summarise(x0 = min(x[child]), x1 = max(x[child]), .groups = "drop") %>%
    transmute(x0, x1, y0 = dep[parent], y1 = dep[parent])
  list(tr = tr, n = n, x = x, dep = dep, tips = tr$tip.label[tip_seq], seg = bind_rows(vseg, hseg))
}

clade_runs <- function(tips, clades) {
  cl <- clades$clade[match(norm_id(tips), clades$ID_norm)]; cl[is.na(cl)] <- ""
  r <- rle(cl); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  tibble(clade = r$values, x0 = s - 0.5, x1 = e + 0.5, len = r$lengths) %>% filter(clade != "") %>%
    group_by(clade) %>% mutate(label = ifelse(len == max(len), clade, "")) %>% ungroup() %>%
    mutate(shade = rep(c(TRUE, FALSE), length.out = n()))
}

tree_panel <- function(TL, clades, bands, xlim, support_min = 0.8) {
  # bootstrap support at each clade's crown node
  sup <- map_dfr(unique(clades$clade), function(cn) {
    tp <- intersect(clades$ID_original[clades$clade == cn], TL$tr$tip.label)
    if (length(tp) < 2) return(NULL)
    nd <- getMRCA(TL$tr, tp); lab <- TL$tr$node.label[nd - TL$n]
    tibble(clade = cn, x = TL$x[nd], y = TL$dep[nd], lab = lab)
  })
  sup <- if (nrow(sup)) sup %>% mutate(v = suppressWarnings(as.numeric(lab))) %>% filter(!is.na(v), v >= support_min)
         else tibble(x = numeric(), y = numeric(), v = numeric())
  maxd <- max(TL$dep)
  ggplot() +
    geom_rect(data = bands %>% filter(shade), aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf),
              fill = "#f0f0f0") +
    geom_segment(data = TL$seg, aes(x = x0, xend = x1, y = y0, yend = y1), linewidth = 0.18, colour = TXT,
                 lineend = "square") +
    geom_point(data = sup, aes(x, y), colour = REF_RED, size = 0.9) +
    geom_text(data = sup, aes(x, y, label = formatC(v, format = "fg", digits = 2)), family = FAM,
              size = 2.3, vjust = 1.6, colour = TXT) +
    annotate("segment", x = 2, xend = 2, y = maxd * 0.55, yend = maxd * 0.55 + 0.5, linewidth = 0.4) +
    annotate("text", x = 4, y = maxd * 0.55 + 0.25, label = "0.5", hjust = 0, family = FAM, size = 2.6) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.01))) + theme_sni()
}

clade_label_panel <- function(bands, xlim, refs) {
  ggplot(bands) +
    geom_rect(data = bands %>% filter(shade), aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf), fill = "#f0f0f0") +
    geom_text(data = bands %>% filter(label != ""), aes(x = (x0 + x1) / 2, y = 1, label = label,
              colour = label %in% refs), family = FAM, angle = 90, hjust = 1, size = 3.3) +
    scale_colour_manual(values = c(`FALSE` = TXT, `TRUE` = REF_RED), guide = "none") +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    theme_sni()
}

# ================================================================== BENCHMARK
bench_path <- rp(opt$benchmark)
if (!is.na(bench_path) && file.exists(bench_path)) {
  cat("== Benchmark\n")
  b <- read_summary(bench_path) %>% filter(ID_norm %in% names(BENCH_LABELS)) %>% clamp() %>% mutate(group = ID_norm)
  G <- group_matrix(b, "group", PARAM_ORDER)
  strip <- setNames(ifelse(rownames(G$raw) %in% BENCH_SENSORS, "NRC-S (sensor)", "NRC (helper)"), rownames(G$raw))
  sc <- c("NRC (helper)" = "#e3e5e3", "NRC-S (sensor)" = "#b5b5b5")
  args <- list(col_labels = BENCH_LABELS, col_strip = strip, strip_cols = sc, ref_labels = BENCH_REFS)
  for (cp in c(FALSE, TRUE)) {
    sfx <- if (cp) "_2way" else ""
    do.call(sni_heatmap, c(list(G$z, G$raw, PARAM_ORDER, file.path(OUT, paste0("Fig1A_benchmark_SNI", sfx)), cluster_params = cp), args))
    do.call(sni_heatmap, c(list(G$z, G$raw, SCI_SET, file.path(OUT, paste0("FigS6A_benchmark_scientist", sfx)), cluster_params = cp), args))
    do.call(sni_heatmap, c(list(G$z, G$raw, COSI_SET, file.path(OUT, paste0("FigS6B_benchmark_coscientist", sfx)), cluster_params = cp), args))
  }
  param_clustering(b, file.path(OUT, "ParamClust_benchmark"), "Benchmark (18 NLRs)")
}

# ================================================================== MAIN SET
main_path <- rp(opt$main)
# NB-ARC tree with bootstrap support (same topology as the iTOL export, which lacks supports)
itol <- list.files(rp(opt$trees), pattern = "_rooted\\.tree$", full.names = TRUE)
itol <- c(itol[!grepl("_ref_", basename(itol))], NA)[1]          # NRC-H tree, not the "+ reference" tree
if (is.na(itol)) itol <- list.files(rp(opt$trees), pattern = "itol\\.txt$", full.names = TRUE)[1]
if (!is.na(main_path) && file.exists(main_path)) {
  cat("== Main\n")
  clades <- read_clades(rp(file.path(opt$trees, "clades")))
  m_all <- read_summary(main_path) %>% clamp()
  m <- m_all %>% inner_join(clades %>% select(ID_norm, clade, ID_original), by = "ID_norm")
  G <- group_matrix(m, "clade", PARAM_ORDER)
  args <- list(col_strip = G$n[rownames(G$raw)], strip_name = "Clade size (n)", ref_labels = MAIN_REFS,
               legend_title = "Clade mean (z-score)")
  for (cp in c(FALSE, TRUE)) {
    sfx <- if (cp) "_2way" else ""
    do.call(sni_heatmap, c(list(G$z, G$raw, PARAM_ORDER, file.path(OUT, paste0("Fig2A_clade_SNI", sfx)), cluster_params = cp), args))
    do.call(sni_heatmap, c(list(G$z, G$raw, SCI_SET, file.path(OUT, paste0("FigS7A_clade_scientist", sfx)), cluster_params = cp), args))
    do.call(sni_heatmap, c(list(G$z, G$raw, COSI_SET, file.path(OUT, paste0("FigS7B_clade_coscientist", sfx)), cluster_params = cp), args))
  }
  param_clustering(m, file.path(OUT, "ParamClust_main"), "Main set (per NLR)")

  # ---------- Fig 2B / S4: tree drawn on the same x as the per-sequence data ----------
  if (!is.na(itol) && file.exists(itol)) {
    TL <- tree_layout(read.tree(itol)); xlim <- c(0.5, TL$n + 0.5)
    bands <- clade_runs(TL$tips, clades)
    seqdf <- m_all %>% left_join(clades %>% select(ID_norm, clade), by = "ID_norm") %>%
      mutate(x = match(ID_norm, norm_id(TL$tips))) %>% filter(!is.na(x))
    cat("sequences on tree:", nrow(seqdf), "of", nrow(m_all), "; tree tips:", TL$n, "\n")
    shade_layer <- geom_rect(data = bands %>% filter(shade), inherit.aes = FALSE,
                             aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf), fill = "#f0f0f0")
    ct <- c("CONTACTS_HBOND", "CONTACTS_SALT_BRIDGE", "CONTACTS_HYDROPHOBIC", "CONTACTS_DISULFIDE", "CONTACTS_VDW")
    ccol <- c(CONTACTS_HBOND = "#360568", CONTACTS_SALT_BRIDGE = "#5B2A86", CONTACTS_HYDROPHOBIC = "#7785AC",
              CONTACTS_DISULFIDE = "#9AC6C5", CONTACTS_VDW = "#A5E6BA")
    cl <- seqdf %>% pivot_longer(all_of(ct), names_to = "type", values_to = "value") %>%
      mutate(type = factor(type, levels = rev(ct)), value = pmax(value, 0))
    axis_y <- theme(axis.text.y = element_text(size = 8, colour = MUTED, hjust = 1, margin = margin(r = 3)),
                    panel.grid.major.y = element_line(colour = "#e6e6e6", linewidth = 0.3),
                    axis.title.y = element_text(size = 9, angle = 90, margin = margin(r = 4)))
    p_bar <- ggplot(cl, aes(x, value, fill = type)) + shade_layer +
      geom_col(width = 1) +
      scale_fill_manual(values = ccol, breaks = ct, name = "Contact type",
                        labels = c("H-bond", "Salt bridge", "Hydrophobic", "Disulfide", "Van der Waals"),
                        guide = guide_legend(nrow = 1, keywidth = unit(4, "mm"), keyheight = unit(3, "mm"))) +
      scale_x_continuous(limits = xlim, expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      labs(y = expression(Sigma[italic(CONTACTS)]~"(LCB per interface)")) +
      theme_sni() + axis_y +
      theme(legend.position = "top", legend.text = element_text(size = 8), legend.title = element_text(size = 9))
    p_tree <- tree_panel(TL, clades, bands, xlim)
    p_lab <- clade_label_panel(bands, xlim, MAIN_REFS)
    f2b <- p_bar / p_tree / p_lab + plot_layout(heights = c(2.1, 1.6, 0.75)) +
      plot_annotation(caption = "Inter-protomer contacts: heavy atoms < 4 Å, PAE < 10 Å; lower confidence bound across interfaces and seeds",
                      theme = theme(plot.caption = element_text(family = FAM, size = 7.5, colour = MUTED)))
    save_fig(f2b, file.path(OUT, "Fig2B_contacts_tree"), 13, 6.2)

    S4 <- c("ipTM_LCB", "Sum_CONTACTS", "BSA_INT_PROTO", "SD_THETA_ROT", "S_PROTO", "D_APEX",
            "THETA_APEX_mean", "L_APEX_mean", "H_ABS_mean", "MU_H_mean", "D_MHD_P_mean")
    panels <- map(S4, function(pn) {
      d <- seqdf %>% select(x, v = all_of(pn)) %>% filter(!is.na(v)) %>% arrange(x)
      ggplot(d, aes(x, v)) + shade_layer +
        geom_line(colour = "#b9b9b9", linewidth = 0.25) +
        geom_smooth(method = "loess", span = 0.06, se = FALSE, formula = y ~ x,
                    colour = CATEGORY_COLOURS[CATEGORY[pn]], linewidth = 0.8) +
        scale_x_continuous(limits = xlim, expand = c(0, 0)) +
        scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0.08, 0.12))) +
        labs(y = lab_expr(pn)) + theme_sni() + axis_y +
        theme(axis.title.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 6)),
              plot.margin = margin(3, 0, 3, 0))
    })
    s4 <- wrap_plots(c(panels, list(p_tree, p_lab)), ncol = 1,
                     heights = c(rep(1, length(S4)), 1.6, 0.75))
    save_fig(s4, file.path(OUT, "FigS4_SNI_along_tree"), 12, 17)
    write_csv(seqdf %>% arrange(x) %>% select(tree_x = x, ID, clade, all_of(S4), all_of(ct)),
              file.path(OUT, "FigS4_Fig2B_values_tree_order.csv"))
  }
}

# ================================================================== v1 vs v2
compare_versions <- function(v1p, v2p, label) {
  if (is.na(v1p) || is.na(v2p) || !file.exists(v1p) || !file.exists(v2p)) return(invisible())
  a <- read_summary(v1p); bb <- read_summary(v2p)
  cols <- intersect(c("Sum_CONTACTS", "BSA_INT_PROTO", "SD_THETA_ROT", "H_ABS_mean", "ipTM_LCB",
                      "S_PROTO", "D_APEX", "THETA_APEX_mean", "D_MHD_P_mean"), intersect(names(a), names(bb)))
  j <- inner_join(a %>% select(ID_norm, all_of(cols)), bb %>% select(ID_norm, all_of(cols)),
                  by = "ID_norm", suffix = c("_v1", "_v2"))
  long <- map_dfr(cols, function(cn) tibble(param = cn, ID = j$ID_norm,
                                            v1 = j[[paste0(cn, "_v1")]], v2 = j[[paste0(cn, "_v2")]]))
  summ <- long %>% group_by(param) %>%
    summarise(n = sum(!is.na(v1) & !is.na(v2)), n_changed = sum(abs(v1 - v2) > 1e-6, na.rm = TRUE),
              max_abs_change = max(abs(v1 - v2), na.rm = TRUE), .groups = "drop")
  write_csv(summ, file.path(OUT, paste0("v1_vs_v2_", label, "_summary.csv")))
  write_csv(long, file.path(OUT, paste0("v1_vs_v2_", label, "_per_NLR.csv")))
  long$param <- factor(long$param, levels = cols, labels = SYM[cols])
  p <- ggplot(long, aes(v1, v2)) + geom_abline(colour = "#bdbdbd", linetype = 2, linewidth = 0.3) +
    geom_point(size = 0.7, alpha = 0.55, colour = TXT, na.rm = TRUE) +
    facet_wrap(~ param, scales = "free", labeller = label_parsed) +
    labs(x = "v1 (published)", y = "v2") +
    theme_minimal(base_family = FAM, base_size = 9) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 10))
  save_fig(p, file.path(OUT, paste0("v1_vs_v2_", label)), 9, 7.5)
}
compare_versions(rp(opt$v1_main), main_path, "main")
compare_versions(rp(opt$v1_benchmark), bench_path, "benchmark")
cat("Done:", OUT, "\n")

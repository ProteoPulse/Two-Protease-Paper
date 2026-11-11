# report_stage.R — pandoc-free, self-contained HTML stage reports.
#
# Every stage/checkpoint emits a standalone HTML report so the user can VALIDATE each decision. At
# every decision point we render the decision plot(s) AND the explicit guiding question(s) inline
# (per the user's requirement). No pandoc/rmarkdown dependency: PNGs are base64-embedded and tables
# are rendered directly, so the HTML opens anywhere.
#
# Usage:
#   blocks <- list(
#     rb_text("Stage 1 qualified 4,201 precursors ..."),
#     rb_question("Which detection/coherence pair holds quality while retaining proteins?"),
#     rb_plot(ggobj_or_pngpath, "Metric distributions by protease"),
#     rb_table(dt, "Threshold sensitivity"))
#   render_stage_report("Stage 1 — Precursor qualification", blocks, out_html)
#
# Depends on: base64enc (for image embedding) if available; else links PNGs relatively.

suppressPackageStartupMessages({ library(data.table) })

rb_text     <- function(x)              list(type = "text", x = x)
rb_question <- function(...)            list(type = "question", x = c(...))
rb_table    <- function(dt, caption="") list(type = "table", dt = dt, caption = caption)
rb_plot     <- function(obj, caption="", width=8, height=5, dpi=110)
  list(type = "plot", obj = obj, caption = caption, width = width, height = height, dpi = dpi)

.html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&","&amp;",x,fixed=TRUE); x <- gsub("<","&lt;",x,fixed=TRUE)
  gsub(">","&gt;",x,fixed=TRUE)
}

.dt_to_html <- function(dt, caption = "", max_rows = 200) {
  dt <- as.data.table(dt)
  n <- nrow(dt); if (n > max_rows) dt <- dt[seq_len(max_rows)]
  fmt <- function(v) if (is.numeric(v)) formatC(v, digits = 4, format = "g") else as.character(v)
  hdr <- paste0("<tr>", paste0("<th>", .html_escape(names(dt)), "</th>", collapse=""), "</tr>")
  rows <- vapply(seq_len(nrow(dt)), function(i)
    paste0("<tr>", paste0("<td>", .html_escape(vapply(dt[i], fmt, character(1))), "</td>", collapse=""), "</tr>"),
    character(1))
  cap <- if (nzchar(caption)) sprintf("<div class='cap'>%s%s</div>", .html_escape(caption),
              if (n > max_rows) sprintf(" (showing %d of %d rows)", max_rows, n) else "") else ""
  paste0(cap, "<table>", hdr, paste(rows, collapse=""), "</table>")
}

.plot_to_img <- function(obj, out_dir, idx, width, height, dpi) {
  png_path <- file.path(out_dir, sprintf("fig_%02d.png", idx))
  if (inherits(obj, "ggplot")) {
    ggplot2::ggsave(png_path, obj, width = width, height = height, dpi = dpi)
  } else if (is.character(obj) && file.exists(obj)) {
    file.copy(obj, png_path, overwrite = TRUE)
  } else return("<div class='cap'>[plot unavailable]</div>")
  if (requireNamespace("base64enc", quietly = TRUE)) {
    uri <- base64enc::dataURI(file = png_path, mime = "image/png")
    sprintf("<img src='%s' />", uri)
  } else sprintf("<img src='%s' />", basename(png_path))
}

render_stage_report <- function(title, blocks, out_html, subtitle = "") {
  out_dir <- dirname(out_html); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  css <- "body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:1000px;margin:24px auto;padding:0 16px;color:#1a1a1a;line-height:1.5}
    h1{border-bottom:2px solid #444;padding-bottom:6px} h2{margin-top:28px;color:#333}
    .q{background:#fff8e1;border-left:4px solid #f5a623;padding:10px 14px;margin:14px 0;border-radius:4px}
    .q b{color:#b26a00} table{border-collapse:collapse;margin:10px 0;font-size:13px}
    th,td{border:1px solid #ddd;padding:4px 8px;text-align:right} th{background:#f4f4f4}
    td:first-child,th:first-child{text-align:left} .cap{color:#666;font-size:13px;margin:6px 0}
    img{max-width:100%;border:1px solid #eee;border-radius:4px;margin:8px 0}
    .sub{color:#777;margin-top:-8px}"
  parts <- c(sprintf("<h1>%s</h1>", .html_escape(title)),
             if (nzchar(subtitle)) sprintf("<p class='sub'>%s</p>", .html_escape(subtitle)) else "")
  fig_i <- 0L
  for (b in blocks) {
    if (b$type == "text")     parts <- c(parts, sprintf("<p>%s</p>", .html_escape(b$x)))
    else if (b$type == "question") parts <- c(parts,
      sprintf("<div class='q'><b>Decision question:</b><ul>%s</ul></div>",
              paste0("<li>", .html_escape(b$x), "</li>", collapse="")))
    else if (b$type == "table") parts <- c(parts, .dt_to_html(b$dt, b$caption))
    else if (b$type == "plot")  { fig_i <- fig_i + 1L
      parts <- c(parts, .plot_to_img(b$obj, out_dir, fig_i, b$width, b$height, b$dpi),
                 if (nzchar(b$caption)) sprintf("<div class='cap'>%s</div>", .html_escape(b$caption)) else "") }
  }
  html <- sprintf("<!doctype html><html><head><meta charset='utf-8'><title>%s</title><style>%s</style></head><body>%s</body></html>",
                  .html_escape(title), css, paste(parts, collapse="\n"))
  writeLines(html, out_html)
  out_html
}

# Build a simple INDEX.html linking a run's stage reports.
render_run_index <- function(run_dir, title, links) {
  parts <- c(sprintf("<h1>%s</h1>", .html_escape(title)), "<ul>")
  for (nm in names(links)) parts <- c(parts, sprintf("<li><a href='%s'>%s</a></li>",
                                                     .html_escape(links[[nm]]), .html_escape(nm)))
  parts <- c(parts, "</ul>")
  writeLines(sprintf("<!doctype html><html><head><meta charset='utf-8'><title>%s</title></head><body>%s</body></html>",
                     .html_escape(title), paste(parts, collapse="\n")),
             file.path(run_dir, "INDEX.html"))
}

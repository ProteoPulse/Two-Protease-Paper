# app.R — Protease Quant Advisor, GUI.
#
# No statistics live in this file. Every number comes from engine/pipeline.R::run_advisor(), which is
# the same function run_advisor.R calls, so the GUI and the CLI cannot drift.
#
# Launch:  Rscript -e 'shiny::runApp("protease-quant-advisor/app", launch.browser = TRUE)'
#
# Inputs:  two pr_matrix files (precursor level — the evidence) and two pg_matrix files (protein level
#          — the quantities being adjudicated). Optionally a FASTA, which unlocks sequence coverage.

suppressPackageStartupMessages({
  library(shiny); library(data.table); library(ggplot2); library(DT); library(yaml)
})

ROOT <- normalizePath(file.path(getwd(), ".."))   # runApp() sets wd to the app folder
source(file.path(ROOT, "engine", "load_engine.R")); load_advisor_engine(ROOT)
CFG0 <- yaml::read_yaml(file.path(ROOT, "config_default.yml"))
options(shiny.maxRequestSize = 1024 * 1024^2)     # pr_matrix files run to tens of MB

TRY_COL <- "#2C6FB3"; HTA_COL <- "#E08214"

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif}
    .well{background:#f7f9fb;border-color:#e3e9ef}
    h4{font-weight:600;margin-top:2px}
    .hint{color:#6b7684;font-size:12px;line-height:1.45;margin:-4px 0 10px}
    .banner{background:#eef4fa;border-left:4px solid #2C6FB3;padding:10px 14px;border-radius:4px;margin-bottom:12px}
    .warn{background:#fff8e1;border-left:4px solid #f5a623;padding:10px 14px;border-radius:4px;margin-bottom:12px}
    .sechd{font-size:11px;letter-spacing:.09em;text-transform:uppercase;color:#7c8794;
           font-weight:600;margin:16px 0 6px;border-bottom:1px solid #e3e9ef;padding-bottom:4px}
  "))),
  titlePanel("Protease Quant Advisor"),
  div(class = "hint", style = "margin-left:15px",
      "Two searches of one injection. For every protein both searches saw, decide which one's protein-group value to trust."),
  sidebarLayout(
    sidebarPanel(width = 3,
      div(class = "sechd", "Precursor matrices (evidence)"),
      fileInput("pr_a", "pr_matrix — channel A", accept = c(".tsv", ".txt")),
      fileInput("pr_b", "pr_matrix — channel B", accept = c(".tsv", ".txt")),
      div(class = "sechd", "Protein matrices (the quant)"),
      fileInput("pg_a", "pg_matrix — channel A", accept = c(".tsv", ".txt")),
      fileInput("pg_b", "pg_matrix — channel B", accept = c(".tsv", ".txt")),
      fluidRow(column(6, textInput("label_a", "Label A", CFG0$label_a)),
               column(6, textInput("label_b", "Label B", CFG0$label_b))),

      div(class = "sechd", "Which runs"),
      div(class = "hint", "Regex on the run name after stripping the path and extension. Leave blank to use every run."),
      fluidRow(column(6, textInput("run_include", "keep runs matching", "")),
               column(6, textInput("run_exclude", "drop runs matching", ""))),
      textInput("run_sample_regex", "run → sample regex (one capture group)", ""),
      div(class = "hint", "Only needed when several runs are replicates of one biological sample."),

      div(class = "sechd", "Stage 1 — qualification"),
      fluidRow(column(6, numericInput("det_min", "min detection", CFG0$qualify$detection_rate_min, 0, 1, .05)),
               column(6, numericInput("coh_min", "min coherence", CFG0$qualify$coherence_min, -1, 1, .05))),
      numericInput("cv_max", "max replicate CV (replicates only)", CFG0$qualify$median_cv_max, 0, 2, .05),

      div(class = "sechd", "Stage 3 — channel selection"),
      selectInput("sel_mode", "how to choose the channel",
                  c("completeness (manuscript method)" = "completeness",
                    "multi-criteria weighted vote" = "multi_criteria"),
                  CFG0$recommend$selection_mode %||% "completeness"),
      div(class = "hint", "Completeness = the search that saw the protein in more samples wins (the paper's method). Multi-criteria adds coherence, reliability, depth, and (with replicates) CV/ICC. Both report the other as vote_channel."),

      div(class = "sechd", "Stage 4 — agreement"),
      fluidRow(column(6, selectInput("cc_method", "statistic",
                        c("Kendall tau-b" = "kendall", "Spearman rho" = "spearman"),
                        CFG0$concordance$method)),
               column(6, numericInput("cc_pct", "null percentile", CFG0$concordance$null_percentile, .8, .999, .01))),
      fluidRow(column(6, numericInput("cc_minpaired", "min paired", CFG0$concordance$min_paired, 3, 100, 1)),
               column(6, numericInput("cc_nperm", "permutations", CFG0$concordance$n_perm, 100, 10000, 100))),

      div(class = "sechd", "QC probe"),
      checkboxInput("probe_on", "Compute the roll-up probe (diagnostic)", CFG0$probe$enabled),
      checkboxInput("probe_score", "Let it vote in the score", CFG0$probe$in_score),
      div(class = "hint", "Off by default: on our data the probe correlates 0.91 with peptide count and did not measurably improve the call."),

      div(class = "sechd", "Stage 5 — gap filling"),
      checkboxInput("imp_on", "Fill gaps from the other protease", CFG0$impute$enabled),
      selectInput("imp_model", "offset model",
                  c("per-protein median (manuscript)" = "protein_median", "sample-kNN" = "sample_knn"),
                  CFG0$impute$offset_model),
      div(class = "hint", "A protein whose channels disagree is never filled. Refusals are reported, not hidden."),

      div(class = "sechd", "Optional"),
      fileInput("fasta", "FASTA (unlocks sequence coverage)", accept = c(".fasta", ".fa", ".txt")),
      tags$hr(),
      actionButton("run", "Run advisor", class = "btn-primary btn-block"),
      tags$br(),
      downloadButton("dl_all", "Download all tables (.zip)", class = "btn-block"),
      downloadButton("dl_report", "Download HTML report", class = "btn-block")
    ),
    mainPanel(width = 9,
      uiOutput("status"),
      tabsetPanel(id = "tabs",
        tabPanel("Overview", br(),
          fluidRow(column(6, tableOutput("headline")), column(6, tableOutput("ingest"))),
          plotOutput("p_detect", height = "300px")),
        tabPanel("1 · Qualification", br(),
          plotOutput("p_sweep", height = "420px"), br(),
          plotOutput("p_qcmet", height = "300px"), br(),
          plotOutput("p_qctrace", height = "360px"), br(), tableOutput("qual_sum")),
        tabPanel("2 · Evidence", br(),
          plotOutput("p_evi", height = "400px"), br(),
          conditionalPanel("input.probe_on", plotOutput("p_probe", height = "340px")),
          br(), DT::dataTableOutput("t_evi")),
        tabPanel("4 · Agreement", br(),
          plotOutput("p_conc", height = "400px"), br(), tableOutput("conc_sum")),
        tabPanel("3 · Recommendation", br(),
          fluidRow(column(6, plotOutput("p_rec", height = "320px")),
                   column(6, plotOutput("p_conseq", height = "320px"))),
          br(), DT::dataTableOutput("t_rec")),
        tabPanel("Protein drill-down", br(),
          uiOutput("pg_pick"),
          fluidRow(column(7, plotOutput("p_drill", height = "380px")),
                   column(5, tableOutput("t_drill")))),
        tabPanel("5 · Gap filling", br(),
          uiOutput("imp_note"),
          plotOutput("p_imp", height = "360px"), br(),
          plotOutput("p_impdb", height = "420px"), br(),
          plotOutput("p_impprot", height = "460px"), br(),
          tableOutput("imp_sum"), br(), DT::dataTableOutput("t_imp")),
        tabPanel("Final matrix", br(),
          div(class = "hint",
              "One value per protein per sample, with the channel it came from and how confident that choice was. Filter value_source to drop imputed cells."),
          DT::dataTableOutput("t_final"))
      )
    )
  )
)

server <- function(input, output, session) {

  cfg_now <- reactive({
    cfg <- CFG0
    cfg$label_a <- input$label_a; cfg$label_b <- input$label_b
    cfg$ingest$run_include <- input$run_include
    cfg$ingest$run_exclude <- input$run_exclude
    cfg$ingest$run_sample_regex <- input$run_sample_regex
    cfg$qualify$detection_rate_min <- input$det_min
    cfg$qualify$coherence_min <- input$coh_min
    cfg$qualify$median_cv_max <- input$cv_max
    cfg$recommend$selection_mode <- input$sel_mode
    cfg$concordance$method <- input$cc_method
    cfg$concordance$null_percentile <- input$cc_pct
    cfg$concordance$min_paired <- input$cc_minpaired
    cfg$concordance$n_perm <- input$cc_nperm
    cfg$probe$enabled <- input$probe_on
    cfg$probe$in_score <- input$probe_score
    cfg$impute$enabled <- input$imp_on
    cfg$impute$offset_model <- input$imp_model
    cfg$fasta <- if (!is.null(input$fasta)) input$fasta$datapath else ""
    cfg
  })

  RES <- eventReactive(input$run, {
    req(input$pr_a, input$pr_b, input$pg_a, input$pg_b)
    withProgress(message = "Running the advisor", value = 0, {
      options(advisor.progress = function(msg, frac) setProgress(value = frac, detail = msg))
      on.exit(options(advisor.progress = NULL), add = TRUE)
      tryCatch(
        run_advisor(input$pr_a$datapath, input$pr_b$datapath,
                    input$pg_a$datapath, input$pg_b$datapath, cfg_now()),
        error = function(e) structure(list(error = conditionMessage(e)), class = "advisor_error"))
    })
  })

  ok <- reactive({ r <- RES(); !inherits(r, "advisor_error") })

  output$status <- renderUI({
    if (input$run == 0)
      return(div(class = "banner", strong("Load four matrices and press Run advisor."),
                 " Two pr_matrix files supply the evidence; two pg_matrix files supply the quantities being chosen between."))
    r <- RES()
    if (inherits(r, "advisor_error"))
      return(div(class = "warn", strong("Could not complete: "), r$error))
    div(class = "banner",
        sprintf("%d samples · %d shared proteins adjudicated (%s) · %d → %s, %d → %s · agreement threshold %.3f",
                length(r$samples), nrow(r$rec),
                r$cfg$recommend$selection_mode %||% "completeness",
                r$rec[recommended_channel == r$cfg$label_a, .N], r$cfg$label_a,
                r$rec[recommended_channel == r$cfg$label_b, .N], r$cfg$label_b, r$threshold),
        tags$div(class = "hint", style = "margin:6px 0 0", r$notes[1]))
  })

  sp <- function(f, ...) renderPlot({ req(ok()); p <- f(RES(), ...); req(!is.null(p)); p })
  output$p_detect <- sp(p_detection);      output$p_sweep  <- sp(p_qual_sweep)
  output$p_evi    <- sp(p_evidence);       output$p_probe  <- sp(p_probe)
  output$p_conc   <- sp(p_concordance);    output$p_rec    <- sp(p_recommendation)
  output$p_conseq <- sp(p_consequential);  output$p_imp    <- sp(p_imputation)
  output$p_qcmet  <- sp(p_qc_metrics);     output$p_qctrace <- sp(p_qc_traces)
  output$p_impdb  <- sp(p_impute_dumbbell); output$p_impprot <- sp(p_impute_protein_dumbbell)

  output$headline <- renderTable({ req(ok()); advisor_headline(RES()) })
  output$ingest   <- renderTable({ req(ok()); RES()$ingest_summary })
  output$qual_sum <- renderTable({ req(ok()); q <- RES()$qual
    data.table(precursors = nrow(q), qualifying = sum(q$qualifies),
               pct = round(100*mean(q$qualifies), 1)) })
  output$conc_sum <- renderTable({ req(ok()); r <- RES()
    c1 <- r$concordance[, .N, by = concordance_class_final]
    rbind(c1, data.table(concordance_class_final = "threshold",
                         N = round(r$threshold, 4)), fill = TRUE) })

  dtopts <- list(pageLength = 15, scrollX = TRUE, dom = "frtip")
  output$t_rec <- DT::renderDataTable({ req(ok())
    DT::datatable(advisor_tables(RES())$recommendation, options = dtopts,
                  filter = "top", rownames = FALSE) })
  output$t_evi <- DT::renderDataTable({ req(ok())
    DT::datatable(RES()$evidence, options = dtopts, filter = "top", rownames = FALSE) })
  output$t_final <- DT::renderDataTable({ req(ok())
    DT::datatable(RES()$final_matrix, options = dtopts, filter = "top", rownames = FALSE) })
  output$t_imp <- DT::renderDataTable({ req(ok()); im <- RES()$imputation; req(!is.null(im))
    d <- if (nrow(im$ledger)) im$ledger else im$refused
    DT::datatable(d, options = dtopts, filter = "top", rownames = FALSE) })
  output$imp_sum <- renderTable({ req(ok()); im <- RES()$imputation; req(!is.null(im)); im$summary })
  output$imp_note <- renderUI({ req(ok())
    if (is.null(RES()$imputation))
      div(class = "warn", "Gap filling is off. Switch it on in the sidebar and re-run.")
    else div(class = "hint",
      "Most gaps are usually refused, and that is the intended behaviour: a donor whose values do not track the recipient is not a valid surrogate.") })

  # ---- drill-down: one protein, both channels side by side ----
  output$pg_pick <- renderUI({ req(ok())
    selectizeInput("pg_sel", "Protein", choices = sort(RES()$rec$pg_canonical),
                   width = "440px", options = list(maxOptions = 2000)) })
  output$p_drill <- renderPlot({ req(ok(), input$pg_sel)
    r <- RES(); d <- r$pg[pg_canonical == input$pg_sel]
    req(nrow(d))
    rr <- r$rec[pg_canonical == input$pg_sel]
    ggplot(d, aes(sample_id, log2q, colour = protease, group = protease)) +
      geom_line(alpha = .8) + geom_point(size = 1.7) +
      scale_colour_manual(values = setNames(c(TRY_COL, HTA_COL), c(r$cfg$label_a, r$cfg$label_b)), name = NULL) +
      labs(title = input$pg_sel,
           subtitle = if (nrow(rr)) sprintf("recommended: %s (%s confidence, margin %.2f)",
                        rr$recommended_channel, rr$recommendation_confidence, rr$recommendation_margin) else "",
           x = NULL, y = "log2 intensity") +
      theme_adv() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  })
  output$t_drill <- renderTable({ req(ok(), input$pg_sel)
    e <- RES()$evidence[pg_canonical == input$pg_sel]
    req(nrow(e))
    cols <- intersect(c("protease","n_prec_qual","n_features","protein_detection_rate",
                        "median_sibling_coherence","peptide_set_jaccard","splithalf_reliability",
                        "lopo_rank_stability"), names(e))
    d <- e[, ..cols]
    for (cn in setdiff(cols, "protease")) if (is.numeric(d[[cn]])) d[[cn]] <- round(d[[cn]], 3)
    d })

  # ---- downloads ----
  output$dl_all <- downloadHandler(
    filename = function() sprintf("protease_advisor_tables_%s.zip", format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      req(ok())
      td <- file.path(tempdir(), paste0("adv_", as.integer(Sys.time())))
      dir.create(file.path(td, "tables"), recursive = TRUE, showWarnings = FALSE)
      tb <- advisor_tables(RES())
      for (nm in names(tb)) if (!is.null(tb[[nm]]) && nrow(tb[[nm]]))
        fwrite(tb[[nm]], file.path(td, "tables", paste0(nm, ".tsv")), sep = "\t")
      owd <- setwd(td); on.exit(setwd(owd), add = TRUE)
      utils::zip(file, "tables")
    })
  output$dl_report <- downloadHandler(
    filename = function() sprintf("protease_advisor_report_%s.html", format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      req(ok())
      td <- file.path(tempdir(), paste0("advrep_", as.integer(Sys.time())))
      dir.create(td, recursive = TRUE, showWarnings = FALSE)
      render_advisor_report(RES(), td, "app run")
      file.copy(file.path(td, "report.html"), file, overwrite = TRUE)
    })
}

shinyApp(ui, server)

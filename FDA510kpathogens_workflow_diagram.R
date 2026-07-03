################################################################################
# pathegenList3_workflow_diagram.R
#
# Reproduces the pathegenList3.R workflow diagram.
# Two rendering options:
#   Option A (default) — DiagrammeR::grViz()  → interactive HTML widget
#   Option B           — DiagrammeR + webshot2 → static PNG / PDF for reports
#
# Required packages:
#   install.packages(c("DiagrammeR", "htmlwidgets"))
#   # For PNG/PDF export (Option B):
#   install.packages(c("webshot2", "chromote"))
################################################################################

library(DiagrammeR)
library(htmlwidgets)

# ── Colour palette (mirrors the SVG diagram) ──────────────────────────────────
# Teal   : ingest / mapping      steps 1-2
# Purple : LLM gate / filter     steps 3-4
# Blue   : NCBI resolution       steps 5-7
# Coral  : LLM NER fallback      step 8
# Amber  : merge and enrich      step 9
# Gray   : final output node

COL <- list(
  teal_fill    = "#E1F5EE",  teal_stroke  = "#0F6E56",  teal_font   = "#085041",
  purple_fill  = "#EEEDFE",  purple_stroke= "#534AB7",  purple_font = "#3C3489",
  blue_fill    = "#E6F1FB",  blue_stroke  = "#185FA5",  blue_font   = "#0C447C",
  coral_fill   = "#FAECE7",  coral_stroke = "#993C1D",  coral_font  = "#712B13",
  amber_fill   = "#FAEEDA",  amber_stroke = "#854F0B",  amber_font  = "#633806",
  gray_fill    = "#F1EFE8",  gray_stroke  = "#5F5E5A",  gray_font   = "#444441"
)

# ── DOT source ────────────────────────────────────────────────────────────────
dot_src <- sprintf('
digraph pathegenList3 {

  node [
    shape     = rectangle,
    style     = "filled,rounded",
    fontname  = "Arial",
    fontsize  = 20,
    penwidth  = 0.8,
    width     = 4.2,
    height    = 0.75,
    margin    = "0.20,0.12"
  ]

  edge [
    fontname  = "Arial",
    fontsize  = 45,
    penwidth  = 1.0,
    color     = "#888780",
    arrowsize = 0.7
  ]

  # ── Step 1 ──────────────────────────────────────────────────────────────────
  S1 [
    label     = "Step 1 — Ingest FDA IVD microbial table\nnucleic-acid-based-tests#microbial → microbial_df",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s"
  ]

  # ── Step 2 ──────────────────────────────────────────────────────────────────
  S2 [
    label     = "Step 2 — Build organism_mapping\nUnique label × Submission ID pairs",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s"
  ]

  # ── Step 3 ──────────────────────────────────────────────────────────────────
  S3 [
    label     = "Step 3 — Technology classification (LLM gate)\nFetch PDF summary for every unique Submission ID\nRun LLM classifier → submission_tech_df\nKeep PCR / Other_Amplification / Hybrid → naat_submissions",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s",
    height    = 1.10
  ]

  # Excluded branch (note node)
  EXCL [
    label     = "Hybridization / Unclear\n→ excluded",
    shape     = note,
    style     = "filled",
    fillcolor = "#FFF8F0",
    color     = "%s",
    fontcolor = "%s",
    fontsize  = 15,
    width     = 1.6,
    height    = 0.55
  ]

  # ── Step 4 ──────────────────────────────────────────────────────────────────
  S4 [
    label     = "Step 4 — Filter to NAAT submissions\nmicrobial_df and organism_mapping reduced to naat_submissions",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s"
  ]

  # ── Step 5 ──────────────────────────────────────────────────────────────────
  S5 [
    label     = "Step 5 — Label partitioning\nSubset 1: binomial / uninomial  ·  Subset 2: noisy / panel / multi-agent",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s"
  ]

  # ── Step 6 ──────────────────────────────────────────────────────────────────
  S6 [
    label     = "Step 6 — Resolve Subset 1\nStop-word filter → spelling patches\nUnique NCBI TaxID → pathogen_set1",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s",
    width     = 3.2,
    height    = 1.0
  ]

  # ── Step 7 ──────────────────────────────────────────────────────────────────
  S7 [
    label     = "Step 7 — Expand Subset 2\nSplit CT/GC  ·  acronym expansion\nUnique NCBI TaxID → pathogen_set2",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s",
    width     = 2.8,
    height    = 0.85
  ]

  # Invisible merge node (combines the two branches cleanly)
  MERGE [
    label     = "",
    shape     = point,
    width     = 0.01,
    height    = 0.01,
    style     = invis
  ]

  # ── Step 8 ──────────────────────────────────────────────────────────────────
  S8 [
    label     = "Step 8 — LLM NER on unmapped submissions\nFetch FDA summary PDF for still-unmapped Submission IDs\nNER extracts organism names → pathogen_set3",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s",
    height    = 0.90
  ]

  # ── Step 9 ──────────────────────────────────────────────────────────────────
  S9 [
    label     = "Step 9 — Merge and enrich\nCombine pathogen_set1/2/3  ·  dedupe by TaxID\nAttach rank / division / genus / species from NCBI",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s",
    height    = 0.90
  ]

  # ── Final output ─────────────────────────────────────────────────────────────
  OUT [
    label     = "FDA510kPathogens\nNAAT-only  ·  TaxID-resolved  ·  taxonomy-enriched",
    fillcolor = "%s",
    color     = "%s",
    fontcolor = "%s",
    width     = 3.4
  ]

  # Force EXCL into the same rank as S3 (same horizontal band) so the
  # edge runs rightward rather than downward.
  { rank = same; S3; EXCL }

  # ── Edges ────────────────────────────────────────────────────────────────────
  S1   -> S2   [color = "%s"]
  S2   -> S3   [color = "%s"]
  # constraint=false prevents EXCL from being pulled into the main vertical
  # spine; minlen=1 keeps the arrow short and pointing right.
  S3   -> EXCL [label = " excluded", style = dashed, arrowhead = open,
                color = "%s", fontcolor = "%s", fontsize = 13,
                constraint = false, minlen = 1]
  S3   -> S4   [color = "%s"]
  S4   -> S5   [color = "%s"]
  S5   -> S6   [label = " Subset 1", fontcolor = "%s",
                color = "%s", fontsize = 13]
  S5   -> S7   [label = "Subset 2 ", fontcolor = "%s",
                color = "%s", fontsize = 13]
  S6   -> MERGE [arrowhead = none, color = "%s"]
  S7   -> MERGE [arrowhead = none, color = "%s"]
  MERGE -> S8  [color = "%s"]
  S8   -> S9   [color = "%s"]
  S9   -> OUT  [color = "%s"]

  # Force steps 6 and 7 into the same rank (side-by-side)
  { rank = same; S6; S7 }
}
',
# Step 1 node
COL$teal_fill,   COL$teal_stroke,   COL$teal_font,
# Step 2 node
COL$teal_fill,   COL$teal_stroke,   COL$teal_font,
# Step 3 node
COL$purple_fill, COL$purple_stroke, COL$purple_font,
# EXCL note
COL$coral_stroke, COL$coral_font,
# Step 4 node
COL$purple_fill, COL$purple_stroke, COL$purple_font,
# Step 5 node
COL$blue_fill,   COL$blue_stroke,   COL$blue_font,
# Step 6 node
COL$blue_fill,   COL$blue_stroke,   COL$blue_font,
# Step 7 node
COL$blue_fill,   COL$blue_stroke,   COL$blue_font,
# Step 8 node
COL$coral_fill,  COL$coral_stroke,  COL$coral_font,
# Step 9 node
COL$amber_fill,  COL$amber_stroke,  COL$amber_font,
# OUT node
COL$gray_fill,   COL$gray_stroke,   COL$gray_font,
# Edge colours
COL$teal_stroke,                          # S1->S2
COL$teal_stroke,                          # S2->S3
COL$coral_stroke, COL$coral_stroke,       # S3->EXCL (dashed)
COL$purple_stroke,                        # S3->S4
COL$purple_stroke,                        # S4->S5
COL$blue_font,  COL$blue_stroke,          # S5->S6 label+line
COL$blue_font,  COL$blue_stroke,          # S5->S7 label+line
COL$blue_stroke,                          # S6->MERGE
COL$blue_stroke,                          # S7->MERGE
COL$blue_stroke,                          # MERGE->S8
COL$coral_stroke,                         # S8->S9
COL$amber_stroke                          # S9->OUT
)

# ── Render ────────────────────────────────────────────────────────────────────

# Option A: interactive HTML widget (opens in RStudio Viewer / browser)
diagram <- DiagrammeR::grViz(dot_src)
print(diagram)

# ── Option B: save to file ─────────────────────────────────────────────────────
# HTML (self-contained, shareable)
htmlwidgets::saveWidget(
  widget   = diagram,
  file     = "pathegenList3_workflow.html",
  selfcontained = TRUE,
  title    = "pathegenList3 Workflow"
)
message("Saved: pathegenList3_workflow.html")

# PNG / PDF  (requires webshot2 + a local Chrome/Chromium install)
# Uncomment to activate:
#
# if (requireNamespace("webshot2", quietly = TRUE)) {
#   webshot2::webshot(
#     url    = "pathegenList3_workflow.html",
#     file   = "pathegenList3_workflow.png",
#     vwidth = 900,
#     vheight = 1200,
#     zoom   = 2          # 2× for retina-quality output
#   )
#   message("Saved: pathegenList3_workflow.png")
#
#   webshot2::webshot(
#     url    = "pathegenList3_workflow.html",
#     file   = "pathegenList3_workflow.pdf",
#     vwidth = 900,
#     vheight = 1200
#   )
#   message("Saved: pathegenList3_workflow.pdf")
# } else {
#   message("Install webshot2 for PNG/PDF export: install.packages('webshot2')")
# }

# ==============================================================================
# SCRIPT:  FDA510kpathogens.R
# PURPOSE: Build a comprehensive resource mapping FDA-cleared/approved
#          Nucleic Acid Amplification Tests (NAATs) to NCBI taxonomy IDs.
#
# PIPELINE OVERVIEW
#   1. Scrape the FDA "Nucleic Acid-Based Tests" microbial table (510(k) records).
#   2. Collapse to one row per unique organism label, with its 510(k) submissions.
#   3. Split labels into "binomial/uninomial" candidates (Subset 1) vs. noisy /
#      multi-agent / panel labels (Subset 2).
#   4. Resolve Subset 1 labels to a unique NCBI Taxonomy ID -> pathogen_set1.
#   5. Expand Subset 2 multi-agent labels (e.g. "CT/GC"), apply acronym
#      expansion, and resolve to NCBI Taxonomy IDs -> pathogen_set2.
#   6. For 510(k) submissions still unmapped, fetch the FDA summary PDF and use
#      an LLM (NER) to extract organism names -> pathogen_set3.
#   7. Merge all three sets, dedupe by TaxID, and attach taxonomic metadata
#      (rank, division, genus, species, ...) from NCBI -> FDA510kPathogens.
#
# CHECKPOINTING
#   Every step that hits an external service (FDA website, FDA PDF server,
#   NCBI Entrez, Gemini LLM) is cached to an .rds file under CACHE_DIR. Re-running
#   the script reuses cached results instead of re-querying, unless
#   FORCE_REFRESH is set TRUE for that step (or the cache file is deleted).
# ==============================================================================

suppressPackageStartupMessages({
  library(rvest)      # HTML table scraping
  library(httr)        # HTTP requests with custom headers
  library(tidyverse)    # dplyr, tidyr, stringr, purrr, etc.
  library(DT)           # interactive table preview (optional, end of script)
  library(rentrez)      # NCBI Entrez E-utilities (taxonomy lookups)
  library(pdftools)     # pdf_text() for FDA 510(k) summary PDFs
  library(jsonlite)     # fromJSON() for parsing LLM NER output
  library(ellmer)       # chat_google_gemini() LLM client
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

CACHE_DIR <- "cache"
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR)

# Set any of these TRUE to force that step to re-hit the external service
# instead of reusing its cache file.
FORCE_REFRESH <- list(
  fda_table        = FALSE,
  taxdb_subset1    = FALSE,
  taxdb_subset2    = FALSE,
  taxdb_subset3    = FALSE,
  fda_pdfs         = FALSE,
  llm_ner          = FALSE,
  taxonomy_summary = FALSE
)

# Minimum delay (seconds) between consecutive NCBI Entrez calls. NCBI's usage
# policy expects no more than ~3 requests/second without an API key; we use a
# conservative default to avoid 429s when looping over hundreds of organisms.
ENTREZ_DELAY <- 0.34

# Acronym / clinical shorthand -> binomial name. Used when expanding
# multi-agent panel labels (e.g. "CT/GC" -> "Chlamydia trachomatis",
# "Neisseria gonorrhoeae").
ACRONYM_MAP <- c(
  "CT"      = "Chlamydia trachomatis",
  "GC"      = "Neisseria gonorrhoeae",
  "NG"      = "Neisseria gonorrhoeae",
  "GBS"     = "Streptococcus agalactiae",
  "GAS"     = "Streptococcus pyogenes",
  "Strep A" = "Streptococcus pyogenes",
  "Strep B" = "Streptococcus agalactiae",
  "TV"      = "Trichomonas vaginalis",
  "HSV1"    = "Human alphaherpesvirus 1",
  "HSV2"    = "Human alphaherpesvirus 2",
  "VZV"     = "Human alphaherpesvirus 3",
  "RSV"     = "Human respiratory syncytial virus",
  "MTB"     = "Mycobacterium tuberculosis"
)

# Manual corrections for FDA regulatory labels that don't resolve cleanly
# against the NCBI Taxonomy database (synonyms, misspellings, deprecated
# names). Keyed by the *original* FDA label text.
NCBI_PATCHES_SUBSET1 <- c(
  "Clostridium difficile" = "Clostridioides difficile",
  "Entercoccus"           = "Enterococcus",
  "Human Papillomvirus"   = "Human papillomavirus",
  "Leishmaniaspecies"     = "Leishmania",
  "Mycobacterium species" = "Mycobacterium",
  "Mycoplasma pneumoniae" = "Mycoplasmoides pneumoniae",
  "Streptococci"          = "Streptococcus",
  "Variola"               = "Variola virus"
)

NCBI_PATCHES_SUBSET2 <- c(
  "Candida parapsilosiis" = "Candida parapsilosis",
  "Candida kruseii"       = "Pichia kudriavzevii"
)

NCBI_PATCHES_SUBSET3 <- c(
  "Streptococcus vini" = "Streptococcus mutans",
  "Candida rugosa"     = "Diutina rugosa",
  "Variola"            = "Variola virus",
  "Enterobacter dispar" = "Enterobacter cloacae"
)

# Regex flagging panel/metadata noise rather than a real organism name
# (e.g. "Carbapenem Resistance", "Vaginitis Panel").
STOP_WORDS_PATTERN <- "(?i)resistance|panel|media|material|toxin|illness|specimen|test|kit|system"

# Regex matching a clean binomial ("Genus species") or uninomial
# ("Adenovirus") organism label: 1-2 alphanumeric words, nothing else.
BINOMIAL_OR_UNINOMIAL_PATTERN <- "^[A-Za-z0-9]+( [A-Za-z0-9]+)?$"

FDA_IVD_MICROBIAL_URL <- "https://www.fda.gov/medical-devices/in-vitro-diagnostics/nucleic-acid-based-tests#microbial"
USER_AGENT <- "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36"

# ------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------

# Loads an .rds cache if present (and refresh not forced); otherwise runs
# `compute_fn`, caches the result, and returns it. Centralizes the
# checkpointing pattern used throughout the pipeline.
with_cache <- function(cache_name, force_refresh, compute_fn) {
  cache_path <- file.path(CACHE_DIR, paste0(cache_name, ".rds"))
  if (!force_refresh && file.exists(cache_path)) {
    message("Using cached result for '", cache_name, "' (", cache_path, ")")
    return(readRDS(cache_path))
  }
  message("Computing '", cache_name, "' ...")
  result <- compute_fn()
  saveRDS(result, cache_path)
  result
}

# One NCBI Taxonomy search with a polite delay before the request, used
# everywhere we loop over many organism names. Returns the full search
# result list (so callers can read $count or $ids as needed).
#
# Hardened against two failure modes seen with LLM-extracted / patched
# labels that can be blank, whitespace-only, or otherwise malformed:
#   1. Blank/empty terms are short-circuited to a zero-result list without
#      ever hitting NCBI (an empty query term is not a meaningful search).
#   2. entrez_search() itself can throw "subscript out of bounds" on some
#      malformed/unusual query terms (a known rentrez XML-parsing failure
#      mode, not a network/rate-limit issue) - wrapped in tryCatch so a
#      single bad label doesn't abort the whole batch resolution.
entrez_search_taxonomy <- function(term) {
  if (is.na(term) || !nzchar(trimws(term))) {
    return(list(count = 0, ids = character(0)))
  }
  
  Sys.sleep(ENTREZ_DELAY)
  tryCatch({
    entrez_search(db = "taxonomy", term = term)
  }, error = function(e) {
    message("entrez_search failed for term '", term, "': ", e$message, " - treating as unresolved.")
    list(count = 0, ids = character(0))
  })
}

# Applies a named patch vector (old_label -> new_label) to a character
# vector, leaving unmatched entries unchanged.
apply_patches <- function(labels, patches) {
  ifelse(labels %in% names(patches), patches[labels], labels)
}

# Resolves a vector of organism labels to NCBI Taxonomy IDs, keeping only
# labels with EXACTLY one matching tax id (ambiguous or unmatched labels are
# returned separately so callers can route them elsewhere, e.g. to an LLM
# fallback). label_col names the output column holding the original label.
resolve_taxids_unique_only <- function(labels, label_col) {
  labels <- unique(labels)
  counts <- vapply(labels, function(x) as.integer(entrez_search_taxonomy(x)$count), integer(1))
  
  unique_labels <- names(which(counts == 1))
  ambiguous_or_unmatched <- names(which(counts != 1))
  
  resolved <- do.call(rbind, lapply(unique_labels, function(x) {
    res <- entrez_search_taxonomy(x)
    data.frame(
      label = x,
      TaxID = if (length(res$ids) > 0) unlist(res$ids)[1] else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(resolved)) resolved <- data.frame(label = character(), TaxID = character())
  names(resolved)[names(resolved) == "label"] <- label_col
  
  list(resolved = resolved, unresolved_labels = ambiguous_or_unmatched)
}

# Aggregates a data frame of label/submission-id rows to one row per unique
# label, with all of that label's submission IDs joined into one string.
# value_col names the input column holding submission IDs (it may be called
# "Submission" or "Submissions" depending on which earlier step produced the
# data frame); the output column is always named "Submissions" (plural).
collapse_submissions_by_label <- function(df, label_col, value_col = "Submissions", sep = ",") {
  df %>%
    group_by(.data[[label_col]]) %>%
    dplyr::summarize(Submissions = paste(unique(.data[[value_col]]), collapse = sep), .groups = "drop")
}

# ==============================================================================
# STEP 1: Scrape the FDA microbial NAAT table (510(k) records)
# ==============================================================================

MicrobialTests <- with_cache("01_fda_microbial_table", FORCE_REFRESH$fda_table, function() {
  response <- GET(FDA_IVD_MICROBIAL_URL, add_headers(`User-Agent` = USER_AGENT))
  if (status_code(response) != 200) {
    stop("Failed to connect to FDA site. Status code: ", status_code(response))
  }
  
  page_content <- content(response, as = "text")
  tables <- read_html(page_content) %>% html_table(fill = TRUE)
  
  # The "Microbial" table is the second table on the page.
  microbial_tests <- tables[[2]]
  colnames(microbial_tests) <- make.names(colnames(microbial_tests))
  
  microbial_tests %>%
    mutate(Submission = str_extract(Submission, "^[A-Za-z]+\\d+"))
})

# ==============================================================================
# STEP 2: Collapse to one row per unique organism label
# ==============================================================================

organism_mapping <- MicrobialTests %>%
  # Squishing whitespace first so e.g. "Adenovirus " and "Adenovirus" match.
  mutate(FDAregulatoryLabels = str_squish(Microorganism)) %>%
  group_by(FDAregulatoryLabels) %>%
  dplyr::summarize(Submissions = paste(unique(Submission), collapse = ","), .groups = "drop")

# ==============================================================================
# STEP 3-4: Partition into binomial/uninomial candidates (Subset 1) vs.
#           noisy / panel / multi-agent labels (Subset 2)
# ==============================================================================

is_binomial_or_uninomial <- str_detect(organism_mapping$FDAregulatoryLabels, BINOMIAL_OR_UNINOMIAL_PATTERN)
subset_1 <- organism_mapping[is_binomial_or_uninomial, ]
subset_2 <- organism_mapping[!is_binomial_or_uninomial, ]

# Move stop-word matches (panel/metadata noise) from Subset 1 into Subset 2.
is_stop_word_subset1 <- str_detect(subset_1$FDAregulatoryLabels, STOP_WORDS_PATTERN)
subset_1_1 <- subset_1[!is_stop_word_subset1, ]
subset_2 <- rbind(subset_2, subset_1[is_stop_word_subset1, ])

# Apply manual NCBI corrections before querying the taxonomy database.
subset_1_1$FDAregulatoryLabels <- apply_patches(subset_1_1$FDAregulatoryLabels, NCBI_PATCHES_SUBSET1)

# ==============================================================================
# STEP 5: Resolve Subset 1 labels to NCBI Taxonomy IDs -> pathogen_set1
# ==============================================================================

subset1_taxonomy <- with_cache("02_subset1_taxonomy", FORCE_REFRESH$taxdb_subset1, function() {
  resolve_taxids_unique_only(subset_1_1$FDAregulatoryLabels, label_col = "FDAregulatoryLabels")
})

# Labels with zero or multiple taxonomy matches aren't reliable single-organism
# entries; route them to Subset 2 instead of dropping them.
subset_2 <- rbind(
  subset_2,
  subset_1_1[subset_1_1$FDAregulatoryLabels %in% subset1_taxonomy$unresolved_labels, ]
)
subset_1_1 <- subset_1_1[!subset_1_1$FDAregulatoryLabels %in% subset1_taxonomy$unresolved_labels, ]

pathogen_set1 <- merge(subset_1_1, subset1_taxonomy$resolved, by = "FDAregulatoryLabels")

# ==============================================================================
# STEP 6: Expand Subset 2 multi-agent labels and resolve to TaxIDs -> pathogen_set2
# ==============================================================================

# Rows that bundle multiple organisms in one label, e.g. "CT/GC", "GAS and GBS".
subset_2_multi_agent <- subset_2 %>%
  dplyr::filter(str_detect(FDAregulatoryLabels, "/|,|\\band\\b"))

subset_2_expanded <- subset_2_multi_agent %>%
  separate_rows(FDAregulatoryLabels, sep = "/|,|\\band\\b") %>%
  mutate(FDAregulatoryLabels = str_squish(FDAregulatoryLabels))

# Drop any newly-exposed panel/metadata noise from the expanded rows.
subset_2_expanded <- subset_2_expanded[!str_detect(subset_2_expanded$FDAregulatoryLabels, STOP_WORDS_PATTERN), ]
subset_2_expanded <- collapse_submissions_by_label(subset_2_expanded, "FDAregulatoryLabels")

subset_2_expanded$FDAregulatoryLabels <- apply_patches(subset_2_expanded$FDAregulatoryLabels, NCBI_PATCHES_SUBSET2)

# Expand acronyms (CT, GC, ...) and strip clinical/diagnostic qualifiers
# (DNA, RNA, Antigen, ...) so labels are ready for taxonomy lookup.
clean_for_taxonomy_lookup <- function(df, acronym_map) {
  df %>%
    mutate(
      Microorganism_Expanded = apply_patches(FDAregulatoryLabels, acronym_map),
      Microorganism_NCBItaxDb_ready = Microorganism_Expanded %>%
        str_remove_all("(?i)\\b(DNA|RNA|rRNA|Culture|Antigen|Toxin|Detection|Nuclear -1)\\b") %>%
        str_remove_all("\\(.*?\\)") %>%
        str_squish()
    )
}

subset_2_cleaned <- clean_for_taxonomy_lookup(subset_2_expanded, ACRONYM_MAP)

# Skip organisms already resolved in Subset 1.
subset_2_cleaned <- subset_2_cleaned[!subset_2_cleaned$Microorganism_NCBItaxDb_ready %in% pathogen_set1$FDAregulatoryLabels, ]
subset_2_cleaned <- collapse_submissions_by_label(subset_2_cleaned, "Microorganism_NCBItaxDb_ready")

subset2_taxonomy <- with_cache("03_subset2_taxonomy", FORCE_REFRESH$taxdb_subset2, function() {
  resolve_taxids_unique_only(subset_2_cleaned$Microorganism_NCBItaxDb_ready, label_col = "Microorganism_NCBItaxDb_ready")
})

subset_2_cleaned <- subset_2_cleaned[!subset_2_cleaned$Microorganism_NCBItaxDb_ready %in% subset2_taxonomy$unresolved_labels, ]

pathogen_set2 <- merge(subset_2_cleaned, subset2_taxonomy$resolved, by = "Microorganism_NCBItaxDb_ready") %>%
  dplyr::rename(FDAregulatoryLabels = Microorganism_NCBItaxDb_ready)

# ==============================================================================
# Combine Subset 1 + Subset 2 results
# ==============================================================================

pathogen_set1_2 <- rbind(
  pathogen_set1[, c("FDAregulatoryLabels", "Submissions", "TaxID")],
  pathogen_set2[, c("FDAregulatoryLabels", "Submissions", "TaxID")]
)

# "Sepsis" is a clinical syndrome, not an organism; drop it if present.
pathogen_set1_2 <- pathogen_set1_2[pathogen_set1_2$FDAregulatoryLabels != "Sepsis", ]

# ==============================================================================
# STEP 7: For 510(k) submissions still unmapped to any organism, fetch the FDA
#         summary PDF and run LLM-based NER to extract organism names -> pathogen_set3
# ==============================================================================

mapped_submissions <- pathogen_set1_2$Submissions %>%
  str_split(",") %>%
  unlist() %>%
  str_squish() %>%
  unique()

all_submissions <- unique(MicrobialTests$Submission)
unmapped_submissions <- setdiff(all_submissions, mapped_submissions)

# FDA URL convention for 510(k) summary PDFs: a "K" number like K220870 lives
# at .../pdf22/K220870.pdf (the year prefix is embedded in the submission ID).
# De Novo (DEN) and PMA (P) numbers use different URL schemes, not handled here.
get_fda_summary_text <- function(submission_id) {
  if (!str_starts(submission_id, "K")) return(NA_character_)
  
  yr_prefix <- substr(submission_id, 2, 3)
  url <- paste0("https://www.accessdata.fda.gov/cdrh_docs/pdf", yr_prefix, "/", submission_id, ".pdf")
  
  tryCatch({
    raw_text <- pdf_text(url)
    paste(raw_text, collapse = " ") %>% str_squish()
  }, error = function(e) {
    message("Could not fetch PDF for: ", submission_id)
    NA_character_
  })
}

fda510kSubmissions <- with_cache("04_fda_510k_pdf_text", FORCE_REFRESH$fda_pdfs, function() {
  texts <- lapply(unmapped_submissions, get_fda_summary_text)
  names(texts) <- unmapped_submissions  # BUGFIX: previously named by content, not by submission ID
  texts
})

fda510kSubmissions_with_text <- fda510kSubmissions[!is.na(fda510kSubmissions)]

# --- LLM-based Named Entity Recognition over the PDF text ---

ner_chat <- chat_google_gemini(
  model = "gemini-2.5-flash",
  system_prompt = paste0(
    "Perform Named Entity Recognition (NER). ",
    "Extract all MICROORGANISM entities from the text. ",
    "Return them as a strictly JSON formatted array of strings. ",
    "Standardize abbreviations to binomial nomenclature."
  ),
  api_args = list(generationConfig = list(temperature = 0))
)

# Extracts organism names from one document's text via the LLM, retrying on
# rate-limit (429) errors and failing fast (after logging) on any other
# error. BUGFIX: the original version's `return()` inside the tryCatch error
# handler only exited that handler's own closure, not get_organisms_safely,
# so a non-429 error would silently retry up to max_retries times instead of
# returning immediately. This version uses an explicit `give_up` flag (set
# via `<<-` since the error handler IS its own closure) to break out of the
# retry loop on non-retryable errors.
get_organisms_safely <- function(doc_id, text_data, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    parsed <- NULL
    give_up <- FALSE
    
    tryCatch({
      raw_response <- ner_chat$chat(text_data)
      clean_json <- str_remove_all(raw_response, "```json\n|```")
      parsed <- fromJSON(clean_json)
    }, error = function(e) {
      if (grepl("429", e$message)) {
        message("Token limit hit on ", doc_id, ". Sleeping 45s to reset quota...")
        Sys.sleep(45)
      } else {
        message("Non-rate-limit error on ", doc_id, ": ", e$message)
        give_up <<- TRUE
      }
    })
    
    if (!is.null(parsed)) {
      # parsed may be character(0) if the LLM found no organism entities in
      # this document (a valid "no match" result, not an error) - in that
      # case return a single NA row rather than a 0-row data.frame, since
      # data.frame(Submission = doc_id, FDAregulatoryLabels = character(0))
      # errors on the row-count mismatch (1 vs 0).
      labels_out <- if (length(parsed) == 0) NA_character_ else parsed
      return(data.frame(Submission = doc_id, FDAregulatoryLabels = labels_out, stringsAsFactors = FALSE))
    }
    if (give_up) break
  }
  # Failed after retries (or a non-retryable error): return NA row so the
  # caller can identify and skip/inspect this document.
  data.frame(Submission = doc_id, FDAregulatoryLabels = NA_character_, stringsAsFactors = FALSE)
}

subset_3 <- with_cache("05_llm_ner_results", FORCE_REFRESH$llm_ner, function() {
  do.call(rbind, lapply(names(fda510kSubmissions_with_text), function(doc_id) {
    get_organisms_safely(doc_id, fda510kSubmissions_with_text[[doc_id]])
  }))
})

subset_3 <- subset_3[!is.na(subset_3$FDAregulatoryLabels), ]
subset_3_mapping <- collapse_submissions_by_label(subset_3, "FDAregulatoryLabels", value_col = "Submission")

# Merge submissions found via NER into pathogen_set1_2 wherever the label
# already matches an organism resolved in Steps 5-6.
pathogen_set1_2 <- pathogen_set1_2 %>%
  left_join(subset_3_mapping, by = "FDAregulatoryLabels", suffix = c("_existing", "_ner")) %>%
  mutate(
    Submissions = case_when(
      is.na(Submissions_ner) ~ Submissions_existing,
      is.na(Submissions_existing) ~ Submissions_ner,
      TRUE ~ paste(Submissions_existing, Submissions_ner, sep = ",")
    )
  ) %>%
  dplyr::select(FDAregulatoryLabels, Submissions, TaxID)

# NER-discovered labels not already covered above need their own taxonomy
# resolution -> pathogen_set3.
subset_3_unmapped <- subset_3_mapping[!subset_3_mapping$FDAregulatoryLabels %in% pathogen_set1_2$FDAregulatoryLabels, ]
subset_3_unmapped$FDAregulatoryLabels <- apply_patches(subset_3_unmapped$FDAregulatoryLabels, NCBI_PATCHES_SUBSET3)
subset_3_unmapped <- collapse_submissions_by_label(subset_3_unmapped, "FDAregulatoryLabels")

subset3_taxonomy <- with_cache("06_subset3_taxonomy", FORCE_REFRESH$taxdb_subset3, function() {
  resolve_taxids_unique_only(subset_3_unmapped$FDAregulatoryLabels, label_col = "FDAregulatoryLabels")
})

subset_3_resolved <- subset_3_unmapped[subset_3_unmapped$FDAregulatoryLabels %in% subset3_taxonomy$resolved$FDAregulatoryLabels, ]
pathogen_set3 <- merge(subset_3_resolved, subset3_taxonomy$resolved, by = "FDAregulatoryLabels")

# ==============================================================================
# Final combination, dedup, and taxonomic annotation
# ==============================================================================

FDA510kPathogens <- rbind(pathogen_set1_2, pathogen_set3)

# Different regulatory documents sometimes use inconsistent capitalization or
# slight typos for the same organism, producing duplicate TaxIDs. Collapse to
# one row per TaxID, keeping the first label seen and merging all submissions.
FDA510kPathogens <- FDA510kPathogens %>%
  group_by(TaxID) %>%
  dplyr::summarise(
    FDAregulatoryLabels = first(FDAregulatoryLabels),
    Submissions = paste(unique(Submissions), collapse = ", "),
    .groups = "drop"
  )

# Fetch taxonomic classification (rank, division, genus, species, ...) for
# each resolved TaxID.
fetch_taxonomy_summary <- function(taxid) {
  Sys.sleep(ENTREZ_DELAY)
  tax_rec <- entrez_summary(db = "taxonomy", taxid)
  fields <- extract_from_esummary(
    tax_rec,
    c("rank", "division", "scientificname", "commonname", "taxid", "genus", "species", "subsp"),
    simplify = TRUE
  )
  as.data.frame(t(fields), stringsAsFactors = FALSE)
}

FDA510kPathogensDivisions <- with_cache("07_taxonomy_summaries", FORCE_REFRESH$taxonomy_summary, function() {
  do.call(rbind, lapply(FDA510kPathogens$TaxID, fetch_taxonomy_summary))
})

FDA510kPathogens <- merge(FDA510kPathogens, FDA510kPathogensDivisions, by.x = "TaxID", by.y = "taxid")

# ==============================================================================
# Done. FDA510kPathogens now holds one row per unique organism (by TaxID),
# its FDA regulatory label, all associated 510(k) submission IDs, and its
# NCBI taxonomic classification.
# ==============================================================================

# Uncomment to preview interactively:
# DT::datatable(FDA510kPathogens)

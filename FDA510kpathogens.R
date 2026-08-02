################################################################################
# SCRIPT:  FDA510kpathogens.R
# PURPOSE: Build a comprehensive resource mapping FDA-cleared/approved
#          Nucleic Acid Amplification Tests (NAATs) to NCBI taxonomy IDs.
#
# Workflow:
#   Step 1  — Ingest FDA IVD microbial table → microbial_df
#   Step 2  — Build organism_mapping (unique label × Submission ID)
#   Step 3  — Fetch PDF summary for every unique Submission ID
#             Run LLM technology classifier → submission_tech_df
#             Filter: select NAAT (amplification-based) submissions only
#             → naat_submissions
#   Step 4  — Filter microbial_df and organism_mapping to naat_submissions
#   Step 5  - Split labels into "binomial/uninomial" candidates (Subset 1) vs. 
#             noisy /multi-agent / panel labels (Subset 2).
#   Step 6  - Resolve Subset 1 labels to a unique NCBI Taxonomy ID -> pathogen_set1.
#   Step 7  - Expand Subset 2 multi-agent labels (e.g. "CT/GC"), apply acronym
#             expansion, and resolve to NCBI Taxonomy IDs -> pathogen_set2.
#   Step 8  - For 510(k) submissions still unmapped, fetch the FDA summary PDF and use
#             an LLM (NER) to extract organism names -> pathogen_set3.
#   Step 9  - Merge all three sets, dedupe by TaxID, and attach taxonomic metadata
#             (rank, division, genus, species, ...) from NCBI -> FDA510kPathogens.
################################################################################

# ── 0.  PACKAGES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(rvest)
  library(httr)
  library(tidyverse)
  library(rentrez)
  library(stringr)
  library(tidyr)
  library(dplyr)
  library(pdftools)      
  library(ellmer)        
  library(jsonlite)      
  library(cachem)       
})


# ── 0.1  CACHE SETUP ──────────────────────────────────────────────────────────
cache_dir  <- "cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir)
disk_cache <- cachem::cache_disk(cache_dir)

# forcing refresh flags — setting any to TRUE to re-run that step from scratch
FORCE_REFRESH <- list(
  fda_table             = FALSE,
  organism_mapping      = FALSE,
  fda_pdfs_all          = FALSE,  
  tech_classify         = FALSE,   
  taxdb_subset1         = FALSE,
  taxdb_subset2         = FALSE,
  llm_ner_resp_viruses  = FALSE,   # Step 7b: NER on "Respiratory Viruses" submissions
  llm_ner               = FALSE,
  taxonomy_summary      = FALSE
)

# ── 0.2 resource  ────────────────────────────────────────────────────────────
FDA_IVD_URL  <- "https://www.fda.gov/medical-devices/in-vitro-diagnostics/nucleic-acid-based-tests#microbial"
USER_AGENT   <- "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36"
GEMINI_MODEL <- "gemini-2.5-flash"

# validating technology categories returned by the classifier
VALID_TECH <- c(
  "PCR",                      
  "Other_Amplification",      
  "Hybridization",            
  "Hybrid_PCR_Hybridization", 
  "Unclear"                  
)

# Technology classifier system prompt (defined once, passed into every fresh chat)
TECH_SYSTEM_PROMPT <- paste0(
  "You are classifying an FDA 510(k) summary document for an in vitro diagnostic device. ",
  "Classify the PRIMARY detection technology as exactly ONE of these five labels:\n",
  "  'PCR'                      — polymerase chain reaction, real-time PCR, RT-PCR, qPCR.\n",
  "  'Other_Amplification'      — TMA, LAMP, SDA, NASBA, or any other isothermal or ",
  "                               target-amplification method that is NOT PCR.\n",
  "  'Hybridization'            — nucleic acid probe hybridization WITHOUT target amplification, ",
  "                               including Hybrid Capture (HC2), branched-DNA (bDNA) signal ",
  "                               amplification, FISH, dot-blot, or in situ hybridization. ",
  "                               Key distinction: the TARGET nucleic acid is NOT amplified.\n",
  "  'Hybrid_PCR_Hybridization' — PCR or other target amplification that ALSO uses hybridization ",
  "                               probes (e.g. TaqMan, molecular beacons) for fluorescent readout. ",
  "                               The target IS amplified; hybridisation is only the detection step.\n",
  "  'Unclear'                  — technology cannot be determined from the provided text.\n\n",
  "Base your decision on whether the TARGET NUCLEIC ACID itself is amplified before detection, ",
  "not merely on whether the word 'hybridization' appears anywhere in the text.\n\n",
  "Return ONLY valid JSON — no markdown fences, no prose — in this exact schema:\n",
  '{"technology": "<label>", "evidence_quote": "<verbatim excerpt from text, max 120 chars>", ',
  '"confidence": "<high|medium|low>"}'
)

NER_RESP_VIRUS_PROMPT <- paste0(
  "You are a regulatory affairs expert analyzing an FDA 510(k) summary for a NAAT assay. ",
  "List every specific respiratory virus that this assay is explicitly designed to detect. ",
  "Use full scientific or accepted common names (e.g. 'Influenza A virus', 'SARS-CoV-2'). ",
  "CRUCIAL: You MUST EXCLUDE any 'off-panel', 'cross-reactivity', or interference viruses. ",
  "Extract the microorganisms and separate them into two categories. ",
  "Return ONLY a valid JSON object matching this exact schema, with no markdown fences:\n",
  '{"targets": ["Virus 1", "Virus 2"], "off_panel": ["Virus 3", "Virus 4"]}'
)

# Named Entity Recognition (NER) system prompt for targeted pathogens using the "decoy bin" strategy
NER_SYSTEM_PROMPT <- paste0(
  "You are a regulatory affairs expert analyzing an FDA 510(k) summary for a NAAT assay. ",
  "510(k) documents contain the targeted pathogens (usually in the Intended Use section) ",
  "and long lists of 'off-panel', 'cross-reactivity', or 'analytical specificity' organisms. ",
  "Extract the microorganisms and separate them into two categories. ",
  "Standardize abbreviations to full binomial nomenclature. ",
  "Return ONLY a valid JSON object matching this exact schema, with no markdown fences:\n",
  '{"targets": ["Organism 1", "Organism 2"], "off_panel": ["Organism 3", "Organism 4"]}'
)

# checking disk_cache, running expr() on miss or force_refresh.
with_cache <- function(key, force_refresh = FALSE, expr) {
  if (!force_refresh && disk_cache$exists(key)) {
    message(sprintf("[cache HIT]  %s", key))
    return(disk_cache$get(key))
  }
  message(sprintf("[cache MISS] %s — computing...", key))
  result <- expr()
  disk_cache$set(key, result)
  result
}

# get_fda_summary_text 
get_fda_summary_text <- function(submission_id) {
  if (str_starts(submission_id, "K")) {
    yr_prefix <- substr(submission_id, 2, 3)
    url <- paste0(
      "https://www.accessdata.fda.gov/cdrh_docs/pdf",
      yr_prefix, "/", submission_id, ".pdf"
    )
  } else {
    return(NA_character_)   # DEN / PMA submissions: not handled yet
  }
  tryCatch({
    raw_text   <- pdftools::pdf_text(url)
    clean_text <- paste(raw_text, collapse = " ") |> str_squish()
    clean_text
  }, error = function(e) {
    message(sprintf("Could not fetch PDF for: %s", submission_id))
    NA_character_
  })
}

################################################################################
# STEP 1 — Ingest FDA IVD microbial table → microbial_df
################################################################################
microbial_df <- with_cache("01_fda_table", FORCE_REFRESH$fda_table, function() {
  response <- GET(FDA_IVD_URL, add_headers(`User-Agent` = USER_AGENT))
  if (status_code(response) != 200)
    stop(sprintf("Failed to fetch FDA page. Status: %d", status_code(response)))
  
  tables <- content(response, as = "text") |>
    read_html() |>
    html_table(fill = TRUE)
  
  # Microbial table is the second table on the page
  tbl <- tables[[2]]
  colnames(tbl) <- make.names(colnames(tbl))
  
  tbl |>
    mutate(
      Submission = str_extract(Submission, "^[A-Za-z]+\\d+")
    ) |>
    rename(
      Microorganism = Microorganism,
      Trade.Name    = Trade.Name,
      Manufacturer  = Manufacturer
    )
})
################################################################################
# STEP 2 — Build organism_mapping: unique FDAregulatoryLabel × Submissions
################################################################################
organism_mapping <- with_cache("02_organism_mapping", FORCE_REFRESH$organism_mapping, function() {
  microbial_df |>
    mutate(FDAregulatoryLabels = str_squish(Microorganism)) |>
    group_by(FDAregulatoryLabels) |>
    summarise(
      Submissions = paste(unique(Submission), collapse = ","),
      .groups     = "drop"
    )
})

################################################################################
# STEP 3 — Technology classification
#
#   3a. collecting every unique Submission ID across all organism_mapping rows
#   3b. fetching PDF summary text for each submission (cached per submission)
#   3c. running LLM technology classifier on each PDF text: submission_tech_df  
#   (columns: Submission, Technology, EvidenceQuote, Confidence)
#   3d. defining naat_submissions: Submission IDs where Technology ∈
#       {PCR, Other_Amplification, Hybrid_PCR_Hybridization}
#       (i.e. any assay where the target nucleic acid is amplified)
################################################################################

# ── 3a. Unique Submission IDs ─────────────────────────────────────────────────
all_submission_ids <- organism_mapping$Submissions |>
  str_split(",") |>
  unlist() |>
  str_squish() |>
  unique() |>
  sort()

# ── 3b. fetching PDF texts for all submissions ────────────────────────────────────
all_pdf_texts <- with_cache("03_fda_pdfs_all", FORCE_REFRESH$fda_pdfs_all, function() {
  texts <- lapply(all_submission_ids, get_fda_summary_text)
  names(texts) <- all_submission_ids
  texts
})

# partitioning into fetchable vs. failed
pdf_na    <- names(which(is.na(all_pdf_texts)))
pdf_valid <- all_pdf_texts[!names(all_pdf_texts) %in% pdf_na]

message(sprintf(
  "Step 3b: %d PDFs fetched successfully, %d failed (NA)",
  length(pdf_valid), length(pdf_na)
))

# ── 3c. LLM Technology Classifier ────────────────────────────────────────────
#   Key design decisions (informed by conversation history):
#     • Fresh chat_google_gemini() object per document → no history contamination
#     • responseMimeType = "application/json" → API-level JSON enforcement
#     • Schema validated before accepting result
#     • give_up flag checked immediately after tryCatch to avoid ghost retries

is_valid_tech_response <- function(p) {
  is.list(p) &&
    isTRUE(p[["technology"]] %in% VALID_TECH) &&
    length(p[["technology"]]) == 1
}

get_technology <- function(doc_id, text_data,
                           max_retries   = 3,
                           model         = GEMINI_MODEL,
                           system_prompt = TECH_SYSTEM_PROMPT) {
  empty_row <- data.frame(
    Submission    = doc_id,
    Technology    = NA_character_,
    EvidenceQuote = NA_character_,
    Confidence    = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (is.na(text_data) || !nzchar(text_data)) {
    message(sprintf("[%s] No text available — skipping classification", doc_id))
    return(empty_row)
  }
  
  for (attempt in seq_len(max_retries)) {
    parsed  <- NULL
    give_up <- FALSE
    
    tryCatch({
      # Fresh object per call — zero cross-document history contamination
      chat_obj <- chat_google_gemini(
        model         = model,
        system_prompt = system_prompt,
        api_args      = list(
          generationConfig = list(
            temperature      = 0,
            responseMimeType = "application/json"
          )
        )
      )
      raw    <- chat_obj$chat(text_data)
      # Strip markdown fences defensively even though responseMimeType should prevent them
      json   <- raw |>
        str_replace("(?s)^\\s*```(?:json)?\\s*", "") |>
        str_replace("```\\s*$", "") |>
        str_trim()
      parsed <- fromJSON(json)
    }, error = function(e) {
      if (grepl("429", conditionMessage(e))) {
        message(sprintf("[%s] Rate-limit (attempt %d) — %s. Sleeping 45s",
                        doc_id, attempt, conditionMessage(e)))
        Sys.sleep(45)
      } else {
        message(sprintf("[%s] Fatal error: %s", doc_id, conditionMessage(e)))
        give_up <<- TRUE
      }
    })
    
    # Exit retry loop as soon as we have a valid result or a fatal error
    if (give_up || is_valid_tech_response(parsed)) break
  }
  
  if (!is_valid_tech_response(parsed)) return(empty_row)
  
  data.frame(
    Submission    = doc_id,
    Technology    = parsed[["technology"]],
    EvidenceQuote = parsed[["evidence_quote"]] %||% NA_character_,
    Confidence    = parsed[["confidence"]]     %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

# running classifier over all fetchable PDFs (cached)
submission_tech_df <- with_cache("03_tech_classify", FORCE_REFRESH$tech_classify, function() {
  results <- Map(
    get_technology,
    doc_id    = names(pdf_valid),
    text_data = pdf_valid
  )
  dplyr::bind_rows(results)
})

# appending NA rows for submissions whose PDF could not be fetched
if (length(pdf_na) > 0) {
  na_rows <- data.frame(
    Submission    = pdf_na,
    Technology    = NA_character_,
    EvidenceQuote = NA_character_,
    Confidence    = NA_character_,
    stringsAsFactors = FALSE
  )
  submission_tech_df <- dplyr::bind_rows(submission_tech_df, na_rows)
}

message(sprintf("Step 3c complete: technology classified for %d submissions", nrow(submission_tech_df)))
message("Technology distribution:")
print(table(submission_tech_df$Technology, useNA = "ifany"))

# ── 3d. defining the NAAT gate ─────────────────────────────────────────────────
# Amplification-based technologies: target nucleic acid is amplified.
# Hybridization-only and Unclear are excluded.
AMPLIFICATION_TECH <- c("PCR", "Other_Amplification", "Hybrid_PCR_Hybridization")

naat_submissions <- submission_tech_df |>
  filter(Technology %in% AMPLIFICATION_TECH) |>
  pull(Submission)

excluded_submissions <- submission_tech_df |>
  filter(!Technology %in% AMPLIFICATION_TECH | is.na(Technology))

message(sprintf(
  "Step 3d: %d NAAT submissions retained | %d excluded (Hybridization / Unclear / NA)",
  length(naat_submissions),
  nrow(excluded_submissions)
))

################################################################################
# STEP 4 — Filter microbial_df and organism_mapping to NAAT submissions only
#
#   A Submission is NAAT if ANY of its comma-separated Submission IDs is in
#   naat_submissions.  We keep the organism_mapping row when at least one of
#   the Submissions linked to that label is a NAAT submission.
################################################################################

# check: does at least one token in the comma-separated Submissions string
# appear in the NAAT set?
has_naat_submission <- function(submissions_str, naat_set) {
  ids <- str_split(submissions_str, ",")[[1]] |> str_squish()
  any(ids %in% naat_set)
}

organism_mapping_naat <- organism_mapping |>
  filter(
    vapply(Submissions,
           has_naat_submission,
           logical(1),
           naat_set = naat_submissions)
  ) |>
  # trimming Submissions field to only keep NAAT IDs in the string
  mutate(
    Submissions = vapply(Submissions, function(s) {
      ids <- str_split(s, ",")[[1]] |> str_squish()
      paste(ids[ids %in% naat_submissions], collapse = ",")
    }, character(1))
  )

# ==============================================================================
# STEP 5: Partition into binomial/uninomial candidates (Subset 1) vs.
#           noisy / panel / multi-agent labels (Subset 2)
# ==============================================================================

# selecting binomial names (Genus + species): filtering for exactly two words 
# consisting only of alphanumeric characters and dealing with the "Uninomial" 
# virus instances: some entries are a single word (e.g. "Adenovirus", 
# "Enterovirus", "Norovirus").
# defining the binomial or uninomial patterns
BiOrUninomial_pattern <- "^[A-Za-z0-9]+( [A-Za-z0-9]+)?$"
stop_words            <- "(?i)resistance|panel|media|material|toxin|illness|specimen|test|kit|system"

ncbi_patches_1 <- c(
  "Clostridium difficile" = "Clostridioides difficile",
  "Entercoccus"           = "Enterococcus",
  "Human Papillomvirus"   = "Human papillomavirus",
  "Leishmaniaspecies"     = "Leishmania",
  "Mycobacterium species" = "Mycobacterium",
  "Mycoplasma pneumoniae" = "Mycoplasmoides pneumoniae",
  "Streptococci"          = "Streptococcus",
  "Variola"               = "Variola virus"
)

ncbi_patches_2 <- c(
  "Candida parapsilosiis" = "Candida parapsilosis",
  "Candida kruseii"       = "Pichia kudriavzevii"
)

taxonomic_patches <- c(
  "Streptococcus vini"  = "Streptococcus mutans",
  "Candida rugosa"      = "Diutina rugosa",
  "Variola"             = "Variola virus",
  "Enterobacter dispar" = "Enterobacter cloacae"
)

acronym_map <- c(
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

# partitionning into subset_1 (bi/uninomial) and subset_2 (other) ──────
subset_1 <- organism_mapping_naat |>
  filter(str_detect(FDAregulatoryLabels, BiOrUninomial_pattern))

subset_2 <- organism_mapping_naat |>
  filter(!str_detect(FDAregulatoryLabels, BiOrUninomial_pattern))

# filtering subset_1 entries based on top-words
subset_1_1 <- subset_1 |>
  filter(!str_detect(FDAregulatoryLabels, stop_words))

subset_2 <- bind_rows(
  subset_2,
  subset_1 |> filter(str_detect(FDAregulatoryLabels, stop_words))
)

# applying spelling patches
subset_1_1 <- subset_1_1 |>
  mutate(
    FDAregulatoryLabels = ifelse(
      FDAregulatoryLabels %in% names(ncbi_patches_1),
      ncbi_patches_1[FDAregulatoryLabels],
      FDAregulatoryLabels
    )
  )

# resolving subset 1 labels to NCBI Taxonomy IDs -> pathogen_set1
taxdbquery1 <- with_cache("taxdbquery1", FORCE_REFRESH$taxdb_subset1, function() {
  sapply(subset_1_1$FDAregulatoryLabels, function(x) {
    entrez_search(db = "taxonomy", term = x)$count
  })
})

subset_2 <- bind_rows(
  subset_2,
  subset_1_1 |> filter(FDAregulatoryLabels %in% names(which(taxdbquery1 != 1)))
)
subset_1_1 <- subset_1_1 |>
  filter(FDAregulatoryLabels %in% names(which(taxdbquery1 == 1)))

pathogen_set1 <- do.call(rbind, lapply(subset_1_1$FDAregulatoryLabels, function(x) {
  res <- entrez_search(db = "taxonomy", term = x)
  data.frame(
    Microorganism = x,
    TaxID         = if (length(res$ids) > 0) unlist(res$ids)[1] else NA_character_,
    stringsAsFactors = FALSE
  )
}))
pathogen_set1 <- merge(subset_1_1, pathogen_set1,
                       by.x = "FDAregulatoryLabels", by.y = "Microorganism")

# ==============================================================================
# Step 7  - Expand Subset 2 multi-agent labels (e.g. "CT/GC"), apply acronym
#             expansion, and resolve to NCBI Taxonomy IDs -> pathogen_set2.
# ==============================================================================

process_subset2_agents <- function(df, mapping_vec) {
  df |>
    mutate(
      Microorganism_Expanded = ifelse(
        FDAregulatoryLabels %in% names(mapping_vec),
        mapping_vec[FDAregulatoryLabels],
        FDAregulatoryLabels
      ),
      Microorganism_NCBItaxDb_ready = Microorganism_Expanded |>
        str_remove_all("(?i)\\b(DNA|RNA|rRNA|Culture|Antigen|Toxin|Detection|Nuclear -1)\\b") |>
        str_remove_all("\\(.*?\\)") |>
        str_squish()
    )
}

subset_2_expanded <- subset_2 |>
  filter(str_detect(FDAregulatoryLabels, "/|,|\\band\\b")) |>
  separate_rows(FDAregulatoryLabels, sep = "/|,|\\band\\b") |>
  mutate(FDAregulatoryLabels = str_squish(FDAregulatoryLabels))

subset_2_expanded_1 <- subset_2_expanded |>
  filter(!str_detect(FDAregulatoryLabels, stop_words)) |>
  group_by(FDAregulatoryLabels) |>
  summarise(Submissions = paste(unique(Submissions), collapse = ", "), .groups = "drop")

subset_2_expanded_1 <- subset_2_expanded_1 |>
  mutate(
    FDAregulatoryLabels = ifelse(
      FDAregulatoryLabels %in% names(ncbi_patches_2),
      ncbi_patches_2[FDAregulatoryLabels],
      FDAregulatoryLabels
    )
  )

subset_2_cleaned <- process_subset2_agents(subset_2_expanded_1, acronym_map) |>
  filter(!Microorganism_NCBItaxDb_ready %in% pathogen_set1$FDAregulatoryLabels) |>
  group_by(Microorganism_NCBItaxDb_ready) |>
  summarise(Submissions = paste(unique(Submissions), collapse = ", "), .groups = "drop")

taxdbQuery2 <- with_cache("taxdbquery2", FORCE_REFRESH$taxdb_subset2, function() {
  sapply(subset_2_cleaned$Microorganism_NCBItaxDb_ready, function(x) {
    entrez_search(db = "taxonomy", term = x)$count
  })
})

subset_2_cleaned <- subset_2_cleaned |>
  filter(Microorganism_NCBItaxDb_ready %in% names(which(taxdbQuery2 == 1)))

pathogen_set2 <- do.call(rbind, lapply(subset_2_cleaned$Microorganism_NCBItaxDb_ready, function(x) {
  res <- entrez_search(db = "taxonomy", term = x)
  data.frame(
    Microorganism = x,
    TaxID         = if (length(res$ids) > 0) unlist(res$ids)[1] else NA_character_,
    stringsAsFactors = FALSE
  )
}))
pathogen_set2 <- merge(subset_2_cleaned, pathogen_set2,
                       by.x = "Microorganism_NCBItaxDb_ready", by.y = "Microorganism") |>
  dplyr::rename(FDAregulatoryLabels = Microorganism_NCBItaxDb_ready)

pathogen_set1_2 <- bind_rows(pathogen_set1, pathogen_set2) |>
  filter(FDAregulatoryLabels != "Sepsis")

# ==============================================================================
# Step 7b — Targeted NER for "Respiratory Viruses" submissions
#
# "Respiratory Viruses" passed the binomial/uninomial filter and mapped to an
# NCBI TaxID, but that TaxID covers all viruses — far too broad to be useful.
# This step:
#   1. Identifies all Submission IDs linked to "Respiratory Viruses".
#   2. Fetches their FDA summary PDFs (reusing all_pdf_texts where available).
#   3. Runs a targeted LLM NER prompt asking for specific respiratory virus names.
#   4. Resolves the extracted names to unique NCBI TaxIDs.
#   5. Removes the "Respiratory Viruses" row from pathogen_set1_2 and appends
#      the resolved granular entries as pathogen_set_resp_viruses.
# ==============================================================================


# 1. Extract Submission IDs linked to "Respiratory Viruses"
resp_virus_row <- pathogen_set1_2 |>
  filter(FDAregulatoryLabels == "Respiratory Viruses")

resp_virus_submissions <- character(0)
if (nrow(resp_virus_row) > 0) {
  resp_virus_submissions <- resp_virus_row$Submissions |>
    str_split(",\\s*") |>
    unlist() |>
    str_squish() |>
    unique()
}

# 2. Fetch PDF texts — reuse all_pdf_texts cache, fetch any missing ones
resp_virus_texts <- lapply(resp_virus_submissions, function(sid) {
  if (sid %in% names(all_pdf_texts) && !is.na(all_pdf_texts[[sid]])) {
    return(all_pdf_texts[[sid]])
  }
  get_fda_summary_text(sid)
})
names(resp_virus_texts) <- resp_virus_submissions
resp_virus_texts_valid  <- resp_virus_texts[!is.na(resp_virus_texts)]

message(sprintf(
  "Step 7b: %d PDFs available for NER (%d could not be fetched)",
  length(resp_virus_texts_valid),
  length(resp_virus_submissions) - length(resp_virus_texts_valid)
))

# 3. Run targeted NER — fresh Chat object per document, no history contamination
resp_virus_ner_raw <- with_cache(
  "07b_llm_ner_resp_viruses",
  FORCE_REFRESH$llm_ner_resp_viruses,
  function() {
    results <- Map(
      get_organisms_safely,          # reuse the existing NER function
      doc_id        = names(resp_virus_texts_valid),
      text_data     = resp_virus_texts_valid,
      MoreArgs      = list(
        max_retries   = 3,
        model         = GEMINI_MODEL,
        system_prompt = NER_RESP_VIRUS_PROMPT   # targeted prompt
      )
    )
    dplyr::bind_rows(results)
  }
)

# 4. Build a label × Submissions mapping from NER results
resp_virus_mapping <- resp_virus_ner_raw |>
  filter(!is.na(FDAregulatoryLabels), nzchar(FDAregulatoryLabels)) |>
  # Apply the same taxonomic patches used elsewhere in the pipeline
  mutate(
    FDAregulatoryLabels = ifelse(
      FDAregulatoryLabels %in% names(taxonomic_patches),
      taxonomic_patches[FDAregulatoryLabels],
      FDAregulatoryLabels
    )
  ) |>
  group_by(FDAregulatoryLabels) |>
  summarise(
    Submissions = paste(unique(Submission), collapse = ", "),
    .groups     = "drop"
  )

# Resolve extracted names to unique NCBI TaxIDs
taxdbQuery_resp <- sapply(resp_virus_mapping$FDAregulatoryLabels, function(x) {
  entrez_search(db = "taxonomy", term = x)$count
})

resp_virus_mapping_1to1 <- resp_virus_mapping |>
  filter(FDAregulatoryLabels %in% names(which(taxdbQuery_resp == 1)))

message(sprintf(
  "Step 7b: %d unique respiratory virus labels → %d resolved to a unique NCBI TaxID",
  nrow(resp_virus_mapping),
  nrow(resp_virus_mapping_1to1)
))

pathogen_set_resp_viruses <- do.call(rbind, lapply(
  resp_virus_mapping_1to1$FDAregulatoryLabels,
  function(x) {
    res <- entrez_search(db = "taxonomy", term = x)
    data.frame(
      FDAregulatoryLabels = x,
      TaxID = if (length(res$ids) > 0) unlist(res$ids)[1] else NA_character_,
      stringsAsFactors = FALSE
    )
  }
))
pathogen_set_resp_viruses <- merge(
  resp_virus_mapping_1to1,
  pathogen_set_resp_viruses,
  by = "FDAregulatoryLabels"
)

# 5. Remove the overly-broad "Respiratory Viruses" row and append granular results
pathogen_set1_2 <- pathogen_set1_2 |>
  filter(FDAregulatoryLabels != "Respiratory Viruses") |>
  bind_rows(pathogen_set_resp_viruses)

message(sprintf(
  "Step 7b complete: 'Respiratory Viruses' replaced by %d specific virus entries",
  nrow(pathogen_set_resp_viruses)
))


# ==============================================================================
# Step 8  - For 510(k) submissions still unmapped, fetching the FDA summary PDF 
#           and using an LLM (NER) to extract organism names -> pathogen_set3.
# ==============================================================================

mapped_submissions <- pathogen_set1_2$Submissions |>
  str_split(",") |>
  unlist() |>
  str_squish() |>
  unique()

# applying the NER on unmapped submissions that are in the NAAT set
unmapped_naat <- setdiff(naat_submissions, mapped_submissions)
message(sprintf("Step 9: %d NAAT submissions still unmapped → running LLM NER", length(unmapped_naat)))

# reusing already-fetched PDF texts where available; fetch any remaining ones
get_text_for_submission <- function(sid) {
  if (sid %in% names(all_pdf_texts) && !is.na(all_pdf_texts[[sid]])) {
    return(all_pdf_texts[[sid]])
  }
  get_fda_summary_text(sid)
}

unmapped_texts <- lapply(unmapped_naat, get_text_for_submission)
names(unmapped_texts) <- unmapped_naat
unmapped_texts_valid  <- unmapped_texts[!is.na(unmapped_texts)]

# NER function — fresh Chat object per document (no history contamination)
get_organisms_safely <- function(doc_id, text_data,
                                 max_retries   = 3,
                                 model         = GEMINI_MODEL,
                                 system_prompt = NER_SYSTEM_PROMPT) {
  empty_row <- data.frame(
    Submission          = doc_id,
    FDAregulatoryLabels = NA_character_,
    stringsAsFactors    = FALSE
  )
  
  if (is.na(text_data) || !nzchar(text_data)) return(empty_row)
  
  for (attempt in seq_len(max_retries)) {
    parsed  <- NULL
    give_up <- FALSE
    
    tryCatch({
      # Fresh Chat object per call — stateless, no history
      ner_chat <- chat_google_gemini(
        model         = model,
        system_prompt = system_prompt,
        api_args      = list(
          generationConfig = list(
            temperature      = 0,
            responseMimeType = "application/json"
          )
        )
      )
      raw    <- ner_chat$chat(text_data)
      json   <- raw |>
        str_replace("(?s)^\\s*```(?:json)?\\s*", "") |>
        str_replace("```\\s*$", "") |>
        str_trim()
      parsed <- fromJSON(json)
    }, error = function(e) {
      if (grepl("429", conditionMessage(e))) {
        message(sprintf("[NER %s] Rate-limit (attempt %d) — sleeping 45s", doc_id, attempt))
        Sys.sleep(45)
      } else {
        message(sprintf("[NER %s] Fatal error: %s", doc_id, conditionMessage(e)))
        give_up <<- TRUE
      }
    })
    
    if (give_up || !is.null(parsed)) break
  }
  
  # if (is.null(parsed) || !is.character(parsed)) return(empty_row)
  # if (length(parsed) == 0) return(empty_row)
  # 
  # data.frame(
  #   Submission          = doc_id,
  #   FDAregulatoryLabels = parsed,
  #   stringsAsFactors    = FALSE
  # Extract only the "targets" array from the JSON object
  targets <- parsed[["targets"]]
  
  if (is.null(targets) || length(targets) == 0) return(empty_row)
  
  data.frame(
    Submission          = doc_id,
    FDAregulatoryLabels = targets,  # Only keep the targeted pathogens
    stringsAsFactors    = FALSE
  )
}

subset_3_raw <- with_cache("07_llm_ner", FORCE_REFRESH$llm_ner, function() {
  results <- Map(
    get_organisms_safely,
    doc_id    = names(unmapped_texts_valid),
    text_data = unmapped_texts_valid
  )
  dplyr::bind_rows(results)
})

subset_3 <- subset_3_raw |>
  filter(!is.na(FDAregulatoryLabels))

subset_3_mapping <- subset_3 |>
  group_by(FDAregulatoryLabels) |>
  summarise(Submissions = paste(unique(Submission), collapse = ","), .groups = "drop")

# merging subset_3 organisms into pathogen_set1_2 where labels already exist
pathogen_set1_2 <- pathogen_set1_2 |>
  left_join(subset_3_mapping, by = "FDAregulatoryLabels", suffix = c("_old", "_new")) |>
  mutate(
    Submissions = case_when(
      is.na(Submissions_new) ~ Submissions_old,
      is.na(Submissions_old) ~ Submissions_new,
      TRUE ~ paste(Submissions_old, Submissions_new, sep = ",")
    )
  ) |>
  select(FDAregulatoryLabels, Submissions, TaxID)

# New organism labels from subset_3 not yet in pathogen_set1_2
subset_3_new <- subset_3_mapping |>
  filter(!FDAregulatoryLabels %in% pathogen_set1_2$FDAregulatoryLabels) |>
  mutate(
    FDAregulatoryLabels = ifelse(
      FDAregulatoryLabels %in% names(taxonomic_patches),
      taxonomic_patches[FDAregulatoryLabels],
      FDAregulatoryLabels
    )
  ) |>
  group_by(FDAregulatoryLabels) |>
  summarise(Submissions = paste(unique(Submissions), collapse = ","), .groups = "drop")

taxdbQuery3 <- sapply(subset_3_new$FDAregulatoryLabels, function(x) {
  entrez_search(db = "taxonomy", term = x)$count
})

subset_3_new <- subset_3_new |>
  filter(FDAregulatoryLabels %in% names(which(taxdbQuery3 == 1)))

pathogen_set3 <- do.call(rbind, lapply(subset_3_new$FDAregulatoryLabels, function(x) {
  res <- entrez_search(db = "taxonomy", term = x)
  data.frame(
    FDAregulatoryLabels = x,
    TaxID               = if (length(res$ids) > 0) unlist(res$ids)[1] else NA_character_,
    stringsAsFactors    = FALSE
  )
}))
pathogen_set3 <- merge(subset_3_new, pathogen_set3, by = "FDAregulatoryLabels")

# ==============================================================================
#   Step 9  - Merge all three sets, dedupe by TaxID, and attach taxonomic metadata
#             (rank, division, genus, species, ...) from NCBI -> FDA510kPathogens.
# ==============================================================================

FDA510kPathogens_NAAT <- bind_rows(pathogen_set1_2, pathogen_set3) |>
  group_by(TaxID) |>
  summarise(
    FDAregulatoryLabels = first(FDAregulatoryLabels),
    Submissions         = paste(unique(unlist(str_split(Submissions, ",\\s*"))),
                                collapse = ", "),
    .groups             = "drop"
  )

message(sprintf("Step 10: %d pathogen entries after combining sets", nrow(FDA510kPathogens_NAAT)))

# getting taxonomic rank level 
FDA510kPathogensDivisions <- with_cache("08_taxonomy_summary", FORCE_REFRESH$taxonomy_summary, function() {
  do.call(rbind, lapply(FDA510kPathogens_NAAT$TaxID, function(uid) {
    tryCatch({
      tax_rec <- entrez_summary(db = "taxonomy", uid)
      as.data.frame(
        t(extract_from_esummary(
          tax_rec,
          c("rank", "division", "scientificname", "commonname",
            "taxid", "genus", "species", "subsp"),
          simplify = TRUE
        ))
      ) |>
        unnest(rank, division, scientificname, commonname,
               taxid, genus, species, subsp)
    }, error = function(e) {
      message(sprintf("entrez_summary failed for TaxID %s: %s", uid, e$message))
      data.frame(taxid = uid, stringsAsFactors = FALSE)
    })
  }))
})

FDA510kPathogens_NAAT <- merge(
  FDA510kPathogens_NAAT,
  FDA510kPathogensDivisions,
  by.x = "TaxID",
  by.y = "taxid"
)

# attaching technology classification to final output 
# expanding Submissions in the tech table to one row per Submission ID,
# then joining the technology back to each pathogen entry via its Submissions.

# building a long-format lookup: SubmissionID → Technology
tech_lookup <- submission_tech_df |>
  select(Submission, Technology, Confidence) |>
  filter(!is.na(Technology))

# gathering the set of technologies from all its linked submissions for each 
# pathogen and collapsing to a single representative value.
# using that technology if all linked submissions agree.
# flagging as "Mixed" if mixed (e.g. one PCR and one Hybridization for the same 
# label.  
resolve_technology <- function(submissions_str, lookup_df) {
  ids   <- str_split(submissions_str, ",\\s*")[[1]] |> str_squish()
  techs <- lookup_df$Technology[lookup_df$Submission %in% ids] |>
    na.omit() |>
    unique()
  if (length(techs) == 0)  return(NA_character_)
  if (length(techs) == 1)  return(techs)
  paste(sort(techs), collapse = " / ")
}

FDA510kPathogens_NAAT <- FDA510kPathogens_NAAT |>
  mutate(
    Technology = vapply(
      Submissions,
      resolve_technology,
      character(1),
      lookup_df = tech_lookup
    )
  )

# getting the total number of submissions per pathogen
FDA510kPathogens_NAAT$Submission_Count <- sapply(1:dim(FDA510kPathogens_NAAT)[1], \(i) stringr::str_count(FDA510kPathogens_NAAT$Submissions[i], ",") + 1) 

# saving outputs 
saveRDS(FDA510kPathogens_NAAT, "FDA510kPathogens_NAAT.rds")
write.csv(FDA510kPathogens_NAAT, "FDA510kPathogens_NAAT.csv", row.names = FALSE)

# also saving excluded (hybridization) submissions
#saveRDS(excluded_submissions, "FDA510k_excluded_hybridization.rds")
write.csv(excluded_submissions, "FDA510k_excluded_hybridization.csv", row.names = FALSE)

message("Done. Outputs: FDA510kPathogens_NAAT.rds / .csv")
###############################################################################

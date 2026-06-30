# ampRegistry
Workflow dedicated to building a comprehensive repository of FDA-approved Nucleic Acid Amplification Tests (NAATs).
 
## Data source:  
The FDA https://www.fda.gov/medical-devices/in-vitro-diagnostics/nucleic-acid-based-tests#microbial which reports approved NAATs applied to detect pathogens.  
The data features 4 attributes: Microorganism, Trade Name, Manufacturer and a Submission ID.    
## Workflow:  
###   Step 1  — Ingest FDA IVD microbial table → microbial_df  
###   Step 2  — Build organism_mapping (unique label × Submission ID)  
###   Step 3  — Fetch PDF summary for every unique Submission ID  
* Run LLM technology classifier → submission_tech_df    
* Filter: select NAAT (amplification-based) submissions only → naat_submissions  
###   Step 4  — Filter microbial_df and organism_mapping to naat_submissions  
###   Step 5  - Split labels into "binomial/uninomial" candidates (Subset 1) vs. noisy /multi-agent / panel labels (Subset 2). 
###   Step 6  - Resolve Subset 1 labels to a unique NCBI Taxonomy ID -> pathogen_set1.  
###   Step 7  - Expand Subset 2 multi-agent labels (e.g. "CT/GC"), apply acronym expansion, and resolve to NCBI Taxonomy IDs -> pathogen_set2.   
###   Step 8  - For 510(k) submissions still unmapped, fetch the FDA summary PDF and use an LLM (NER) to extract organism names -> pathogen_set3.  
###   Step 9  - Merge all three sets, dedupe by TaxID, and attach taxonomic metadata (rank, division, genus, species, ...) from NCBI -> FDA510kPathogens.  

## =====================================================================
## Chlorella vulgaris — High Light (HL) vs Low Light (LL) DIA proteomics
## Full workflow mirroring course lessons 02–08, adapted to this dataset.
##
## Design: 3 HL vs 3 LL, unpaired, single two-group contrast.
##   -> maps directly onto the taught `~0 + group` + makeContrasts model.
##
## Main DIA-NN report to import:
##   DIANNout/Chlorella_light_proteomics.parquet
##   (NOT the *_generated_speclib.parquet — that lacks per-run quantities)
## =====================================================================

## ---- Libraries (same stack as the course) --------------------------
library("arrow")
library("QFeatures")
library("limma")
library("tidyr")
library("dplyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("factoextra")
library("patchwork")

data_dir <- "C:/Users/Raiyan Subedar/Downloads/Kaust BioSummer/Cambridge/New/Week2_Expression_Proteomics/Project/project_proteomics/chlorella_light_proteomics/DIANNout"
report   <- file.path(data_dir, "Chlorella_light_proteomics.parquet")

## =====================================================================
## Lesson 02 — Import & infrastructure
## =====================================================================
diann_df <- read_parquet(report)
## The HL/LL design is encoded in the Run names, so we build the sample
## metadata directly rather than reading an external sheet.
unique(diann_df$Run)
sample_metadata <- diann_df %>%
  distinct(Run) %>%
  mutate(
    runCol = Run,
    group  = if_else(grepl("HIGH-LIGHT", Run), "HL", "LL"),
    replicate = sub(".*_(\\d)$", "\\1", Run)
  )
sample_metadata

## =====================================================================
## Lesson 03
## =====================================================================

## Read the DIA-NN report into a QFeatures object (one set per run).
dia_qf <- readQFeaturesFromDIANN(
  assayData = diann_df,
  colData   = sample_metadata,
  quantCols = "Precursor.Quantity",
  runCol    = "Run",
  fnames    = "Precursor.Id"
)


dia_qf <- dia_qf %>%
  filterFeatures(~ Q.Value      <= 0.01) %>%   # run-level precursor FDR
  filterFeatures(~ PG.Q.Value   <= 0.01) %>%   # run-level protein-group FDR
  filterFeatures(~ Lib.Q.Value  <= 0.01) %>%   # library precursor FDR
  filterFeatures(~ Lib.PG.Q.Value <= 0.01)     # library protein-group FDR

hist(
  nrows(dia_qf),
  xlab = "Precursors",
  main = "Precursor counts per sample\nafter Q-value filtering"
)

## Merge the per-run sets into one precursor matrix (NA where not detected)
dia_qf <- joinAssays(dia_qf, i = names(dia_qf), name = "precursors")
dia_qf

## Drop the individual per-run sets, keep the joined "precursors" set
sets_to_rm <- which(names(dia_qf) != "precursors")
if (length(sets_to_rm)) dia_qf <- removeAssay(dia_qf, sets_to_rm)


rd <- rowData(dia_qf[['precursors']]) %>%
  data.frame()
table(grepl('Cont_', rd$Protein.Ids))
dia_qf <- addAssay(dia_qf, dia_qf[["precursors"]], name = "precursors_no_cont")
dia_qf <- addAssayLink(dia_qf, from = "precursors", to = "precursors_no_cont")
plot(dia_qf)

dia_qf <- filterFeatures(
  dia_qf,
  ~ !grepl("Cont_", Protein.Ids),
  i = "precursors_no_cont"
)
dia_qf


dia_qf <- addAssay(
  dia_qf,
  dia_qf[["precursors_no_cont"]],
  name = "precursors_filtered_missing"
)

dia_qf <- addAssayLink(
  dia_qf,
  from = "precursors_no_cont",
  to = "precursors_filtered_missing"
)
plot(dia_qf)


nNA(dia_qf, i = "precursors_no_cont")
## Zeros are missing values in DIA-NN output
dia_qf <- dia_qf %>%
  filterNA(pNA = 0.5, i = "precursors_filtered_missing")
dia_qf  
## =====================================================================
## Lesson 04 — Log2 transform, aggregate to protein, normalise
## =====================================================================
dia_qf <- logTransform(dia_qf, i = "precursors_filtered_missing",
                        name = "precursors_log")
dia_qf


dia_qf[["precursors_log"]] %>%
  assay() %>%
  longForm() %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  theme_bw() +
  xlab("Abundance (raw)")

## robustSummary: precursor -> protein-group abundance (Protein.Group as fcol)
dia_qf <- aggregateFeatures(
  dia_qf,
  i    = "precursors_log",
  fcol = "Protein.Group",
  name = "proteins",
  fun  = MsCoreUtils::robustSummary,
  na.rm = TRUE
)
dia_qf

# dia_qf[,,'proteins'] subsets the QFeatures object to the "proteins" set
longForm(dia_qf[,, 'proteins'], colvars = 'group') %>%
  data.frame() %>%
  head()

## Median normalisation across samples
dia_qf <- normalize(dia_qf, i = "proteins",
                     name = "norm_proteins", method = "center.median")

pre_norm <- longForm(dia_qf[,, 'proteins'], colvars = 'group') %>%
  ggplot(aes(x = value, colour = group, group = colname)) +
  geom_density() +
  theme_classic() +
  xlab("log2 (Abundance)") +
  ggtitle('Pre-normalisation')

post_norm <- longForm(dia_qf[,, 'norm_proteins'], colvars = 'group') %>%
  ggplot(aes(x = value, colour = group, group = colname)) +
  geom_density() +
  theme_classic() +
  xlab("log2 (Abundance)") +
  ggtitle('Post-normalisation')

pre_norm + post_norm

## =====================================================================
## Lesson 05 — Exploration: PCA of HL vs LL
## =====================================================================
protein_pca <- dia_qf[["norm_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp(scale = TRUE, center = TRUE)

summary(protein_pca)

fviz_screeplot(protein_pca)


p <- protein_pca$x %>%
  merge(colData(dia_qf), by = 'row.names') %>%
  data.frame() %>%
  ggplot(aes(x = PC1, y = PC2, colour = group, shape = replicate)) +
  geom_point(size = 3) +
  theme_bw()

print(p)
print(p + aes(x = PC3, y = PC4))
## =====================================================================
## Lesson 06 — Differential abundance with limma
##   Two groups, unpaired -> single HL vs LL contrast.
## =====================================================================
se <- getWithColData(dia_qf, "norm_proteins")
se

## Keep proteins quantified in >= 2 of 3 replicates in BOTH groups
hl <- se[, se$group == "HL"]; ll <- se[, se$group == "LL"]
tokeep <- (rowSums(!is.na(assay(hl))) >= 2) &
          (rowSums(!is.na(assay(ll))) >= 2)
dia_qf <- addAssay(dia_qf, se[tokeep, ], name = "norm_proteins_replicated")
dia_qf <- addAssayLink(dia_qf, from = "norm_proteins",
                       to = "norm_proteins_replicated")
print(table(tokeep))
plot(dia_qf)
dia_qf

## No-intercept design (estimates each group mean directly)
group  <- factor(dia_qf[["norm_proteins_replicated"]]$group, levels = c("LL", "HL"))
design <- model.matrix(formula(~ 0 + group))
colnames(design) <- levels(group)
design

contrasts_matrix <- makeContrasts(HL_vs_LL = HL - LL, levels = design)
contrasts_matrix

data <- assay(dia_qf[["norm_proteins_replicated"]])
fit  <- lmFit(data, design)
fit  <- contrasts.fit(fit, contrasts_matrix)
fit  <- eBayes(fit, trend = TRUE, robust = TRUE)   # variance moderation

plotSA(fit = fit, cex = 0.5, xlab = "Average log2 abundance")

## Format results
limma_results_F <- topTable(
  fit = fit,
  coef = NULL,
  adjust.method = "BH", # Method for multiple hypothesis testing
  number = Inf
) %>% # Print results for all proteins
  rownames_to_column("UniprotID")

head(limma_results_F)

# Check distribution of p-values
limma_results_F %>%
  ggplot(aes(x = P.Value)) +
  geom_histogram() +
  theme_classic() +
  ggtitle("Distribution of p-values for F-test")

## Descriptive protein names live only in the pg_matrix (the long report
## carries UniProt IDs), so pull them in for annotation by Protein.Group.
annot <- read.delim(file.path(data_dir, "Chlorella_light_proteomics.pg_matrix.tsv"),
                    check.names = FALSE) %>%
  select(Protein.Group, Genes, First.Protein.Description) %>%
  distinct()

res <- topTable(fit, coef = "HL_vs_LL", number = Inf, sort.by = "p") %>%
  data.frame() %>% rownames_to_column("Protein.Group") %>%
  merge(annot, by = "Protein.Group", all.x = TRUE) %>%
  arrange(P.Value)


sig <- topTable(fit, coef="HL_vs_LL", number=Inf)
sum(sig$adj.P.Val < 0.05 & abs(sig$logFC) > 1)          # the one to quote
table(sign(sig$logFC[sig$adj.P.Val < 0.05 & abs(sig$logFC) > 1]))  # up vs down

## =====================================================================
## Lesson 07 — Volcano plot
## =====================================================================
top30 <- res %>% filter(adj.P.Val < 0.05) %>% slice_min(P.Value, n = 30)
view(res)
ggplot(res, aes(logFC, -log10(P.Value),
                colour = adj.P.Val < 0.05, alpha = adj.P.Val < 0.05)) +
  geom_point(pch = 20, size = 2.5) +
  geom_text_repel(data = top30, aes(label = Genes), size = 2, colour = "grey20") +
  scale_colour_manual(values = c("grey", "orangered3")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  labs(x = "log2 FC (HL / LL)", y = "-log10 P") +
  theme_classic() + theme(legend.position = "none")

## =====================================================================
## Lesson 08 — GO over-representation analysis (adapted for Chlorella)
##   The course uses enrichGO(OrgDb = org.Hs.eg.db) — HUMAN only. Chlorella is
##   non-model with no Bioconductor OrgDb, so we do what the lesson-08 note
##   recommends: fetch GO annotations DIRECTLY for the UniProt IDs, then run
##   the generic clusterProfiler::enricher() with a custom TERM2GENE table.
##   Split up- and down-sets rather than merging them.
##   Coverage caveat: many locus-tag proteins (D9Q98_*) lack GO in UniProt,
##   so enrichment is biased toward well-annotated processes.
## =====================================================================
library("httr")
library("readr")
library("clusterProfiler")
library("GO.db")          # GO term names + ontology (BP/MF/CC)

## 1. Protein groups -> a single UniProt accession (first ID of each group)
res <- res %>% mutate(acc = sub(";.*", "", Protein.Group))
up_HL <- res %>% filter(adj.P.Val < 0.05, logFC > 0) %>% pull(acc)  # up in high light
up_LL <- res %>% filter(adj.P.Val < 0.05, logFC < 0) %>% pull(acc)  # up in low light
bg    <- unique(res$acc)                                            # background

## 2. Fetch GO from UniProt REST, in batches of 100 accessions
fetch_go <- function(accs) {
  q <- paste0("accession:", accs, collapse = " OR ")
  r <- GET("https://rest.uniprot.org/uniprotkb/search",
           query = list(query = q, fields = "accession,go_id",
                        format = "tsv", size = 500))
  stop_for_status(r)
  read_tsv(content(r, "text", encoding = "UTF-8"), show_col_types = FALSE)
}
chunks <- split(bg, ceiling(seq_along(bg) / 100))
go_raw <- bind_rows(lapply(chunks, function(x) { Sys.sleep(0.3); fetch_go(x) }))

## 3. Build TERM2GENE (GO id -> accession) and TERM2NAME (GO id -> term)
term2gene <- go_raw %>%
  rename(acc = Entry, go = `Gene Ontology IDs`) %>%
  filter(!is.na(go), go != "") %>%
  separate_rows(go, sep = ";\\s*") %>%
  transmute(go = trimws(go), acc) %>%
  distinct()
go_info <- AnnotationDbi::select(GO.db, keys = unique(term2gene$go),
                                 columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
term2name <- go_info[, c("GOID", "TERM")]
cat("GO-annotated proteins:", n_distinct(term2gene$acc), "of", length(bg), "background\n")


## 4. Over-representation, split by direction
ego_up <- enricher(up_HL, universe = bg, TERM2GENE = term2gene,
                   TERM2NAME = term2name, pvalueCutoff = 0.05)
ego_dn <- enricher(up_LL, universe = bg, TERM2GENE = term2gene,
                   TERM2NAME = term2name, pvalueCutoff = 0.05)

## 5. Visualise. (simplify() needs GOSemSim data from an OrgDb, unavailable for
##    Chlorella, so we just cap the number of terms to keep the plot readable.)
dotplot(ego_up, showCategory = 15) + ggtitle("Enriched in HIGH light")
dotplot(ego_dn, showCategory = 10) + ggtitle("Enriched in LOW light")

write.csv(as.data.frame(ego_up), "chlorella_results/GO_enriched_up_in_HL.csv", row.names = FALSE)
write.csv(as.data.frame(ego_dn), "chlorella_results/GO_enriched_up_in_LL.csv", row.names = FALSE)

## =====================================================================
## UNIQUE ANGLE — Coordinated photo-acclimation as a mechanism,
##   not a gene list. Tag functional groups and show they move as blocks.
## =====================================================================
panel <- res %>%
  mutate(desc = tolower(First.Protein.Description),
         cat = case_when(
           grepl("light-harvest|chlorophyll a-b|lhc", desc) ~ "LHC antenna",
           grepl("photosystem", desc)                        ~ "Photosystem core",
           grepl("carbonic anhydrase", desc)                 ~ "Carbonic anhydrase (CCM)",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(cat))

## order: LHC on top, CCM at bottom (first factor level plots at bottom)
lev <- c("Carbonic anhydrase (CCM)", "Photosystem core", "LHC antenna")
panel$cat <- factor(panel$cat, levels = lev)

## n per group, for the right-hand labels
counts <- panel %>% count(cat) %>% mutate(lab = paste0("n=", n))

cols <- c("LHC antenna" = "#6a8f1f",
          "Photosystem core" = "#1b9aa8",
          "Carbonic anhydrase (CCM)" = "#b06fb0")

ggplot(panel, aes(logFC, cat, colour = cat, fill = cat)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.7) +
  geom_boxplot(alpha = 0.30, outlier.shape = NA, width = 0.55,
               colour = "black") +
  geom_jitter(height = 0.14, size = 2.2, stroke = 0.3,
              shape = 21, colour = "black") +
  geom_text(data = counts, aes(x = 1.75, y = cat, label = lab),
            inherit.aes = FALSE, hjust = 0, size = 4.5, colour = "grey20") +
  annotate("text", x = -2.9, y = 0.35, label = "\u2190 reduced in high light",
           hjust = 0, fontface = "italic", size = 4, colour = "#b03030") +
  annotate("text", x = 1.9, y = 0.35, label = "increased \u2192",
           hjust = 1, fontface = "italic", size = 4, colour = "#3060b0") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(-3.2, 2.0)) +
  labs(title = "Antenna, photosystems and carbon-capture all down-regulate in high light",
       x = "log2 fold-change  (high light / low light)", y = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 14))
## =====================================================================
## END. Save results table.
## =====================================================================
write.csv(res, "chlorella_results/chlorella_HL_vs_LL_differential_abundance.csv",
          row.names = FALSE)
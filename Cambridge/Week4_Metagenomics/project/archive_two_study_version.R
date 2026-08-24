# =============================================================================
#  Week 4 mini project - 16S amplicon analysis
#
#  QUESTION
#    Is the "C. difficile infection" gut microbiome signature specific to
#    C. difficile, or is it largely a signature of having diarrhoea?
#
#  DATA
#    Duvallet et al. 2017, Nat Commun 8:1784 (s41467-017-01973-8)
#      PART 1  cdi_schubert   - Schubert et al. 2014 mBio.  H / nonCDI / CDI
#      PART 2  cdi_youngster  - Youngster et al. 2014.      H / CDI / postFMT
#
#  HOW TO RUN
#    Open in RStudio, set the working directory to the folder holding the
#    *_results directories, then run top to bottom (Ctrl+Enter line by line,
#    or Ctrl+Shift+S to source the whole file).
#      setwd("~/Downloads/Week 4 mini project")
#
#  RELATION TO THE COURSE MATERIAL ("13-16S sequencing")
#    The published data are already denoised - 100% identity de novo OTUs are
#    equivalent to DADA2 ASVs - so the DADA2 pre-processing steps
#    (filterAndTrim / learnErrors / dada / mergePairs / removeBimeraDenovo)
#    cannot be re-run: there are no FASTQ files and no quality scores.
#    We join the course workflow at assignTaxonomy() and continue through
#    phyloseq exactly as in the practical.
# =============================================================================


# =============================================================================
#  PART 0 - SETUP
# =============================================================================

library("dada2")
library("phyloseq")
library("Biostrings")
library("tidyverse")
library("vegan")
theme_set(theme_bw())

OUT_DIR       <- "results"    # figures and tables are written here
MIN_DEPTH     <- 1000         # drop samples with fewer reads than this
RAREFY_DEPTH  <- 2000         # even depth, used for alpha diversity only
MIN_PREV      <- 0.10         # a genus must appear in >=10% of samples to test
MIN_OTU_READS <- 1            # drop OTUs seen fewer times than this in total
SEED          <- 42

# Taxonomic ranks kept for every OTU. Counts are always summed at GENUS level
# (see 1.6), but the full lineage is carried through so that plots can be drawn
# at any of these levels, exactly as in the course practical:
#     plot_bar(ps, fill = "Phylum")   /   "Class"   /   "Order"   /   "Family"
RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# SILVA's "wSpecies" training set also returns a Species column. Set
# USE_SPECIES to FALSE to stop at Genus.
#
# Expect most Species entries to be NA regardless: a 150-200 bp fragment
# of the 16S gene usually cannot pin down a species, and dada2 deliberately
# leaves it blank rather than guessing.
USE_SPECIES <- TRUE
ALL_RANKS   <- if (USE_SPECIES) c(RANKS, "Species") else RANKS

PLOT_LEVEL <- "Genus"         # level used for the composition bar plots

# About MIN_OTU_READS. Naming the sequences (assignTaxonomy, below) is by far
# the slowest step, and most OTUs are near-singletons. Raising this to 10 makes
# taxonomy roughly 5x faster:
#     cdi_schubert   19,314 -> 3,976 OTUs, still 96.8% of all reads
#     cdi_youngster  70,819 -> 12,033 OTUs, still 91.6% of all reads
#
# But it is NOT free. Those rare OTUs still contribute to genus totals, and
# dropping them changes the headline result:
#     MIN_OTU_READS = 1   ->  197 genera, 53 CDI-associated, 92% shared
#     MIN_OTU_READS = 10  ->  143 genera, 51 CDI-associated, 86% shared
# So it is off by default. Set it to 10 only if you want a quick first pass,
# and quote the numbers from a full run.

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# To also capture the console output to a file, uncomment these two lines and
# the matching sink() at the very bottom of the script:
# sink(file.path(OUT_DIR, "analysis_log.txt"), split = TRUE)


# ---- Two small helpers ------------------------------------------------------
# Everything else in this script is written out in sequence. These two are
# functions only because they are each used several times further down.

# (1) Cliff's delta: a non-parametric effect size between -1 and +1.
#     0 means the two groups overlap completely.
cliffs_delta <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
}

# (2) Wilcoxon test of every genus between two groups, with BH correction.
#     Used three times in Part 1.
test_genera <- function(rel_mat, groups, g1, g2, min_prev = MIN_PREV) {
  i1 <- which(groups == g1)
  i2 <- which(groups == g2)
  prev <- rowMeans(rel_mat[, c(i1, i2), drop = FALSE] > 0)
  m    <- rel_mat[prev >= min_prev, , drop = FALSE]

  res <- data.frame(
    genus = rownames(m),
    p     = apply(m, 1, function(v)
              suppressWarnings(wilcox.test(v[i1], v[i2])$p.value)),
    delta = apply(m, 1, function(v) cliffs_delta(v[i1], v[i2])),
    row.names = NULL
  )
  res$q         <- p.adjust(res$p, method = "BH")
  res$direction <- sign(res$delta)
  res[order(res$q), ]
}


# ---- The SILVA reference ----------------------------------------------------
# assignTaxonomy() needs the SILVA training set. It is ~130 MB, so it is
# downloaded once and reused. A failed download leaves a zero-byte stub behind,
# which would silently be treated as a valid reference on the next run - so the
# file size is checked before trusting it.

SILVA_FILE <- "silva_nr99_v138.1_wSpecies_train_set.fa.gz"
SILVA_URL  <- paste0("https://zenodo.org/record/4587955/files/",
                     "silva_nr99_v138.1_wSpecies_train_set.fa.gz")

silva_ok <- file.exists(SILVA_FILE) && file.info(SILVA_FILE)$size >= 50e6

if (!silva_ok) {
  # Only ever delete a file that actually fails the size check, so that setting
  # silva_ok by hand cannot destroy a good 130 MB download.
  if (file.exists(SILVA_FILE) && file.info(SILVA_FILE)$size < 50e6) {
    cat("Removing truncated SILVA file (", file.info(SILVA_FILE)$size, " bytes)\n")
    unlink(SILVA_FILE)
  }
  cat("Downloading the SILVA training set (~130 MB), one time only ...\n")
  try(download.file(SILVA_URL, SILVA_FILE, mode = "wb"), silent = TRUE)
  silva_ok <- file.exists(SILVA_FILE) && file.info(SILVA_FILE)$size >= 50e6
  if (!silva_ok && file.exists(SILVA_FILE)) unlink(SILVA_FILE)
}

# SILVA is the only taxonomy source in this script, so stop rather than carry on
# with something half-built. If the download keeps failing, fetch the file
# manually from SILVA_URL and drop it next to this script.
if (!silva_ok)
  stop("SILVA reference not available. Download it manually from:\n  ", SILVA_URL,
       "\nand place it in: ", normalizePath("."), call. = FALSE)

cat("SILVA reference found (",
    round(file.info(SILVA_FILE)$size / 1e6), "MB ) - assigning taxonomy with assignTaxonomy()\n")


# #############################################################################
#  PART 1 - cdi_schubert
#  Is the CDI signature specific to C. difficile, or is it diarrhoea?
# #############################################################################

cat("\n===== PART 1: cdi_schubert =====\n")

# ---- 1.1 Read the OTU count table ------------------------------------------
# Rows are OTUs (denovo1, denovo2, ...), columns are samples.

otu_sch <- as.matrix(read.delim(
  "cdi_schubert_results/cdi_schubert.otu_table.100.denovo",
  row.names = 1, check.names = FALSE))

cat("OTU table:", nrow(otu_sch), "OTUs x", ncol(otu_sch), "samples\n")


# ---- 1.2 Read the sample metadata ------------------------------------------
# quote = "" and comment.char = "" stop R from choking on apostrophes and #
# characters inside the free-text columns.

md_sch <- read.delim("cdi_schubert_results/cdi_schubert.metadata.txt",
                     check.names = FALSE, colClasses = "character",
                     quote = "", comment.char = "")
md_sch <- md_sch[!duplicated(md_sch[[1]]), ]
rownames(md_sch) <- md_sch[[1]]

cat("Metadata:", nrow(md_sch), "samples,", ncol(md_sch), "columns\n")
table(md_sch$DiseaseState)


  # ---- 1.3 Drop very rare OTUs (speeds up taxonomy) --------------------------

# abundant <- rowSums(otu_sch) >= MIN_OTU_READS
# cat("Keeping", sum(abundant), "of", length(abundant), "OTUs with >=",
#     MIN_OTU_READS, "reads (",
#     sprintf("%.1f%%", 100 * sum(otu_sch[abundant, ]) / sum(otu_sch)),
#     "of reads )\n")
# otu_sch <- otu_sch[abundant, ]


# ---- 1.4 Assign taxonomy ----------------------------------------------------
# This is the course step:
#     taxa <- assignTaxonomy(seqtab.nochim, "silva_..._train_set.fa.gz")
# Our equivalent input is the representative sequence for each OTU.
#
# IMPORTANT - SEQUENCE ORIENTATION. The two studies store their sequences in
# OPPOSITE directions relative to the SILVA reference. Probing 30-mers from 60
# random OTUs against 150,000 SILVA reference sequences gives:
#     cdi_schubert    0 forward matches, 100,229 reverse-complement matches
#     cdi_youngster   6,181 forward matches, 0 reverse-complement matches
# So cdi_schubert is stored reverse-complemented and cdi_youngster is not.
# tryRC = TRUE makes assignTaxonomy test both directions, which is why it is
# essential here - without it, cdi_schubert would come back almost entirely
# unassigned. (Do not "optimise" this by flipping the sequences yourself: that
# would fix one study and break the other.)
#
# Expect this step to take tens of minutes the first time. The result is
# cached to an .rds file, so later runs are instant.

seqs_sch <- readDNAStringSet("cdi_schubert_results/cdi_schubert.otu_seqs.100.fasta")
seqs_sch <- seqs_sch[names(seqs_sch) %in% rownames(otu_sch)]

cache_sch <- file.path(OUT_DIR, "cdi_schubert.silva_taxonomy.rds")
if (file.exists(cache_sch)) {
  tax_sch <- readRDS(cache_sch)
  cat("Loaded cached SILVA taxonomy\n")
} else {
  cat("Running assignTaxonomy() on", length(seqs_sch),
      "sequences - go and make a coffee ...\n")
  tax_sch <- assignTaxonomy(as.character(seqs_sch), SILVA_FILE,
                            multithread = TRUE, tryRC = TRUE, minBoot = 50)
  rownames(tax_sch) <- names(seqs_sch)
  saveRDS(tax_sch, cache_sch)
}

taxa_sch  <- tax_sch[, ALL_RANKS]
genus_sch <- taxa_sch[, "Genus"]

# Sanity check on orientation. If tryRC had failed, almost nothing would be
# named and everything downstream would be quietly empty - so fail loudly here
# instead of producing an near-empty genus table 200 lines later.
if (mean(!is.na(genus_sch)) < 0.20)
  stop("Only ", round(100 * mean(!is.na(genus_sch))), "% of OTUs got a genus. ",
       "That usually means the sequences are in the wrong orientation - ",
       "check tryRC = TRUE.", call. = FALSE)

# How much of the data survives at each level? Genus is the strictest.
cat("Assignment rate by rank:",
    paste(sprintf("%s %.0f%%", colnames(taxa_sch),
                  100 * colMeans(!is.na(taxa_sch))),
          collapse = ", "), "\n")


# ---- 1.5 Match the tables up and drop shallow samples -----------------------

keep_otu <- intersect(rownames(otu_sch), rownames(taxa_sch))
otu_sch   <- otu_sch[keep_otu, ]
taxa_sch  <- taxa_sch[keep_otu, , drop = FALSE]
genus_sch <- taxa_sch[, "Genus"]

keep_samp <- intersect(colnames(otu_sch), rownames(md_sch))
cat("Dropped", ncol(otu_sch) - length(keep_samp), "samples with no metadata\n")
otu_sch <- otu_sch[, keep_samp]
md_sch  <- md_sch[keep_samp, ]

deep <- colSums(otu_sch) >= MIN_DEPTH
cat("Dropped", sum(!deep), "samples with <", MIN_DEPTH, "reads\n")
otu_sch <- otu_sch[, deep]
md_sch  <- md_sch[deep, ]

# keep only the three groups we are comparing
in_groups <- md_sch$DiseaseState %in% c("H", "nonCDI", "CDI")
otu_sch <- otu_sch[, in_groups]
md_sch  <- md_sch[in_groups, ]
md_sch$DiseaseState <- factor(md_sch$DiseaseState, levels = c("H", "nonCDI", "CDI"))

cat("Final:", nrow(otu_sch), "OTUs x", ncol(otu_sch), "samples\n")
print(table(md_sch$DiseaseState))


# ---- 1.6 Collapse to genus level -------------------------------------------
# ~19,000 individual sequences is too fine-grained to interpret. Summing the
# counts of every OTU that shares a genus gives ~200 named groups.
# (Equivalent to phyloseq's tax_glom(ps, "Genus"), but much faster.)

named   <- !is.na(genus_sch)
gen_sch <- rowsum(otu_sch[named, ], group = genus_sch[named])

cat("Genus level:", nrow(gen_sch), "genera;",
    sprintf("%.1f%%", 100 * sum(gen_sch) / sum(otu_sch)), "of reads kept\n")

# Carry the full lineage across to the collapsed table, so every genus still
# knows its Phylum / Class / Order / Family. Each genus name maps to exactly
# one lineage in these data (checked - no ambiguous cases), so taking the
# first OTU per genus is safe.
first_of_genus <- taxa_sch[named, , drop = FALSE][!duplicated(genus_sch[named]), , drop = FALSE]
rownames(first_of_genus) <- genus_sch[named][!duplicated(genus_sch[named])]
lineage_sch <- first_of_genus[rownames(gen_sch), , drop = FALSE]


# ---- 1.7 Build the phyloseq object -----------------------------------------
# Same three ingredients as the course practical: counts, sample metadata,
# and taxonomy.

ps_sch <- phyloseq(otu_table(gen_sch, taxa_are_rows = TRUE),
                   sample_data(md_sch),
                   tax_table(as.matrix(lineage_sch)))
ps_sch

# The taxonomy table now holds all six ranks, so any of them can be used for
# plotting or subsetting:
head(tax_table(ps_sch))
# e.g. subset_taxa(ps_sch, Phylum == "Firmicutes")


# ---- 1.8 Alpha diversity ----------------------------------------------------
# Diversity depends on how deeply a sample was sequenced, so all samples are
# first subsampled to the same depth.

counts_sch <- t(otu_sch)                                  # samples x OTUs
deep_enough <- rowSums(counts_sch) >= RAREFY_DEPTH
set.seed(SEED)
rare_sch <- rrarefy(counts_sch[deep_enough, ], RAREFY_DEPTH)

alpha_sch <- data.frame(
  Observed     = specnumber(rare_sch),
  Shannon      = diversity(rare_sch, index = "shannon"),
  DiseaseState = md_sch$DiseaseState[deep_enough]
)
cat("Rarefied to", RAREFY_DEPTH, "reads:", nrow(alpha_sch), "samples retained\n")

cat(sprintf("Shannon, Kruskal-Wallis across the 3 groups: p = %.3g\n",
            kruskal.test(Shannon ~ DiseaseState, data = alpha_sch)$p.value))

for (pair in list(c("H", "nonCDI"), c("H", "CDI"), c("nonCDI", "CDI"))) {
  x <- alpha_sch$Shannon[alpha_sch$DiseaseState == pair[1]]
  y <- alpha_sch$Shannon[alpha_sch$DiseaseState == pair[2]]
  cat(sprintf("   %-7s vs %-7s : p = %-10.3g Cliff's delta = %+.2f\n",
              pair[1], pair[2],
              suppressWarnings(wilcox.test(x, y)$p.value), cliffs_delta(x, y)))
}

p_alpha_sch <- alpha_sch %>%
  pivot_longer(c(Observed, Shannon), names_to = "metric", values_to = "value") %>%
  ggplot(aes(DiseaseState, value, fill = DiseaseState)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = NULL, y = NULL, title = "Alpha diversity (rarefied)") +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "schubert_alpha_diversity.png"), p_alpha_sch,
       width = 7, height = 4, dpi = 150)


# ---- 1.9 Beta diversity, ordination and PERMANOVA --------------------------
# Convert to relative abundance, then measure how different every pair of
# samples is (Bray-Curtis) and squash that down to two dimensions (NMDS).

rel_sch  <- transform_sample_counts(ps_sch, function(x) x / sum(x))
bray_sch <- phyloseq::distance(rel_sch, method = "bray")

ord_sch <- ordinate(rel_sch, method = "NMDS", distance = "bray",
                    trymax = 20, trace = 0)
cat(sprintf("NMDS stress = %.3f%s\n", ord_sch$stress,
            if (ord_sch$stress > 0.2) "  <- POOR (>0.2), read the plot with caution" else ""))

p_ord_sch <- plot_ordination(rel_sch, ord_sch, color = "DiseaseState") +
  geom_point(size = 2) + stat_ellipse(level = 0.68) +
  labs(title = "Bray-Curtis NMDS, genus level")
ggsave(file.path(OUT_DIR, "schubert_nmds.png"), p_ord_sch,
       width = 6.5, height = 5, dpi = 150)

adonis_sch <- adonis2(bray_sch ~ DiseaseState, data = md_sch, permutations = 999)
cat(sprintf("PERMANOVA DiseaseState: R2 = %.3f, p = %.3g\n",
            adonis_sch$R2[1], adonis_sch$`Pr(>F)`[1]))

# Is non-CDI diarrhoea intermediate between healthy and CDI, or not?
bm_sch  <- as.matrix(bray_sch)
healthy <- rownames(md_sch)[md_sch$DiseaseState == "H"]
d2h_sch <- rowMeans(bm_sch[, healthy])

cat("Mean Bray-Curtis distance to the healthy group:\n")
for (g in levels(md_sch$DiseaseState))
  cat(sprintf("   %-7s %.3f\n", g, mean(d2h_sch[md_sch$DiseaseState == g])))


# ---- 1.10 THE KEY TEST -------------------------------------------------------
# How much of the "CDI signature" is shared with diarrhoea that is NOT CDI?

rel_mat_sch <- as(otu_table(rel_sch), "matrix")
grp_sch     <- as.character(md_sch$DiseaseState)

cdi_vs_h    <- test_genera(rel_mat_sch, grp_sch, "CDI",    "H")
noncdi_vs_h <- test_genera(rel_mat_sch, grp_sch, "nonCDI", "H")
cdi_vs_non  <- test_genera(rel_mat_sch, grp_sch, "CDI",    "nonCDI")

sig_cdi <- cdi_vs_h %>% filter(q < 0.05)

compare <- sig_cdi %>%
  select(genus, delta_cdi = delta, dir_cdi = direction) %>%
  left_join(noncdi_vs_h %>% select(genus, q_non = q, delta_non = delta,
                                   dir_non = direction),
            by = "genus")

shared <- compare %>% filter(!is.na(q_non), q_non < 0.05, dir_non == dir_cdi)

cat("\n--- HEADLINE NUMBERS ---\n")
cat("Genera tested (present in >=", 100 * MIN_PREV, "% of samples): ",
    nrow(cdi_vs_h), "\n", sep = "")
cat("Genera differing between CDI and healthy (q<0.05): ", nrow(sig_cdi), "\n", sep = "")
cat("   ...of those, also shifted the SAME WAY in non-CDI diarrhoea: ",
    nrow(shared), sprintf(" (%.0f%%)", 100 * nrow(shared) / nrow(sig_cdi)), "\n", sep = "")
cat("Genera differing between CDI and non-CDI diarrhoea: ",
    sum(cdi_vs_non$q < 0.05, na.rm = TRUE), "\n\n", sep = "")

write.csv(cdi_vs_h,    file.path(OUT_DIR, "schubert_CDI_vs_healthy.csv"),    row.names = FALSE)
write.csv(noncdi_vs_h, file.path(OUT_DIR, "schubert_nonCDI_vs_healthy.csv"), row.names = FALSE)
write.csv(cdi_vs_non,  file.path(OUT_DIR, "schubert_CDI_vs_nonCDI.csv"),     row.names = FALSE)

p_cmp_sch <- compare %>%
  filter(!is.na(delta_non)) %>%
  ggplot(aes(delta_non, delta_cdi)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.7) +
  labs(x = "Effect size: non-CDI diarrhoea vs healthy",
       y = "Effect size: CDI vs healthy",
       title = "Genera altered in CDI are altered the same way in diarrhoea alone")
ggsave(file.path(OUT_DIR, "schubert_effect_size_comparison.png"), p_cmp_sch,
       width = 6.5, height = 5.5, dpi = 150)


# ---- 1.11 The obvious confounder: antibiotics -------------------------------
# These patients had also taken far more antibiotics, which damage gut
# bacteria on their own. Note the column is literally named "antibiotics >3mo",
# and as.data.frame() would rewrite that name via make.names().

abx_sch <- tolower(trimws(md_sch[["antibiotics >3mo"]]))
has_abx <- abx_sch %in% c("yes", "no")

cat("--- Antibiotic exposure ---\n")
abx_tab <- table(abx_sch[has_abx], md_sch$DiseaseState[has_abx])
for (g in colnames(abx_tab))
  cat(sprintf("   %-7s %3d/%3d exposed (%.0f%%)\n", g, abx_tab["yes", g],
              sum(abx_tab[, g]), 100 * abx_tab["yes", g] / sum(abx_tab[, g])))

# Antibiotics entered FIRST, so disease is only credited with variance that
# antibiotic exposure cannot already explain.
sub_sch <- prune_samples(has_abx, rel_sch)
adonis_abx <- adonis2(
  phyloseq::distance(sub_sch, "bray") ~ antibiotics + DiseaseState,
  data = data.frame(antibiotics  = factor(abx_sch[has_abx]),
                    DiseaseState = md_sch$DiseaseState[has_abx]),
  permutations = 999, by = "terms")
cat(sprintf("   PERMANOVA, antibiotics first: antibiotics R2 = %.3f (p = %.3g); disease R2 = %.3f (p = %.3g)\n",
            adonis_abx$R2[1], adonis_abx$`Pr(>F)`[1],
            adonis_abx$R2[2], adonis_abx$`Pr(>F)`[2]))

# Strongest version of the question: among antibiotic-exposed patients only,
# can CDI still be told apart from non-CDI diarrhoea?
expd <- has_abx & abx_sch == "yes" & md_sch$DiseaseState %in% c("CDI", "nonCDI")
ps_expd <- prune_samples(expd, rel_sch)
adonis_expd <- adonis2(
  phyloseq::distance(ps_expd, "bray") ~ DiseaseState,
  data = data.frame(DiseaseState = factor(as.character(md_sch$DiseaseState[expd]))),
  permutations = 999)
cat(sprintf("   Within antibiotic-exposed patients only, CDI vs nonCDI (n = %d): R2 = %.3f, p = %.3g\n\n",
            sum(expd), adonis_expd$R2[1], adonis_expd$`Pr(>F)`[1]))


# ---- 1.12 Composition bar plot (course step) --------------------------------

# This is the course step:
#     plot_bar(ps.top20, fill = "Class") + facet_grid(~Facility, ...)
# PLOT_LEVEL is set in Part 0 - change it to "Phylum", "Class", "Order",
# "Family" or "Genus" and re-run this block. Nothing else needs to change.
#
# Higher levels are also better covered, because fewer OTUs go unnamed - the
# script prints the exact rates per rank in 1.4. A Phylum plot therefore uses
# more of your data than a Genus one, and Species least of all.

# GOTCHA: phyloseq's psmelt() silently drops taxonomy columns that are entirely
# NA, so asking to fill by a rank that was never assigned (Species is the likely
# one) fails with a confusing "object 'Species' not found" rather than an empty
# plot. Check before plotting and drop back to Genus.
plot_level_sch <- PLOT_LEVEL
if (all(is.na(tax_table(ps_sch)[, plot_level_sch]))) {
  cat("NOTE:", plot_level_sch,
      "was never assigned in this dataset - plotting at Genus instead\n")
  plot_level_sch <- "Genus"
}

top15_sch <- names(sort(taxa_sums(rel_sch), decreasing = TRUE))[1:15]
p_bar_sch <- prune_taxa(top15_sch, rel_sch) %>%
  plot_bar(fill = plot_level_sch) +
  facet_grid(~DiseaseState, scales = "free_x", space = "free") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, title = paste("Top 15 genera, coloured by", plot_level_sch))
ggsave(file.path(OUT_DIR, paste0("schubert_composition_", plot_level_sch, ".png")),
       p_bar_sch, width = 11, height = 5, dpi = 150)

# To see every level at once, as the practical does with Class / Order:
for (lvl in c("Phylum", "Class", "Order", "Family")) {
  p <- prune_taxa(top15_sch, rel_sch) %>%
    plot_bar(fill = lvl) +
    facet_grid(~DiseaseState, scales = "free_x", space = "free") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = NULL, title = paste("Top 15 genera, coloured by", lvl))
  ggsave(file.path(OUT_DIR, paste0("schubert_composition_", lvl, ".png")),
         p, width = 11, height = 5, dpi = 150)
}


# #############################################################################
#  PART 2 - cdi_youngster
#  Does a faecal transplant reverse it?
# #############################################################################

cat("\n===== PART 2: cdi_youngster =====\n")

# ---- 2.1 Read the OTU table and metadata -----------------------------------

otu_yng <- as.matrix(read.delim(
  "cdi_youngster_results/cdi_youngster.otu_table.100.denovo",
  row.names = 1, check.names = FALSE))

md_yng <- read.delim("cdi_youngster_results/cdi_youngster.metadata.txt",
                     check.names = FALSE, colClasses = "character",
                     quote = "", comment.char = "")
md_yng <- md_yng[!duplicated(md_yng[[1]]), ]
rownames(md_yng) <- md_yng[[1]]

cat("OTU table:", nrow(otu_yng), "OTUs x", ncol(otu_yng), "samples\n")


# ---- 2.2 Assign taxonomy ----------------------------------------------------

seqs_yng <- readDNAStringSet("cdi_youngster_results/cdi_youngster.otu_seqs.100.fasta")
seqs_yng <- seqs_yng[names(seqs_yng) %in% rownames(otu_yng)]

cache_yng <- file.path(OUT_DIR, "cdi_youngster.silva_taxonomy.rds")
if (file.exists(cache_yng)) {
  tax_yng <- readRDS(cache_yng)
  cat("Loaded cached SILVA taxonomy\n")
} else {
  cat("Running assignTaxonomy() on", length(seqs_yng),
      "sequences - this is the bigger of the two studies ...\n")
  tax_yng <- assignTaxonomy(as.character(seqs_yng), SILVA_FILE,
                            multithread = TRUE, tryRC = TRUE, minBoot = 50)
  rownames(tax_yng) <- names(seqs_yng)
  saveRDS(tax_yng, cache_yng)
}

taxa_yng  <- tax_yng[, ALL_RANKS]
genus_yng <- taxa_yng[, "Genus"]

cat("Assignment rate by rank:",
    paste(sprintf("%s %.0f%%", colnames(taxa_yng),
                  100 * colMeans(!is.na(taxa_yng))), collapse = ", "), "\n")

if (mean(!is.na(genus_yng)) < 0.20)
  stop("Only ", round(100 * mean(!is.na(genus_yng))), "% of OTUs got a genus. ",
       "Check tryRC = TRUE.", call. = FALSE)


# ---- 2.3 Match up, filter, and label the groups ----------------------------

keep_otu <- intersect(rownames(otu_yng), rownames(taxa_yng))
otu_yng   <- otu_yng[keep_otu, ]
taxa_yng  <- taxa_yng[keep_otu, , drop = FALSE]
genus_yng <- taxa_yng[, "Genus"]

keep_samp <- intersect(colnames(otu_yng), rownames(md_yng))
cat("Dropped", ncol(otu_yng) - length(keep_samp), "samples with no metadata\n")
otu_yng <- otu_yng[, keep_samp]
md_yng  <- md_yng[keep_samp, ]

deep <- colSums(otu_yng) >= MIN_DEPTH
cat("Dropped", sum(!deep), "samples with <", MIN_DEPTH, "reads\n")
otu_yng <- otu_yng[, deep]
md_yng  <- md_yng[deep, ]

md_yng$Group <- factor(
  recode(md_yng$DiseaseState,
         "CDI" = "CDI pre-FMT", "postFMT_CDI" = "Post-FMT", "H" = "Donor (healthy)"),
  levels = c("CDI pre-FMT", "Post-FMT", "Donor (healthy)"))

print(table(md_yng$Group))


# ---- 2.4 Genus level and phyloseq object -----------------------------------

named   <- !is.na(genus_yng)
gen_yng <- rowsum(otu_yng[named, ], group = genus_yng[named])
cat("Genus level:", nrow(gen_yng), "genera;",
    sprintf("%.1f%%", 100 * sum(gen_yng) / sum(otu_yng)), "of reads kept\n")

first_of_genus <- taxa_yng[named, , drop = FALSE][!duplicated(genus_yng[named]), , drop = FALSE]
rownames(first_of_genus) <- genus_yng[named][!duplicated(genus_yng[named])]
lineage_yng <- first_of_genus[rownames(gen_yng), , drop = FALSE]

ps_yng  <- phyloseq(otu_table(gen_yng, taxa_are_rows = TRUE),
                    sample_data(md_yng),
                    tax_table(as.matrix(lineage_yng)))
rel_yng <- transform_sample_counts(ps_yng, function(x) x / sum(x))

# Composition plot for this study too, at whichever level PLOT_LEVEL is set to
plot_level_yng <- PLOT_LEVEL
if (all(is.na(tax_table(ps_yng)[, plot_level_yng]))) {
  cat("NOTE:", plot_level_yng,
      "was never assigned in this dataset - plotting at Genus instead\n")
  plot_level_yng <- "Genus"
}

top15_yng <- names(sort(taxa_sums(rel_yng), decreasing = TRUE))[1:15]
p_bar_yng <- prune_taxa(top15_yng, rel_yng) %>%
  plot_bar(fill = plot_level_yng) +
  facet_grid(~Group, scales = "free_x", space = "free") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, title = paste("Top 15 genera, coloured by", plot_level_yng))
ggsave(file.path(OUT_DIR, paste0("youngster_composition_", plot_level_yng, ".png")),
       p_bar_yng, width = 10, height = 5, dpi = 150)


# ---- 2.5 THE KEY PLAYERS: does FMT bring them back? ------------------------
#
#  This is the heart of Part 2, and the basis of slide 2.
#
#  Part 1 found which genera collapse (or bloom) in CDI. Here we take that
#  SAME list and ask one simple question of each one:
#
#      after a faecal transplant, does it return to the level seen in donors?
#
#  Part 1 must have been run first, because the list comes from there.

if (!exists("cdi_vs_h"))
  stop("Run Part 1 first - the key player list is built from its results.",
       call. = FALSE)

# --- (a) the list of key players, most heavily depleted first ---------------

key <- cdi_vs_h[cdi_vs_h$q < 0.05, c("genus", "delta")]
key <- key[order(key$delta), ]
key <- key[key$genus %in% taxa_names(rel_yng), ]   # must also exist in this study

cat("\nKey players carried over from Part 1:", nrow(key), "of",
    sum(cdi_vs_h$q < 0.05), "\n")

# --- (b) average abundance of each one, in each of the three groups ---------

rel_mat_yng <- as(otu_table(rel_yng), "matrix")
grp_yng     <- as.character(md_yng$Group)

recovery <- data.frame(
  genus  = key$genus,
  status = ifelse(key$delta < 0, "lost in CDI", "blooms in CDI"),
  cdi    = NA_real_,   # mean abundance before FMT
  post   = NA_real_,   # mean abundance after FMT
  donor  = NA_real_,   # mean abundance in the healthy donors
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(recovery))) {
  abundances <- rel_mat_yng[recovery$genus[i], ]
  recovery$cdi[i]   <- mean(abundances[grp_yng == "CDI pre-FMT"])
  recovery$post[i]  <- mean(abundances[grp_yng == "Post-FMT"])
  recovery$donor[i] <- mean(abundances[grp_yng == "Donor (healthy)"])
}

# --- (c) how far back towards the donor level did it get? ------------------
#     0.0 = did not move at all
#     1.0 = landed exactly on the donor level
#    >1.0 = overshot past the donor level
#
# CAUTION: this is a ratio, so it blows up when the CDI and donor levels are
# nearly the same - dividing by a gap of almost nothing gives nonsense like
# "962% recovered". Only compute it when there is a real gap to close.

MIN_GAP <- 0.001                      # 0.1% relative abundance
gap <- recovery$donor - recovery$cdi
recovery$recovered <- ifelse(abs(gap) < MIN_GAP, NA,
                             (recovery$post - recovery$cdi) / gap)
cat("(", sum(is.na(recovery$recovered)),
    "genera have too small a CDI-to-donor gap for a recovery % to mean anything)\n")

# --- (d) two simple tests for each genus -----------------------------------
#     did it move?           post-FMT vs CDI pre-FMT   (small q = yes, it moved)
#     is it back to normal?  post-FMT vs donors        (LARGE q = yes, it matches)

recovery$p_moved    <- NA_real_
recovery$p_vs_donor <- NA_real_

for (i in seq_len(nrow(recovery))) {
  abundances <- rel_mat_yng[recovery$genus[i], ]
  before <- abundances[grp_yng == "CDI pre-FMT"]
  after  <- abundances[grp_yng == "Post-FMT"]
  donor  <- abundances[grp_yng == "Donor (healthy)"]
  recovery$p_moved[i]    <- suppressWarnings(wilcox.test(before, after)$p.value)
  recovery$p_vs_donor[i] <- suppressWarnings(wilcox.test(after,  donor)$p.value)
}

recovery$q_moved    <- p.adjust(recovery$p_moved,    method = "BH")
recovery$q_vs_donor <- p.adjust(recovery$p_vs_donor, method = "BH")

# A genus counts as RESTORED if it moved after FMT and can no longer be told
# apart from the donors.
#
# BE HONEST ABOUT THIS ONE. "No longer significantly different from donors" is
# not the same as "the same as donors" - with only 18 donor samples the test has
# limited power, so some genera pass simply because there is not enough data to
# show a difference. Read `restored` alongside `recovered`: if a genus is
# flagged restored but only closed 20% of the gap, the honest reading is
# "moved in the right direction, not demonstrably back to normal".
recovery$restored <- recovery$q_moved < 0.05 & recovery$q_vs_donor > 0.05

# --- (e) the headline numbers for slide 2 ----------------------------------

cat("\n--- KEY PLAYER RECOVERY AFTER FMT ---\n")
cat("Key genera tracked:                    ", nrow(recovery), "\n")
cat("Moved significantly after FMT:         ", sum(recovery$q_moved < 0.05, na.rm = TRUE), "\n")
cat("Restored (moved, and now match donors):", sum(recovery$restored, na.rm = TRUE), "\n")
cat("Still differ from donors after FMT:    ", sum(recovery$q_vs_donor < 0.05, na.rm = TRUE), "\n\n")

print(head(recovery[, c("genus", "status", "cdi", "post", "donor",
                        "recovered", "restored")], 12), digits = 2)

write.csv(recovery, file.path(OUT_DIR, "youngster_key_player_recovery.csv"),
          row.names = FALSE)

# --- (f) the slide 2 figure -------------------------------------------------
# One row per genus, three dots: where it sits before FMT, after FMT, and in
# the donors. If FMT works, the middle dot sits near the donor dot.

show   <- head(recovery, 12)        # the 12 most depleted key players
PSEUDO <- 1e-5                      # so that zeros can be shown on a log axis

plot_dat <- data.frame(
  genus     = rep(show$genus, 3),
  group     = rep(c("CDI pre-FMT", "Post-FMT", "Donor (healthy)"),
                  each = nrow(show)),
  abundance = c(show$cdi, show$post, show$donor)
)
plot_dat$group <- factor(plot_dat$group,
                         levels = c("CDI pre-FMT", "Post-FMT", "Donor (healthy)"))
plot_dat$genus <- factor(plot_dat$genus, levels = rev(show$genus))

p_key <- ggplot(plot_dat, aes(abundance + PSEUDO, genus, colour = group)) +
  geom_line(aes(group = genus), colour = "grey75") +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(x = "Mean relative abundance (log scale)", y = NULL, colour = NULL,
       title = "Key players lost in CDI, and where they sit after FMT")
ggsave(file.path(OUT_DIR, "youngster_key_players.png"), p_key,
       width = 8, height = 6, dpi = 150)


# ---- 2.6 How quickly does it happen? ---------------------------------------
# One summary number per sample: how far its whole community sits from the
# healthy donors. Then plot that against days since the transplant.

bray_yng <- phyloseq::distance(rel_yng, method = "bray")
bm_yng   <- as.matrix(bray_yng)
donors   <- rownames(md_yng)[md_yng$Group == "Donor (healthy)"]
d2d_yng  <- rowMeans(bm_yng[, donors])

post_yng <- rownames(md_yng)[md_yng$Group == "Post-FMT"]
traj <- data.frame(
  sample = post_yng,
  days   = suppressWarnings(as.numeric(md_yng[post_yng, "days_since_fmt"])),
  dist   = d2d_yng[post_yng]
) %>% filter(!is.na(days))

ct <- suppressWarnings(cor.test(traj$days, traj$dist, method = "spearman"))
cat(sprintf("Distance to donors vs days since FMT: rho = %.2f, p = %.3g (n = %d)\n",
            ct$estimate, ct$p.value, nrow(traj)))

p_traj <- ggplot(traj, aes(days, dist)) +
  geom_point(size = 2) + geom_smooth(method = "loess", se = TRUE, span = 1) +
  scale_x_log10() +
  labs(x = "Days since FMT (log scale)",
       y = "Bray-Curtis distance to donor pool",
       title = "Recovery trajectory after FMT")
ggsave(file.path(OUT_DIR, "youngster_trajectory.png"), p_traj,
       width = 6, height = 4.5, dpi = 150)



# #############################################################################
#  PART 3 - SUPPORTING ANALYSES
#
#  Everything below backs up Parts 1 and 2 but is NOT slide material for a
#  5-minute talk. Keep it for backup slides and for answering questions.
# #############################################################################

cat("\n===== PART 3: supporting analyses =====\n")

# ---- 3.1 Ordination and PERMANOVA ------------------------------------------
bray_yng <- phyloseq::distance(rel_yng, method = "bray")
ord_yng  <- ordinate(rel_yng, method = "NMDS", distance = "bray",
                     trymax = 20, trace = 0)
cat(sprintf("NMDS stress = %.3f%s\n", ord_yng$stress,
            if (ord_yng$stress > 0.2) "  <- POOR (>0.2), read the plot with caution" else ""))

p_ord_yng <- plot_ordination(rel_yng, ord_yng, color = "Group") +
  geom_point(size = 2.5) + stat_ellipse(level = 0.68) +
  labs(title = "FMT moves patients toward the healthy donor range")
ggsave(file.path(OUT_DIR, "youngster_nmds.png"), p_ord_yng,
       width = 6.5, height = 5, dpi = 150)

adonis_yng <- adonis2(bray_yng ~ Group, data = md_yng, permutations = 999)
cat(sprintf("PERMANOVA Group: R2 = %.3f, p = %.3g\n",
            adonis_yng$R2[1], adonis_yng$`Pr(>F)`[1]))



# ---- 3.2 Distance to the donors, by group ----------------------------------
bm_yng  <- as.matrix(bray_yng)
donors  <- rownames(md_yng)[md_yng$Group == "Donor (healthy)"]
d2d_yng <- rowMeans(bm_yng[, donors])

cat("Mean Bray-Curtis distance to the healthy donors:\n")
for (g in levels(md_yng$Group))
  cat(sprintf("   %-16s %.3f\n", g, mean(d2d_yng[md_yng$Group == g])))
cat(sprintf("   CDI pre-FMT vs Post-FMT: p = %.3g\n",
            suppressWarnings(wilcox.test(
              d2d_yng[md_yng$Group == "CDI pre-FMT"],
              d2d_yng[md_yng$Group == "Post-FMT"])$p.value)))



# ---- 3.3 Alpha diversity ---------------------------------------------------
counts_yng  <- t(otu_yng)
deep_enough <- rowSums(counts_yng) >= RAREFY_DEPTH
set.seed(SEED)
rare_yng <- rrarefy(counts_yng[deep_enough, ], RAREFY_DEPTH)

alpha_yng <- data.frame(
  Observed = specnumber(rare_yng),
  Shannon  = diversity(rare_yng, index = "shannon"),
  Group    = md_yng$Group[deep_enough]
)

cat("Mean Shannon diversity by group (rarefied, n =", nrow(alpha_yng), "):\n")
for (g in levels(alpha_yng$Group))
  cat(sprintf("   %-16s %.2f\n", g, mean(alpha_yng$Shannon[alpha_yng$Group == g])))

p_alpha_yng <- ggplot(alpha_yng, aes(Group, Shannon, fill = Group)) +
  geom_boxplot(outlier.size = 0.6) +
  labs(x = NULL, title = "Diversity recovers after FMT") +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "youngster_alpha_diversity.png"), p_alpha_yng,
       width = 5.5, height = 4, dpi = 150)



# ---- 3.4 Is engraftment donor-specific? ------------------------------------
# Do patients end up resembling THEIR OWN donor, or just healthy people
# in general?
#
# The donor sample IDs in the metadata (FMT1, FMT2) are written without the
# zero padding or the "A" suffix used in the OTU table (FMT01, FMT1A), so they
# have to be matched up by hand.

own_donor <- vapply(md_yng[post_yng, "donor_sample"], function(id) {
  id <- trimws(id)
  if (!nzchar(id)) return(NA_character_)
  if (id %in% donors) return(id)
  n <- sub("^FMT0*", "", id)
  cand <- donors[donors %in% c(paste0("FMT", n),  paste0("FMT0", n),
                               paste0("FMT", n, "A"), paste0("FMT0", n, "A"))]
  if (length(cand)) cand[1] else NA_character_
}, character(1))

ds <- data.frame(sample = post_yng, own_donor = own_donor) %>% filter(!is.na(own_donor))
ds$d_own   <- bm_yng[cbind(ds$sample, ds$own_donor)]
ds$d_other <- vapply(seq_len(nrow(ds)), function(i)
  mean(bm_yng[ds$sample[i], setdiff(donors, ds$own_donor[i])]), numeric(1))

observed_gap <- mean(ds$d_other - ds$d_own)   # positive = own donor is closer

# A paired test alone is not enough. Most recipients share only a handful of
# donors, so a donor that happens to sit centrally in the data would look
# "closer" to everyone by construction. Shuffling which donor is assigned to
# which recipient gives the null distribution this needs.
#
# With k donors, the gap for recipient i assigned donor j is
#     (S_i - D_ij)/(k-1) - D_ij  =  S_i/(k-1) - D_ij * k/(k-1)
# where D is the recipient x donor distance matrix and S_i its row sum.

D <- bm_yng[ds$sample, donors, drop = FALSE]
k <- length(donors)
S <- rowSums(D)
gap_for <- function(j) mean(S / (k - 1) - D[cbind(seq_len(nrow(D)), j)] * k / (k - 1))

set.seed(SEED)
null_gaps <- replicate(999, gap_for(sample(match(ds$own_donor, donors))))
p_perm    <- (sum(null_gaps >= observed_gap) + 1) / 1000

cat("\n--- DONOR SPECIFICITY ---\n")
cat("Recipients matched to their own donor sample: ", nrow(ds), " of ",
    length(post_yng), "\n", sep = "")
cat("Distinct donors actually used:                ", length(unique(ds$own_donor)), "\n", sep = "")
cat(sprintf("Mean distance to OWN donor:    %.3f\n", mean(ds$d_own)))
cat(sprintf("Mean distance to OTHER donors: %.3f\n", mean(ds$d_other)))
cat(sprintf("Paired Wilcoxon:               p = %.3g\n",
            suppressWarnings(wilcox.test(ds$d_own, ds$d_other, paired = TRUE)$p.value)))
cat(sprintf("Donor-label permutation test:  observed gap = %+.3f, null mean = %+.3f, p = %.3g\n\n",
            observed_gap, mean(null_gaps), p_perm))

write.csv(ds, file.path(OUT_DIR, "youngster_donor_specificity.csv"), row.names = FALSE)



# =============================================================================
#  DONE
# =============================================================================

cat("Figures and tables written to:", normalizePath(OUT_DIR), "\n")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

# sink()   # uncomment if you enabled the sink() at the top

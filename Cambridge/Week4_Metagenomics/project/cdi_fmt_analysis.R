# =============================================================================
#  Week 4 mini project - 16S amplicon analysis
#
#  QUESTION
#    Which gut bacteria are lost when a patient has a C. difficile infection,
#    and does a faecal transplant bring them back?
#
#  DATA - one study, three groups
#    Duvallet et al. 2017, Nat Commun 8:1784   (s41467-017-01973-8)
#    cdi_youngster (Youngster et al. 2014, Clin Infect Dis 58:1515)
#
#    100 samples, but only 23 PEOPLE - patients were sampled repeatedly:
#        CDI pre-FMT      27 samples from 19 patients
#        Post-FMT         55 samples from the same 19 patients
#        Donor (healthy)  18 samples from just 4 donors
#    The paper says "20 patients" because it counts people. See 1.5 - the two
#    slide analyses average each person first, so nobody is counted twice.
#
#  THE TWO SLIDES
#    PART 2   slide 1 - does FMT restore gut diversity?
#    PART 3   slide 2 - which "key players" are lost, and do they come back?
#    PART 4   supporting analyses - backup slides only
#
#  HOW TO RUN
#    Open in RStudio, point the working directory at the folder holding the
#    cdi_youngster_results directory, then run from top to bottom.
#        setwd("~/Downloads/Week 4 mini project")
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

STUDY        <- "cdi_youngster"
OUT_DIR      <- "results"
MIN_DEPTH    <- 1000      # drop samples with fewer reads than this
RAREFY_DEPTH <- 2000      # even depth, used for diversity only
MIN_PREV     <- 0.10      # a genus must appear in >=10% of samples to be tested
SEED         <- 42

RANKS      <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
PLOT_LEVEL <- "Genus"     # "Phylum" / "Class" / "Order" / "Family" / "Genus"

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ---- One helper -------------------------------------------------------------
# Everything else is written out in sequence. This one is a function because it
# is used three times below.

# Compare every genus between two groups.
# Returns one row per genus: p value, effect size, and BH-adjusted q value.
compare_groups <- function(mat, groups, group1, group2, min_prev = MIN_PREV) {
  i1 <- which(groups == group1)
  i2 <- which(groups == group2)

  # only test genera that actually show up in a reasonable number of samples
  prevalence <- rowMeans(mat[, c(i1, i2), drop = FALSE] > 0)
  mat <- mat[prevalence >= min_prev, , drop = FALSE]

  out <- data.frame(genus = rownames(mat), p = NA_real_, delta = NA_real_)

  for (i in seq_len(nrow(mat))) {
    x <- mat[i, i1]
    y <- mat[i, i2]
    out$p[i] <- suppressWarnings(wilcox.test(x, y)$p.value)
    # Cliff's delta: effect size from -1 to +1. Positive means higher in group1.
    out$delta[i] <- (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) /
                    (length(x) * length(y))
  }

  out$q <- p.adjust(out$p, method = "BH")   # correct for testing many genera
  out[order(out$q), ]
}


# ---- The SILVA reference ----------------------------------------------------
# assignTaxonomy() needs the SILVA training set (~130 MB). It is downloaded
# once, and the classified result is cached, so this only ever happens once.
# If the cache already exists the reference is never opened - which means you
# can run the slow step on one machine, copy the .rds file, and work elsewhere.

SILVA_FILE <- "silva_nr99_v138.1_wSpecies_train_set.fa.gz"
SILVA_URL  <- paste0("https://zenodo.org/record/4587955/files/",
                     "silva_nr99_v138.1_wSpecies_train_set.fa.gz")
CACHE_FILE <- file.path(OUT_DIR, paste0(STUDY, ".silva_taxonomy.rds"))

cache_ready <- file.exists(CACHE_FILE)
silva_ok    <- file.exists(SILVA_FILE) && file.info(SILVA_FILE)$size >= 50e6

if (cache_ready) {
  cat("Taxonomy cache found - the SILVA reference is not needed\n")
} else if (!silva_ok) {
  # a failed download leaves a stub behind; delete it rather than trust it
  if (file.exists(SILVA_FILE)) unlink(SILVA_FILE)
  cat("Downloading the SILVA training set (~130 MB), one time only ...\n")
  try(download.file(SILVA_URL, SILVA_FILE, mode = "wb"), silent = TRUE)
  silva_ok <- file.exists(SILVA_FILE) && file.info(SILVA_FILE)$size >= 50e6
  if (!silva_ok && file.exists(SILVA_FILE)) unlink(SILVA_FILE)
}

if (!cache_ready && !silva_ok)
  stop("SILVA reference not available. Download it from:\n  ", SILVA_URL,
       "\nand put it in: ", normalizePath("."), call. = FALSE)


# =============================================================================
#  PART 1 - LOAD THE DATA
# =============================================================================

cat("\n===== PART 1: loading", STUDY, "=====\n")

# ---- 1.1 Read the count table and the sample information -------------------
# The OTU table has one row per unique sequence and one column per sample.

otu <- as.matrix(read.delim(
  file.path(paste0(STUDY, "_results"), paste0(STUDY, ".otu_table.100.denovo")),
  row.names = 1, check.names = FALSE))

# quote = "" and comment.char = "" stop R choking on apostrophes and # symbols
# inside the free-text metadata columns
md <- read.delim(
  file.path(paste0(STUDY, "_results"), paste0(STUDY, ".metadata.txt")),
  check.names = FALSE, colClasses = "character", quote = "", comment.char = "")
md <- md[!duplicated(md[[1]]), ]
rownames(md) <- md[[1]]

cat("OTU table:", nrow(otu), "OTUs x", ncol(otu), "samples\n")


# ---- 1.2 Give the sequences names (the course assignTaxonomy step) ---------
# The course practical runs:
#     taxa <- assignTaxonomy(seqtab.nochim, "silva_..._train_set.fa.gz")
# Our equivalent input is the representative sequence for each OTU.
#
# tryRC = TRUE is essential. Sequences in these datasets are not always stored
# in the same direction as the reference, and without it almost nothing gets
# named. The check further down will catch it if that ever goes wrong.
#
# The first run takes a long time. The result is cached, so later runs are
# instant.

seqs <- readDNAStringSet(
  file.path(paste0(STUDY, "_results"), paste0(STUDY, ".otu_seqs.100.fasta")))
seqs <- seqs[names(seqs) %in% rownames(otu)]

if (file.exists(CACHE_FILE)) {
  taxa_raw <- readRDS(CACHE_FILE)
  cat("Loaded cached taxonomy\n")
} else {
  cat("Running assignTaxonomy() on", length(seqs), "sequences - this is slow ...\n")
  taxa_raw <- assignTaxonomy(as.character(seqs), SILVA_FILE,
                             multithread = TRUE, tryRC = TRUE, minBoot = 50)
  rownames(taxa_raw) <- names(seqs)
  saveRDS(taxa_raw, CACHE_FILE)
}

taxa  <- taxa_raw[, RANKS]
genus <- taxa[, "Genus"]

cat("Named at each level:",
    paste(sprintf("%s %.0f%%", RANKS, 100 * colMeans(!is.na(taxa))),
          collapse = ", "), "\n")

# If hardly anything got a name, something is wrong - stop now rather than
# produce an almost-empty table 100 lines later.
if (mean(!is.na(genus)) < 0.20)
  stop("Only ", round(100 * mean(!is.na(genus))), "% of OTUs got a genus name. ",
       "Check that tryRC = TRUE.", call. = FALSE)


# ---- 1.3 Line the tables up, and label the three groups --------------------

keep_otu <- intersect(rownames(otu), rownames(taxa))
otu   <- otu[keep_otu, ]
taxa  <- taxa[keep_otu, , drop = FALSE]
genus <- taxa[, "Genus"]

# Some columns of the OTU table are laboratory controls with no metadata.
# Matching on sample name removes them, and keeps both tables in the same order.
keep_samples <- intersect(colnames(otu), rownames(md))
cat("Dropped", ncol(otu) - length(keep_samples), "samples with no metadata\n")
otu <- otu[, keep_samples]
md  <- md[keep_samples, ]

deep <- colSums(otu) >= MIN_DEPTH
cat("Dropped", sum(!deep), "samples with fewer than", MIN_DEPTH, "reads\n")
otu <- otu[, deep]
md  <- md[deep, ]

md$Group <- factor(
  recode(md$DiseaseState,
         "CDI"         = "CDI pre-FMT",
         "postFMT_CDI" = "Post-FMT",
         "H"           = "Donor (healthy)"),
  levels = c("CDI pre-FMT", "Post-FMT", "Donor (healthy)"))

print(table(md$Group))


# ---- 1.4 Collapse to genus level, and build the phyloseq object ------------
# ~70,000 individual sequences is far too fine-grained to interpret, and a
# single species is split across hundreds of them. Adding up every sequence
# that shares a genus puts each organism back together.
# (Same as tax_glom(ps, "Genus"), but much faster.)

named    <- !is.na(genus)
gen_tab  <- rowsum(otu[named, ], group = genus[named])

cat("Genus level:", nrow(gen_tab), "genera;",
    sprintf("%.1f%%", 100 * sum(gen_tab) / sum(otu)), "of reads kept\n")

# keep the full lineage for each genus, so plots can use any rank
first_row <- taxa[named, , drop = FALSE][!duplicated(genus[named]), , drop = FALSE]
rownames(first_row) <- genus[named][!duplicated(genus[named])]
lineage <- first_row[rownames(gen_tab), , drop = FALSE]

ps  <- phyloseq(otu_table(gen_tab, taxa_are_rows = TRUE),
                sample_data(md),
                tax_table(as.matrix(lineage)))
rel <- transform_sample_counts(ps, function(x) x / sum(x))   # relative abundance

rel_mat <- as(otu_table(rel), "matrix")
groups  <- as.character(md$Group)

ps


# ---- 1.5 Composition overview - the course practical's bar plot ------------
# This is the figure from "13-16S sequencing", adapted from Facility to Group:
#     top20    <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
#     ps.top20 <- prune_taxa(top20, transform_sample_counts(ps, function(x) x/sum(x)))
#     plot_bar(ps.top20, fill = "Class") + facet_grid(~Facility, ...)
#
# It is purely descriptive - no test, no filter that affects any statistic.
# The top 20 selection is only about what fits legibly in a legend.

#  plot_bar()'s defaults do not survive 20 categories - the palette repeats
#  hues and it outlines every segment in black, which turns the bars into
#  static. So we melt the object and draw it ourselves. Same data, same idea,
#  three changes: no segment borders, a distinguishable palette, and samples
#  sorted within each group so the pattern reads left to right.

# 15 rather than the practical's 20: at this figure width a 21-entry legend
# runs off the canvas and gets clipped. 15 still covers 59% of the average
# sample (20 covers 68%), and the rest is pooled into the grey "Other" band.
top_n   <- 15
top_gen <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:top_n]
ps_top  <- prune_taxa(top_gen, rel)

dat <- psmelt(ps_top)

# everything outside the top genera becomes one grey band, so bars reach 100%
other <- data.frame(Sample = sample_names(rel),
                    Abundance = 1 - colSums(as(otu_table(ps_top), "matrix")),
                    Genus = paste0("Other (not in top ", top_n, ")"),
                    Group = get_variable(rel, "Group"))
dat <- rbind(dat[, c("Sample", "Abundance", "Genus", "Group")], other)

# order genera by overall abundance, with Other last; order samples by their
# most abundant genus so similar samples sit together
lev <- c(names(sort(taxa_sums(ps_top), decreasing = TRUE)), paste0("Other (not in top ", top_n, ")"))
dat$Genus <- factor(dat$Genus, levels = rev(lev))
top_genus <- lev[1]
ord <- dat[dat$Genus == top_genus, ]
dat$Sample <- factor(dat$Sample, levels = ord$Sample[order(ord$Group, -ord$Abundance)])

pal <- c(setNames(colorRampPalette(
           c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
             "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#e377c2",
             "#7f7f7f","#bcbd22","#17becf"))(top_n), lev[1:top_n]),
         setNames("grey78", paste0("Other (not in top ", top_n, ")")))

# No plot title - the slide carries it. Legend along the bottom so the panels
# get the full width, and a wide/short aspect so the figure fills a slide.
p_comp <- ggplot(dat, aes(Sample, Abundance, fill = Genus)) +
  geom_col(width = 1, colour = NA) +
  facet_grid(~Group, scales = "free_x", space = "free") +
  scale_fill_manual(values = pal, breaks = lev) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        legend.position = "bottom", legend.key.size = unit(0.38, "cm"),
        legend.text = element_text(size = 10),
        legend.margin = margin(t = -2)) +
  labs(x = NULL, y = "Relative abundance", fill = NULL)

# Sized to the slide slot (9.2 x 3.1 in) so the fonts above are the sizes the
# audience actually sees. Saving bigger and letting PowerPoint shrink it is
# what made the text small.
ggsave(file.path(OUT_DIR, "composition_Genus.png"), p_comp,
       width = 9.6, height = 3.35, dpi = 200)
cat("Composition bar plot written (top 20 genera + Other)\n")


# ---- 1.6 One row per PERSON, not per sample --------------------------------
#
#  IMPORTANT. The 100 samples do NOT come from 100 people. They come from 23:
#      19 patients  (a pre-FMT sample, plus up to 11 post-FMT samples each)
#       4 donors    (2-8 stool samples each)
#
#  The original paper (Youngster et al. 2014, Clin Infect Dis) reports 20
#  patients, and Duvallet's table lists "CDI 19, H 4" - they are counting
#  PEOPLE. Our 100 is counting SAMPLES.
#
#  Treating repeated samples from one person as independent (pseudoreplication)
#  makes p values look far better than they should - the donor group in
#  particular is really n = 4, not n = 18. So for the two slides we average
#  each person's samples within each group, and test people against people.
#  Part 4 keeps the per-sample data, because the time course needs it.

person <- paste(md$subject, md$Group, sep = "|")

# average that person's samples, one column per person-and-group
rel_person <- sapply(unique(person),
                     function(p) rowMeans(rel_mat[, person == p, drop = FALSE]))
groups_person <- sub(".*\\|", "", colnames(rel_person))

cat("\nSamples:", ncol(rel_mat), " -> people:", ncol(rel_person), "\n")
for (g in levels(md$Group))
  cat(sprintf("   %-16s %2d samples from %2d people\n", g,
              sum(groups == g), sum(groups_person == g)))


# =============================================================================
#  PART 2 - SLIDE 1: does FMT restore gut diversity?
# =============================================================================
#  Shannon diversity is one number per sample: how many different bacteria are
#  present, and how evenly. Low = a few species dominate. High = a rich, mixed
#  community.
#
#  Deeper sequencing finds more bacteria, so every sample is first subsampled
#  to the same number of reads. Then, because 19 patients were each sampled
#  before AND after their transplant, we can use a PAIRED test - every patient
#  is their own control, which is the strongest comparison available here.

cat("\n===== PART 2: diversity =====\n")

counts <- t(otu)                                  # samples x OTUs
deep_enough <- rowSums(counts) >= RAREFY_DEPTH
set.seed(SEED)
rarefied <- rrarefy(counts[deep_enough, ], RAREFY_DEPTH)

alpha <- data.frame(
  sample  = rownames(rarefied),
  Shannon = diversity(rarefied, index = "shannon"),
  Group   = md$Group[deep_enough],
  subject = md$subject[deep_enough]
)

cat("Mean Shannon diversity (", nrow(alpha), "samples ):\n")
for (g in levels(alpha$Group))
  cat(sprintf("   %-16s %.2f\n", g, mean(alpha$Shannon[alpha$Group == g])))

# ---- 2.1 The paired before-and-after comparison ----------------------------
# One value per patient per group, then match each patient's before to their
# own after.

per_patient <- aggregate(Shannon ~ subject + Group, data = alpha, FUN = mean)
before <- per_patient[per_patient$Group == "CDI pre-FMT", ]
after  <- per_patient[per_patient$Group == "Post-FMT", ]
paired <- merge(before, after, by = "subject", suffixes = c("_before", "_after"))

improved  <- sum(paired$Shannon_after > paired$Shannon_before)
paired_p  <- suppressWarnings(
  wilcox.test(paired$Shannon_before, paired$Shannon_after, paired = TRUE)$p.value)

cat("\n--- SLIDE 1 HEADLINE ---\n")
cat("Patients with a before AND an after sample:", nrow(paired), "\n")
cat(sprintf("Mean Shannon before FMT: %.2f\n", mean(paired$Shannon_before)))
cat(sprintf("Mean Shannon after FMT : %.2f\n", mean(paired$Shannon_after)))
cat("Diversity improved in", improved, "of", nrow(paired), "patients\n")
cat(sprintf("Paired Wilcoxon test: p = %.2g\n\n", paired_p))

# is the post-FMT level now the same as the donors? (large p = yes)
donor_div <- per_patient$Shannon[per_patient$Group == "Donor (healthy)"]
cat(sprintf("Post-FMT vs donors: %.2f vs %.2f, p = %.2g (large p = indistinguishable)\n",
            mean(after$Shannon), mean(donor_div),
            suppressWarnings(wilcox.test(after$Shannon, donor_div)$p.value)))

write.csv(per_patient, file.path(OUT_DIR, "diversity_per_patient.csv"), row.names = FALSE)

# ---- 2.2 The slide 1 figure: one line per patient --------------------------
# Each grey line is one patient, before -> after. Almost all of them go up.

paired_long <- data.frame(
  subject = rep(paired$subject, 2),
  when    = factor(rep(c("Before FMT", "After FMT"), each = nrow(paired)),
                   levels = c("Before FMT", "After FMT")),
  Shannon = c(paired$Shannon_before, paired$Shannon_after)
)

p_paired <- ggplot(paired_long, aes(when, Shannon)) +
  geom_line(aes(group = subject), colour = "grey70") +
  geom_point(aes(colour = when), size = 2.5) +
  geom_hline(yintercept = mean(donor_div), linetype = 2) +
  annotate("text", x = 0.7, y = mean(donor_div), vjust = -0.6,
           label = "healthy donors", size = 3.5) +
  labs(x = NULL, y = "Shannon diversity", colour = NULL,
       title = sprintf("Diversity rose in %d of %d patients after FMT",
                       improved, nrow(paired))) +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "slide1_diversity_paired.png"), p_paired,
       width = 5.5, height = 5, dpi = 150)

# ---- 2.3 The same thing as a boxplot, if you prefer it ---------------------

p_box <- ggplot(alpha, aes(Group, Shannon, fill = Group)) +
  geom_boxplot(outlier.size = 0.8) +
  labs(x = NULL, y = "Shannon diversity",
       title = "Gut diversity: before FMT, after FMT, and in healthy donors") +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "slide1_diversity_boxplot.png"), p_box,
       width = 5.5, height = 4, dpi = 150)


# =============================================================================
#  PART 3a - SLIDE 2: which key players are lost in CDI?
# =============================================================================
#  Compare patients before their transplant against the healthy donors, one
#  genus at a time.

cat("\n===== PART 3a: which genera differ in CDI? =====\n")

cdi_vs_donor <- compare_groups(rel_person, groups_person,
                               "CDI pre-FMT", "Donor (healthy)")

# delta is positive when a genus is HIGHER in CDI patients
key_players <- cdi_vs_donor[cdi_vs_donor$q < 0.05, ]
key_players$status <- ifelse(key_players$delta < 0, "lost in CDI", "blooms in CDI")
key_players <- key_players[order(key_players$delta), ]   # most depleted first

cat("Genera tested:               ", nrow(cdi_vs_donor), "\n")
cat("Significantly different:     ", nrow(key_players), "\n")
cat("  lost in CDI:               ", sum(key_players$status == "lost in CDI"), "\n")
cat("  bloom in CDI:              ", sum(key_players$status == "blooms in CDI"), "\n\n")

cat("Most depleted:", paste(head(key_players$genus, 8), collapse = ", "), "\n")
cat("Most bloomed :", paste(head(rev(key_players$genus), 5), collapse = ", "), "\n")

write.csv(key_players, file.path(OUT_DIR, "key_players_CDI_vs_donor.csv"),
          row.names = FALSE)

# ---- the slide 1 figure ----------------------------------------------------
# The 10 most depleted and 10 most bloomed genera, drawn as effect sizes.

slide1_dat <- rbind(head(key_players, 10), tail(key_players, 10))
slide1_dat$genus <- factor(slide1_dat$genus, levels = slide1_dat$genus)

p_slide1 <- ggplot(slide1_dat, aes(delta, genus, fill = status)) +
  geom_col() +
  geom_vline(xintercept = 0) +
  labs(x = "Effect size  (negative = lost in CDI, positive = blooms in CDI)",
       y = NULL, fill = NULL,
       title = "Which gut bacteria change during C. difficile infection")
ggsave(file.path(OUT_DIR, "slide1_key_players.png"), p_slide1,
       width = 8, height = 6, dpi = 150)


# =============================================================================
#  PART 3b - SLIDE 2: do the key players come back after FMT?
# =============================================================================
#  Take the same genera from Part 2 and ask, for each one:
#      after the transplant, does it return to the level seen in the donors?

cat("\n===== PART 3b: recovery after FMT =====\n")

recovery <- data.frame(
  genus  = key_players$genus,
  status = key_players$status,
  cdi    = NA_real_,   # mean abundance before FMT
  post   = NA_real_,   # mean abundance after FMT
  donor  = NA_real_,   # mean abundance in the donors
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(recovery))) {
  abundance <- rel_person[recovery$genus[i], ]
  recovery$cdi[i]   <- mean(abundance[groups_person == "CDI pre-FMT"])
  recovery$post[i]  <- mean(abundance[groups_person == "Post-FMT"])
  recovery$donor[i] <- mean(abundance[groups_person == "Donor (healthy)"])
}

# How far back towards the donor level did it get?
#     0 = did not move,  1 = exactly the donor level,  >1 = overshot
# CAUTION: this is a ratio, so it becomes meaningless when the CDI and donor
# levels are already almost the same. Only compute it when there is a real gap.
MIN_GAP <- 0.001
gap <- recovery$donor - recovery$cdi
recovery$recovered <- ifelse(abs(gap) < MIN_GAP, NA,
                             (recovery$post - recovery$cdi) / gap)

# Two simple questions per genus:
#   did it move?          post-FMT vs pre-FMT   (small q = yes, it moved)
#   is it back to normal? post-FMT vs donors    (LARGE q = looks like a donor)
moved     <- compare_groups(rel_person, groups_person, "Post-FMT", "CDI pre-FMT")
vs_donor  <- compare_groups(rel_person, groups_person, "Post-FMT", "Donor (healthy)")

recovery$q_moved    <- moved$q[match(recovery$genus, moved$genus)]
recovery$q_vs_donor <- vs_donor$q[match(recovery$genus, vs_donor$genus)]

# BE HONEST ABOUT THIS. "No longer significantly different from the donors" is
# not the same as "identical to the donors" - with only 18 donor samples the
# test has limited power, so some genera pass simply because there is not
# enough data to show a difference. Always read `restored` next to `recovered`:
# a genus flagged restored that only closed 20% of the gap is better described
# as "moved in the right direction".
recovery$restored <- recovery$q_moved < 0.05 & recovery$q_vs_donor > 0.05

cat("Key players tracked:                    ", nrow(recovery), "\n")
cat("Moved significantly after FMT:          ", sum(recovery$q_moved < 0.05, na.rm = TRUE), "\n")
cat("Restored (moved, and now match donors): ", sum(recovery$restored, na.rm = TRUE), "\n")
cat("Still differ from donors after FMT:     ", sum(recovery$q_vs_donor < 0.05, na.rm = TRUE), "\n\n")

print(head(recovery[, c("genus", "status", "cdi", "post", "donor",
                        "recovered", "restored")], 12), digits = 2)

write.csv(recovery, file.path(OUT_DIR, "key_player_recovery.csv"), row.names = FALSE)

# ---- the slide 2 figure ----------------------------------------------------
# One row per genus, three dots: before FMT, after FMT, and the donor level.
# If the transplant worked, the middle dot sits close to the donor dot.

show   <- head(recovery, 12)     # the 12 most depleted key players
PSEUDO <- 1e-5                   # so zeros can be shown on a log axis

slide2_dat <- data.frame(
  genus     = rep(show$genus, 3),
  group     = rep(c("CDI pre-FMT", "Post-FMT", "Donor (healthy)"), each = nrow(show)),
  abundance = c(show$cdi, show$post, show$donor)
)
slide2_dat$group <- factor(slide2_dat$group,
                           levels = c("CDI pre-FMT", "Post-FMT", "Donor (healthy)"))
slide2_dat$genus <- factor(slide2_dat$genus, levels = rev(show$genus))

p_slide2 <- ggplot(slide2_dat, aes(abundance + PSEUDO, genus, colour = group)) +
  geom_line(aes(group = genus), colour = "grey75") +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(x = "Mean relative abundance (log scale)", y = NULL, colour = NULL,
       title = "Key players lost in CDI, and where they sit after FMT")
ggsave(file.path(OUT_DIR, "slide2_recovery.png"), p_slide2,
       width = 8, height = 6, dpi = 150)


# =============================================================================
#  PART 4 - SUPPORTING ANALYSES
#
#  None of this belongs on your two slides. Keep it for backup slides and for
#  answering questions.
# =============================================================================

cat("\n===== PART 4: supporting analyses =====\n")

# ---- 4.1 Richness (how many genera, ignoring evenness) ---------------------
# Shannon is already covered in Part 2. This is the simpler count-based view.

alpha$Observed <- specnumber(rarefied)
cat("Mean number of OTUs observed:\n")
for (g in levels(alpha$Group))
  cat(sprintf("   %-16s %.0f\n", g, mean(alpha$Observed[alpha$Group == g])))

p_rich <- ggplot(alpha, aes(Group, Observed, fill = Group)) +
  geom_boxplot() +
  labs(x = NULL, y = "OTUs observed", title = "Richness by group") +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "supp_richness.png"), p_rich,
       width = 5.5, height = 4, dpi = 150)


# ---- 4.2 Whole-community comparison ----------------------------------------
# Bray-Curtis measures how different two samples are overall. NMDS squashes
# that down to two dimensions so it can be plotted.

bray <- phyloseq::distance(rel, method = "bray")
ord  <- ordinate(rel, method = "NMDS", distance = "bray", trymax = 20, trace = 0)
cat(sprintf("NMDS stress = %.3f%s\n", ord$stress,
            if (ord$stress > 0.2) "  <- POOR (>0.2), read the plot with caution" else ""))

p_ord <- plot_ordination(rel, ord, color = "Group") +
  geom_point(size = 2.5) + stat_ellipse(level = 0.68) +
  labs(title = "FMT moves patients toward the donor range")
ggsave(file.path(OUT_DIR, "supp_ordination.png"), p_ord,
       width = 6.5, height = 5, dpi = 150)

perm <- adonis2(bray ~ Group, data = md, permutations = 999)
cat(sprintf("PERMANOVA Group: R2 = %.3f, p = %.3g\n", perm$R2[1], perm$`Pr(>F)`[1]))

# distance from every sample to the healthy donors
bray_mat <- as.matrix(bray)
donors   <- rownames(md)[md$Group == "Donor (healthy)"]
dist_to_donor <- rowMeans(bray_mat[, donors])

cat("Mean distance to the donors:\n")
for (g in levels(md$Group))
  cat(sprintf("   %-16s %.3f\n", g, mean(dist_to_donor[md$Group == g])))


# ---- 4.3 How quickly does recovery happen? ---------------------------------

post_samples <- rownames(md)[md$Group == "Post-FMT"]
traj <- data.frame(
  days = suppressWarnings(as.numeric(md[post_samples, "days_since_fmt"])),
  dist = dist_to_donor[post_samples]
)
traj <- traj[!is.na(traj$days), ]

time_test <- suppressWarnings(cor.test(traj$days, traj$dist, method = "spearman"))
cat(sprintf("Distance to donors vs days since FMT: rho = %.2f, p = %.3g (n = %d)\n",
            time_test$estimate, time_test$p.value, nrow(traj)))

p_traj <- ggplot(traj, aes(days, dist)) +
  geom_point(size = 2) + geom_smooth(method = "loess", se = TRUE, span = 1) +
  scale_x_log10() +
  labs(x = "Days since FMT (log scale)", y = "Distance to the donors",
       title = "Recovery happens quickly, then levels off")
ggsave(file.path(OUT_DIR, "supp_trajectory.png"), p_traj,
       width = 6, height = 4.5, dpi = 150)


# ---- 4.4 Do patients end up resembling their OWN donor? --------------------
# Donor IDs in the metadata (FMT1, FMT2) are written without the zero padding
# or "A" suffix used in the OTU table (FMT01, FMT1A), so they must be matched
# up by hand.

own_donor <- vapply(md[post_samples, "donor_sample"], function(id) {
  id <- trimws(id)
  if (!nzchar(id)) return(NA_character_)
  if (id %in% donors) return(id)
  n <- sub("^FMT0*", "", id)
  hit <- donors[donors %in% c(paste0("FMT", n),  paste0("FMT0", n),
                              paste0("FMT", n, "A"), paste0("FMT0", n, "A"))]
  if (length(hit)) hit[1] else NA_character_
}, character(1))

pairs <- data.frame(sample = post_samples, own_donor = own_donor)
pairs <- pairs[!is.na(pairs$own_donor), ]
pairs$d_own   <- bray_mat[cbind(pairs$sample, pairs$own_donor)]
pairs$d_other <- vapply(seq_len(nrow(pairs)), function(i)
  mean(bray_mat[pairs$sample[i], setdiff(donors, pairs$own_donor[i])]), numeric(1))

observed_gap <- mean(pairs$d_other - pairs$d_own)  # positive = own donor closer

# A paired test alone is not enough: most recipients share only a few donors,
# so a donor sitting centrally in the data would look "closer" to everyone by
# construction. Shuffling the donor labels gives the null this needs.
D <- bray_mat[pairs$sample, donors, drop = FALSE]
k <- length(donors)
S <- rowSums(D)
gap_for <- function(j) mean(S / (k - 1) - D[cbind(seq_len(nrow(D)), j)] * k / (k - 1))

set.seed(SEED)
null_gaps <- replicate(999, gap_for(sample(match(pairs$own_donor, donors))))
p_perm    <- (sum(null_gaps >= observed_gap) + 1) / 1000

cat("\nDonor specificity:\n")
cat(sprintf("   distance to OWN donor:    %.3f\n", mean(pairs$d_own)))
cat(sprintf("   distance to OTHER donors: %.3f\n", mean(pairs$d_other)))
cat(sprintf("   permutation test: observed gap %+.3f vs null %+.3f, p = %.3g\n",
            observed_gap, mean(null_gaps), p_perm))

write.csv(pairs, file.path(OUT_DIR, "supp_donor_specificity.csv"), row.names = FALSE)


# =============================================================================
#  PART 5 - OPTIONAL: DESeq2 cross-check
# =============================================================================
#  Nothing here is needed for the slides. It re-asks the Part 3a question -
#  which genera differ between patients and donors - using a completely
#  different statistical model, so you can say the answer does not depend on
#  the method you happened to pick.
#
#  Install once:   if (!require("BiocManager")) install.packages("BiocManager")
#                  BiocManager::install("DESeq2")

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  cat("\nDESeq2 not installed - skipping the cross-check.\n",
      '  BiocManager::install("DESeq2")\n')
} else {
library(DESeq2)
cat("\n===== PART 5: DESeq2 cross-check =====\n")

# ---- 5.1 One column per PERSON, counts SUMMED ------------------------------
# DESeq2 models whole-number counts, so each person's samples are summed rather
# than averaged. That makes people with many samples look more deeply
# sequenced - which is precisely what size factors exist to correct.

person_counts <- t(rowsum(t(gen_tab), group = person))
person_group  <- sub(".*\\|", "", colnames(person_counts))

sel <- person_group %in% c("CDI pre-FMT", "Donor (healthy)")
cts <- person_counts[, sel, drop = FALSE]
coldata <- data.frame(
  Group = factor(person_group[sel], levels = c("Donor (healthy)", "CDI pre-FMT")),
  row.names = colnames(cts))

cat("Input:", nrow(cts), "genera x", ncol(cts), "people (",
    sum(coldata$Group == "CDI pre-FMT"), "CDI vs",
    sum(coldata$Group == "Donor (healthy)"), "donors )\n")

# same prevalence filter as the main analysis, so the two are comparable
cts <- cts[rowMeans(cts > 0) >= MIN_PREV, , drop = FALSE]
cat("After the", 100 * MIN_PREV, "% prevalence filter:", nrow(cts), "genera\n")

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Group)

# ---- 5.2 Size factors - the microbiome gotcha ------------------------------
# DESeq2's default size factors use a geometric mean ACROSS samples, which
# requires at least one feature present in EVERY sample. In this table only
# 1 genus out of 169 qualifies, because 65% of the entries are zero. The
# default call therefore either errors ("every gene contains at least one
# zero") or silently rests on that single genus.
#
# type = "poscounts" computes the geometric mean over the non-zero entries
# only. This is the standard fix for sparse microbiome counts - do not remove it.

dds <- estimateSizeFactors(dds, type = "poscounts")

# fitType = "local": the default parametric dispersion fit frequently fails to
# converge on sparse microbiome data. "local" is the usual fallback.
dds <- DESeq(dds, fitType = "local", quiet = TRUE)

# ---- 5.3 Results, with fold-change shrinkage -------------------------------
#
# IMPORTANT. Since DESeq2 v1.16, DESeq() does NOT shrink fold changes - you get
# raw maximum-likelihood estimates. On sparse data those are wild: a genus seen
# in 6 patients and 0 donors came out at log2FC = +24, i.e. "17 million fold",
# with a standard error of 3.9. That is noise wearing a decimal point.
#
# lfcShrink() pulls poorly-supported estimates toward zero while leaving
# well-supported ones alone. Always plot the shrunken values.
# apeglm is the recommended estimator; ashr and "normal" are fallbacks.

shrink_type <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else
               if (requireNamespace("ashr",   quietly = TRUE)) "ashr"   else "normal"
cat("Shrinking fold changes with type =", shrink_type,
    if (shrink_type == "normal")
      '  (install apeglm for the better estimator: BiocManager::install("apeglm"))'
    else "", "\n")

res <- as.data.frame(lfcShrink(dds, coef = 2, type = shrink_type, quiet = TRUE))
res$genus <- rownames(res)
res <- res[!is.na(res$padj), ]

# ---- Flag genera that are absent from an entire group ----------------------
# A genus with zero counts in all 4 donors has no meaningful fold change - you
# are dividing by a pseudocount, not by data. Shrinkage tames these but cannot
# make them interpretable, so they are excluded from the bar plot and marked
# in the CSV.

donor_cols <- coldata$Group == "Donor (healthy)"
cdi_cols   <- coldata$Group == "CDI pre-FMT"
absent <- (rowSums(cts[, donor_cols, drop = FALSE]) == 0) |
          (rowSums(cts[, cdi_cols,   drop = FALSE]) == 0)
res$absent_in_a_group <- absent[res$genus]

cat("Genera absent from one group entirely (fold change not interpretable):",
    sum(res$absent_in_a_group, na.rm = TRUE), "\n")

up <- res[res$padj < 0.05 & res$log2FoldChange > 0, ]
dn <- res[res$padj < 0.05 & res$log2FoldChange < 0, ]

cat("Significant (padj < 0.05):", nrow(up) + nrow(dn),
    "  enriched in CDI:", nrow(up), "  depleted:", nrow(dn), "\n")
cat("  top enriched:",
    paste(head(up$genus[order(-up$log2FoldChange)], 8), collapse = ", "), "\n")
cat("  top depleted:",
    paste(head(dn$genus[order(dn$log2FoldChange)], 8), collapse = ", "), "\n")

write.csv(res[order(res$padj), ],
          file.path(OUT_DIR, "deseq2_CDI_vs_donor.csv"), row.names = FALSE)

# ---- 5.4 Plots -------------------------------------------------------------

res$direction <- ifelse(res$padj >= 0.05, "Not significant",
                 ifelse(res$log2FoldChange > 0, "Enriched in CDI", "Depleted in CDI"))
res$direction <- factor(res$direction,
  levels = c("Depleted in CDI", "Enriched in CDI", "Not significant"))
PAL <- c("Depleted in CDI" = "#3A7CA5",
         "Enriched in CDI" = "#E4572E",
         "Not significant" = "grey78")

# (a) VOLCANO - every genus at once.
#     x = how big the difference is, y = how confident we are.
#     Far left  = strongly depleted in patients
#     Far right = strongly enriched in patients
#     Above the dashed line = passes the 5% false-discovery cutoff
res$negLogP <- -log10(pmax(res$padj, 1e-300))
lab <- res[res$padj < 0.05 & !res$absent_in_a_group, ]
lab <- lab[order(lab$padj), ][seq_len(min(10, nrow(lab))), ]

# hollow points = absent from one group entirely, so the x position is not
# a real fold change. Shown, but visually set apart.
p_volcano <- ggplot(res, aes(log2FoldChange, negLogP, colour = direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey55") +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_point(aes(shape = absent_in_a_group), size = 2.2, alpha = 0.85) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c(`FALSE` = "present in both groups",
                                `TRUE`  = "absent from one group"), name = NULL) +
  { if (nrow(lab)) geom_text(data = lab, aes(label = genus), size = 2.9,
                             vjust = -0.9, show.legend = FALSE, check_overlap = TRUE) } +
  scale_colour_manual(values = PAL) +
  labs(x = "log2 fold change   (CDI patients vs healthy donors)",
       y = "-log10 adjusted p", colour = NULL,
       title = "DESeq2: which genera differ between patients and donors") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "deseq2_volcano.png"), p_volcano,
       width = 8, height = 5.5, dpi = 150)

# (b) RANKED BARS - the presentable version.
#     Shrunken fold changes only, and only genera present in BOTH groups, so
#     every bar is a number you can defend. Error bars are +/- 1 standard error.
# Also require the effect to be bigger than its own uncertainty. apeglm shrinks
# the fold change but the p value still comes from the UNSHRUNK Wald test, so a
# genus can be "significant" while its shrunken estimate sits on zero - three do
# here. Requiring |log2FC| > lfcSE removes bars that would otherwise be drawn at
# zero with a confidence interval straddling it.
sig <- res[res$padj < 0.05 & !res$absent_in_a_group &
           abs(res$log2FoldChange) > res$lfcSE, ]
if (nrow(sig)) {
  top <- rbind(head(sig[order(sig$log2FoldChange), ], 10),
               head(sig[order(-sig$log2FoldChange), ], 10))
  top <- top[!duplicated(top$genus), ]
  top$genus <- factor(top$genus, levels = top$genus[order(top$log2FoldChange)])

  p_bars <- ggplot(top, aes(log2FoldChange, genus, fill = direction)) +
    geom_col() +
    geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE),
                   height = 0.3, colour = "grey25", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30") +
    scale_fill_manual(values = PAL) +
    labs(x = "log2 fold change, shrunken   (negative = lost in CDI, positive = blooms in CDI)",
         y = NULL, fill = NULL,
         title = "DESeq2: genera that change most during infection",
         subtitle = paste0("shrunken estimates; ",
                           sum(res$padj < 0.05 & res$absent_in_a_group, na.rm = TRUE),
                           " genera absent from one group entirely are not shown")) +
    theme(legend.position = "bottom")
  ggsave(file.path(OUT_DIR, "deseq2_top_genera.png"), p_bars,
         width = 8, height = 6, dpi = 150)
  cat("  bar plot drawn from", nrow(sig), "interpretable genera\n")
} else {
  cat("  no interpretable genera passed padj < 0.05 - skipping the bar plot\n")
}

# (c) DISPERSION FIT - a DESeq2 diagnostic, not a result.
#     Points should scatter around the fitted red line. If the fit looks wrong,
#     try fitType = "mean" in the DESeq() call above.
png(file.path(OUT_DIR, "deseq2_dispersion.png"), width = 1200, height = 900, res = 150)
plotDispEsts(dds)
dev.off()

# ---- 5.5 Does it agree with the main analysis? -----------------------------
if (exists("key_players")) {
  main_dn <- key_players$genus[key_players$delta < 0]
  main_up <- key_players$genus[key_players$delta > 0]
  cat("Agreement with Part 3a:\n")
  cat("   depleted in both:", length(intersect(dn$genus, main_dn)),
      "of", length(main_dn), "called depleted there\n")
  cat("   enriched in both:", length(intersect(up$genus, main_up)),
      "of", length(main_up), "called enriched there\n")
  cat("   NOTE: DESeq2 normalises by size factors rather than by sample total,\n")
  cat("   so it is expected to find MORE enriched genera than the relative-\n")
  cat("   abundance analysis did. Disagreement on that side is the point.\n")
}
}


# =============================================================================
#  DONE
# =============================================================================

cat("\nFigures and tables written to:", normalizePath(OUT_DIR), "\n")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

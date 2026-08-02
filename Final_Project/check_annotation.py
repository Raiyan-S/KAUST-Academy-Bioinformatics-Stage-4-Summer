import anndata as ad

a = ad.read_h5ad("packer2019.h5ad")
print(a)
print("shape:", a.shape)
print("obs cols:", list(a.obs.columns))
print("var cols:", list(a.var.columns))

# don't assume orientation — find whichever axis carries the cell metadata
if "lineage" in a.var.columns:
    cells = a.var
elif "lineage" in a.obs.columns:
    cells = a.obs
else:
    raise SystemExit("no 'lineage' column on either axis")

print("n cells:", len(cells))

# missing values in these files are often strings, not NaN
lin = cells["lineage"].astype(str).str.strip()
missing = {"nan", "NA", "N/A", "", "None", "unknown"}
usable = ~lin.str.lower().isin({m.lower() for m in missing})

print("lineage usable:", usable.sum(), "of", len(cells))
print(lin[usable].value_counts().head(30))

ct = cells["cell_type"].astype(str).str.strip()
print("cell_type usable:", (~ct.str.lower().isin({m.lower() for m in missing})).sum())

print(cells["passed_qc"].value_counts(dropna=False))




cells = a.obs[a.obs["passed_qc"] == True].copy()
lin = cells["lineage"].astype(str).str.strip()
ok = ~lin.str.lower().isin({"nan", "na", "n/a", "", "none", "unknown"})
sub = cells[ok].copy()
sub["lineage"] = lin[ok].str.split(":").str[0]      # drop ':pseudotime_bin_N' suffixes

sub["n_alt"] = sub["lineage"].str.count("/") + 1
print(sub["n_alt"].value_counts().sort_index())

single = sub[sub["n_alt"] == 1]
print("unambiguous cells:", len(single))
print("distinct lineages:", single["lineage"].nunique())
print("depth distribution:")
print(single["lineage"].str.len().value_counts().sort_index())
print("founders:")
print(single["lineage"].str[:2].value_counts())



def label_at_depth(s, d):
    parts = [p[:d] for p in s.split("/")]
    if len(set(parts)) == 1 and len(parts[0]) == d:
        return parts[0]
    return None

print(f"{'depth':>5} {'cells':>7} {'classes':>8} {'cells/class':>12}")
for d in range(2, 12):
    lab = sub["lineage"].apply(lambda s: label_at_depth(s, d))
    n, k = lab.notna().sum(), lab.nunique()
    print(f"{d:>5} {n:>7} {k:>8} {n/max(k,1):>12.0f}")
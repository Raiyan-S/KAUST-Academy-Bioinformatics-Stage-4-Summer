I didn't keep good notes of what I did, but here's what I can tell from what's in S3.

Looks like I downloaded the fastq's from SRA and split spots using --split-files.

I then used USEARCH with default settings to merge fastqs.
	merged fastqs and merge logs are in merged_fastqs/

Looks like there's a bunch of samples that have IDs that start with MGH.
I'm not sure what these are, and it seems that I only have the metadata for
the samples that start with FMT.

There are also some metagenomics runs (WGS).

I processed only the 16S samples with metadata, i.e. those with sampleIDs that
start with FMT and which are 'AMPLICON' assay type.
These samples are indicated in fq2sid_amplicon_fmt.txt

metadata_Ilan.xlsx is the metadata that Jay got from Ilan Youngster.

I added a DiseaseState column to this metadata. 
This is cdi_youngster.metadata.txt

DUPLICATE SAMPLES
-----------------

Some patients seem to have duplicate samples (within the same pre or post condition).
Samples are also dated (though some samples are taken on the same date for the same
patient...).
For downstream analyses, I want to pick only one sample per patient/case pair. I added
an 'n_samples' column to the metadata to indicate how many samples there are for
each subject-type pair combination.
Downstream analyses should pick just the first (or a random) sample for each pair.

     meta = pd.read_csv('cdi_youngster.metadata.txt', sep='\t', index_col=0)
     meta = meta.sort_values(by='sample_date', ascending=True)
     meta = meta.replace('post ', 'post')
     for grp, subdf in meta.groupby(['type', 'subject']):
     	    meta.loc[subdf.index, 'n_sample'] = range(subdf.shape[0])

I also converted the sample_date serial number column to human-readable dates
using Excel (by manually copy-pasting and changing the data format).

This is all now in cdi_youngster.metadata.txt        

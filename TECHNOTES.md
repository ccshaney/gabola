## GABOLA Technical Notes
#### § Preprocess Module:

#### Step I. Produce non-duplicate split fastqs sorted by barcodes from raw linked reads

**Input:**

Ø   Path to the directory storing raw linked reads (see: https://github.com/10XGenomics/longranger for read file format)
Output:
Ø  Non-duplicate, trimmed and filtered FASTQ files for Read 1 and Read 2 separately (NonDupR1.fq and NonDupR2.fq)
Ø   Directories for every barcode containing their respective reads that are trimmed and filtered (OUTPUT_FOLDER/nonDupFq/split). For example,
Barcode name: AAAACCCCGTGTGTGT
Format: OUTPUT_FOLDER/nonDupFq/split/AAAACCCCGTGTGTGT_1

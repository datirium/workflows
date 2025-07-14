cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-genelists:v3.0.0
inputs:
  threads:
    type: int
    inputBinding:
      prefix: -t
    doc: |
      Number of threads for parallel processing
  genelist_names:
    type: string[]
    inputBinding:
      prefix: -a
      itemSeparator: ','
    doc: |
      Array of genelist aliases/sample names.
  feature_files:
    type: File[]
    inputBinding:
      prefix: -b
      itemSeparator: ','
    doc: |
      Array of TSV files with differential genes from DESeq or diffbind pipelines.
  filtered_files:
    type: File[]
    inputBinding:
      prefix: -c
      itemSeparator: ','
    doc: |
      Array of filtered differential genelists from DESeq or diffbind pipelines.
  sample_names_nabinding:
    type: string[]?
    inputBinding:
      prefix: -d
      itemSeparator: ','
    doc: |
      Array of aliases for ChIP/ATAC/CRT-Seq experiments for row metadata.
  sample_names_rnaseq:
    type: string[]?
    inputBinding:
      prefix: -e
      itemSeparator: ','
    doc: |
      Array of aliases for RNA-Seq experiments for column metadata.
  bam_files:
    type: File[]?
    inputBinding:
      prefix: -f
      itemSeparator: ','
    doc: |
      Array of sample coordinate sorted BAM alignment and BAI index files.
  expression_files:
    type: File[]?
    inputBinding:
      prefix: -g
      itemSeparator: ','
    doc: |
      Array of sample TSV files containing gene annotations with associated TotalReads and Rpkm counts.
outputs:
  master_samplesheet:
    type: File
    outputBinding:
      glob: master_samplesheet.tsv
    doc: |
      contains formatted information of the input data and files
  output_row_metadata:
    type: File
    outputBinding:
      glob: output_row_metadata.tsv
    doc: |
      row metadata for GCT formatter
  output_col_metadata:
    type: File
    outputBinding:
      glob: output_col_metadata.tsv
    doc: |
      column metadata for GCT formatter
  output_counts:
    type: File
    outputBinding:
      glob: output_counts.tsv
    doc: |
      peak average read depth per TSS window and gene expression counts matrix
  heatmap_gct:
    type: File
    outputBinding:
      glob: heatmap.gct
    doc: |
      GCT formatted peak and expression data for morpheus viewer
  heatmap_html:
    type: File
    outputBinding:
      glob: heatmap.html
    doc: |
      HTML preconfigured morpheus heatmap
  heatmap_peaknorm95_html:
    type: File
    outputBinding:
      glob: heatmap_peaknorm95.html
    doc: |
      HTML preconfigured morpheus heatmap scaled to 95th percentile
  heatmap_peaknorm99_html:
    type: File
    outputBinding:
      glob: heatmap_peaknorm99.html
    doc: |
      HTML preconfigured morpheus heatmap scaled to 99th percentile
  log_file_stdout:
    type: stdout
  log_file_stderr:
    type: stderr
baseCommand:
- run_genelists.sh
stdout: genelists_stdout.log
stderr: genelists_stderr.log
doc: "A CWL tool for producing a GCT data file for the morpheus heatmap, and an html heatmap.\nUses both ATAC/ChIP/CRT (NA [nucleic acid] binding) and RNA-Seq data to derive visualization data.\nNA binding data in the form of BAM files per sample is processed to output an average read depth per window +/-5Kbp of each gene's TSS (transcription start site).\nRNA-Seq data in the form of gene expression count matrices are processed to output TotalReads and Rpkm values per gene.\nThese data are then integrated into a single count matrix, a row, and a column metadata file as input to an Rscript that will format the 3 files into GCT format for morpheus heatmap viewer.\nThe HTML heatmap is then produced with preconfigured sorting and grouping settings.\n\n\nPrimary Output files:\n- heatmap.gct, GCT formatted peak and expression data for morpheus viewer\n- heatmap.html, html of morpheus heatmap with preconfigured settings, peak data scaled among all samples\n- heatmap_peaknorm95.html, html of morpheus heatmap with preconfigured settings, peak data scaled per individual sample to 95th percentile\n- heatmap_peaknorm99.html, html of morpheus heatmap with preconfigured settings, peak data scaled per individual sample to 99th percentile\n\nSecondary Output files:\n- master_samplesheet.tsv, contains formatted information of the input data and files\n- output_row_metadata.tsv, row metadata for GCT formatter\n- output_col_metadata.tsv, column metadata for GCT formatter\n- output_counts.tsv, peak average read depth per TSS window and gene expression counts matrix\n\nPARAMS:\n  SECTION 1: general\n  -h\thelp\t\t    show this message\n  -t  INT\t\t\t    number of threads\n  -a\tARRAY\t\t    array of genelist sample names (no commas in names)\n  -b  FILE ARRAY\tarray of associated annotation files for each gene list from (-c), with header\n  -c  FILE ARRAY\tarray of filtered gene list TSVs (must be headerless, columns are: chr, txStart, txEnd, geneID, L2FC, Strand)\n  -d\tARRAY\t\t    array of sample names from NA binding experiments (no commas in names)\n  -e\tARARY\t\t    array of sample names from RNA-Seq experiments (no commas in names)\n  -f\tFILE ARRAY\tarray of BAM files from NA binding experiments\n  -g\tFILE ARARY\tarray of expression table files from RNA-Seq experiments\t\n\n\n____________________________________________________________________________________________________\nReferences:\n- Morpheus, https://software.broadinstitute.org/morpheus\n    "
label: genelists-deseq-diffbind

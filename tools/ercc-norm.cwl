cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-erccnorm:v2.0.0
inputs:
  threads_count:
    type: int
    inputBinding:
      prefix: -t
    doc: |
      Number of threads for parallel processing
  unaligned_R1_fastq:
    type: File
    inputBinding:
      prefix: -u
    doc: |
      unaligned R1 reads post-primary alignment
  unaligned_R2_fastq:
    type: File
    inputBinding:
      prefix: -v
    doc: |
      unaligned R2 reads post-primary alignment
  dilution_factor:
    type: float
    inputBinding:
      prefix: -d
    doc: |
      dilution factor used for ERCC ExFold mix 1 before spike-in
  uL_per_M_cells:
    type: float
    inputBinding:
      prefix: -m
    doc: |
      volume of ERCC ExFold mix 1 spike-in to sample per million cells
  rnaseq_counts:
    type: File
    inputBinding:
      prefix: -c
    doc: |
      csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)
outputs:
  ercc_sam:
    type: File
    outputBinding:
      glob: unaligned_pairs-to-ERCC.sam
    doc: |
      unaligned input reads (against primary reference) aligned to ERCC sequences sam file
  ercc_counts:
    type: File
    outputBinding:
      glob: ercc_counts.tsv
    doc: |
      row metadata for GCT formatter
  ercc_pdf_plot:
    type: File
    outputBinding:
      glob: ercc_expected_v_actual_count_plot.pdf
    doc: |
      ERCC molecules per cell counts (log10) expected vs observed
  rpkm_isoforms_ercc_norm:
    type: File
    outputBinding:
      glob: isoforms.ercc_norm_rpkm.csv
    doc: |
      isoform RPKM counts normalized to ERCC ExFold mix 1 spike-in
  log_file_stdout:
    type: stdout
  log_file_stderr:
    type: stderr
baseCommand:
- run_ercc_norm.sh
stdout: ercc_norm_stdout.log
stderr: ercc_norm_stderr.log
doc: "\nTool for building linear regression function from ERCC ExFold mix 1 RPKM (molecule per cell vs RPKM), and applying this for normalization of RNA-Seq RPKM count data.\n\n  Primary Output files:\n   - unaligned_pairs-to-ERCC.sam\n   - ercc_counts.tsv\n   - isoforms.ercc_norm_rpkm.csv\n   - ercc_expected_v_actual_count_plot.pdf\n\n  PARAMS:\n  -h  help\tshow this message\n  -t  INT\tnumber of threads\n  -u  FILE   array of unaligned \"R1,R2\" reads post-primary alignment\n  -d  FLOAT  dilution factor used for ERCC ExFold mix 1 before spike-in\n  -m  FLOAT  volume of ERCC ExFold mix 1 spike-in to sample per million cells\n  -c  FILE   csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)\n\n____________________________________________________________________________________________________\nReferences:\n- Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.\n- Twelve years of SAMtools and BCFtools. Danecek et al. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008\n    "
label: genelists-deseq-diffbind

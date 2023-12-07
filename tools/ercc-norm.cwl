cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-erccnorm:dev


inputs:

  threads:
    type: int
    inputBinding:
      prefix: "-t"
    doc: |
      Number of threads for parallel processing

  unaligned_fastq_files:
    type: File[]
    inputBinding:
      prefix: "-u"
      itemSeparator: ","
    doc: |
      unaligned R1 reads post-primary alignment

  dilution_factor:
    type: float
    inputBinding:
      prefix: "-d"
    doc: |
      dilution factor used for ERCC ExFold mix 1 before spike-in

  uL_per_M_cells:
    type: float
    inputBinding:
      prefix: "-m"
    doc: |
      volume of ERCC ExFold mix 1 spike-in to sample per million cells

  rnaseq_counts:
    type: File
    inputBinding:
      prefix: "-c"
    doc: |
      csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)


outputs:

  ercc_counts_tsv:
    type: File
    outputBinding:
      glob: ercc_counts.tsv
    doc: |
      row metadata for GCT formatter

  ercc_plot:
    type: File
    outputBinding:
      glob: ercc_expected_v_actual_count_plot.pdf
    doc: |
      expected molecules per cell vs actual ERCC molecule counts (log10)

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


baseCommand: ["run_ercc_norm.sh"]
stdout: ercc_norm_stdout.log
stderr: ercc_norm_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "genelists-deseq-diffbind"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/genelists-deseq-diffbind.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  A CWL tool for producing a GCT data file for the morpheus heatmap, and an html heatmap.
  Uses both ATAC/ChIP/CRT (NA [nucleic acid] binding) and RNA-Seq data to derive visualization data.
  NA binding data in the form of BAM files per sample is processed to output an average read depth per window +/-5Kbp of each gene's TSS (transcription start site).
  RNA-Seq data in the form of gene expression count matrices are processed to output TotalReads and Rpkm values per gene.
  These data are then integrated into a single count matrix, a row, and a column metadata file as input to an Rscript that will format the 3 files into GCT format for morpheus heatmap viewer.
  The HTML heatmap is then produced with preconfigured sorting and grouping settings.


  Primary Output files:
  - heatmap.gct, GCT formatted peak and expression data for morpheus viewer
  - heatmap.html, html of morpheus heatmap with preconfigured settings, peak data scaled among all samples
  - heatmap_peaknorm95.html, html of morpheus heatmap with preconfigured settings, peak data scaled per individual sample to 95th percentile
  - heatmap_peaknorm99.html, html of morpheus heatmap with preconfigured settings, peak data scaled per individual sample to 99th percentile

  Secondary Output files:
  - master_samplesheet.tsv, contains formatted information of the input data and files
  - output_row_metadata.tsv, row metadata for GCT formatter
  - output_col_metadata.tsv, column metadata for GCT formatter
  - output_counts.tsv, peak average read depth per TSS window and gene expression counts matrix

  PARAMS:
    SECTION 1: general
    -h	help		    show this message
    -t  INT			    number of threads
    -a	ARRAY		    array of genelist sample names (no commas in names)
    -b  FILE ARRAY	array of associated annotation files for each gene list from (-c), with header
    -c  FILE ARRAY	array of filtered gene list TSVs (must be headerless, columns are: chr, txStart, txEnd, geneID, L2FC, Strand)
    -d	ARRAY		    array of sample names from NA binding experiments (no commas in names)
    -e	ARARY		    array of sample names from RNA-Seq experiments (no commas in names)
    -f	FILE ARRAY	array of BAM files from NA binding experiments
    -g	FILE ARARY	array of expression table files from RNA-Seq experiments	


  ____________________________________________________________________________________________________
  References:
  - Morpheus, https://software.broadinstitute.org/morpheus
      
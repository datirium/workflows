cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-genelists:dev


inputs:

  threads:
    type: int
    inputBinding:
      prefix: "-t"
    doc: |
      Number of threads for parallel processing

  genelist_names:
    type: string[]
    inputBinding:
      prefix: "-a"
      itemSeparator: ","
    doc: |
      Array of genelist aliases/sample names.

  feature_files:
    type: File[]
    inputBinding:
      prefix: "-b"
      itemSeparator: ","
    doc: |
      Array of TSV files with differential genes from DESeq or diffbind pipelines.

  filtered_files:
    type: File[]
    inputBinding:
      prefix: "-c"
      itemSeparator: ","
    doc: |
      Array of filtered differential genelists from DESeq or diffbind pipelines.

  sample_names_rnaseq:
    type: string[]?
    inputBinding:
      prefix: "-e"
      itemSeparator: ","
    doc: |
      Array of aliases for RNA-Seq experiments for column metadata.

  expression_files:
    type: File[]?
    inputBinding:
      prefix: "-g"
      itemSeparator: ","
    doc: |
      Array of sample TSV files containing gene annotations with associated TotalReads and Rpkm counts.


outputs:

  master_samplesheet_scaled:
    type: File
    outputBinding:
      glob: master_samplesheet_scaled.tsv
    doc: |
      contains formatted information of the input data and files for scaled heatmaps

  master_samplesheet_vst:
    type: File
    outputBinding:
      glob: master_samplesheet_vst.tsv
    doc: |
      contains formatted information of the input data and files for vst normalized heatmap

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
      scaled RPKM (0-99) gene expression counts matrix

  output_counts_vst:
    type: File
    outputBinding:
      glob: output_counts_vst.tsv
    doc: |
      VST normalized gene expression counts matrix

  heatmap_gct:
    type: File
    outputBinding:
      glob: heatmap.gct
    doc: |
      GCT formatted scaled RPKM values (0-99) expression data for morpheus viewer

  heatmap_vst_gct:
    type: File
    outputBinding:
      glob: heatmap_vst.gct
    doc: |
      GCT formatted VST expression data for morpheus viewer

  heatmap_html:
    type: File
    outputBinding:
      glob: heatmap_nonorm.html
    doc: |
      HTML preconfigured morpheus heatmap with scaled RPKM values (0-99)

  heatmap_norm95_html:
    type: File
    outputBinding:
      glob: heatmap_norm95.html
    doc: |
      HTML preconfigured morpheus heatmap with scaled RPKM values (0-99) and 95th percentile max cutoff

  heatmap_norm99_html:
    type: File
    outputBinding:
      glob: heatmap_norm99.html
    doc: |
      HTML preconfigured morpheus heatmap with scaled RPKM values (0-99) and 99th percentile max cutoff

  heatmap_vst_html:
    type: File
    outputBinding:
      glob: heatmap_vst.html
    doc: |
      HTML preconfigured morpheus heatmap with VST values

  log_file_stdout:
    type: stdout

  log_file_stderr:
    type: stderr


baseCommand: ["run_genelists_rnaseq.sh"]
stdout: genelists_stdout.log
stderr: genelists_stderr.log


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
  Uses RNA-Seq data to derive visualization data.
  RNA-Seq data in the form of gene expression count matrices are processed to output TotalReads and Rpkm values per gene.
  The data are then integrated into a single count matrix, a row, and a column metadata file as input to an Rscript that will format the 3 files into GCT format for morpheus heatmap viewer.
  The HTML heatmap is then produced with preconfigured sorting and grouping settings.


  Primary Output files:
  - heatmap.gct, GCT formatted RPKM expression data for morpheus viewer
  - heatmap_vst.gct, GCT formatted VST expression data for morpheus viewer (transformed TotalReads counts)
  - heatmap_vst.html, html of morpheus heatmap with preconfigured settings, VST, no data scaling
  - heatmap_nonorm.html, html of morpheus heatmap with preconfigured settings, RPKM, no data scaling
  - heatmap_norm95.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled per individual sample to 95th percentile
  - heatmap_norm99.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled per individual sample to 99th percentile

  Secondary Output files:
  - master_samplesheet_scaled.tsv, contains formatted information of the input data and files for scaled heatmaps
  - master_samplesheet_vst.tsv, contains formatted information of the input data and files for vst normalized heatmaps
  - output_row_metadata.tsv, row metadata for GCT formatter
  - output_col_metadata.tsv, column metadata for GCT formatter
  - output_counts.tsv, gene expression counts matrix (TotalReads and RPKM)
  - output_counts_vst.tsv, gene expression counts matrix (VST values)

  PARAMS:
    SECTION 1: general
    -h	help		    show this message
    -t  INT			    number of threads
    -a	ARRAY		    array of genelist sample names (no commas in names)
    -b  FILE ARRAY	array of associated annotation files for each gene list from (-c), with header
    -c  FILE ARRAY	array of filtered gene list TSVs (must be headerless, columns are: chr, txStart, txEnd, geneID, L2FC, Strand)
    -e	ARARY		    array of sample names from RNA-Seq experiments (no commas in names)
    -g	FILE ARARY	array of expression table files from RNA-Seq experiments	


  ____________________________________________________________________________________________________
  References:
  - Morpheus, https://software.broadinstitute.org/morpheus
      
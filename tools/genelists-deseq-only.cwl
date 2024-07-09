cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-genelists:v5.0.0


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
      contains formatted information of the input data and files for raw and scaled heatmaps

  master_samplesheet_vst:
    type: File
    outputBinding:
      glob: master_samplesheet_vst.tsv
    doc: |
      contains formatted information of the input data and files for vst normalized heatmap

  master_samplesheet_vst_zscore:
    type: File
    outputBinding:
      glob: master_samplesheet_vst_zscore.tsv
    doc: |
      contains formatted information of the input data and files for vst z-score heatmap

  heatmap_TotalReads_html:
    type: File
    outputBinding:
      glob: heatmap_TotalReads.html
    doc: |
      html of morpheus heatmap with preconfigured settings, TotalReads, no data scaling

  heatmap_vst_html:
    type: File
    outputBinding:
      glob: heatmap_vst.html
    doc: |
      html of morpheus heatmap with preconfigured settings, VST values, no data scaling

  heatmap_vst_zscore_html:
    type: File
    outputBinding:
      glob: heatmap_vst_zscore.html
    doc: |
      html of morpheus heatmap with preconfigured settings, Z-scores of VST values, no data scaling

  heatmap_Rpkm_html:
    type: File
    outputBinding:
      glob: heatmap_Rpkm.html
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, no data scaling

  heatmap_scaled100_html:
    type: File
    outputBinding:
      glob: heatmap_scaled100.html
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, no percentile cutoff

  heatmap_scaled99_html:
    type: File
    outputBinding:
      glob: heatmap_scaled99.html
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, max value set to 99th percentile RPKM value

  heatmap_scaled95_html:
    type: File
    outputBinding:
      glob: heatmap_scaled95.html
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, data scaled  0-99, max value set to 95th percentile RPKM value

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
  Uses only RNA-Seq data to derive visualization data.
  RNA-Seq data in the form of gene expression count matrices are processed to output TotalReads and Rpkm values per gene.
  Additionally, TotalReads are used to produce a separate count matrix of VST values, and VST z-scores.
  This data is then integrated into a single count matrix (one for scaled, a separate matrix for VST, and another for VST z-scores), a row, and a column metadata file.
  These are used as input to an Rscript that will format the 3 files into GCT format for morpheus heatmap generation.


  Primary Output files:
  - heatmap_TotalReads.html, html of morpheus heatmap with preconfigured settings, TotalReads, no data scaling
  - heatmap_vst.html, html of morpheus heatmap with preconfigured settings, VST values, no data scaling
  - heatmap_vst_zscore.html, html of morpheus heatmap with preconfigured settings, Z-scores of VST values, no data scaling
  - heatmap_Rpkm.html, html of morpheus heatmap with preconfigured settings, RPKM, no data scaling
  - heatmap_scaled100.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, no percentile cutoff
  - heatmap_scaled99.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, max value set to 99th percentile RPKM value
  - heatmap_scaled95.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled  0-99, max value set to 95th percentile RPKM value

  Secondary Output files:
  - master_samplesheet_scaled.tsv, contains formatted information of the input data and files for raw and scaled heatmaps
  - master_samplesheet_vst.tsv, contains formatted information of the input data and files for vst normalized heatmaps
  - master_samplesheet_vst_zscore.tsv, contains formatted information of the input data and files for vst normalized heatmaps

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
      
cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genelists:
    - "filter-deseq-for-heatmap.cwl"
    - "filter-diffbind-for-heatmap.cwl"
  samples_nabinding:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "cutandrun-macs2-pe.cwl"
    - "cutandrun-seacr-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"
  samples_rnaseq:
    - "mirna-mirdeep2-se.cwl"
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "trim-quantseq-mrnaseq-se.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  genelist_names:
    type:
    - "null"
    - string[]
    default: null
    label: "Genelist sample names"
    doc: "Array of genelist aliases/sample names, param (-a)"
    'sd:upstreamSource': "genelists/alias"
    'sd:localLabel': true

  genelist_feature_files:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3475"
    label: "feature files from genelist inputs"
    doc: "Array of TSV files with differential genes from DESeq or diffbind pipelines, param (-b)"
    'sd:upstreamSource': "genelists/feature_file"

  genelist_filtered_files:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Filtered differential genes"
    doc: "Array of filtered differential genelists from DESeq or diffbind pipelines, param (-c)"
    'sd:upstreamSource': "genelists/filtered_file"

  sample_names_nabinding:
    type:
      - "null"
      - string[]?
    default: null
    label: "Sample names for ChIP/ATAC/CRT-Seq experiments"
    doc: "Array of aliases for ChIP/ATAC/CRT-Seq experiments for row metadata, param (-d)"
    'sd:upstreamSource': "samples_nabinding/alias"

  sample_names_rnaseq:
    type:
      - "null"
      - string[]?
    default: null
    label: "Sample names for RNA-Seq experiments"
    doc: "Array of aliases for RNA-Seq experiments for column metadata, param (-e)"
    'sd:upstreamSource': "samples_rnaseq/alias"

  datafiles_nabinding:
    type:
    - "null"
    - File[]?
    default: null
    format: "http://edamontology.org/format_2572"
    label: "Array of sample coordinate sorted BAM alignment and BAI index files"
    doc: "Array of sample coordinate sorted BAM alignment and BAI index files, param (-f)"
    'sd:upstreamSource': "samples_nabinding/bambai_pair"

  datafiles_rnaseq:
    type:
    - "null"
    - File[]?
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Array of sample TSV files containing gene annotations with associated TotalReads and Rpkm counts"
    doc: "Array of sample TSV files containing gene annotations with associated TotalReads and Rpkm counts, param (-g)"
    'sd:upstreamSource': "samples_rnaseq/rpkm_genes"

  threads:
    type: int?
    default: 1
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  master_samplesheet:
    type: File
    label: "contains formatted information of the input data and files"
    format: "http://edamontology.org/format_3475"
    doc: "contains formatted information of the input data and files"
    outputSource: data_integration/master_samplesheet

  output_row_metadata:
    type: File
    label: "row metadata for GCT formatter"
    format: "http://edamontology.org/format_3475"
    doc: "row metadata for GCT formatter"
    outputSource: data_integration/output_row_metadata

  output_col_metadata:
    type: File
    label: "column metadata for GCT formatter"
    format: "http://edamontology.org/format_3475"
    doc: "column metadata for GCT formatter"
    outputSource: data_integration/output_col_metadata

  output_counts:
    type: File
    label: "peak average read depth per TSS window and gene expression counts matrix"
    format: "http://edamontology.org/format_3475"
    doc: "peak average read depth per TSS window and gene expression counts matrix"
    outputSource: data_integration/output_counts

  heatmap_gct:
    type: File
    label: "GCT formatted peak and expression data for morpheus viewer"
    format: "http://edamontology.org/format_3709"
    doc: "GCT formatted peak and expression data for morpheus viewer"
    outputSource: data_integration/heatmap_gct

  heatmap_html:
    type: File
    outputSource: morpheus_heatmap/heatmap_html
    label: "Heatmap of peak and expression data"
    doc: |
      Morpheus heatmap in HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  overview:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "combined sample and log files"
    doc: "combined sample and log files"
    outputSource: aggregate_logs_for_overview/overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'


steps:

  data_integration:
    run: ../tools/genelists-deseq-diffbind.cwl
    in:
      threads: threads
      genelist_names: genelist_names
      feature_files: genelist_feature_files
      filtered_files: genelist_filtered_files
      sample_names_nabinding: sample_names_nabinding
      sample_names_rnaseq: sample_names_rnaseq
      bam_files: datafiles_nabinding
      expression_files: datafiles_rnaseq
    out:
      - master_samplesheet
      - output_row_metadata
      - output_col_metadata
      - output_counts
      - heatmap_gct
      - log_file_stdout
      - log_file_stderr

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: data_integration/heatmap_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log

  aggregate_logs_for_overview:
    in:
      master_samplesheet: data_integration/master_samplesheet
      genelists_stdout: data_integration/log_file_stdout
      genelists_stderr: data_integration/log_file_stderr
      heatmap_stdout: morpheus_heatmap/stdout_log
      heatmap_stderr: morpheus_heatmap/stderr_log
    out:
      - overview
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: robertplayer/scidap-genelists:dev
      - class: InitialWorkDirRequirement
        listing:
          - entryname: aggregate.sh
            entry: |
              #!/bin/bash
              printf "# $1\n" > overview.txt
              while read x; do
                printf "\n"
                printf "$x\n"
              done < $1 >> overview.txt
              printf "# $2\n" >> overview.txt
              while read x; do
                printf "\n"
                printf "$x\n"
              done < $2 >> overview.txt
              printf "# $3\n" >> overview.txt
              while read x; do
                printf "\n"
                printf "$x\n"
              done < $3 >> overview.txt
              printf "# $4\n" >> overview.txt
              while read x; do
                printf "\n"
                printf "$x\n"
              done < $4 >> overview.txt
              printf "# $5\n" >> overview.txt
              while read x; do
                printf "\n"
                printf "$x\n"
              done < $5 >> overview.txt
      inputs:
        master_samplesheet:
          type: File
          inputBinding:
            position: 1
        genelists_stdout:
          type: File
          inputBinding:
            position: 2
        genelists_stderr:
          type: File
          inputBinding:
            position: 3
        heatmap_stdout:
          type: File
          inputBinding:
            position: 4
        heatmap_stderr:
          type: File
          inputBinding:
            position: 5
      outputs:
        overview:
          type: File
          outputBinding:
            glob: "overview.txt"
      baseCommand: ["bash", "-c", "aggregate.sh"]     


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Genelists heatmap - peak and expression data visualized together"
label: "Genelists heatmap - peak and expression data visualized together"
s:alternateName: "Genelists heatmap with ChIP/ATAC-Seq peak and RNA-Seq expression data visualized together"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/genelists-deseq-diffbind.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
  - class: s:Organization
    s:legalName: "Datirium, LLC"
    s:member:
      - class: s:Person
        s:name: Artem BArski
        s:email: mailto:Artem.Barski@datirum.com
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@datirium.com
        s:sameAs:
          - id: http://orcid.org/0000-0001-9102-5681


doc: |
  # Genelists heatmap - peak and expression data visualized together

  This visualization workflow takes as input 1 or more genelists derived from the DESeq and/or diffbind workflows along with user-selected samples and visualizes the ChIP/ATAC-Seq peak and/or RNA-Seq expression data visualized together in a single morpheus heatmap.

  ### __References__
    - Morpheus, https://software.broadinstitute.org/morpheus
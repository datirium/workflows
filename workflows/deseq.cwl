cwlVersion: v1.0
class: Workflow

'sd:upstream':
  rnaseq_sample_untreated:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
  rnaseq_sample_treated:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  untreated_files:
    type: File[]
    format:
     - "http://edamontology.org/format_3752"
     - "http://edamontology.org/format_3475"
    label: "Untreated input CSV/TSV files"
    doc: "Untreated input CSV/TSV files"
    'sd:upstreamSource': "rnaseq_sample_untreated/rpkm_common_tss"
    'sd:localLabel': true


  treated_files:
    type: File[]
    format:
     - "http://edamontology.org/format_3752"
     - "http://edamontology.org/format_3475"
    label: "Treated input CSV/TSV files"
    doc: "Treated input CSV/TSV files"
    'sd:upstreamSource': "rnaseq_sample_treated/rpkm_common_tss"
    'sd:localLabel': true

  untreated_col_suffix:
    type: string?
    label: "Untreated RPKM column suffix"
    doc: "Suffix for untreated RPKM column name"
    'sd:layout':
      advanced: true

  treated_col_suffix:
    type: string?
    label: "Treated RPKM column suffix"
    doc: "Suffix for treated RPKM column name"
    'sd:layout':
      advanced: true

  output_filename:
    type: string?
    label: "Output TSV filename"
    doc: "Output TSV filename"
    default: "deseq_results.tsv"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    default: 2
    'sd:layout':
      advanced: true


outputs:

  diff_expr_file:
    type: File
    label: "DESeq resutls, TSV"
    format: "http://edamontology.org/format_3475"
    doc: "DESeq generated list of differentially expressed items grouped by isoforms, genes or common TSS"
    outputSource: deseq/diff_expr_file


steps:

  deseq:
    in:
      untreated_files: untreated_files
      treated_files: treated_files
      untreated_col_suffix: untreated_col_suffix
      treated_col_suffix: treated_col_suffix
      output_filename: output_filename
      threads: threads
    out: [diff_expr_file]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: InlineJavascriptRequirement
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap-deseq:v0.0.5
      inputs:
        untreated_files:
          type: File[]
          inputBinding:
            position: 5
            prefix: "-u"
          doc: "Untreated input CSV/TSV files"
        treated_files:
          type: File[]
          inputBinding:
            position: 6
            prefix: "-t"
          doc: "Treated input CSV/TSV files"
        untreated_col_suffix:
          type: string?
          inputBinding:
            position: 7
            prefix: "-un"
          doc: "Suffix for untreated RPKM column name"
        treated_col_suffix:
          type: string?
          inputBinding:
            position: 8
            prefix: "-tn"
          doc: "Suffix for treated RPKM column name"
        output_filename:
          type: string
          inputBinding:
            position: 9
            prefix: "-o"
          doc: "Output TSV filename"
        threads:
          type: int?
          inputBinding:
            position: 10
            prefix: '-p'
          doc: "Run script using multiple threads"
      outputs:
        diff_expr_file:
          type: File
          outputBinding:
            glob: $(inputs.output_filename)
      baseCommand: [run_deseq.R]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "DESeq - differential gene expression analysis"
label: "DESeq - differential gene expression analysis"
s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/deseq.cwl
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
  Differential gene expression analysis based on the negative binomial distribution

s:about: |
  Differential gene expression analysis
  =====================================

  Differential gene expression analysis based on the negative binomial distribution

  Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.

  DESeq1
  ------

  High-throughput sequencing assays such as RNA-Seq, ChIP-Seq or barcode counting provide quantitative readouts
  in the form of count data. To infer differential signal in such data correctly and with good statistical power,
  estimation of data variability throughout the dynamic range and a suitable error model are required.
  Simon Anders and Wolfgang Huber propose a method based on the negative binomial distribution, with variance and mean
  linked by local regression and present an implementation, [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html),
  as an R/Bioconductor package

  DESeq2
  ------

  In comparative high-throughput sequencing assays, a fundamental task is the analysis of count data,
  such as read counts per gene in RNA-seq, for evidence of systematic changes across experimental conditions.
  Small replicate numbers, discreteness, large dynamic range and the presence of outliers require a
  suitable statistical approach. [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html),
  a method for differential analysis of count data,
  using shrinkage estimation for dispersions and fold changes to improve stability and interpretability of estimates.
  This enables a more quantitative analysis focused on the strength rather than the mere presence of differential expression.

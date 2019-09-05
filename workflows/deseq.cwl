cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:metadata':
- "../metadata/advanced-header.cwl"


'sd:upstream':
  rnaseq_cond_1:
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
  rnaseq_cond_2:
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"


inputs:

  rpkm_isoforms_cond_1:
    type: File[]
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 1, aka 'untreated')"
    doc: "CSV/TSV input files grouped by isoforms (condition 1, aka 'untreated')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_isoforms"
    'sd:localLabel': true

  rpkm_isoforms_cond_2:
    type: File[]
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, aka 'treated')"
    doc: "CSV/TSV input files grouped by isoforms (condition 2, aka 'treated')"
    'sd:upstreamSource': "rnaseq_cond_2/rpkm_isoforms"
    'sd:localLabel': true

  group_by:
    type:
      - "null"
      - type: enum
        symbols: ["isoforms", "genes", "common tss"]
    default: "genes"
    label: "Group by"
    doc: "Grouping method for features: isoforms, genes or common tss"

  threads:
    type: int?
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    default: 1
    'sd:layout':
      advanced: true


outputs:

  diff_expr_file:
    type: File
    label: "Differentially expressed features grouped by isoforms, genes or common TSS"
    format: "http://edamontology.org/format_3475"
    doc: "DESeq generated file of differentially expressed features grouped by isoforms, genes or common TSS in TSV format"
    outputSource: deseq/diff_expr_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Differential Expression Analysis'
        Title: 'Combined DESeq results'
    - scatter:
        tab: 'Volcano Plot'
        Title: 'Volcano'
        xAxisTitle: 'log fold change'
        yAxisTitle: '-log10(pAdj)'
        colors: ["#b3de69"]
        height: 600
        data: [$9, $13]

  plot_lfc_vs_mean:
    type: File
    label: "Plot of normalised mean versus log2 fold change"
    format: "http://edamontology.org/format_3603"
    doc: "Plot of the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples"
    outputSource: deseq/plot_lfc_vs_mean
    'sd:visualPlugins':
    - image:
        tab: 'Other Plots'
        Caption: 'LFC vs mean'

  gene_expr_heatmap:
    type: File
    label: "Heatmap of the 30 most highly expressed features"
    format: "http://edamontology.org/format_3603"
    doc: "Heatmap showing the expression data of the 30 most highly expressed features grouped by isoforms, genes or common TSS, based on the variance stabilisation transformed data"
    outputSource: deseq/gene_expr_heatmap
    'sd:visualPlugins':
    - image:
        tab: 'Other Plots'
        Caption: 'The 30 most highly expressed features'

steps:

  group_isoforms_cond_1:
    run: ../tools/group-isoforms-batch.cwl
    in:
      isoforms_file: rpkm_isoforms_cond_1
    out:
      - genes_file
      - common_tss_file

  group_isoforms_cond_2:
    run: ../tools/group-isoforms-batch.cwl
    in:
      isoforms_file: rpkm_isoforms_cond_2
    out:
      - genes_file
      - common_tss_file

  deseq:
    run: ../tools/deseq-advanced.cwl
    in:
      untreated_files:
        source: [group_by, rpkm_isoforms_cond_1, group_isoforms_cond_1/genes_file, group_isoforms_cond_1/common_tss_file]
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      treated_files:
        source: [group_by, rpkm_isoforms_cond_2, group_isoforms_cond_2/genes_file, group_isoforms_cond_2/common_tss_file]
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      output_filename:
        default: "deseq_results.tsv"
      threads: threads
    out: [diff_expr_file, plot_lfc_vs_mean, gene_expr_heatmap]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "DESeq - differential gene expression analysis"
label: "DESeq - differential gene expression analysis"
s:alternateName: "Differential gene expression analysis based on the negative binomial distribution"

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

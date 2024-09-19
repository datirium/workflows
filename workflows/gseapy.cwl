cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  deseq_experiment:
    - "deseq.cwl"
    - "deseq-for-spikein.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  read_counts_file:
    type: File
    format: "http://edamontology.org/format_3709"
    label: "DESeq experiment"
    doc: "Input gene expression dataset file in txt or gct format. Same with GSEA"
    'sd:upstreamSource': "deseq_experiment/read_counts_file"
    'sd:localLabel': true

  phenotypes_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DESeq experiment"
    doc: "Input class vector (phenotype) file in CLS format. Same with GSEA"
    'sd:upstreamSource': "deseq_experiment/phenotypes_file"
    'sd:localLabel': true

  gene_set_database:
    type:
    - "null"
    - type: enum
      name: "genesetdatabase"
      symbols:
      - H_hallmark_gene_sets
      - C1_positional_gene_sets
      - C2_curated_gene_sets
      - C3_regulatory_target_gene_sets
      - C4_computational_gene_sets
      - C5_ontology_gene_sets
      - C6_oncogenic_signature_gene_sets
      - C7_immunologic_signature_gene_sets
      - C8_cell_type_signature_gene_sets
      - KEGG_2021_Human
      - Reactome_2022
      - WikiPathway_2023_Human
    default: "H_hallmark_gene_sets"
    label: "Gene set database. Ignored if GMT file is privided"
    doc: "Gene set database"

  gene_set_database_file:
    type: File?
    format: "http://edamontology.org/format_2330"
    default: null
    label: "Gene set database file in GMT format"
    doc: "Gene set database file in GMT (Gene Matrix Transposed) format"

  permutation_type:
    type:
    - "null"
    - type: enum
      name: "permutationtype"
      symbols:
      - gene_set
      - phenotype
    default: "gene_set"
    label: "Permutation type"
    doc: "Permutation type. Default: gene_set"

  permutation_count:
    type: int?
    default: 1000
    label: "Number of random permutations"
    doc: "Number of random permutations. For calculating esnulls. Default: 1000"

  min_gene_set_size:
    type: int?
    default: 15
    label: "Min size of input genes presented in Gene Sets"
    doc: "Min size of input genes presented in Gene Sets. Default: 15"
    'sd:layout':
      advanced: true

  max_gene_set_size:
    type: int?
    default: 500
    label: "Max size of input genes presented in Gene Sets"
    doc: "Max size of input genes presented in Gene Sets. Default: 500"
    'sd:layout':
      advanced: true

  ranking_metrics:
    type:
    - "null"
    - type: enum
      name: "rankingmetrics"
      symbols:
      - signal_to_noise
      - t_test
      - ratio_of_classes
      - diff_of_classes
      - log2_ratio_of_classes
    default: "signal_to_noise"
    label: "Methods to calculate correlations of ranking metrics"
    doc: "Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes"

  graphs_pvalue:
    type: float?
    default: 0.05
    label: "Set p-value for enrichment plot and heatmap output"
    doc: "Output only graphs from gene sets with p-value less than this set value. Default: 0.05"
    'sd:layout':
      advanced: true

  seed:
    type: int?
    default: 123
    label: "Number of random seed. Default: None"
    doc: "Number of random seed. Default: None"
    'sd:layout':
      advanced: true

  threads_count:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  gseapy_enrichment_report_csv:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Enrichment report"
    doc: "Enrichment report in original csv format"
    outputSource: run_gseapy/enrichment_report

  gseapy_enrichment_report_tsv:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Enrichment report"
    doc: "Enrichment report converted to tsv format for scidap table"
    outputSource: convert_to_tsv/output_file
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Gene Set Enrichment Table"
        Title: "Gene Set Enrichment Table"

  gseapy_report_summary:
    type: File?
    label: "Markdown formatted table with summary stats"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted table with summary stats"
    outputSource: run_gseapy/report_summary
    'sd:visualPlugins':
    - markdownView:
        tab: "Overview"

  gseapy_enrichment_plots:
    type: File
    label: "Compressed TAR with enrichment plots"
    doc: "Compressed TAR with enrichment plots"
    outputSource: run_gseapy/enrichment_plots

  gseapy_enrichment_heatmaps:
    type: File
    label: "Compressed TAR with enrichment heatmaps"
    doc: "Compressed TAR with enrichment heatmaps"
    outputSource: run_gseapy/enrichment_heatmaps

  gseapy_filtered_enrichment_plots:
    type: File
    label: "Compressed TAR with filtered enrichment plots"
    doc: "Compressed TAR with enrichment plots having a p-value less than graphs_pvalue input value"
    outputSource: run_gseapy/filtered_enrichment_plots

  gseapy_filtered_enrichment_heatmaps:
    type: File
    label: "Compressed TAR with filtered enrichment heatmaps"
    doc: "Compressed TAR with enrichment heatmaps having a p-value less than graphs_pvalue input value"
    outputSource: run_gseapy/filtered_enrichment_heatmaps

  gseapy_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "GSEApy stdout log"
    doc: "GSEApy stdout log"
    outputSource: run_gseapy/stdout_log

  gseapy_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "GSEApy stderr log"
    doc: "GSEApy stderr log"
    outputSource: run_gseapy/stderr_log


steps:

  run_gseapy:
    run: ../tools/gseapy.cwl
    in:
      read_counts_file: read_counts_file
      phenotypes_file: phenotypes_file
      gene_set_database:
        source: [gene_set_database, gene_set_database_file]
        valueFrom: $(self[1]?self[1]:self[0])
      permutation_type: permutation_type
      permutation_count: permutation_count
      min_gene_set_size: min_gene_set_size
      max_gene_set_size: max_gene_set_size
      ranking_metrics: ranking_metrics
      graphs_pvalue: graphs_pvalue
      seed: seed
      threads: threads_count
    out:
      - enrichment_report
      - report_summary
      - enrichment_plots
      - enrichment_heatmaps
      - filtered_enrichment_plots
      - filtered_enrichment_heatmaps
      - stdout_log
      - stderr_log

  convert_to_tsv:
    run: ../tools/custom-bash.cwl
    in:
      input_file: run_gseapy/enrichment_report
      script:
        default: |
          cat "$0" | tr "," "\t" > `basename $0 csv`tsv
    out: [output_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "GSEApy - Gene Set Enrichment Analysis in Python"
label: "GSEApy - Gene Set Enrichment Analysis in Python"
s:alternateName: "GSEApy - Gene Set Enrichment Analysis in Python"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/gseapy.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  GSEAPY: Gene Set Enrichment Analysis in Python
  ==============================================

  Gene Set Enrichment Analysis is a computational method that determines whether an a priori
  defined set of genes shows statistically significant, concordant differences between two
  biological states (e.g. phenotypes).

  GSEA requires as input an expression dataset, which contains expression profiles for multiple samples.
  While the software supports multiple input file formats for these datasets, the tab-delimited GCT format
  is the most common. The first column of the GCT file contains feature identifiers (gene ids or symbols in
  the case of data derived from RNA-Seq experiments). The second column contains a description of the feature;
  this column is ignored by GSEA and may be filled with “NA”s. Subsequent columns contain the expression
  values for each feature, with one sample's expression value per column. It is important to note that there
  are no hard and fast rules regarding how a GCT file's expression values are derived. The important point is
  that they are comparable to one another across features within a sample and comparable to one another
  across samples. Tools such as DESeq2 can be made to produce properly normalized data (normalized counts)
  which are compatible with GSEA.

  Documents
  ==============================================
  - GSEA Home Page: https://www.gsea-msigdb.org/gsea/index.jsp
  - Results Interpretation: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results
  - GSEA User Guide: https://gseapy.readthedocs.io/en/latest/faq.html
  - GSEAPY Docs: https://gseapy.readthedocs.io/en/latest/introduction.html

  References
  ==============================================
  - Subramanian, Tamayo, et al. (2005, PNAS), https://www.pnas.org/content/102/43/15545
  - Mootha, Lindgren, et al. (2003, Nature Genetics), http://www.nature.com/ng/journal/v34/n3/abs/ng1180.html
  - Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14).
  - Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377 .
  - Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A. Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90
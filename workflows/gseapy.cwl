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
      - WikiPathways_2019_Human
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

  ascending_rank_sorting:
    type: boolean?
    default: false
    label: "Ascending rank metric sorting order"
    doc: "Ascending rank metric sorting order. Default: False"

  graphs_count:
    type: int?
    default: 20
    label: "Numbers of top graphs produced"
    doc: "Numbers of top graphs produced. Default: 20"
    'sd:layout':
      advanced: true

  seed:
    type: int?
    default: 123
    label: "Number of random seed. Default: None"
    doc: "Number of random seed. Default: None"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  gseapy_enrichment_report:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Enrichment report"
    doc: "Enrichment report"
    outputSource: convert_to_tsv/output_file
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Gene Set Enrichment"
        Title: "Gene Set Enrichment"

  gseapy_enrichment_plots:
    type: File
    label: "Compressed TAR with enrichment plots"
    doc: "Compressed TAR with enrichment plots"
    outputSource: rename_enrichment_plots/target_file

  gseapy_enrichment_heatmaps:
    type: File
    label: "Compressed TAR with enrichment heatmaps"
    doc: "Compressed TAR with enrichment heatmaps"
    outputSource: rename_enrichment_heatmaps/target_file

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

  summary_report:
    type: File
    format: "http://edamontology.org/format_3835"
    label: "Enrichment report"
    doc: "Enrichment report"
    outputSource: report_summary/summary_file
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  summary_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr log"
    doc: "stderr log"
    outputSource: report_summary/log_file_stderr

  summary_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout log"
    doc: "stdout log"
    outputSource: report_summary/log_file_stdout


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
      ascending_rank_sorting: ascending_rank_sorting
      graphs_count: graphs_count
      seed: seed
      threads: threads
    out:
      - enrichment_report
      - enrichment_plots
      - enrichment_heatmaps
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

  enrichment_plots_to_folder:
    run: ../tools/files-to-folder.cwl
    in:
      input_files: run_gseapy/enrichment_plots
    out: [folder]

  compress_enrichment_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: enrichment_plots_to_folder/folder
    out: [compressed_folder]

  rename_enrichment_plots:
    run: ../tools/rename.cwl
    in:
      source_file: compress_enrichment_plots/compressed_folder
      target_filename:
        default: "enrichment_plots.tar.gz"
    out: [target_file]

  enrichment_heatmaps_to_folder:
    run: ../tools/files-to-folder.cwl
    in:
      input_files: run_gseapy/enrichment_heatmaps
    out: [folder]

  compress_enrichment_heatmaps:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: enrichment_heatmaps_to_folder/folder
    out: [compressed_folder]

  rename_enrichment_heatmaps:
    run: ../tools/rename.cwl
    in:
      source_file: compress_enrichment_heatmaps/compressed_folder
      target_filename:
        default: "enrichment_heatmaps.tar.gz"
    out: [target_file]

  report_summary:
    run: ../tools/gseapy-reportsummary.cwl
    in:
      read_counts_file: read_counts_file
      phenotypes_file: phenotypes_file
      enrichment_report: convert_to_tsv/output_file
    out:
      - summary_file
      - log_file_stderr
      - log_file_stdout


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


# doc:
#   $include: ../descriptions/gseapy.md


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

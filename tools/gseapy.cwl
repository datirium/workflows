cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-gseapy:v1.0.0
inputs:
  read_counts_file:
    type: File
    inputBinding:
      prefix: -d
    doc: Input gene expression dataset file in txt or gct format (from DESeq workflow)
  phenotypes_file:
    type: File
    inputBinding:
      prefix: -c
    doc: Input class vector (phenotype) file in CLS format (from DESeq workflow)
  gene_set_database:
    type:
    - File
    - type: enum
      name: genesetdatabase
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
    inputBinding:
      prefix: -g
    doc: Gene set database
  permutation_type:
    type:
    - 'null'
    - type: enum
      name: permutationtype
      symbols:
      - gene_set
      - phenotype
    inputBinding:
      prefix: -r
    doc: 'Permutation type. Default: gene_set'
  permutation_count:
    type: int?
    inputBinding:
      prefix: -n
    doc: 'Number of random permutations. For calculating esnulls. Default: 1000'
  min_gene_set_size:
    type: int?
    inputBinding:
      prefix: -w
    doc: 'Min size of input genes presented in Gene Sets. Default: 15'
  max_gene_set_size:
    type: int?
    inputBinding:
      prefix: -x
    doc: 'Max size of input genes presented in Gene Sets. Default: 500'
  ranking_metrics:
    type:
    - 'null'
    - type: enum
      name: rankingmetrics
      symbols:
      - signal_to_noise
      - t_test
      - ratio_of_classes
      - diff_of_classes
      - log2_ratio_of_classes
    inputBinding:
      prefix: -m
    doc: 'Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes'
  graphs_pvalue:
    type: float?
    inputBinding:
      prefix: -p
    doc: 'Output only graphs from gene sets with less than thsi set p-value. Default: 0.05'
  seed:
    type: int?
    inputBinding:
      prefix: -s
    doc: 'Number of random seed. Default: None'
  threads:
    type: int?
    inputBinding:
      prefix: -t
    doc: Threads number
outputs:
  enrichment_report:
    type: File?
    outputBinding:
      glob: GSEApy_reports/*.report.csv
  report_summary:
    type: File?
    outputBinding:
      glob: reportsummary.md
  enrichment_plots:
    type: File
    outputBinding:
      glob: all_gseapy_enrichment_plots.tar.gz
  enrichment_heatmaps:
    type: File
    outputBinding:
      glob: all_gseapy_heatmap_plots.tar.gz
  filtered_enrichment_plots:
    type: File
    outputBinding:
      glob: filtered_enrichment_plots.tar.gz
  filtered_enrichment_heatmaps:
    type: File
    outputBinding:
      glob: filtered_heatmap_plots.tar.gz
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- /usr/local/bin/run_gseapy.sh
stdout: gseapy_stdout.log
stderr: gseapy_stderr.log
doc: "GSEAPY: Gene Set Enrichment Analysis using Python\n==============================================\n\nGene Set Enrichment Analysis is a computational method that determines whether an a priori\ndefined set of genes shows statistically significant, concordant differences between two\nbiological states (e.g. phenotypes).\n\nGSEA requires as input an expression dataset, which contains expression profiles for multiple samples.\nWhile the software supports multiple input file formats for these datasets, the tab-delimited GCT format\nis the most common. The first column of the GCT file contains feature identifiers (gene ids or symbols in\nthe case of data derived from RNA-Seq experiments). The second column contains a description of the feature;\nthis column is ignored by GSEA and may be filled with “NA”s. Subsequent columns contain the expression\nvalues for each feature, with one sample's expression value per column. It is important to note that there\nare no hard and fast rules regarding how a GCT file's expression values are derived. The important point is\nthat they are comparable to one another across features within a sample and comparable to one another\nacross samples. Tools such as DESeq2 can be made to produce properly normalized data (normalized counts)\nwhich are compatible with GSEA.\n\nPrimary Output files:\n- gseapy.gsea.gene_set.report.csv, table of gene set enrichment scores and pvalues in CSV format\n- enrichment plots and heatmaps per gene set, depending on user-selected dataset\n\nSecondary Output files:\n- reportsummary.md, markdown with general filtering and analysis statistics\n\nPARAMS:\n    SECTION 1: general\n  -h\thelp\t\tshow this message\n  -t  INT\t\t\tnumber of threads\n  -c  FILE\t\tInput class vector (phenotype) file in CLS format (from DESeq workflow)\n  -d\tFILE\t\tInput gene expression dataset file in txt or gct format (from DESeq workflow)\n  -g\tFILE\t\tGene set database, from prefilled dropdown or use-provided .gmt file\n  -r\tSTRING\t\tPermutation type. Default: \"gene_set\"\n  -n\tINT\t\t\tNumber of random permutations. For calculating esnulls. Default: 1000\n  -w\tINT\t\t\tMin size of input genes presented in Gene Sets. Default: 15\n  -x\tINT\t\t\tMax size of input genes presented in Gene Sets. Default: 500\n  -m\tSTRING\t\tMethods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes\n  -s\tINT\t\t\tNumber of random seed. Default: None\n  -p\tFLOAT\t\tOutput only graphs from gene sets with less than thsi set p-value. Default: 1.00\n"
label: gseapy

cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/gseapy:v0.0.6


inputs:

  read_counts_file:
    type: File
    inputBinding:
      prefix: "-d"
    doc: "Input gene expression dataset file in txt or gct format. Same with GSEA"

  phenotypes_file:
    type: File
    inputBinding:
      prefix: "-c"
    doc: "Input class vector (phenotype) file in CLS format. Same with GSEA"

  gene_set_database:
    type:
    - File
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
    inputBinding:
      prefix: "-g"
    doc: "Gene set database"

  permutation_type:
    type:
    - "null"
    - type: enum
      name: "permutationtype"
      symbols:
      - gene_set
      - phenotype
    inputBinding:
      prefix: "-t"
    doc: "Permutation type. Default: gene_set"

  permutation_count:
    type: int?
    inputBinding:
      prefix: "-n"
    doc: "Number of random permutations. For calculating esnulls. Default: 1000"

  min_gene_set_size:
    type: int?
    inputBinding:
      prefix: "--min-size"
    doc: "Min size of input genes presented in Gene Sets. Default: 15"

  max_gene_set_size:
    type: int?
    inputBinding:
      prefix: "--max-size"
    doc: "Max size of input genes presented in Gene Sets. Default: 500"

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
    inputBinding:
      prefix: "-m"
    doc: "Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes"

  ascending_rank_sorting:
    type: boolean?
    inputBinding:
      prefix: "-a"
    doc: "Ascending rank metric sorting order. Default: False"

  graphs_count:
    type: int?
    inputBinding:
      prefix: "--graph"
    doc: "Numbers of top graphs produced. Default: 20"

  seed:
    type: int?
    inputBinding:
      prefix: "-s"
    doc: "Number of random seed. Default: None"

  threads:
    type: int?
    inputBinding:
      prefix: "-p"
    doc: "Threads number"


outputs:

  enrichment_report:
    type: File?
    outputBinding:
      glob: "GSEApy_reports/*.report.csv"

  enrichment_plots:
    type: File[]
    outputBinding:
      glob: "GSEApy_reports/*.gsea.pdf"

  enrichment_heatmaps:
    type: File[]
    outputBinding:
      glob: "GSEApy_reports/*.heatmap.pdf"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["gseapy", "gsea"]

stdout: gseapy_stdout.log
stderr: gseapy_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "gseapy"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/gseapy.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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

s:about: |
  usage: gseapy gsea [-h] -d DATA -c CLS -g GMT [-t perType] [-o] [-f]
                    [--fs width height] [--graph int] [--no-plot] [-v]
                    [-n nperm] [--min-size int] [--max-size int] [-w float]
                    [-m] [-a] [-s] [-p procs]

  optional arguments:
    -h, --help            show this help message and exit

  Input files arguments:
    -d DATA, --data DATA  Input gene expression dataset file in txt format.Same
                          with GSEA.
    -c CLS, --cls CLS     Input class vector (phenotype) file in CLS format.
                          Same with GSEA.
    -g GMT, --gmt GMT     Gene set database in GMT format. Same with GSEA.
    -t perType, --permu-type perType
                          Permutation type. Same with GSEA, choose from
                          {'gene_set', 'phenotype'}

  Output arguments:
    -o , --outdir         The GSEApy output directory. Default: the current
                          working directory
    -f , --format         File extensions supported by Matplotlib active
                          backend, choose from {'pdf', 'png', 'jpeg','ps',
                          'eps','svg'}. Default: 'pdf'.
    --fs width height, --figsize width height
                          The figsize keyword argument need two parameters to
                          define. Default: (6.5, 6)
    --graph int           Numbers of top graphs produced. Default: 20
    --no-plot             Speed up computing by suppressing the plot output.This
                          is useful only if data are interested. Default: False.
    -v, --verbose         Increase output verbosity, print out progress of your
                          job

  GSEA advanced arguments:
    -n nperm, --permu-num nperm
                          Number of random permutations. For calculating
                          esnulls. Default: 1000
    --min-size int        Min size of input genes presented in Gene Sets.
                          Default: 15
    --max-size int        Max size of input genes presented in Gene Sets.
                          Default: 500
    -w float, --weight float
                          Weighted_score of rank_metrics. For weighting input
                          genes. Choose from {0, 1, 1.5, 2}. Default: 1
    -m , --method         Methods to calculate correlations of ranking metrics.
                          Choose from {'signal_to_noise', 't_test',
                          'ratio_of_classes',
                          'diff_of_classes','log2_ratio_of_classes'}. Default:
                          'log2_ratio_of_classes'
    -a, --ascending       Rank metric sorting order. If the -a flag was chosen,
                          then ascending equals to True. Default: False.
    -s , --seed           Number of random seed. Default: None
    -p procs, --threads procs
                          Number of Processes you are going to use. Default: 1
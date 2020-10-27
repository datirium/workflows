cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/soupx:v0.0.1


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--counts"
    doc: |
      Path to the output folder produced by 'cellranger count' command

  genelist_file:
    type: File?
    inputBinding:
      prefix: "--genes"
    doc: |
      Path to the file with target genes (headerless, one gene per line)

  expression_threshold:
    type: float?
    inputBinding:
      prefix: "--threshold"
    doc: |
      Expression threshold for displaying target genes on a plot (expression > threshold)
      Default: 0.0

  fdr:
    type: float?
    inputBinding:
      prefix: "--fdr"
    doc: |
      FDR cutoff for expression ratio plots
      Default: 0.05

  round_counts:
    type: boolean?
    inputBinding:
      prefix: "--round"
    doc: |
      Round adjusted counts to integers
      Default: False

  matrix_format_version:
    type:
    - "null"
    - type: enum
      name: "matrix_format_version"
      symbols: ["2","3"]
    inputBinding:
      prefix: "--format"
    doc: |
      Output matrix format version. Corresponds to the latest Cell Ranger matrix format
      Default: 3

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix
      Default: soupx


outputs:

  adjusted_feature_bc_matrices_folder:
    type: Directory
    outputBinding:
      glob: "*adjusted_counts"
    doc: |
      Folder with adjusted feature-barcode matrices in MEX format

  adjusted_feature_bc_matrices_h5:
    type: File
    outputBinding:
      glob: "*adjusted_counts.h5"
    doc: |
      Adjusted feature-barcode matrices in HDF5 format

  raw_gene_expression_plots:
    type: File?
    outputBinding:
      glob: "*raw_gene_expression_plots.pdf"
    doc: |
      Raw gene expression plots

  adjusted_gene_expression_plots:
    type: File?
    outputBinding:
      glob: "*adjusted_gene_expression_plots.pdf"
    doc: |
      Adjusted gene expression plots

  raw_gene_expression_to_pure_soup_ratio_plots:
    type: File?
    outputBinding:
      glob: "*raw_gene_expression_to_pure_soup_ratio_plots.pdf"
    doc: |
      Raw gene expression to pure soup ratio plots
  
  raw_to_adjusted_gene_expression_ratio_plots:
    type: File?
    outputBinding:
      glob: "*raw_to_adjusted_gene_expression_ratio_plots.pdf"
    doc: |
      Raw to adjusted gene expression ratio plots

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_soupx.R"]


stdout: soupx_stdout.log
stderr: soupx_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "SoupX - an R package for the estimation and removal of cell free mRNA contamination"
label: "SoupX - an R package for the estimation and removal of cell free mRNA contamination"
s:alternateName: "SoupX - an R package for the estimation and removal of cell free mRNA contamination"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/soupx.cwl
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
  In droplet based, single cell RNA-seq experiments, there is always a certain amount of background
  mRNAs present in the dilution that gets distributed into the droplets with cells and sequenced
  along with them. The net effect of this is to produce a background contamination that represents
  expression not from the cell contained within a droplet, but the solution that contained the cells.

  This collection of cell free mRNAs floating in the input solution (henceforth referred to as
  “the soup”) is created from cells in the input solution being lysed. Because of this, the soup looks
  different for each input solution and strongly resembles the expression pattern obtained by summing
  all the individual cells.

  The aim of this package is to provide a way to estimate the composition of this soup, what fraction
  of UMIs are derived from the soup in each droplet and produce a corrected count table with the soup
  based expression removed.

  The method to do this consists of three parts:

  - Calculate the profile of the soup.
  - Estimate the cell specific contamination fraction.
  - Infer a corrected expression matrix.  


s:about: |
  usage: run_soupx.R [-h] --counts COUNTS [--genes GENES]
                                    [--threshold THRESHOLD] [--fdr FDR]
                                    [--round] [--format {2,3}] [--output OUTPUT]

  Runs SoupX - an R package for the estimation and removal of cell free mRNA
  contamination

  optional arguments:
    -h, --help            show this help message and exit
    --counts COUNTS       Path to the output folder produced by 'cellranger
                          count' command
    --genes GENES         Path to the file with target genes (headerless, one
                          gene per line)
    --threshold THRESHOLD
                          Expression threshold for displaying target genes on a
                          plot (expression > threshold). Default: 0.0
    --fdr FDR             FDR cutoff for expression ratio plots. Default: 0.05
    --round               Round adjusted counts to integers
    --format {2,3}        Output matrix format version. Default: 3 (corresponds
                          to the latest Cell Ranger matrix format)
    --output OUTPUT       Output prefix. Default: soupx
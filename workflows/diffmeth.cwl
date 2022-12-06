cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement

'sd:serviceTag': "Analysis"

'sd:upstream':
  biological_condition1:
   - "bismark-methylation-pe.cwl"
   - "bismark-methylation-se.cwl"
  biological_condition2:
   - "bismark-methylation-pe.cwl"
   - "bismark-methylation-se.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  c1_name:
    type: string
    default: "condition1"
    label: "Condition 1 name, single word with letters and numbers only"
    doc: "Condition 1 name, single word with letters and numbers only"

  c2_name:
    type: string
    default: "condition2"
    label: "Condition 2 name, single word with letters and numbers only"
    doc: "Condition 2 name, single word with letters and numbers only"

  c1_files:
    type: File[]
    format: "http://edamontology.org/format_3003"
    'sd:upstreamSource': "biological_condition1/bismark_coverage_file"
    label: "Bismark coverage files for condition1 (minumum of 2):"
    'sd:localLabel': true
    doc: "Array of input coverage files from Bismark workflow."

  c2_files:
    type: File[]
    format: "http://edamontology.org/format_3003"
    'sd:upstreamSource': "biological_condition2/bismark_coverage_file"
    label: "Bismark coverage files for condition2 (minumum of 2):"
    'sd:localLabel': true
    doc: "Array of input coverage files from Bismark workflow."

  c1_aliases:
    type: string[]
    'sd:upstreamSource': "biological_condition1/alias"
    doc: "Array of input coverage files from Bismark workflow."

  c2_aliases:
    type: string[]
    'sd:upstreamSource': "biological_condition2/alias"
    doc: "Array of input coverage files from Bismark workflow."

  genome:
    type:
    - "null"
    - type: enum
      name: "genomes"
      symbols:
      - hg19
      - hg38
      - mm9
      - mm10
      - rn5
    default: "hg38"
    label: "Sample genome:"
    'sd:localLabel': true
    doc: "Sample genomes available: hg19, hg38, mm9, mm10, rn5"

  threads:
    type: int
    default: 2
    label: "Number of threads for parallel processing:"
    'sd:localLabel': true
    doc: "Number of threads for parallel processing"


outputs:

  samplesheet_md:
    type: File
    label: "Samplesheet with condition labels"
    doc: "Samplesheet with condition labels"
    outputSource: run_rnbeads_diff/samplesheet

  report_tar:
    type: File
    label: "Compressed TAR with RnBeads reports"
    doc: "Compressed TAR with RnBeads reports"
    outputSource: run_rnbeads_diff/report_tar

  report_directory:
    type: Directory
    label: "Output directory for all rnbeads results, preserves structure for html references"
    doc: "Output directory for all rnbeads results, preserves structure for html references"
    outputSource: run_rnbeads_diff/report_directory

  stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout log"
    doc: "stdout log"
    outputSource: run_rnbeads_diff/stdout_log

  stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr log"
    doc: "stderr log"
    outputSource: run_rnbeads_diff/stderr_log

  samplesheet_overview:
    type: File
    label: "Samplesheet with condition labels"
    doc: "Samplesheet with condition labels"
    outputSource: run_rnbeads_diff/samplesheet_overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  report_data_import_html:
    type: File
    label: "Data Import HTML report"
    doc: "Data Import HTML report"
    outputSource: run_rnbeads_diff/report_data_import_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  report_qc_html:
    type: File
    label: "Quality Control HTML report"
    doc: "Quality Control HTML report"
    outputSource: run_rnbeads_diff/report_qc_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  report_preprocessing_html:
    type: File
    label: "Preprocessing HTML report"
    doc: "Preprocessing HTML report"
    outputSource: run_rnbeads_diff/report_preprocessing_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  report_differential_methylation_html:
    type: File
    label: "Differential methylation HTML report"
    doc: "Differential methylation HTML report"
    outputSource: run_rnbeads_diff/report_differential_methylation_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  dm_sites_stats_table:
    type: File
    label: "DM sites"
    doc: "DM sites"
    outputSource: run_rnbeads_diff/dm_sites_stats

  dm_cpg_stats_table:
    type: File
    label: "DM CpGs"
    doc: "DM CpGs"
    outputSource: run_rnbeads_diff/dm_cpg_stats

  dm_tiling_stats_table:
    type: File
    label: "DM tiling"
    doc: "DM tiling"
    outputSource: run_rnbeads_diff/dm_tiling_stats

  dm_genes_stats_table:
    type: File
    label: "DM genes"
    doc: "DM genes"
    outputSource: run_rnbeads_diff/dm_genes_stats
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Differentially Methylated Genes'
        Title: 'Table of differentially methylated genes.'

  dm_sites_group1_igv:
    type: File
    label: "DM sites"
    doc: "DM sites fdr<0.10"
    outputSource: run_rnbeads_diff/dm_sites_group1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm sites grp1"
        height: 60

  dm_cpg_group1_igv:
    type: File
    label: "DM CpGs"
    doc: "DM CpGs fdr<0.10"
    outputSource: run_rnbeads_diff/dm_cpg_group1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm CpGs grp1"
        height: 60

  dm_tiling_group1_igv:
    type: File
    label: "DM tiling"
    doc: "DM tiling fdr<0.10"
    outputSource: run_rnbeads_diff/dm_tiling_group1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm tiling grp1"
        height: 60

  dm_genes_group1_igv:
    type: File
    label: "DM genes"
    doc: "DM genes fdr<0.10"
    outputSource: run_rnbeads_diff/dm_genes_group1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm genes grp1"
        height: 60

  dm_sites_group2_igv:
    type: File
    label: "DM sites"
    doc: "DM sites fdr<0.10"
    outputSource: run_rnbeads_diff/dm_sites_group2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm sites grp2"
        height: 60

  dm_cpg_group2_igv:
    type: File
    label: "DM CpGs"
    doc: "DM CpGs fdr<0.10"
    outputSource: run_rnbeads_diff/dm_cpg_group2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm CpGs grp2"
        height: 60

  dm_tiling_group2_igv:
    type: File
    label: "DM tiling"
    doc: "DM tiling fdr<0.10"
    outputSource: run_rnbeads_diff/dm_tiling_group2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm tiling grp2"
        height: 60

  dm_genes_group2_igv:
    type: File
    label: "DM genes"
    doc: "DM genes fdr<0.10"
    outputSource: run_rnbeads_diff/dm_genes_group2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "dm genes grp2"
        height: 60


steps:

  run_rnbeads_diff:
    label: "Run RnBeads wrapper for differential methylation analysis between groups"
    doc: |
      Runs differential methylation analysis using 2 or more BismarkCov formatted bed files per
      condition. Outputs a compressed tar of standard RnBeads reports directory.
    run: ../tools/rnbeads-diff.cwl
    in:
      threads: threads
      genome:
        source: genome
        valueFrom: $(self)
      condition1_name: c1_name
      condition2_name: c2_name
      condition1_filepaths: c1_files
      condition2_filepaths: c2_files
      condition1_aliases: c1_aliases
      condition2_aliases: c2_aliases
    out:
      - samplesheet
      - samplesheet_overview
      - report_tar
      - report_directory
      - report_data_import_html
      - report_qc_html
      - report_preprocessing_html
      - report_differential_methylation_html
      - dm_sites_stats
      - dm_cpg_stats
      - dm_tiling_stats
      - dm_genes_stats
      - dm_sites_group1
      - dm_cpg_group1
      - dm_tiling_group1
      - dm_genes_group1
      - dm_sites_group2
      - dm_cpg_group2
      - dm_tiling_group2
      - dm_genes_group2
      - stdout_log
      - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Differential Methylation Workflow"
label: "Differential Methylation Workflow"
s:alternateName: "Differential Methylation using RnBeads with BismarkCov outputs"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/diffmeth.cwl
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
  A basic differential methylation analysis workflow using BismarkCov formatted bed files as input. Analysis is conducted on site and region level according to the sample groups specified by user (limited to 2 conditions in this workflow implementation). See report html files for detailed descriptions of analyses and results interpretation.

  ### __Inputs__
  *General Info (required\*):*
  - Experiment short name/Alias* - a unique name for the sample (e.g. what was used on tubes while processing it)
  - Condition 1 name - name defining condition/group 1
  - Condition 2 name - name defining condition/group 2
  - Bismark coverage files* for condition1 - minumum of 2 is required for analysis
  - Bismark coverage files* for condition2 - minumum of 2 is required for analysis
  - Sample genome - available options: hg19, hg38, mm9, mm10, rn5
  - Number of threads for steps that support multithreading - default set to `2`

  *BismarkCov formatted bed:
    https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
    The genome-wide cytosine report (optional) is tab-delimited in the following format (1-based coords):
    <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

  ### __Outputs__
  Intermediate and final downloadable outputs include:
  - stdout and stderr log files
  - reports directory containing:
      reports/
      ├── configuration
      ├── data_import.html
      ├── data_import_data
      ├── data_import_images
      ├── data_import_pdfs
      ├── differential_methylation.html
      ├── differential_methylation_data
      ├── differential_methylation_images
      ├── differential_methylation_pdfs
      ├── preprocessing.html
      ├── preprocessing_data
      ├── preprocessing_images
      ├── preprocessing_pdfs
      ├── quality_control.html
      ├── quality_control_data
      ├── quality_control_images
      ├── quality_control_pdfs
      ├── tracks_and_tables.html
      ├── tracks_and_tables_data
      ├── tracks_and_tables_images
      └── tracks_and_tables_pdfs

  Other outputs include tables and bed files for IGV for DM sites, CpG, tiling, and genes:
  
  https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/RnBeads/inst/doc/RnBeads_Annotations.pdf
  2.1 Sites ($sites)
    Currently, every data package examines cytosines in the context of CpG and contains an
    annnotation table of all CpGs in the respective genome. CpG density and GC content are
    also computed for the neighborhood of length 100 base pairs centered on each locus. The
    total number of dinucleotides annotated in HG19 is 28,217,009 represented both on the
    forward and reverse DNA strands.
  2.4 Regions ($tiling, $cpg, $genes)
    Every data package defines the following sets of regions for the dedicated assembly:
    - GpG islands The CpG island track is downloaded from the dedicated FTP directory of
      the UCSC Genome Browser.
    - Tiling regions Tiling regions with a window size of 5 kilobases are defined over the
      whole genome.
    - Genes and promoters Ensembl3 gene definitions are downloaded using the biomaRt package.
      A promoter isdefined as the region spanning 1,500 bases upstream and 500 bases
      downstream of the transcription start site of the corresponding gene.
    CpG density and GC content are computed for all region types listed above.

  ### __Data Analysis Steps__
  1. generate sample sheet with associated conditions for testing in RnBeads
  2. setup rnbeads analyses in R, and run differential methylation analysis
  3. save outputs

  ### __References__
    - https://rnbeads.org/materials/example_3/differential_methylation.html
    - Makambi, K. (2003) Weighted inverse chi-square method for correlated significance tests. Journal of Applied Statistics, 30(2), 225234
    - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4216143/
    - Assenov Y, Müller F, Lutsik P, Walter J, Lengauer T, Bock C. Comprehensive analysis of DNA methylation data with RnBeads. Nat Methods. 2014 Nov;11(11):1138-1140. doi: 10.1038/nmeth.3115. Epub 2014 Sep 28. PMID: 25262207; PMCID: PMC4216143.

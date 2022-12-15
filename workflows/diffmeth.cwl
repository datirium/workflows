cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

'sd:serviceTag': "Analysis"

'sd:upstream':
  biological_condition1:
   - "bismark-methylation-pe.cwl"
   - "bismark-methylation-se.cwl"
  biological_condition2:
   - "bismark-methylation-pe.cwl"
   - "bismark-methylation-se.cwl"
  genome_indices:
   - "bismark-index.cwl"


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

  c1_methperc_bigwigs:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Methylation statuses bigWig coverage file"
    doc: "bigWig file summarising cytosine methylation as a percent in bigWig format"
    'sd:upstreamSource': "biological_condition1/bigwig_coverage_file"

  c2_methperc_bigwigs:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Methylation statuses bigWig coverage file"
    doc: "bigWig file summarising cytosine methylation as a percent in bigWig format"
    'sd:upstreamSource': "biological_condition2/bigwig_coverage_file"

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
    default: 4
    label: "Number of threads for parallel processing:"
    'sd:localLabel': true
    doc: "Number of threads for parallel processing"

  annotation_file:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"
    'sd:upstreamSource': "genome_indices/annotation"


outputs:

  c1_bismark_meth_perc:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "percent methylation per site for biological condition 1"
    doc: "percent methylation bigWig file(s) for biological condition 1"
    outputSource: run_rnbeads_diff/c1_bigwigs
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "cond1:%meth"
        height: 120

  c2_bismark_meth_perc:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "percent methylation per site for biological condition 2"
    doc: "percent methylation bigWig file(s) for biological condition 2"
    outputSource: run_rnbeads_diff/c2_bigwigs
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "cond2:%meth"
        height: 120

  samplesheet_csv:
    type: File
    label: "Samplesheet with condition labels"
    doc: "Samplesheet with condition labels"
    outputSource: run_rnbeads_diff/samplesheet

  samplesheet_md:
    type: File
    label: "Samplesheet with condition labels"
    doc: "Samplesheet with condition labels"
    outputSource: run_rnbeads_diff/samplesheet_overview
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

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

  meth_stats_genes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "DM genes"
    doc: "DM genes"
    outputSource: run_rnbeads_diff/meth_stats_genes

  meth_stats_promoters:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "DM genes"
    doc: "DM genes"
    outputSource: run_rnbeads_diff/meth_stats_promoters

  meth_stats_cpg:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "DM CpGs"
    doc: "DM CpGs"
    outputSource: run_rnbeads_diff/meth_stats_cpg

  meth_stats_tiling:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "DM tiling"
    doc: "DM tiling"
    outputSource: run_rnbeads_diff/meth_stats_tiling

  meth_stats_sites:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "DM sites"
    doc: "DM sites"
    outputSource: run_rnbeads_diff/meth_stats_sites

  sig_dm_sites_annotated_table:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "gene annotation for sig diff meth sites"
    doc: "gene annotation for sig diff meth sites"
    outputSource: run_rnbeads_diff/sig_dm_sites_annotated
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Diff Meth Sites'
        Title: 'Table of differentially methylated sites with closest gene annotations.'

  sig_dm_sites_igvtrack:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "bed track for sig diffmeth sites"
    doc: "bed track for sig diffmeth sites, see workflow docs for format details"
    outputSource: run_rnbeads_diff/sig_dm_sites_igvtrack
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "sigDiffMeth Sites"
        height: 40

  sig_dm_cpg_annotated_table:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "gene annotation for sig diff meth cpgs"
    doc: "gene annotation for sig diff meth cpgs"
    outputSource: run_rnbeads_diff/sig_dm_cpg_annotated
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Diff Meth CpG Islands'
        Title: 'Table of differentially methylated CpG Islands with closest gene annotations.'

  sig_dm_cpg_igvtrack_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "bed track for sig diffmeth CpG Islands"
    doc: "bed track for sig diffmeth CpG Islands, see workflow docs for format details"
    outputSource: run_rnbeads_diff/sig_dm_cpg_igvtrack
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "sigDiffMeth CpG Islands"
        height: 40


steps:

  run_rnbeads_diff:
    label: "Run wrapper for RnBeads differential methylation pipeline, with downstream processing for tables and IGV"
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
      condition1_methperc_bigwigs: c1_methperc_bigwigs
      condition2_methperc_bigwigs: c2_methperc_bigwigs
      condition1_aliases: c1_aliases
      condition2_aliases: c2_aliases
      refgene_annotations: annotation_file
    out:
      - c1_bigwigs
      - c2_bigwigs
      - samplesheet
      - samplesheet_overview
      - report_tar
      - report_directory
      - report_data_import_html
      - report_qc_html
      - report_preprocessing_html
      - report_differential_methylation_html
      - meth_stats_genes
      - meth_stats_promoters
      - meth_stats_cpg
      - meth_stats_tiling
      - meth_stats_sites
      - sig_dm_sites_igvtrack
      - sig_dm_sites_annotated
      - sig_dm_cpg_igvtrack
      - sig_dm_cpg_annotated
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
  A basic differential methylation analysis workflow using BismarkCov formatted bed files as input to the RnBeads tool. Analysis is conducted on region and sites levels according to the sample groups specified by user (limited to 2 conditions in this workflow implementation). See report html files for detailed descriptions of analyses and results interpretation.

  ### __Inputs__
  *General Info:*
  - Experiment short name/Alias* - a unique name for the sample (e.g. what was used on tubes while processing it)
  - Condition 1 name - name defining condition/group 1
  - Condition 2 name - name defining condition/group 2
  - Bismark coverage files* for condition1 - minumum of 2 is required for analysis
  - Bismark coverage files* for condition2 - minumum of 2 is required for analysis
  - Sample genome - available options: hg19, hg38, mm9, mm10, rn5
  - Genome type - indicate mismark index used for upstream samples (input for conditions 1 and 2)

  *Advanced:*
  - Number of threads for steps that support multithreading - default set to `4`

  *[BismarkCov formatted bed](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf):
    The genome-wide cytosine report (optional) is tab-delimited in the following format (1-based coords):
    <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

  ### __Outputs__
  Intermediate and final downloadable outputs include:
  - sig_dm_sites.bed ([bed for IGV](https://genome.ucsc.edu/FAQ/FAQformat.html#format1); sig diff meth sites)
  - sig_dm_sites_annotated.tsv (tsv for TABLE; for each site above, closest single gene annotation)
    - Site_id, unique indentifer per methylated site
    - Site_Chr, chromosome of methylated site
    - Site_position, 1-based position in chr of methylated site
    - Site_strand, strand of methylated site
    - Log2_Meth_Quotient, log2 of the quotient in methylation: log2((mean.g1+epsilon)/(mean.g2+epsilon)), where epsilon:=0.01. In case of paired analysis, it is the mean of the pairwise quotients.
    - FDR, adjusted p-values, all <0.10 assumed to be significant
    - Coverage_score, value between 0-1000 reflects strength of mean coverage difference between conditions and equals [1000-(1000/(meancov_g1-meancov_g2)^2](https://www.wolframalpha.com/input?i=solve+1000-%281000%2F%28x%5E2%29%29), if meancov_g1-meancov_g2==0, score=0, elif score<1==1, else score
    - meancov_g1, mean coverage of condition1
    - meancov_g2, mean coverage of condition2
    - refSeq_id, RefSeq gene id
    - Gene_id, gene symbol
    - Chr, gene chromosome
    - txStart, gene transcription start position
    - tsEnd, gene transcription end position
    - txStrand, gene strand
  - stdout and stderr log files
  - Packaged RnBeads reports directory (reports.tar.gz) contains:
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

  Reported methylation is in the form of regions (genes, promoters, cpg, tiling) and specific sites:
   - genes - Ensembl gene definitions are downloaded using the biomaRt package.
   - promoters - A promoter is defined as the region spanning 1,500 bases upstream and 500 bases downstream of the transcription start site of the corresponding gene
   - cpg - the CpG islands from the UCSC Genome Browser
   - tiling - a window size of 5 kilobases are defined over the whole genome
   - sites - all cytosines in the context of CpGs in the respective genome

  ### __Data Analysis Steps__
  1. generate sample sheet with associated conditions for testing in RnBeads
  2. setup rnbeads analyses in R, and run differential methylation analysis
  3. process output diffmeth files for regions and sites
  4. find single closest gene annotations for all significantly diffmeth sites
  5. package and save rnbeads report directory
  6. clean up report dir for html outputs

  ### __References__
    - https://rnbeads.org/materials/example_3/differential_methylation.html
    - Makambi, K. (2003) Weighted inverse chi-square method for correlated significance tests. Journal of Applied Statistics, 30(2), 225234
    - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4216143/
    - Assenov Y, Müller F, Lutsik P, Walter J, Lengauer T, Bock C. Comprehensive analysis of DNA methylation data with RnBeads. Nat Methods. 2014 Nov;11(11):1138-1140. doi: 10.1038/nmeth.3115. Epub 2014 Sep 28. PMID: 25262207; PMCID: PMC4216143.
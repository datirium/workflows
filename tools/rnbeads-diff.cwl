cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-rnbeads:stable


inputs:

  threads:
    type: int
    inputBinding:
      prefix: "-t"
    doc: "Number of threads for parallel processing"

  genome:
    type: string
    inputBinding:
      prefix: "-g"
    doc: "Sample genome, available options: hg19, hg38, mm9, mm10, rn5"

  condition1_name:
    type: string
    inputBinding:
      prefix: "-a"
    doc: "Name for group 1, default: condition1"

  condition2_name:
    type: string
    inputBinding:
      prefix: "-b"
    doc: "Name for group 2, default: condition2"

  condition1_filepaths:
    type: File[]
    inputBinding:
      prefix: "-c"
      itemSeparator: ","
    doc: "List of absolute filepaths to BismarkCov formatted bed files for group 1"

  condition2_filepaths:
    type: File[]
    inputBinding:
      prefix: "-d"
      itemSeparator: ","
    doc: "List of absolute filepaths to BismarkCov formatted bed files for group 2"

  condition1_aliases:
    type: string[]
    inputBinding:
      prefix: "-j"
      itemSeparator: ","
    doc: "List of alias sample names for group 1"

  condition2_aliases:
    type: string[]
    inputBinding:
      prefix: "-k"
      itemSeparator: ","
    doc: "List of alias sample names for group 2"

outputs:

  samplesheet:
    type: File
    outputBinding:
      glob: sample_annotation.csv
    doc: |
      Prepared samplesheet from user inputs, used for condition labeling by RnBeads.

  samplesheet_overview:
    type: File
    outputBinding:
      glob: sample_annotation.md
    doc: |
      Prepared samplesheet from user inputs, used for display on Overview tab.

  report_tar:
    type: File
    outputBinding:
      glob: reports.tar.gz
    doc: |
      Compressed reports dir tarball containing html summary files and additional files.

  report_data_import:
    type: Directory
    outputBinding:
      glob: reports/data_import_images
    doc: |
      image directory for data_import report html

  report_data_import_html:
    type: File
    outputBinding:
      glob: reports/data_import.html
    doc: |
      Differential

  report_qc:
    type: Directory
    outputBinding:
      glob: reports/quality_control_images
    doc: |
      image directory for quality_control report html

  report_qc_html:
    type: File
    outputBinding:
      glob: reports/quality_control.html
    doc: |
      Differential

  report_preprocessing:
    type: Directory
    outputBinding:
      glob: reports/preprocessing_images
    doc: |
      image directory for preprocessing report html

  report_preprocessing_html:
    type: File
    outputBinding:
      glob: reports/preprocessing.html
    doc: |
      preprocessing report

  report_differential_methylation:
    type: Directory
    outputBinding:
      glob: reports/differential_methylation_images
    doc: |
      image directory for differential_methylation report html

  report_differential_methylation_html:
    type: File
    outputBinding:
      glob: reports/differential_methylation.html
    doc: |
      Differential methylation HTML report for visualization.  

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_rnbeads_diff.sh"]
stdout: rnbeads_diff_stdout.log
stderr: rnbeads_diff_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "rnbeads-diff"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/rnbeads-diff.cwl
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
  Wrapper for RnBeads differential methylation pipeline.
  Output reports directory in container at '/tmp/reports/', includes:
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

  OPTIONS:
  -h  help	show this message
  -t  INT	number of threads
  -g  STRING   Sample genome, available options: hg19, hg38, mm9, mm10, rn5
  -a  STRING     name of condition1
  -b  STRING     name of condition2
  -c  LIST	comma separated list of absolute filepaths to all condition1 bed files (BismarkCov format)
  -d  LIST	comma separated list of absolute filepaths to all condition2 bed files (BismarkCov format)
  OPTIONAL:
  -o  DIR        absolute path to output directory for reports, default '/tmp/[reports]'

  BismarkCov formatted bed:
      https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
      The genome-wide cytosine report (optional) is tab-delimited in the following format (1-based coords):
      <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

  ____________________________________________________________________________________________________
  References:
      https://rnbeads.org/materials/example_3/differential_methylation.html
          Makambi, K. (2003) Weighted inverse chi-square method for correlated significance tests.
          Journal of Applied Statistics, 30(2), 225234
      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4216143/
          Assenov Y, Müller F, Lutsik P, Walter J, Lengauer T, Bock C. Comprehensive analysis of DNA
          methylation data with RnBeads. Nat Methods. 2014 Nov;11(11):1138-1140. doi: 10.1038/nmeth.3115.
          Epub 2014 Sep 28. PMID: 25262207; PMCID: PMC4216143.

cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/diffbind:v0.0.3


inputs:

  read_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-r1"
    doc: "Read files for condition 1"

  read_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-r2"
    doc: "Read files for condition 2"

  peak_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-p1"
    doc: "Peak files for condition 1. Format corresponds to -pf"

  peak_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-p2"
    doc: "Peak files for condition 2. Format corresponds to -pf"

  sample_names_cond_1:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "-n1"
    doc: "Sample names for condition 1. Default: basenames of -r1 without extensions"

  sample_names_cond_2:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "-n2"
    doc: "Sample names for condition 2. Default: basenames of -r2 without extensions"

  peakformat:
    type:
      - "null"
      - type: enum
        name: "peakformat"
        symbols: ["raw","bed","narrow","macs","bayes","tpic","sicer","fp4","swembl","csv","report"]
    inputBinding:
      prefix: "-pf"
    doc: "Peak files format. One of [raw, bed, narrow, macs, bayes, tpic, sicer, fp4, swembl, csv, report]. Default: macs"

  name_cond_1:
    type: string?
    inputBinding:
      prefix: "-c1"
    doc: "Condition 1 name, single word with letters and numbers only. Default: condition_1"

  name_cond_2:
    type: string?
    inputBinding:
      prefix: "-c2"
    doc: "Condition 2 name, single word with letters and numbers only. Default: condition_2"

  fragmentsize:
    type: int?
    inputBinding:
      prefix: "-fs"
    doc: "Extended each read from its endpoint along the appropriate strand. Default: 125bp"

  remove_duplicates:
    type: boolean?
    inputBinding:
      prefix: "-rd"
    doc: "Remove reads that map to exactly the same genomic position. Default: false"

  analysis_method:
    type:
      - "null"
      - type: enum
        name: "method"
        symbols: ["deseq2", "edger"]
    inputBinding:
      prefix: "-me"
    doc: "Method by which to analyze differential binding affinity. Default: deseq2"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output prefix. Default: diffbind"

  threads:
    type: int?
    inputBinding:
      prefix: "-th"
    doc: "Threads number. Default: 1"


outputs:

  diffbind_report_file:
    type: File
    outputBinding:
      glob: "*_diffpeaks.tsv"
    doc: "Differential binding analysis results exported as TSV"

  peak_correlation_heatmap:
    type: File
    outputBinding:
      glob: "*_peak_overlap_correlation_heatmap.png"
    doc: "Peak overlap correlation heatmap"

  counts_correlation_heatmap:
    type: File
    outputBinding:
      glob: "*_counts_correlation_heatmap.png"
    doc: "Counts correlation heatmap"

  all_data_correlation_heatmap:
    type: File
    outputBinding:
      glob: "*_correlation_heatmap_based_on_all_normalized_data.png"
    doc: "Correlation heatmap based on all normalized data"

  db_sites_correlation_heatmap:
    type: File
    outputBinding:
      glob: "*_correlation_heatmap_based_on_db_sites_only.png"
    doc: "Correlation heatmap based on DB sites only"

  db_sites_binding_heatmap:
    type: File
    outputBinding:
      glob: "*_binding_heatmap_based_on_db_sites.png"
    doc: "Binding heatmap based on DB sites"

  pca_plot:
    type: File
    outputBinding:
      glob: "*_pca.png"
    doc: "PCA plot using affinity data for only differentially bound sites"

  ma_plot:
    type: File
    outputBinding:
      glob: "*_ma.png"
    doc: "MA plot for tested conditions"

  volcano_plot:
    type: File
    outputBinding:
      glob: "*_volcano.png"
    doc: "Volcano plot for tested conditions"

  boxplot_plot:
    type: File
    outputBinding:
      glob: "*_boxplot.png"
    doc: "Box plots of read distributions for significantly differentially bound (DB) sites"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_diffbind.R"]
stdout: diffbind_stdout.log
stderr: diffbind_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

label: "DiffBind - Differential Binding Analysis of ChIP-Seq Peak Data"
s:name: "DiffBind - Differential Binding Analysis of ChIP-Seq Peak Data"
s:alternateName: "Compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data."

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/diffbind.cwl
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
  Runs R script to compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data.

s:about: |
  usage: /Users/kot4or/workspaces/cwl_ws/workflows/tools/dockerfiles/scripts/run_diffbind.R
        [-h] -r1 READ1 [READ1 ...] -r2 READ2 [READ2 ...] -p1 PEAK1 [PEAK1 ...]
        -p2 PEAK2 [PEAK2 ...] [-n1 [NAME1 [NAME1 ...]]]
        [-n2 [NAME2 [NAME2 ...]]]
        [-pf {raw,bed,narrow,macs,bayes,tpic,sicer,fp4,swembl,csv,report}]
        [-c1 CONDITION1] [-c2 CONDITION2] [-fs FRAGMENTSIZE] [-rd]
        [-me {edger,deseq2}] [-th THREADS] [-o OUTPUT]

  Differential binding analysis of ChIP-Seq experiments using affinity (read
  count) data

  optional arguments:
    -h, --help            show this help message and exit
    -r1 READ1 [READ1 ...], --read1 READ1 [READ1 ...]
                          Read files for condition 1. Minimim 2 files in BAM
                          format
    -r2 READ2 [READ2 ...], --read2 READ2 [READ2 ...]
                          Read files for condition 2. Minimim 2 files in BAM
                          format
    -p1 PEAK1 [PEAK1 ...], --peak1 PEAK1 [PEAK1 ...]
                          Peak files for condition 1. Minimim 2 files in format
                          set with -pf
    -p2 PEAK2 [PEAK2 ...], --peak2 PEAK2 [PEAK2 ...]
                          Peak files for condition 2. Minimim 2 files in format
                          set with -pf
    -n1 [NAME1 [NAME1 ...]], --name1 [NAME1 [NAME1 ...]]
                          Sample names for condition 1. Default: basenames of
                          -r1 without extensions
    -n2 [NAME2 [NAME2 ...]], --name2 [NAME2 [NAME2 ...]]
                          Sample names for condition 2. Default: basenames of
                          -r2 without extensions
    -pf {raw,bed,narrow,macs,bayes,tpic,sicer,fp4,swembl,csv,report}, --peakformat {raw,bed,narrow,macs,bayes,tpic,sicer,fp4,swembl,csv,report}
                          Peak files format. One of [raw, bed, narrow, macs,
                          bayes, tpic, sicer, fp4, swembl, csv, report].
                          Default: macs
    -c1 CONDITION1, --condition1 CONDITION1
                          Condition 1 name, single word with letters and numbers
                          only. Default: condition_1
    -c2 CONDITION2, --condition2 CONDITION2
                          Condition 2 name, single word with letters and numbers
                          only. Default: condition_2
    -fs FRAGMENTSIZE, --fragmentsize FRAGMENTSIZE
                          Extended each read from its endpoint along the
                          appropriate strand. Default: 125bp
    -rd, --removedup      Remove reads that map to exactly the same genomic
                          position. Default: false
    -me {edger,deseq2}, --method {edger,deseq2}
                          Method by which to analyze differential binding
                          affinity. Default: deseq2
    -th THREADS, --threads THREADS
                          Threads to use
    -o OUTPUT, --output OUTPUT
                          Output prefix. Default: diffbind
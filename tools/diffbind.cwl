cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/diffbind:v0.0.1


inputs:

  read_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-r1"
    doc: "Read files for condition 1. Minimim 2 files in BAM format"

  read_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-r2"
    doc: "Read files for condition 2. Minimim 2 files in BAM format"

  peak_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-p1"
    doc: "Peak files for condition 1. Minimim 2 files in format set with -pf"

  peak_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-p2"
    doc: "Peak files for condition 2. Minimim 2 files in format set with -pf"

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

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output prefix. Default: diffbind"


outputs:

  diffbind_report_file:
    type: File
    outputBinding:
      glob: "*_diffpeaks.tsv"
    doc: "Differential binding analysis results exported as TSV"

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
        [-c1 CONDITION1] [-c2 CONDITION2] [-o OUTPUT]

  Differential binding analysis of ChIP-Seq peak data

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
    -o OUTPUT, --output OUTPUT
                          Output prefix. Default: diffbind
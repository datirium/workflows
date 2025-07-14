cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    http_proxy: $(inputs.http_proxy?inputs.http_proxy:"")
    https_proxy: $(inputs.https_proxy?inputs.https_proxy:"")
- class: ResourceRequirement
  ramMin: 7024
  coresMin: 1
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastqdwnld:v0.0.5
inputs:
  srr_id:
    type:
    - string
    - type: array
      items: string
    inputBinding:
      position: 60
    doc: |
      SRR identifiers
  split_files:
    type: boolean?
    inputBinding:
      position: 10
      prefix: --split-files
    doc: "Write reads into separate files. Read \nnumber will be suffixed to the file name.  \nNOTE! The `--split-3` option is recommended. \nIn cases where not all spots have the same \nnumber of reads, this option will produce \nfiles that WILL CAUSE ERRORS in most programs \nwhich process split pair fastq files. \n"
  split_3:
    type: boolean?
    inputBinding:
      position: 11
      prefix: --split-3
    doc: "3-way splitting for mate-pairs. For each \nspot, if there are two biological reads \nsatisfying filter conditions, the first is \nplaced in the `*_1.fastq` file, and the \nsecond is placed in the `*_2.fastq` file. If \nthere is only one biological read \nsatisfying the filter conditions, it is \nplaced in the `*.fastq` file.All other \nreads in the spot are ignored.\n"
  http_proxy:
    type: string?
    doc: |
      Optional HTTP proxy settings
  https_proxy:
    type: string?
    doc: |
      Optional HTTPS proxy settings
outputs:
  fastq_files:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*.gz'
  metadata_xml:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*.xml'
  report_md:
    type: File
    outputBinding:
      glob: report.md
  collected_metadata:
    type: File
    outputBinding:
      glob: collected_metadata.tsv
  run_acc:
    type:
    - 'null'
    - type: array
      items: string
    outputBinding:
      loadContents: true
      glob: collected_metadata.tsv
      outputEval: |
        ${
          var pattern = /run_acc\:.*/;
          var splitted_line = self[0].contents.match(pattern)[0].trim().split(" ").slice(1);
          return (!!splitted_line.length)?splitted_line:null;
        }
  experiment_acc:
    type:
    - 'null'
    - type: array
      items: string
    outputBinding:
      loadContents: true
      glob: collected_metadata.tsv
      outputEval: |
        ${
          var pattern = /experiment_acc\:.*/;
          var splitted_line = self[0].contents.match(pattern)[0].trim().split(" ").slice(1);
          return (!!splitted_line.length)?splitted_line:null;
        }
  study_acc:
    type:
    - 'null'
    - type: array
      items: string
    outputBinding:
      loadContents: true
      glob: collected_metadata.tsv
      outputEval: |
        ${
          var pattern = /study_acc\:.*/;
          var splitted_line = self[0].contents.match(pattern)[0].trim().split(" ").slice(1);
          return (!!splitted_line.length)?splitted_line:null;
        }
  biosample:
    type:
    - 'null'
    - type: array
      items: string
    outputBinding:
      loadContents: true
      glob: collected_metadata.tsv
      outputEval: |
        ${
          var pattern = /biosample\:.*/;
          var splitted_line = self[0].contents.match(pattern)[0].trim().split(" ").slice(1);
          return (!!splitted_line.length)?splitted_line:null;
        }
  bioproject:
    type:
    - 'null'
    - type: array
      items: string
    outputBinding:
      loadContents: true
      glob: collected_metadata.tsv
      outputEval: |
        ${
          var pattern = /bioproject\:.*/;
          var splitted_line = self[0].contents.match(pattern)[0].trim().split(" ").slice(1);
          return (!!splitted_line.length)?splitted_line:null;
        }
  log_stdout:
    type: stdout
  log_stderr:
    type: stderr
baseCommand:
- sra_download.sh
stdout: sra_download_stdout.log
stderr: sra_download_stderr.log
label: FASTQ Download
doc: |
  FASTQ Download

  Assists in downloading problematic single-cell sequencing
  data from Sequence Read Archive (SRA)

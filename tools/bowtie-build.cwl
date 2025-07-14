cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 15250
  coresMin: 1
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {
          "class": "Directory",
          "basename": inputs.index_base_name,
          "listing": [],
          "writable": true}
      ]
    }
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bowtie:v1.2.0
inputs:
  fasta_file:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      itemSeparator: ','
      position: 25
    doc: |
      comma-separated list of files with ref sequences
  index_base_name:
    type: string
    inputBinding:
      position: 26
      valueFrom: $("./"+self+"/"+self)
    doc: |
      write Ebwt data to files with this dir/basename
  force_large_index:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 3
      prefix: --large-index
    doc: |
      force generated index to be 'large', even if ref has fewer than 4 billion nucleotides
  color:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 4
      prefix: --color
    doc: |
      build a colorspace index
  noauto:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 5
      prefix: --noauto
    doc: |
      disable automatic -p/--bmax/--dcv memory-fitting
  packed:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 6
      prefix: --packed
    doc: |
      use packed strings internally; slower, less memory
  bmax:
    type:
    - 'null'
    - int
    inputBinding:
      position: 7
      prefix: --bmax
    doc: |
      max bucket sz for blockwise suffix-array builder
  bmaxdivn:
    type:
    - 'null'
    - int
    inputBinding:
      position: 8
      prefix: --bmaxdivn
    doc: |
      max bucket sz as divisor of ref len (default: 4)
  dcv:
    type:
    - 'null'
    - int
    inputBinding:
      position: 9
      prefix: --dcv
    doc: |
      diff-cover period for blockwise (default: 1024)
  nodc:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 10
      prefix: --nodc
    doc: |
      disable diff-cover (algorithm becomes quadratic)
  noref:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 11
      prefix: --noref
    doc: |
      don't build .3/.4 index files
  justref:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 12
      prefix: --justref
    doc: |
      just build .3/.4 index files
  offrate:
    type:
    - 'null'
    - int
    inputBinding:
      position: 13
      prefix: --offrate
    doc: |
      SA is sampled every 2^<int> BWT chars (default: 5)
  ftabchars:
    type:
    - 'null'
    - int
    inputBinding:
      position: 14
      prefix: --ftabchars
    doc: |
      # of chars consumed in initial lookup (default: 10)
  ntoa:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 15
      prefix: --ntoa
    doc: |
      convert Ns in reference to As
  seed:
    type:
    - 'null'
    - int
    inputBinding:
      position: 16
      prefix: --seed
    doc: |
      seed for random number generator
  quiet:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 17
      prefix: --quiet
    doc: |
      verbose output (for debugging)
outputs:
  indices_folder:
    type: Directory
    outputBinding:
      glob: $(inputs.index_base_name)
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bowtie-build
stderr: bowtie_stderr.log
stdout: bowtie_stdout.log
doc: |
  Tool runs bowtie-build
  Not supported parameters:
    -c  -  reference sequences given on cmd line (as <seq_in>)
label: bowtie-build

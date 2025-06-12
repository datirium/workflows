cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_folder_name = function() { if (inputs.output_folder_name == ""){ var root = inputs.genome_fasta_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.genome_fasta_file.basename:root; } else { return inputs.output_folder_name; } };
hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:8.0.1
inputs:
  genome_fasta_file:
    type: File
    inputBinding:
      position: 5
      prefix: --fasta
    doc: |
      Genome FASTA file
  annotation_gtf_file:
    type: File
    inputBinding:
      position: 6
      prefix: --genes
    doc: |
      GTF annotation file
  output_folder_name:
    type: string?
    inputBinding:
      position: 7
      prefix: --genome
      valueFrom: $(get_output_folder_name())
    default: ''
    doc: |
      Unique genome name, used to name output folder
  threads:
    type: int?
    inputBinding:
      valueFrom: $(["--nthreads", self, "--localcores", self])
      position: 8
    doc: |
      Number of threads used during STAR
      genome index. And the max cores the
      pipeline may request at one time.
      Default: 1 for --nthreads and all
      available for --localcores
  memory_limit:
    type: int?
    inputBinding:
      position: 9
      prefix: --memgb
    doc: |
      Maximum memory (GB) used when aligning reads with STAR
      Defaults: 16
outputs:
  indices_folder:
    type: Directory
    outputBinding:
      glob: $(get_output_folder_name())
    doc: |
      Cellranger-compatible reference folder that includes
      STAR indices and some additional files
  chrom_length_file:
    type: File
    outputBinding:
      glob: $(get_output_folder_name() + "/star/chrNameLength.txt")
    doc: |
      Chromosome length file in TSV format
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- cellranger
- mkref
stdout: cellranger_mkref_stdout.log
stderr: cellranger_mkref_stderr.log
label: Cell Ranger Reference (RNA)
doc: |
  Cell Ranger Reference (RNA)

  Builds a reference genome of a selected species
  for quantifying gene expression.

  Both --nthreads and --localcores parameters are
  configured through "threads" input.

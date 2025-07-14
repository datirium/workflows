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
      Genome FASTA file. Hard/soft-masked
      files are not allowed.
  annotation_gtf_file:
    type: File
    inputBinding:
      position: 6
      prefix: --genes
    doc: |
      GTF annotation file. Should include
      gene_biotype/transcript_biotype fields
  output_folder_name:
    type: string?
    inputBinding:
      position: 7
      prefix: --genome
      valueFrom: $(get_output_folder_name())
    default: ''
    doc: |
      Unique genome name, used
      to name the output folder
  threads:
    type: int?
    inputBinding:
      position: 8
      prefix: --localcores
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available
  memory_limit:
    type: int?
    inputBinding:
      position: 9
      prefix: --memgb
    doc: |
      Maximum memory (GB) used.
      Defaults: 16
outputs:
  indices_folder:
    type: Directory
    outputBinding:
      glob: $(get_output_folder_name())
    doc: |
      Cell Ranger V(D)J-compatible reference
      folder. This folder will include V(D)J
      segment FASTA file.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- cellranger
- mkvdjref
stdout: cellranger_mkvdjref_stdout.log
stderr: cellranger_mkvdjref_stderr.log
label: Cell Ranger Reference (VDJ)
doc: |
  Cell Ranger Reference (VDJ)

  Builds a reference genome of a selected species for V(D)J
  contigs assembly and clonotype calling.

  Input --seqs is not implemented.

  Chromosome names in GTF file should correspond to the chromosome
  names in FASTA file.

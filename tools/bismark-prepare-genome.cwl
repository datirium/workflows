cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.genome_folder,
                  "writable": true
                }
              ]
    }
inputs:
  genome_folder:
    type: Directory
    inputBinding:
      position: 2
    label: Genome folder
    doc: Genome folder with FASTA files
outputs:
  indices_folder:
    type: Directory
    label: Bismark indices folder
    doc: Bismark generated indices folder
    outputBinding:
      glob: $(inputs.genome_folder.basename)
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bismark_genome_preparation
stderr: bismark_build_stderr.log
stdout: bismark_build_stdout.log
doc: |
  bismark_genome_preparation script generates indices using Bowtie2 aligners by default.
label: bismark-prepare-genome

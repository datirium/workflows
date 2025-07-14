cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2
inputs:
  indices_folder:
    type: Directory
    inputBinding:
      position: 3
    label: Bismark indices folder
    doc: Path to Bismark generated indices folder
  fastq_file_r1:
    type: File
    inputBinding:
      position: 5
    label: FASTQ read 1 file
    doc: Uncompressed or gzipped FASTQ read 1 file
  fastq_file_r2:
    type: File?
    inputBinding:
      position: 7
    label: Optional FASTQ read 2 file
    doc: Optional uncompressed or gzipped FASTQ read 2 file
  processes:
    type: int?
    inputBinding:
      position: 1
      prefix: --multicore
    label: Number of Bismark instances to run
    doc: Set the number of parallel Bismark instances to run concurrently. Each Bismark instance runs four Bowtie2 aligners
  threads:
    type: int?
    inputBinding:
      position: 2
      prefix: -p
    label: Number of Bowtie2 threads to use
    doc: Set the number of threads for each Bowtie2 aligner
outputs:
  bam_file:
    type: File
    label: BAM alignment file
    doc: Bismark generated BAM alignment file
    outputBinding:
      glob: '*.bam'
  alignment_report:
    type: File
    label: Bismark alignment and methylation report
    doc: Bismark generated alignment and methylation summary report
    outputBinding:
      glob: '*.txt'
baseCommand:
- bismark
- --non_directional
arguments:
- valueFrom: |
    ${
      if (inputs.fastq_file_r1 && inputs.fastq_file_r2){
        return "-1";
      }
      return null;
    }
  position: 4
- valueFrom: |
    ${
      if (inputs.fastq_file_r1 && inputs.fastq_file_r2){
        return "-2";
      }
      return null;
    }
  position: 6
doc: |
  Default aligner - Bowtie2.
  Parameters used:
  --non_directional
    The sequencing library was constructed in a non strand-specific manner, alignments to all four bisulfite strands
    will be reported. (The current Illumina protocol for BS-Seq is directional, in which case the strands complementary
    to the original strands are merely theoretical and should not exist in reality. Specifying directional alignments
    (which is the default) will only run 2 alignment threads to the original top (OT) or bottom (OB) strands in parallel
    and report these alignments. This is the recommended option for strand-specific libraries).
label: bismark-align

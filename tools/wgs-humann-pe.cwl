cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: ResourceRequirement
  ramMin: 30510
  coresMin: 8
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-metaphlan:v1.0.0
inputs:
  script_command:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nStdout log file for wgs-humann-pe.cwl tool:\n"
      R1=$0
      R2=$1
      printf "EXECUTION:\n"
      #   commands start
      printf "\trun humann analysis for PE reads\n"
      echo "concat paired reads together (humann does not take paired information into account)"
      cp $R1 combined.fq
      cat $R2 >> combined.fq
      humann --threads 8 --input combined.fq --output humann_outdir
      echo "normalize RPKs to relative abundance"
      humann_renorm_table --input humann_outdir/*_genefamilies.tsv \
        --output humann_outdir/normalized_genefamilies-cpm.tsv --units cpm --update-snames
      echo "regroup genes to other functional categories"
      humann_regroup_table --input humann_outdir/normalized_genefamilies-cpm.tsv \
        --output humann_outdir/rxn-cpm.tsv --groups uniref90_rxn
      echo "attach human-readable names to features"
      humann_rename_table --input humann_outdir/rxn-cpm.tsv \
          --output humann_outdir/rxn-cpm-named.tsv --names metacyc-rxn
      echo "cleaning up tmp files"
      rm combined.fq
      printf "END OF SCRIPT\n"
    inputBinding:
      position: 1
  read1file:
    type: File
    label: R1 fastq
    inputBinding:
      position: 6
    doc: FASTQ file 1 of paired end read data.
  read2file:
    type: File
    label: R2 fastq
    inputBinding:
      position: 7
    doc: FASTQ file 2 of paired end read data.
outputs:
  genefamilies_rpk:
    type: File
    outputBinding:
      glob: humann_outdir/combined_genefamilies.tsv
  genefamilies_cpm:
    type: File
    outputBinding:
      glob: humann_outdir/normalized_genefamilies-cpm.tsv
  regroup_to_rxn_cpm:
    type: File
    outputBinding:
      glob: humann_outdir/rxn-cpm.tsv
  regroup_to_rxn_cpm_named:
    type: File
    outputBinding:
      glob: humann_outdir/rxn-cpm-named.tsv
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: log.stdout
stderr: log.stderr
doc: |
  Tool runs humann3 for functional assignment of metagenomic data sets. Databases used are the
    chocophlan full nucleotide database and filtered uniref90 protein database.

  HUMAnN (the HMP Unified Metabolic Analysis Network) is a method for efficiently and accurately
    profiling the abundance of microbial metabolic pathways and other molecular functions from
    metagenomic or metatranscriptomic sequencing data. It is appropriate for any type of microbial
    community, not just the human microbiome (the name "HUMAnN" is a historical product of the
    method's origins in the Human Microbiome Project).
label: humann3

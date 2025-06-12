cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-tgif:v1.0.0
inputs:
  input_fastq:
    type: File
    inputBinding:
      prefix: -f
    doc: FASTQ file of WGS or NCATS enriched library (ideally long read, short read single-end only)
  read_type:
    type:
    - 'null'
    - type: enum
      name: read type
      symbols:
      - map-ont
      - map-pb
      - sr
    inputBinding:
      prefix: -x
    doc: 'sequencing read type (minimap2 param: map-ont, map-pb, sr)'
  reference_fasta:
    type: File
    inputBinding:
      prefix: -r
    doc: reference fasta file for mapping
  plasmid_fasta:
    type: File
    inputBinding:
      prefix: -i
    doc: plasmid fasta file for mapping (only a single sequence, entire sequence should be on a single line in the file, i.e. max lines in the file is 2)
  yn_igv_output:
    type:
    - 'null'
    - type: enum
      symbols:
      - y
      - n
    inputBinding:
      prefix: -s
    doc: 'Output sorted bam of downselected reads for IGV? Default: y'
  yn_plot_output:
    type:
    - 'null'
    - type: enum
      symbols:
      - y
      - n
    inputBinding:
      prefix: -p
    doc: 'Output read pileup plot per probable insertion site? Default: y'
  threads_count:
    type: int?
    inputBinding:
      prefix: -t
    doc: Threads number
outputs:
  tgif_insertions_all:
    type: File
    outputBinding:
      glob: tgif_ncats-*/insertions_all.tsv
  tgif_insertions_filtered:
    type: File
    outputBinding:
      glob: tgif_ncats-*/insertions_filtered.tgif
  insertion_site_plots:
    type: File?
    outputBinding:
      glob: insertion_site_plots.tar
  alignment_files:
    type: File?
    outputBinding:
      glob: alignment_files.tar.gz
  script_log_file:
    type: File?
    outputBinding:
      glob: tgif_ncats-*/log
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- /TgIF/tgif_ncats.sh
stdout: tgif-ncats_stdout.log
stderr: tgif-ncats_stderr.log
doc: "TgIF (trans-gene insertion finder)\n==============================================\n\nhttps://github.com/jhuapl-bio/TgIF/tree/main\n\nThe TgIF algorithm returns a list of probable insertion sites in a target organism. It requires the user to provided a FASTQ file of ONT (Oxford Nanopore Technologies) reads (-f), the reference FASTA of the trans-gene (Tg) vector containing the insertion sequence (-i), and the reference FASTA of the target organism (-r). The algorithm is tailored for ONT reads from a modified nCATS[1] (nanopore Cas9-targeted sequencing) enriched library, however the algorithm will also produce informative results from a FASTQ derived from WGS (shotgun) sequencing libraries. The modified nCATS method is described here, and a brief overview can be found below.\n\nThe basic workflow of TgIF is alignment (using minimap2[2]) of reads (-f) to a combined reference of the Tg vector (containing the desired insertion sequence) and target organism (ie. -i and -r are concatenated), and then searching for valleys (or gaps) in the resulting pileup of reads that map to both references at MAPQ>=30. A starting position (ps) of a valley is where the depth (d) at dp=0 and dp-1>0, an ending position (pe) of a valley is where the depth at dp=0 and dp+1>0, and a potential insertion scar is the gap between and including ps and pe.\n\nPseudocode:\n\n1. Due to potential sequence homology between the vector and target genome sequences, reads are aligned to a minimap2 index containing both\n2. Alignments with MAPQ<30 are discarded, and resulting alignments are then down-selected for the best per read alignment based on highest MAPQ score\n3. Use changes in read depth as indication of insertion site:\n  - find gaps in both forward (5' to 3') and reverse (3' to 5') directions\n  - ensure each forward gap location (+strand) has a match in reverse (-strand)\n  - find rows where gap length and flanking depths are both >0\n  - these are the potential insertion sites, likely most will have flanking position depths of 1 and the gap length will be large (>5000 bp)\n4. Filter and apply confidence estimate to detected sites based on flanking depth and strandedness of aligned reads (see output format col11 for details)\n\nPARAMS:\n  HELP/OUTFMT\n    -h      help\tshow this message\n    -o\t\tformat\tformat of tsv output\n  REQUIRED\n    -t\tINT\t  number of threads to GNU parallel over\n    -f\tREADS\tsequencing reads fasta/q file run NCATS enriched library\n    -x  TYPE  sequencing read type (minimap2 param: map-ont, map-pb, sr)\n    -r\tFASTA\tfasta reference of target organism\n    -i\tFASTA\tfasta of plasmid/inserted gene(s)\n  OPTIONAL\n    -s\ty/n\t\toutput sorted bam of downselected reads for IGV use [n]\n    -p\ty/n\t\tgenerate read pileup png with R/ggplot2 per insertion site [n]\n"
label: tgif-ncats

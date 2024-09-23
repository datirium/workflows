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
      prefix: "-f"
    doc: "FASTQ file of WGS or NCATS enriched library (ideally long read, short read single-end only)"

  read_type:
    type:
    - "null"
    - type: enum
      name: "read type"
      symbols:
      - map-ont
      - map-pb
      - sr
    inputBinding:
      prefix: "-x"
    doc: "sequencing read type (minimap2 param: map-ont, map-pb, sr)"

  reference_fasta:
    type: File
    inputBinding:
      prefix: "-r"
    doc: "reference fasta file for mapping"

  plasmid_fasta:
    type: File
    inputBinding:
      prefix: "-i"
    doc: "plasmid fasta file for mapping (only a single sequence, entire sequence should be on a single line in the file, i.e. max lines in the file is 2)"

  yn_igv_output:
    type:
    - "null"
    - type: enum
      symbols:
      - y
      - n
    inputBinding:
      prefix: "-s"
    doc: "Output sorted bam of downselected reads for IGV? Default: y"

  yn_plot_output:
    type:
    - "null"
    - type: enum
      symbols:
      - y
      - n
    inputBinding:
      prefix: "-p"
    doc: "Output read pileup plot per probable insertion site? Default: y"

  threads_count:
    type: int?
    inputBinding:
      prefix: "-t"
    doc: "Threads number"


outputs:

  tgif_insertions_all:
    type: File
    outputBinding:
      glob: "tgif_ncats-*/insertions_all.tsv"

  tgif_insertions_filtered:
    type: File
    outputBinding:
      glob: "tgif_ncats-*/insertions_filtered.tgif"

  insertion_site_plots:
    type: File?
    outputBinding:
      glob: "insertion_site_plots.tar"

  alignment_files:
    type: File?
    outputBinding:
      glob: "alignment_files.tar.gz"

  script_log_file:
    type: File?
    outputBinding:
      glob: "tgif_ncats-*/log"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["/TgIF/tgif_ncats.sh"]

stdout: tgif-ncats_stdout.log
stderr: tgif-ncats_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "tgif-ncats"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/tgif-ncats.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

      
doc: |
  TgIF (trans-gene insertion finder)
  ==============================================

  https://github.com/jhuapl-bio/TgIF/tree/main

  The TgIF algorithm returns a list of probable insertion sites in a target organism. It requires the user to provided a FASTQ file of ONT (Oxford Nanopore Technologies) reads (-f), the reference FASTA of the trans-gene (Tg) vector containing the insertion sequence (-i), and the reference FASTA of the target organism (-r). The algorithm is tailored for ONT reads from a modified nCATS[1] (nanopore Cas9-targeted sequencing) enriched library, however the algorithm will also produce informative results from a FASTQ derived from WGS (shotgun) sequencing libraries. The modified nCATS method is described here, and a brief overview can be found below.

  The basic workflow of TgIF is alignment (using minimap2[2]) of reads (-f) to a combined reference of the Tg vector (containing the desired insertion sequence) and target organism (ie. -i and -r are concatenated), and then searching for valleys (or gaps) in the resulting pileup of reads that map to both references at MAPQ>=30. A starting position (ps) of a valley is where the depth (d) at dp=0 and dp-1>0, an ending position (pe) of a valley is where the depth at dp=0 and dp+1>0, and a potential insertion scar is the gap between and including ps and pe.

  Pseudocode:

  1. Due to potential sequence homology between the vector and target genome sequences, reads are aligned to a minimap2 index containing both
  2. Alignments with MAPQ<30 are discarded, and resulting alignments are then down-selected for the best per read alignment based on highest MAPQ score
  3. Use changes in read depth as indication of insertion site:
    - find gaps in both forward (5' to 3') and reverse (3' to 5') directions
    - ensure each forward gap location (+strand) has a match in reverse (-strand)
    - find rows where gap length and flanking depths are both >0
    - these are the potential insertion sites, likely most will have flanking position depths of 1 and the gap length will be large (>5000 bp)
  4. Filter and apply confidence estimate to detected sites based on flanking depth and strandedness of aligned reads (see output format col11 for details)

  PARAMS:
    HELP/OUTFMT
      -h      help	show this message
      -o		format	format of tsv output
    REQUIRED
      -t	INT	  number of threads to GNU parallel over
      -f	READS	sequencing reads fasta/q file run NCATS enriched library
      -x  TYPE  sequencing read type (minimap2 param: map-ont, map-pb, sr)
      -r	FASTA	fasta reference of target organism
      -i	FASTA	fasta of plasmid/inserted gene(s)
    OPTIONAL
      -s	y/n		output sorted bam of downselected reads for IGV use [n]
      -p	y/n		generate read pileup png with R/ggplot2 per insertion site [n]

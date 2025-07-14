cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() { if (inputs.output_filename == ""){ var ext = ".tar.gz"; var root = inputs.target_regions_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.target_regions_file.basename+ext:root+ext; } else { return inputs.output_filename; } };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.genome_fasta_file,
                  "entryname": inputs.genome_fasta_file.basename,
                  "writable": true
                }
              ]
    }
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      echo findMotifsGenome.pl $0 $1 homer_results ${@:3}
      findMotifsGenome.pl $0 $1 homer_results ${@:3}
      echo tar -czf $2 homer_results
      tar -czf $2 homer_results
    inputBinding:
      position: 5
    doc: |
      Bash script to run findMotifsGenome.pl and to compress results to tar.gz
  target_regions_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Target regions BED file to be scanned for motifs. [chr] [start] [end] [unique ID] [dummy] [strand]
  genome_fasta_file:
    type: File
    inputBinding:
      position: 7
    doc: |
      Reference genome FASTA file. Includes all chromosomes in a single file
  output_filename:
    type: string?
    inputBinding:
      position: 8
      valueFrom: $(get_output_filename())
    default: ''
    doc: |
      Name for compressed output folder to keep all the results
  background_regions_file:
    type: File?
    inputBinding:
      position: 9
      prefix: -bg
    doc: |
      Background regions BED file. [chr] [start] [end] [unique ID] [dummy] [strand]
  chopify_background_regions:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -chopify
    doc: |
      Chop up large background regions to the avg size of target regions
  search_size:
    type: string?
    inputBinding:
      position: 11
      prefix: -size
    doc: |
      Fragment size to use for motif finding.
      <#> - i.e. -size 300 will get sequences from -150 to +150 relative from center
      <#,#> - i.e. -size -100,50 will get sequences from -100 to +50 relative from center
      given - will use the exact regions you give it.
      Default=200
  motif_length:
    type: string?
    inputBinding:
      position: 12
      prefix: -len
    doc: |
      <#>[,<#>,<#>...] - motif length
      Default=8,10,12
  apply_mask_on_genome:
    type: boolean?
    inputBinding:
      position: 13
      prefix: -mask
    doc: |
      Use the repeat-masked sequence (all repeats will be masked by N)
  use_hypergeometric:
    type: boolean?
    inputBinding:
      position: 14
      prefix: -h
    doc: |
      Use hypergeometric for p-values, instead of default binomial.
      Usefull if the number of background sequences is smaller than target sequences
  skip_denovo:
    type: boolean?
    inputBinding:
      position: 15
      prefix: -nomotif
    doc: |
      Don't search for de novo motif enrichment
  skip_known:
    type: boolean?
    inputBinding:
      position: 16
      prefix: -noknown
    doc: |
      Don't search for known motif enrichment
  motifs_db:
    type:
    - 'null'
    - type: enum
      name: motifs
      symbols:
      - vertebrates
      - insects
      - worms
      - plants
      - yeast
      - all
    default: vertebrates
    inputBinding:
      position: 17
      prefix: -mset
    doc: |
      Set motifs DB to check against. Default: vertebrates
  gc_normalization:
    type: boolean?
    inputBinding:
      position: 18
      prefix: -gc
    doc: |
      Use GC% for sequence content normalization, now the default
  cpg_normalization:
    type: boolean?
    inputBinding:
      position: 19
      prefix: -cpg
    doc: |
      Use CpG% instead of GC% for sequence content normalization
  skip_normalization:
    type: boolean?
    inputBinding:
      position: 20
      prefix: -noweight
    doc: |
      Skip CG correction
  threads:
    type: int?
    inputBinding:
      position: 21
      prefix: -p
    doc: |
      Number of threads to use
outputs:
  compressed_results_folder:
    type: File
    outputBinding:
      glob: $(get_output_filename())
    doc: |
      Compressed folder with all generated results
  known_motifs:
    type: File?
    outputBinding:
      glob: homer_results/knownResults.html
    doc: |
      Known motifs html file
  denovo_motifs:
    type: File?
    outputBinding:
      glob: homer_results/homerResults.html
    doc: |
      de novo motifs html file
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: homer_find_motifs_genome_stdout.log
stderr: homer_find_motifs_genome_stderr.log
doc: |
  HOMER contains a novel motif discovery algorithm that was designed for regulatory element analysis
  in genomics applications (DNA only, no protein).  It is a differential motif discovery algorithm,
  which means that it takes two sets of sequences and tries to identify the regulatory elements that
  are specifically enriched in on set relative to the other.  It uses ZOOPS scoring (zero or one
  occurrence per sequence) coupled with the hypergeometric enrichment calculations (or binomial) to
  determine motif enrichment. HOMER also tries its best to account for sequenced bias in the dataset.
  It was designed with ChIP-Seq and promoter analysis in mind, but can be applied to pretty much any
  nucleic acids motif finding problem.

  Only selected parameters are implemented.
label: homer-find-motifs

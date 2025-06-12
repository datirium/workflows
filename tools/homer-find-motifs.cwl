cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() { if (inputs.output_filename == ""){ var ext = ".tar.gz"; var root = inputs.target_fasta_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.target_fasta_file.basename+ext:root+ext; } else { return inputs.output_filename; } };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.target_fasta_file,
                  "entryname": "target.fa",
                  "writable": true
                },
                {
                  "entry": inputs.background_fasta_file,
                  "entryname": "background.fa",
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
      ls -la
      echo findMotifs.pl $0 dummy homer_results ${@:2}
      echo samtools faidx $0
      echo samtools faidx $3
      samtools faidx $0
      samtools faidx $3
      findMotifs.pl $0 dummy homer_results ${@:2}
      echo tar -czf $1 homer_results
      tar -czf $1 homer_results
    inputBinding:
      position: 5
    doc: |
      Bash script to run findMotifs.pl and to compress results to tar.gz
  target_fasta_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Target FASTA file to scan for motifs
  output_filename:
    type: string?
    inputBinding:
      position: 7
      valueFrom: $(get_output_filename())
    default: ''
    doc: |
      Name for compressed output folder to keep all the results
  background_fasta_file:
    type: File
    inputBinding:
      position: 8
      prefix: -fasta
    doc: |
      Background FASTA file suitable for use as a null distribution
  skip_denovo:
    type: boolean?
    inputBinding:
      position: 9
      prefix: -nomotif
    doc: |
      Don't search for de novo motif enrichment
  skip_known:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -noknown
    doc: |
      Don't search for known motif enrichment
  use_binomial:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -b
    doc: |
      Use binomial distribution to calculate p-values (default is hypergeometric)
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
      position: 12
      prefix: -mset
    doc: |
      Set motifs DB to check against. Default: vertebrates
  threads:
    type: int?
    inputBinding:
      position: 13
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
stdout: homer_find_motifs_stdout.log
stderr: homer_find_motifs_stderr.log
doc: |
  HOMER contains a novel motif discovery algorithm that was designed for regulatory element analysis
  in genomics applications (DNA only, no protein).  It is a differential motif discovery algorithm,
  which means that it takes two sets of sequences and tries to identify the regulatory elements that
  are specifically enriched in on set relative to the other.  It uses ZOOPS scoring (zero or one
  occurrence per sequence) coupled with the hypergeometric enrichment calculations (or binomial) to
  determine motif enrichment.  HOMER also tries its best to account for sequenced bias in the dataset.
  It was designed with ChIP-Seq and promoter analysis in mind, but can be applied to pretty much any
  nucleic acids motif finding problem.

  Only selected parameters are implemented.
label: homer-find-motifs

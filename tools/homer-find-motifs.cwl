cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() {
          if (inputs.output_filename == ""){
            var ext = ".tar.gz";
            var root = inputs.target_fasta_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.target_fasta_file.basename+ext:root+ext;
          } else {
            return inputs.output_filename;
          }
        };
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
    default: ""
    doc: |
      Name for compressed output folder to keep all the results

  background_fasta_file:
    type: File
    inputBinding:
      position: 8
      prefix: "-fasta"
    doc: |
      Background FASTA file suitable for use as a null distribution

  skip_denovo:
    type: boolean?
    inputBinding:
      position: 9
      prefix: "-nomotif"
    doc: |
      Don't search for de novo motif enrichment
      
  skip_known:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "-noknown"
    doc: |
      Don't search for known motif enrichment

  use_binomial:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "-b"
    doc: |
      Use binomial distribution to calculate p-values (default is hypergeometric)

  motifs_db:
    type:
      - "null"
      - type: enum
        name: "motifs"
        symbols: ["vertebrates", "insects", "worms", "plants", "yeast", "all"]
    default: "vertebrates"
    inputBinding:
      position: 12
      prefix: "-mset"
    doc: |
      Set motifs DB to check against. Default: vertebrates

  threads:
    type: int?
    inputBinding:
      position: 13
      prefix: "-p"
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
      glob: "homer_results/knownResults.html"
    doc: |
      Known motifs html file

  denovo_motifs:
    type: File?
    outputBinding:
      glob: "homer_results/homerResults.html"
    doc: |
      de novo motifs html file

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]


stdout: homer_find_motifs_stdout.log
stderr: homer_find_motifs_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/homer-metadata.yaml

s:name: "homer-find-motifs"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/homer-find-motifs.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


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

  
s:about: |
  Usage:  findMotifs.pl <input list> <promoter set> <output directory> [additoinal options]

  example: findMotifs.pl genelist.txt mouse motifResults/ -len 10

  FASTA example: findMotifs.pl targets.fa fasta motifResults/ -fasta background.fa

  Available Promoter Sets: Add custom promoters sets with loadPromoters.pl

  Try typing "perl /opt/homer/.//configureHomer.pl -list" to see available promoter sets
  Typing "perl /opt/homer/.//configureHomer.pl -install NNN" to install promoter set NNN

  Basic options:
  -len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the program
    to run out of memmory - in these cases decrease the number of sequences analyzed]
  -bg <background file> (ids to use as background, default: all genes)
  -start <#> (offset from TSS, default=-300) [max=based on Promoter Set]
  -end <#> (offset from TSS, default=50) [max=based on Promoter Set]
  -rna (output RNA motif logos and compare to RNA motif database, automatically sets -norevopp)
  -mask/-nomask (use/don't use repeatmasked files, default: -mask)
  -S <#> (Number of motifs to optimize, default: 25)
  -mis <#> (global optimization: searches for strings with # mismatches, default: 1)
  -noconvert (will not worry about converting input files into unigene ids)
  -norevopp (do not search the reverse strand for motifs)
  -nomotif (don't search for de novo motif enrichment)

  Scanning sequence for motifs
  -find <motif file> (This will cause the program to only scan for motifs)

  Including Enhancers - peak files of enhancer location, peak ID should be gene ID
  -enhancers <peak file> <genome verion>
    (enhancers to include in search space, peaks/sequences should be named with a gene ID
    If multiple enhancers per gene, use the same gene ID, and all will be included)
  -enhancersOnly (do not include promoter sequence in motif search)

  FASTA files: If you prefer to use your own fasta files, place target sequences and 
  background sequences in two separate FASTA formated files (must have unique identifiers)
  Target File - use in place of <input list> (i.e. the first argument)
  Background File - after output directory (with additional options) use the argument:
    -fastaBg <background fasta file> (This is recommended for fasta based analysis)
  In place of the promoter set use "fasta", or any valid set (this parameter is ignored)
  When finding motifs [-find], only the target file with be searched)
    -chopify (chops up background regions to match size of target regions)
      i.e. if background is a full genome or all mRNAs

  Known Motif Options/Visualization:
  -mset <vertebrates|insects|worms|plants|yeast|all> (check against motif collects, default: auto)
  -basic (don't check de novo motifs for similarity to known motifs)
  -bits (scale sequence logos by information content, default: doesn't scale)
  -nocheck (don't check for similarity between novo motif motifs and known motifs)
  -mcheck <motif file> (known motifs to check against de novo motifs,
  -noknown (don't search for known motif enrichment, default: -known)
  -mknown <motif file> (known motifs to check for enrichment,
  -nofacts (omit humor)
  -seqlogo (uses weblogo/seqlogo/ghostscript to visualize motifs, default uses SVG)

  Advanced options:
  -b (use binomial distribution to calculate p-values, hypergeometric is default)
  -nogo (don't search for gene ontology enrichment)
  -humanGO (Convert IDs to human for GO analysis)
  -ontology <ont.genes> [ont.genes] ... (custom ontologies for GO analysis)
  -noweight (no CG correction)
  -noredun (Don't remove predetermined redundant promoters/sequences)
  -g (input file is a group file, i.e. 1st column = id, 2nd = 0 or 1 [1=target,0=back])
  -cpg (use CpG% instead of GC% for sequence normalization)
  -rand (randomize labels for target and backgound sequences)
  -maskMotif <motif file 1> [motif file 2] ... (motifs to mask before motif finding)
  -opt <motif file 1> [motif file 2] ... (motifs to optimize/change length)
  -peaks (will produce peak file of promoters to use with findMotifsGenome.pl)
  -nowarn (no warnings)
  -keepFiles (don't delete temporary files)
  -dumpFasta (create target.fa and background.fa files)
  -min <#> (remove sequences shorter than #, default: 0)
  -max <#> (remove sequences longer than #, default: 1e10)
  -reuse (rerun homer using old seq files etc. with new options
      and ignores input list, organism)
  -fdr <#> (Calculate empirical FDR for de novo discovery #=number of randomizations)

  homer2 specific options:
  -homer2 (use homer2 instead of original homer, default)
  -nlen <#> (length of lower-order oligos to normalize - general sequences, default: 3)
    -nmax <#> (Max normalization iterations, default: 160)
    -neutral (weight sequences to neutral frequencies, i.e. 25%, 6.25%, etc.)
  -olen <#> (lower-order oligo normalization for oligo table, use if -nlen isn't working well)
  -p <#> (Number of processors to use, default: 1)
  -e <#> (Maximum expected motif instance per bp in random sequence, default: 0.01)
  -cache <#> (size in MB for statistics cache, default: 500)
  -quickMask (skip full masking after finding motifs, similar to original homer)
  -homer1 (to force the use of the original homer)
  -minlp <#> (stop looking for motifs when seed logp score gets above #, default: -10)

  Original homer specific options:
  -homer1 (to force the use of the original homer)
  -depth [low|med|high|allnight] (time spent on local optimization default: med)
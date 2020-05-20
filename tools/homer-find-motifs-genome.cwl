cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() {
          if (inputs.output_filename == ""){
            var ext = ".tar.gz";
            var root = inputs.target_regions_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.target_regions_file.basename+ext:root+ext;
          } else {
            return inputs.output_filename;
          }
        };
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
    default: ""
    doc: |
      Name for compressed output folder to keep all the results

  background_regions_file:
    type: File?
    inputBinding:
      position: 9
      prefix: "-bg"
    doc: |
      Background regions BED file. [chr] [start] [end] [unique ID] [dummy] [strand]

  chopify_background_regions:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "-chopify"
    doc: |
      Chop up large background regions to the avg size of target regions

  search_size:
    type: string?
    inputBinding:
      position: 11
      prefix: "-size"
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
      prefix: "-len"
    doc: |
      <#>[,<#>,<#>...] - motif length
      Default=8,10,12

  apply_mask_on_genome:
    type: boolean?
    inputBinding:
      position: 13
      prefix: "-mask"
    doc: |
      Use the repeat-masked sequence (all repeats will be masked by N)

  use_hypergeometric:
    type: boolean?
    inputBinding:
      position: 14
      prefix: "-h"
    doc: |
      Use hypergeometric for p-values, instead of default binomial.
      Usefull if the number of background sequences is smaller than target sequences

  skip_denovo:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "-nomotif"
    doc: |
      Don't search for de novo motif enrichment

  skip_known:
    type: boolean?
    inputBinding:
      position: 16
      prefix: "-noknown"
    doc: |
      Don't search for known motif enrichment

  motifs_db:
    type:
      - "null"
      - type: enum
        name: "motifs"
        symbols: ["vertebrates", "insects", "worms", "plants", "yeast", "all"]
    default: "vertebrates"
    inputBinding:
      position: 17
      prefix: "-mset"
    doc: |
      Set motifs DB to check against. Default: vertebrates

  threads:
    type: int?
    inputBinding:
      position: 18
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


stdout: homer_find_motifs_genome_stdout.log
stderr: homer_find_motifs_genome_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/homer-metadata.yaml

s:name: "homer-find-motifs"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/homer-find-motifs-genome.cwl
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
  determine motif enrichment. HOMER also tries its best to account for sequenced bias in the dataset.
  It was designed with ChIP-Seq and promoter analysis in mind, but can be applied to pretty much any
  nucleic acids motif finding problem.
  
  Only selected parameters are implemented.

  
s:about: |
  Usage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]
  Example: findMotifsGenome.pl peaks.txt mm8r peakAnalysis -size 200 -len 8

  Possible Genomes:
      -- or --
    Custom: provide the path to genome FASTA files (directory or single file)
      Heads up: will create the directory "preparsed/" in same location.

  Basic options:
    -mask (mask repeats/lower case sequence, can also add 'r' to genome, i.e. mm9r)
    -bg <background position file> (genomic positions to be used as background, default=automatic)
      removes background positions overlapping with target positions
      -chopify (chop up large background regions to the avg size of target regions)
    -len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the program
      to run out of memory - in these cases decrease the number of sequences analyzed (-N),
      or try analyzing shorter sequence regions (i.e. -size 100)]
    -size <#> (fragment size to use for motif finding, default=200)
      -size <#,#> (i.e. -size -100,50 will get sequences from -100 to +50 relative from center)
      -size given (uses the exact regions you give it)
    -S <#> (Number of motifs to optimize, default: 25)
    -mis <#> (global optimization: searches for strings with # mismatches, default: 2)
    -norevopp (don't search reverse strand for motifs)
    -nomotif (don't search for de novo motif enrichment)
    -rna (output RNA motif logos and compare to RNA motif database, automatically sets -norevopp)

  Scanning sequence for motifs
    -find <motif file> (This will cause the program to only scan for motifs)

  Known Motif Options/Visualization
    -mset <vertebrates|insects|worms|plants|yeast|all> (check against motif collects, default: auto)
    -basic (just visualize de novo motifs, don't check similarity with known motifs)
    -bits (scale sequence logos by information content, default: doesn't scale)
    -nocheck (don't search for de novo vs. known motif similarity)
    -mcheck <motif file> (known motifs to check against de novo motifs,
    -noknown (don't search for known motif enrichment, default: -known)
    -mknown <motif file> (known motifs to check for enrichment,
    -nofacts (omit humor)
    -seqlogo (use weblogo/seqlogo/ghostscript to generate logos, default uses SVG now)

  Sequence normalization options:
    -gc (use GC% for sequence content normalization, now the default)
    -cpg (use CpG% instead of GC% for sequence content normalization)
    -noweight (no CG correction)
    Also -nlen <#>, -olen <#>, see homer2 section below.

  Advanced options:
    -h (use hypergeometric for p-values, binomial is default)
    -N <#> (Number of sequences to use for motif finding, default=max(50k, 2x input)
    -local <#> (use local background, # of equal size regions around peaks to use i.e. 2)
    -redundant <#> (Remove redundant sequences matching greater than # percent, i.e. -redundant 0.5)
    -maxN <#> (maximum percentage of N's in sequence to consider for motif finding, default: 0.7)
    -maskMotif <motif file1> [motif file 2]... (motifs to mask before motif finding)
    -opt <motif file1> [motif file 2]... (motifs to optimize or change length of)
    -rand (randomize target and background sequences labels)
    -ref <peak file> (use file for target and background - first argument is list of peak ids for targets)
    -oligo (perform analysis of individual oligo enrichment)
    -dumpFasta (Dump fasta files for target and background sequences for use with other programs)
    -preparse (force new background files to be created)
    -preparsedDir <directory> (location to search for preparsed file and/or place new files)
    -keepFiles (keep temporary files)
    -fdr <#> (Calculate empirical FDR for de novo discovery #=number of randomizations)

  homer2 specific options:
    -homer2 (use homer2 instead of original homer, default)
    -nlen <#> (length of lower-order oligos to normalize in background, default: -nlen 3)
      -nmax <#> (Max normalization iterations, default: 160)
      -neutral (weight sequences to neutral frequencies, i.e. 25%, 6.25%, etc.)
    -olen <#> (lower-order oligo normalization for oligo table, use if -nlen isn't working well)
    -p <#> (Number of processors to use, default: 1)
    -e <#> (Maximum expected motif instance per bp in random sequence, default: 0.01)
    -cache <#> (size in MB for statistics cache, default: 500)
    -quickMask (skip full masking after finding motifs, similar to original homer)
    -minlp <#> (stop looking for motifs when seed logp score gets above #, default: -10)

  Original homer specific options:
    -homer1 (to force the use of the original homer)
    -depth [low|med|high|allnight] (time spent on local optimization default: med)

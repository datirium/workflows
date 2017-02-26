#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/bowtie:v1.1.2
  dockerFile: >
    $import: ./dockerfiles/bowtie-Dockerfile

inputs:

  reference_in:
    type:
      - File
      - type: array
        items: File
    doc: >
      comma-separated list of files with ref sequences
    inputBinding:
      itemSeparator: ","
      position: 25

  index_base:
    type: string?
    doc: |
      write Ebwt data to files with this dir/basename
    inputBinding:
      position: 26
      shellQuote: false
    default: bowtie_indices

  f:
    type:
      - "null"
      - boolean
    doc: |
      reference files are Fasta (default)
    inputBinding:
      position: 1
      prefix: '-f'

  c:
    type:
      - "null"
      - boolean
    doc: |
      reference sequences given on cmd line (as <seq_in>)
    inputBinding:
      position: 2
      prefix: '-c'

  large_index:
    type:
      - "null"
      - boolean
    doc: |
      force generated index to be 'large', even if ref has fewer than 4 billion nucleotides
    inputBinding:
      position: 3
      prefix: '--large-index'

  color:
    type:
      - "null"
      - boolean
    doc: |
      build a colorspace index
    inputBinding:
      position: 4
      prefix: '--color'

  noauto:
    type:
      - "null"
      - boolean
    doc: |
      disable automatic -p/--bmax/--dcv memory-fitting
    inputBinding:
      position: 5
      prefix: '--noauto'

  packed:
    type:
      - "null"
      - boolean
    doc: |
      use packed strings internally; slower, less memory
    inputBinding:
      position: 6
      prefix: '--packed'

  bmax:
    type:
      - "null"
      - int
    doc: |
      max bucket sz for blockwise suffix-array builder
    inputBinding:
      position: 7
      prefix: '--bmax'

  bmaxdivn:
    type:
      - "null"
      - int
    doc: |
      max bucket sz as divisor of ref len (default: 4)
    inputBinding:
      position: 8
      prefix: '--bmaxdivn'

  dcv:
    type:
      - "null"
      - int
    doc: |
      diff-cover period for blockwise (default: 1024)
    inputBinding:
      position: 9
      prefix: '--dcv'

  nodc:
    type:
      - "null"
      - boolean
    doc: |
      disable diff-cover (algorithm becomes quadratic)
    inputBinding:
      position: 10
      prefix: '--nodc'

  noref:
    type:
      - "null"
      - boolean
    doc: |
      don't build .3/.4 index files
    inputBinding:
      position: 11
      prefix: '--noref'

  justref:
    type:
      - "null"
      - boolean
    doc: |
      just build .3/.4 index files
    inputBinding:
      position: 12
      prefix: '--justref'

  offrate:
    type:
      - "null"
      - int
    doc: |
      SA is sampled every 2^<int> BWT chars (default: 5)
    inputBinding:
      position: 13
      prefix: '--offrate'

  ftabchars:
    type:
      - "null"
      - int
    doc: |
      # of chars consumed in initial lookup (default: 10)
    inputBinding:
      position: 14
      prefix: '--ftabchars'

  ntoa:
    type:
      - "null"
      - boolean
    doc: |
      convert Ns in reference to As
    inputBinding:
      position: 15
      prefix: '--ntoa'

  seed:
    type:
      - "null"
      - int
    doc: |
      seed for random number generator
    inputBinding:
      position: 16
      prefix: '--seed'

  quiet:
    type:
      - "null"
      - boolean
    doc: |
      verbose output (for debugging)
    inputBinding:
      position: 17
      prefix: '--quiet'

outputs:

  indices:
    type: File
    outputBinding:
      glob: $(inputs.index_base + ".1.ebwt*")
    secondaryFiles: |
      ${
        var ext = self.location.split('/').slice(-1)[0].split('.').slice(-1)[0];
        var basename = self.location.split("/").slice(-1)[0].split(".").slice(0, -2).join (".");
        var dirname = self.location.split("/").slice(0,-1).join("/");
        return [{"class": "File", "location": dirname + "/" + basename + ".2." + ext},
                {"class": "File", "location": dirname + "/" + basename + ".3." + ext},
                {"class": "File", "location": dirname + "/" + basename + ".4." + ext},
                {"class": "File", "location": dirname + "/" + basename + ".rev.1." + ext},
                {"class": "File", "location": dirname + "/" + basename + ".rev.2." + ext}
        ]
      }

  output_log:
    type: File
    outputBinding:
      glob: $(inputs.index_base + ".log")

baseCommand:
  - bowtie-build

arguments:
  - valueFrom: $('> ' + inputs.index_base + '.log')
    position: 100000
    shellQuote: false


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bowtie2-metadata.yaml

s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/tools/bowtie-build.cwl
s:codeRepository: https://github.com/SciDAP/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: Michael Kotliar
  s:email: mailto:michael.kotliar@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0002-6486-3898
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab

doc: |
  Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>
      reference_in            comma-separated list of files with ref sequences
      ebwt_outfile_base       write Ebwt data to files with this dir/basename
  Options:
      -f                      reference files are Fasta (default)
      -c                      reference sequences given on cmd line (as <seq_in>)
      --large-index           force generated index to be 'large', even if ref
                              has fewer than 4 billion nucleotides
      -C/--color              build a colorspace index
      -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
      -p/--packed             use packed strings internally; slower, uses less mem
      --bmax <int>            max bucket sz for blockwise suffix-array builder
      --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
      --dcv <int>             diff-cover period for blockwise (default: 1024)
      --nodc                  disable diff-cover (algorithm becomes quadratic)
      -r/--noref              don't build .3/.4.ebwt (packed reference) portion
      -3/--justref            just build .3/.4.ebwt (packed reference) portion
      -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)
      -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
      --ntoa                  convert Ns in reference to As
      --seed <int>            seed for random number generator
      -q/--quiet              verbose output (for debugging)
      -h/--help               print detailed description of tool and its options
      --usage                 print this usage message
      --version               print version information and quit


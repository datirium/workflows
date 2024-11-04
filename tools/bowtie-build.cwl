cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ResourceRequirement
    ramMin: 15250
    coresMin: 1
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: |
      ${
        return [
          {
            "class": "Directory",
            "basename": inputs.index_base_name,
            "listing": [],
            "writable": true}
        ]
      }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bowtie:v1.2.0


inputs:

  fasta_file:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      itemSeparator: ","
      position: 25
    doc: |
      comma-separated list of files with ref sequences

  index_base_name:
    type: string
    inputBinding:
      position: 26
      valueFrom: $("./"+self+"/"+self)
    doc: |
      write Ebwt data to files with this dir/basename

  force_large_index:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 3
      prefix: '--large-index'
    doc: |
      force generated index to be 'large', even if ref has fewer than 4 billion nucleotides

  color:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 4
      prefix: '--color'
    doc: |
      build a colorspace index

  noauto:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 5
      prefix: '--noauto'
    doc: |
      disable automatic -p/--bmax/--dcv memory-fitting

  packed:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 6
      prefix: '--packed'
    doc: |
      use packed strings internally; slower, less memory

  bmax:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: '--bmax'
    doc: |
      max bucket sz for blockwise suffix-array builder

  bmaxdivn:
    type:
      - "null"
      - int
    inputBinding:
      position: 8
      prefix: '--bmaxdivn'
    doc: |
      max bucket sz as divisor of ref len (default: 4)

  dcv:
    type:
      - "null"
      - int
    inputBinding:
      position: 9
      prefix: '--dcv'
    doc: |
      diff-cover period for blockwise (default: 1024)

  nodc:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 10
      prefix: '--nodc'
    doc: |
      disable diff-cover (algorithm becomes quadratic)

  noref:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 11
      prefix: '--noref'
    doc: |
      don't build .3/.4 index files

  justref:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 12
      prefix: '--justref'
    doc: |
      just build .3/.4 index files

  offrate:
    type:
      - "null"
      - int
    inputBinding:
      position: 13
      prefix: '--offrate'
    doc: |
      SA is sampled every 2^<int> BWT chars (default: 5)

  ftabchars:
    type:
      - "null"
      - int
    inputBinding:
      position: 14
      prefix: '--ftabchars'
    doc: |
      # of chars consumed in initial lookup (default: 10)

  ntoa:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 15
      prefix: '--ntoa'
    doc: |
      convert Ns in reference to As

  seed:
    type:
      - "null"
      - int
    inputBinding:
      position: 16
      prefix: '--seed'
    doc: |
      seed for random number generator

  quiet:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 17
      prefix: '--quiet'
    doc: |
      verbose output (for debugging)


outputs:

  indices_folder:
    type: Directory
    outputBinding:
      glob: $(inputs.index_base_name)

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bowtie-build"]
stderr: bowtie_stderr.log
stdout: bowtie_stdout.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/bowtie-metadata.yaml

s:name: "bowtie-build"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bowtie-build.cwl
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
  Tool runs bowtie-build
  Not supported parameters:
    -c  -  reference sequences given on cmd line (as <seq_in>)

s:about: |
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

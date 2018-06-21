cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_prefix = function() {
        let ext = '.';
        let root = inputs.bam_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.bam_file.basename+ext:root+ext;
    };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.bam_file,
                  "entryname": inputs.bam_file.basename,
                  "writable": true
                }
              ]
    }

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/geep:v0.0.3


inputs:

  bam_file:
    type:
      - File
    inputBinding:
      position: 10
      prefix: --bam
    doc: |
      Set the path to the coordinate sorted BAM file. Required

  annotation_file:
    type:
      - File
    inputBinding:
      position: 11
      prefix: --annotation
    doc: |
      Set the path to the GTF or TAB-delimited file. Required

  log_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 12
      prefix: --log
    doc: |
      Set the path to save LOG file. Default: /dev/null

  output_prefix:
    type:
      - "null"
      - string
    inputBinding:
      position: 13
      prefix: --output
      valueFrom: |
        ${
            if (self == null){
              return default_output_prefix();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Set the prefix to save output files. Default: ""

  rpkm_threshold:
    type:
      - "null"
      - double
    inputBinding:
      position: 14
      prefix: --threshold
    doc: |
      Set rpkm cutoff threshold, below which everything will be changed to value set with --cutoff. Default: 0

  rpkm_cutoff:
    type:
      - "null"
      - double
    inputBinding:
      position: 15
      prefix: --cutoff
    doc: |
      Set rpkm cutoff value to be used for all rpkms below --threshold. Default: 0

  exclude_chr:
    type:
      - "null"
      - string
    inputBinding:
      position: 16
      prefix: --exclude
    doc: |
      Coma separated list of chromosomes to be ignored. Default: ""

  min_interval_length:
    type:
      - "null"
      - int
    inputBinding:
      position: 17
      prefix: --minIntLen
    doc: |
      Set the minimal interval length. All shorter intervals will be discarded. Default: 0

  min_read_length:
    type:
      - "null"
      - int
    inputBinding:
      position: 18
      prefix: --minReadLen
    doc: |
      Set the minimal read length. All parts of spliced reads that intersect with exon in
      less than minReadLen nucleotides will be discarded. Default: 0

  keep_unique:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 19
      prefix: --keepUnique
    doc: |
      Set this flag if you want prevent distributing the isoform unique reads among other isoforms. Default: False

  dutp:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 20
      prefix: --dutp
    doc: |
      Set this dutp flag if strand specific analysys should be made. Default: False

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 21
      prefix: --threads
    doc: |
      Set the number of threads. Default: 1


outputs:

  isoforms_file:
    type: File
    outputBinding:
      glob: "*isoforms.csv"

  genes_file:
    type: File
    outputBinding:
      glob: "*genes.csv"

  raw_file:
    type: File
    outputBinding:
      glob: "*raw.txt"

  log_file:
    type:
     - "null"
     - File
    outputBinding:
      glob: $(inputs.log_filename)

baseCommand: [geep]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/geep-metadata.yaml

s:name: "geep"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/geep.cwl
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
  Tool calculates RPKM values grouped by isoforms or genes.

  `default_output_prefix` function returns default prefix based on `bam_file` basename, if `output_prefix` is not
  provided.

  Before running `baseCommand` `bam_file` is staged into output directory with write permissions (`"writable": true`).
  This allow to automatically generate index file at the same directory as input `bam_file`. In case when index file
  is provided in `secondaryFiles` of `bam_file`, it's not generated twice.

s:about: |
  Usage:
    geep [params]
      -b,--bam         Set the path to the BAM file. Required
      -g,--annotation  Set the path to the GTF or TAB-delimited file. Required
      -l,--log         Set the path to save LOG file. Default: /dev/null
      -o,--output      Set the prefix to save output files. Default: ""
      -t,--threshold   Set rpkm cutoff threshold, below which everything will be changed to value set with --cutoff. Default: 0
      -c,--cutoff      Set rpkm cutoff value to be used for all rpkms below --threshold. Default: 0
      -x,--exclude     Coma separated list of chromosomes to be ignored. Default: ""
      -i,--minIntLen   Set the minimal interval length. All shorter intervals will be discarded. Default: 0
      -r,--minReadLen  Set the minimal read length. All parts of spliced reads that intersect with exon in less than minReadLen nucleotides will be discarded. Default: 0
      -p,--threads     Set the number of threads. Default: 1
      -u,--keepUnique  Set this flag if you want prevent distributing the isoform unique reads among other isoforms. Default: False
      -d,--dutp        Set this dutp flag if strand specific analysys should be made. Default: False
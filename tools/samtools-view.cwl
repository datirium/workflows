cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      if (inputs.iscram && !inputs.isbam){
        return '.cram';
      } else if (inputs.isbam) {
        return '.bam';
      } else {
        return '.sam';
      }
    };
  - var default_output_name = function() {
        return inputs.view_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+ext();
    }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  isbam:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -b
    doc: |
      output in BAM format

  iscram:
    type: boolean?
    inputBinding:
      position: 6
      prefix: -C
    doc: |
      output in CRAM format

  fastcompression:
    type: boolean?
    inputBinding:
      position: 7
      prefix: '-1'
    doc: |
      use fast BAM compression (implies -b)

  uncompressed:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -u
    doc: |
      uncompressed BAM output (implies -b)

  samheader:
    type: boolean?
    inputBinding:
      position: 9
      prefix: -h
    doc: |
      include header in the output

  header_only:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -H
    doc: |
      include header in the output

  count:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -c
    doc: |
      Instead of printing the alignments, only count them and print the total number.
      All filter options, such as -f, -F, and -q, are taken into account.

  output_name:
    type: string
    inputBinding:
      position: 12
      prefix: -o
      valueFrom: |
        ${
            if (self == null){
              return default_output_name();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Output to file

  filtered_out:
    type: string?
    inputBinding:
      position: 13
      prefix: -U
    doc: |
      Write alignments that are not selected by the various filter options to FILE.
      When this option is used, all alignments (or all alignments intersecting the
      regions specified) are written to either the output file or this file, but never both.

  chr_length:
    type: File?
    inputBinding:
      position: 14
      prefix: -t
    doc: |
      A tab-delimited FILE. Each line must contain the reference name in the first column and the length
      of the reference in the second column, with one line for each distinct reference.
      Any additional fields beyond the second column are ignored. This file also defines the order of
      the reference sequences in sorting. If you run: `samtools faidx <ref.fa>',
      the resulting index file <ref.fa>.fai can be used as this FILE.

  reference_fasta:
    type: File?
    inputBinding:
      position: 15
      prefix: -T
    doc: |
      A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed by samtools faidx.
      If an index is not present, one will be generated for you.

  bed_overlap:
    type: File?
    inputBinding:
      position: 16
      prefix: -L
    doc: |
      only include reads overlapping this BED FILE

  reads_in_group:
    type: string?
    inputBinding:
      position: 17
      prefix: -r
    doc: |
      Only output alignments in read group STR [null].

  readsquality:
    type: int?
    inputBinding:
      position: 18
      prefix: -q
    doc: |
      Skip alignments with MAPQ smaller than INT [0].

  reads_in_library:
    type: string?
    inputBinding:
      position: 19
      prefix: -l
    doc: |
      Only output alignments in library STR [null].

  cigar:
    type: int?
    inputBinding:
      position: 20
      prefix: -m
    doc: |
      Only output alignments with number of CIGAR bases consuming query sequence ≥ INT [0]

  reads_with_bits:
    type: int?
    inputBinding:
      position: 21
      prefix: -f
    doc: |
      only include reads with all bits set in INT set in FLAG [0]

  reads_without_bits:
    type: int?
    inputBinding:
      position: 22
      prefix: -F
    doc: |
      only include reads with none of the bits set in INT set in FLAG [0]

  reads_tag:
    type: string?
    inputBinding:
      position: 23
      prefix: -x
    doc: |
      Read tag to exclude from output (repeatable) [null]

  collapse_cigar:
    type: boolean?
    inputBinding:
      position: 24
      prefix: -B
    doc: |
      collapse the backward CIGAR operation

  random_seed:
    type: float?
    inputBinding:
      position: 25
      prefix: -s
    doc: |
      Integer part is used to seed the random number generator [0]. Part after the decimal
      point sets the fraction of templates/pairs to subsample [no subsampling].

  threads:
    type: int?
    inputBinding:
      position: 26
      prefix: -@
    doc: |
      number of BAM compression threads [0]

  reads_in_group_file:
    type: File?
    inputBinding:
      position: 27
      prefix: -R
    doc: |
      only include reads with read group listed in FILE [null]

  view_input:
    type: File
    inputBinding:
      position: 30
    doc: |
      Input SAM, BAM, or CRAM file.

  region:
    type: string?
    inputBinding:
      position: 31
    doc: |
      RNAME[:STARTPOS[-ENDPOS]] and all position coordinates are 1-based.

outputs:
  view_file:
    type: File
    outputBinding:
      glob: $(default_output_name())

  filtered_file:
    type: File?
    outputBinding:
      glob: $(inputs.filtered_out)

baseCommand: [samtools, view]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/tools/samtools-view.cwl
s:codeRepository: https://github.com/SciDAP/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: Andrey Kartashov
  s:email: mailto:Andrey.Kartashov@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab
doc: |
  samtools-view.cwl is developed for CWL consortium
  Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]

  Options:
    -b       output BAM
    -C       output CRAM (requires -T)
    -1       use fast BAM compression (implies -b)
    -u       uncompressed BAM output (implies -b)
    -h       include header in SAM output
    -H       print SAM header only (no alignments)
    -c       print only the count of matching records
    -o FILE  output file name [stdout]
    -U FILE  output reads not selected by filters to FILE [null]
    -t FILE  FILE listing reference names and lengths (see long help) [null]
    -L FILE  only include reads overlapping this BED FILE [null]
    -r STR   only include reads in read group STR [null]
    -R FILE  only include reads with read group listed in FILE [null]
    -q INT   only include reads with mapping quality >= INT [0]
    -l STR   only include reads in library STR [null]
    -m INT   only include reads with number of CIGAR operations consuming
             query sequence >= INT [0]
    -f INT   only include reads with all bits set in INT set in FLAG [0]
    -F INT   only include reads with none of the bits set in INT set in FLAG [0]
    -s FLOAT subsample reads (given INT.FRAC option value, 0.FRAC is the
             fraction of templates/read pairs to keep; INT part sets seed)
    -x STR   read tag to strip (repeatable) [null]
    -B       collapse the backward CIGAR operation
    -?       print long help, including note about region specification
    -S       ignored (input format is auto-detected)
        --input-fmt-option OPT[=VAL]
                 Specify a single input file format option in the form
                 of OPTION or OPTION=VALUE
    -O, --output-fmt FORMAT[,OPT[=VAL]]...
                 Specify output format (SAM, BAM, CRAM)
        --output-fmt-option OPT[=VAL]
                 Specify a single output file format option in the form
                 of OPTION or OPTION=VALUE
    -T, --reference FILE
                 Reference sequence FASTA FILE [null]
    -@, --threads INT
                 Number of additional threads to use [0]

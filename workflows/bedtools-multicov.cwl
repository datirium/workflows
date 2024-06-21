cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  sample:
  - "chipseq-se.cwl"
  - "chipseq-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  - "rnaseq-se.cwl"
  - "rnaseq-pe.cwl"
  - "rnaseq-se-dutp.cwl"
  - "rnaseq-pe-dutp.cwl"
  - "trim-rnaseq-pe.cwl"
  - "trim-rnaseq-se.cwl"
  - "trim-rnaseq-pe-dutp.cwl"
  - "trim-rnaseq-se-dutp.cwl"
  - "trim-rnaseq-pe-smarter-dutp.cwl"
  - "trim-quantseq-mrnaseq-se-strand-specific.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  alignment_files:
    type: File[]
    format: "http://edamontology.org/format_2572"
    secondaryFiles:
    - .bai
    label: "Sample to produce alignment file"
    doc: "Coordinate sorted and indexed alignment BAM file(s)"
    'sd:upstreamSource': "sample/bambai_pair"
    'sd:localLabel': true

  alignment_names:
    type: string[]
    label: "Sample to produce alignment file"
    doc: "Names for input alignment files. Order corresponds to the alignment_files"
    'sd:upstreamSource': "sample/alias"

  intervals_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Regions of interest, headerless BED file"
    doc: "Intervals file defined in a BED/GFF/VCF format"

  split:
    type: boolean?
    default: false
    label: "Treat 'split' BAM or BED12 entries as distinct BED intervals"
    doc: "Treat 'split' BAM or BED12 entries as distinct BED intervals"
    'sd:layout':
      advanced: true

  same_strand:
    type: boolean?
    default: false
    label: "Report only hits that overlap the same strand"
    doc: |
      Require same strandedness. That is, only report hits in B
      that overlap A on the _same_ strand. By default, overlaps
      are reported without respect to strand.
    'sd:layout':
      advanced: true

  opposite_strand:
    type: boolean?
    default: false
    label: "Report only hits that overlap the opposite strand"
    doc: |
      Require different strandedness. That is, only report hits in B
      that overlap A on the _opposite_ strand. By default, overlaps
      are reported without respect to strand.
    'sd:layout':
      advanced: true

  min_overlap_fraction:
    type: float?
    default: 0.000000001
    label: "Minimum overlap required as a fraction. Default value correposponds to 1bp"
    doc: |
      Minimum overlap required as a fraction of each A.
      Default is 1E-9 (i.e., 1bp).
    'sd:layout':
      advanced: true

  reciprocal_overlap:
    type: boolean?
    default: false
    label: "Require that the fraction overlap be reciprocal"
    doc: |
      Require that the fraction overlap be reciprocal for each A and B.
      In other words, if -f is 0.90 and -r is used, this requires that
      B overlap 90% of A and A _also_ overlaps 90% of B.
    'sd:layout':
      advanced: true

  min_mapping_quality:
    type: int?
    default: 0
    label: "Minimum mapping quality allowed"
    doc: |
      Minimum mapping quality allowed. Default is 0.
    'sd:layout':
      advanced: true

  include_duplicate_reads:
    type: boolean?
    default: false
    label: "Include duplicate reads"
    doc: |
      Include duplicate reads.  Default counts non-duplicates only.
    'sd:layout':
      advanced: true

  include_failed_qc_reads:
    type: boolean?
    default: false
    label: "Include failed-QC reads"
    doc: |
      Include failed-QC reads.  Default counts pass-QC reads only.
    'sd:layout':
      advanced: true

  only_count_proper_pairs:
    type: boolean?
    default: false
    label: "Only count proper pairs"
    doc: |
      Only count proper pairs. Default counts all alignments with
      MAPQ > -q argument, regardless of the BAM FLAG field.
    'sd:layout':
      advanced: true


outputs:

  overlap_report_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Interval overlapping alignments counts"
    doc: "Interval overlapping alignments counts"
    outputSource: add_columns_names/output_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Overlapping counts'
        Title: 'Interval overlapping alignments counts'

  overlap_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "bedtools multicov stderr log"
    doc: "bedtools multicov stderr log"
    outputSource: overlap/stderr_log


steps:

  overlap:
    run: ../tools/bedtools-multicov.cwl
    in:
      alignment_files: alignment_files
      intervals_file: intervals_file
      split: split
      same_strand: same_strand
      opposite_strand: opposite_strand
      min_overlap_fraction: min_overlap_fraction
      reciprocal_overlap: reciprocal_overlap
      min_mapping_quality: min_mapping_quality
      include_duplicate_reads: include_duplicate_reads
      include_failed_qc_reads: include_failed_qc_reads
      only_count_proper_pairs: only_count_proper_pairs
    out:
    - report_file
    - stderr_log

  add_columns_names:
    run: ../tools/custom-bash.cwl
    in:
      input_file: overlap/report_file
      param: alignment_names
      script:
        default: |
          TOTAL=`awk '{print NF}' "$0" | sort -nu | tail -n 1`
          ALIASES=$#
          INDICES=`expr $TOTAL - $ALIASES`
          echo "Total number of columns in the file is $TOTAL" 
          echo "Among which $ALIASES last columns should have names $@"
          echo "And first $INDICES columns should have names as indices"
          for (( i=1; i<=$INDICES; ++i)); do
            echo "Add $i index"
            echo -n -e $i'\t' >> `basename $0`;
          done
          count=1
          for i in "$@"; do
            echo "Add $i name"
            echo -n -e $i >> `basename $0`;
            if [ $count != $ALIASES ]; then
              echo -n -e "\t" >> `basename $0`;
            fi
            (( count++ ))
          done;
          echo "" >> `basename $0`
          cat "$0" >> `basename $0`
    out:
    - output_file



$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Interval overlapping alignments counts"
label: "Interval overlapping alignments counts"
s:alternateName: "Interval overlapping alignments counts"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/bedtools-multicov.cwl
s:codeRepository: https://github.com/datirium/workflows
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
  Interval overlapping alignments counts
  ======================================

  Reports the count of alignments from multiple samples that overlap specific intervals.
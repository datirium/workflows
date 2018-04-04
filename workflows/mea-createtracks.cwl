cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


inputs:

  bam_file:
    type: File
    label: "Input BAM file"
    doc: "Input BAM file, sorted by coordinates"

  reference_uniquely_mapped_reads_number:
    type: int?
    label: "Uniquely mapped to reference genome reads number"
    doc: "Uniquely mapped to reference genome reads number from STAR Log.final.out"

  reference_chrom_length_file:
    type: File
    label: "Reference genome chromosome length file"
    doc: "Tab delimited chromosome length file: <chromName><TAB><chromSize>"

  refmap_file:
    type: File
    label: "Strain specific refmap file"
    doc: "Refmap file generated while making strain genome"

outputs:

  bigwig_file:
    type: File
    outputSource: sorted_projected_bedgraph_to_bigwig/bigwig_file
    label: "Strain specific bigWig file"
    doc: "Generated bigWig file for the specific strain"

steps:

  bam_to_bedgraph:
    run: ../tools/bedtools-genomecov.cwl
    in:
      input_file: bam_file
      depth:
        default: "-bg"
      split:
        default: True
      mapped_reads_number: reference_uniquely_mapped_reads_number
    out: [genome_coverage_file]

  sort_bedgraph:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: bam_to_bedgraph/genome_coverage_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  bedgraph_to_wig:
    run: ../tools/mea-bedgraphtowig.cwl
    in:
      bedgraph_file: sort_bedgraph/sorted_file
    out: [wig_file]

  mea_alea_project:
    run: ../tools/mea-alea-project.cwl
    in:
      wig_file: bedgraph_to_wig/wig_file
      refmap_file: refmap_file
    out: [bedgraph_file]

  bedgraph_header_filter:
    run: ../tools/mea-filter.cwl
    in:
      input_file: mea_alea_project/bedgraph_file
      script:
        default: 'cat "$0" | grep -v "track type=bedGraph" > `basename $0`'
    out: [filtered_file]

  sort_projected_bedgraph:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: bedgraph_header_filter/filtered_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  sorted_projected_bedgraph_to_bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      bedgraph_file: sort_projected_bedgraph/sorted_file
      chrom_length_file: reference_chrom_length_file
    out: [bigwig_file]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "mea-createtracks"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/mea-createtracks.cwl
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
  Workflow converts input `bam_file` into bedGraph, then Wig and run ALEA project command to project it to the
  reference genome using `refmap_file` file. Output is sorted and converted to bigWig format.

s:about: |
  Workflow corresponds to MEA createTracks command from
  https://github.com/julienrichardalbert/MEA/blob/e3de228734bafd957cc2072dd8a6a0e84d554724/src/scripts/createTracks.sh
cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


'sd:upstream':
  chipseq_sample:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  alignment_file:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "ChIP-Seq experiment(s)"
    doc: "Array of alignment files in BAM format"
    'sd:upstreamSource': "chipseq_sample/bambai_pair"
    'sd:localLabel': true

  genelist_file:
    type: File
    format: "http://edamontology.org/format_3585"
    label: "Genelist file, BED6 with unique ID in the 4th column"
    doc: "Genelist file in BED6 format, column 4 must be a unique peak ID, column 5 is not used"

  fragment_size:
    type: int[]
    label: "ChIP-Seq experiment(s)"
    doc: "Array of fragment sizes for input BAM files, order corresponds to the alignment_file"
    'sd:upstreamSource': "chipseq_sample/estimated_fragment_size"

  mapped_reads_number:
    type: int[]
    label: "ChIP-Seq experiment(s)"
    doc: "Array of mapped reads numners for input BAM files, order corresponds to the alignment_file"
    'sd:upstreamSource': "chipseq_sample/mapped_reads_number"

  hist_width:
    type: int?
    default: 10000
    label: "Histogram / heatmap width, bp"
    doc: "Histogram / heatmap width, bp"
    'sd:layout':
      advanced: true

  hist_bin_size:
    type: int?
    default: 50
    label: "Histogram / heatmap bin size, bp"
    doc: "Histogram / heatmap bin size, bp"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  heatmap_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "TSS centered heatmap"
    doc: "TSS centered heatmap"
    outputSource: make_tss_heatmap/histogram_file

  histogram_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "TSS centered histogram"
    doc: "TSS centered histogram"
    outputSource: make_tss_histogram/histogram_file


steps:

  make_tag_folders:
    run: ../tools/heatmap-prepare.cwl
    in:
      bam_file: alignment_file
      fragment_size: fragment_size
      total_reads: mapped_reads_number
    out: [tag_folder]

  center_genelist_on_tss:
    run: ../tools/custom-bash.cwl
    in:
      input_file: genelist_file
      script:
        default: cat "$0" | awk '{tss=$2; if ($6=="-") tss=$3; print $1"\t"tss"\t"tss"\t"$4"\t"$5"\t"$6}' > `basename $0`
    out: [output_file]

  make_tss_heatmap:
    run: ../tools/homer-annotate-peaks-hist.cwl
    in:
      peak_file: center_genelist_on_tss/output_file
      tag_folders: make_tag_folders/tag_folder
      hist_width: hist_width
      hist_bin_size: hist_bin_size
      export_heatmap:
        default: True
      threads: threads
      histogram_filename:
        source: genelist_file
        valueFrom: $(get_root(self.basename)+"_heatmap.tsv")
    out: [histogram_file]

  make_tss_histogram:
    run: ../tools/homer-annotate-peaks-hist.cwl
    in:
      peak_file: center_genelist_on_tss/output_file
      tag_folders: make_tag_folders/tag_folder
      hist_width: hist_width
      hist_bin_size: hist_bin_size
      export_heatmap:
        default: False
      threads: threads
      histogram_filename:
        source: genelist_file
        valueFrom: $(get_root(self.basename)+"_histogram.tsv")
    out: [histogram_file]



$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "Tag density profile around TSS centered gene list"
label: "Tag density profile around TSS centered gene list"
s:alternateName: "Generate tag density heatmap and histogram around centered on TSS genes from genelist TSV file"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/heatmap.cwl
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
  Generates tag density heatmap and histogram around centered by TSS features from input genelist file

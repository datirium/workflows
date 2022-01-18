cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_output_filename = function(input_file) {
          if (inputs.estimates_filename == "") {
            var ext = "_preseq_estimates.tsv";
            var root = input_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.input_file.basename+ext:root+ext;
          } else {
            return inputs.estimates_filename;
          }
      };


hints:
  - class: DockerRequirement
    dockerPull: stevetsa/preseq:2.0


inputs:

  confidence_level:
    type: float?
    inputBinding:
      position: 5
      prefix: "-cval"
    doc: "Level for confidence intervals, default: 0.95"

  extrapolation:
    type: float?
    inputBinding:
      position: 6
      prefix: "-extrap"
    doc: "Maximum extrapolation, default: 1e+10"

  max_fragment_size:
    type: int?
    inputBinding:
      position: 7
      prefix: "-seg_len"
    doc: "Maximum segment length when merging paired end bam reads, default: 5000"

  bootstraps:
    type: int?
    inputBinding:
      position: 8
      prefix: "-bootstraps"
    doc: "Number of bootstraps, default: 100"

  extrapolations_step:
    type: float?
    inputBinding:
      position: 9
      prefix: "-step"
    doc: "Step size in extrapolations, default: 1e+06"

  terms:
    type: int?
    inputBinding:
      position: 10
      prefix: "-terms"
    doc: "Maximum number of terms"

  defects_mode:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "-defects"
    doc: "Defects mode to extrapolate without testing for defects"

  quick mode:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "-quick"
    doc: "Quick mode, estimate yield without bootstrapping for confidence intervals"

  verbose_mode:
    type: boolean?
    inputBinding:
      position: 13
      prefix: "-verbose"
    doc: "Verbose mode"

  estimates_filename:
    type: string?
    inputBinding:
      position: 14
      prefix: "-output"
      valueFrom: $(get_output_filename(inputs.bam_file))
    default: ""
    doc: "Output filename"

  pe_mode:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "-pe"
    doc: "Input is paired end read file"

  bam_file:
    type: File
    inputBinding:
      position: 16
    doc: "Coordinate sorted BAM file"

outputs:

  estimates_file:
    type: File?
    outputBinding:
      glob: $(get_output_filename(inputs.bam_file))


baseCommand: ["preseq", "lc_extrap", "-bam"]


successCodes: [1, 139]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/atdp-metadata.yaml

s:name: "preseq-lc-extrap"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/preseq-lc-extrap.cwl
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
  Tool runs preseq lc_extrap. Only BAM input file is supported (-B option is used by default)
  successCodes: [1] - is used to pass this tool as a step in a workflow in case the BAM file was not correct for Preseq
  Discarded arguments:
    -V, -vals        input is a text file containing only the observed counts
    -H, -hist        input is a text file containing the observed histogram


s:about: |
  Usage: lc_extrap [OPTIONS] <sorted-bed-file>

  Options:
    -o, -output      yield output file (default: stdout)
    -e, -extrap      maximum extrapolation (default: 1e+10)
    -s, -step        step size in extrapolations (default: 1e+06)
    -n, -bootstraps  number of bootstraps (default: 100),
    -c, -cval        level for confidence intervals (default: 0.95)
    -x, -terms       maximum number of terms
    -v, -verbose     print more information
    -B, -bam         input is in BAM format
    -l, -seg_len     maximum segment length when merging paired end bam reads
                     (default: 5000)
    -P, -pe          input is paired end read file
    -V, -vals        input is a text file containing only the observed counts
    -H, -hist        input is a text file containing the observed histogram
    -Q, -quick       quick mode, estimate yield without bootstrapping for
                     confidence intervals
    -D, -defects     defects mode to extrapolate without testing for defects

  Help options:
    -?, -help        print this help message
        -about       print about message

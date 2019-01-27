cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/clip_toolkit:v1.1.3

inputs:

  big:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-big'
    label: "big input file"


  separate_strands:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-ss'
    label: "separate the two strands"

  dbkey:
    type: string?
    inputBinding:
      position: 1
      prefix: '--dbkey'
    label: "[string]: species to retrieve the default gene bed file (mm10|hg19)"

  gene:
    type: File?
    inputBinding:
      position: 1
      prefix: '--gene'
    label: "[file]: custom gene bed file for scan statistics (will override --dbkey)"

  valley_seeking:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '--valley-seeking'
    label: "find candidate peaks by valley seeking"

  infile:
    type: File
    inputBinding:
      position: 4
    label: "<tag.bed> : BED file of unique CLIP tags, input"

  outfile:
    type: string?
    inputBinding:
      position: 5
      valueFrom: |
        ${
          return inputs.outfile == "" ? inputs.infile.nameroot + ".peaks.bed" : inputs.outfile;
        }
    default: ""
    label: "<peak.bed>: BED file of called peaks, output"

outputs:
  peaks_bed:
    type: File
    outputBinding:
      glob: |
        ${
          return inputs.outfile == "" ? inputs.infile.nameroot + ".peaks.bed" : inputs.outfile;
        }
    doc: "<peak.bed>: BED file of called peaks, output"

baseCommand: [tag2peak.pl]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html


s:name: "clip-toolkit-tag2peak"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/clip-toolkit-tag2peak.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium, LLC"
  s:logo: "https://datirium.com/assets/images/datirium_llc.svg"
  s:email: mailto:support@datirium.com
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45226"
    s:streetAddress: "3559 Kroger Ave"
  s:member:
  - class: s:Person
    s:name: Artem BArski
    s:email: mailto:Artem.Barski@datirum.com
  - class: s:Person
    s:name: Andrey Kartashov
    s:email: mailto:Andrey.Kartashov@datirium.com
    s:sameAs:
    - id: http://orcid.org/0000-0001-9102-5681

doc: |
    detecting peaks from CLIP data
    Usage: tag2peak.pl [options] <tag.bed> <peak.bed>
     <tag.bed> : BED file of unique CLIP tags, input
     <peak.bed>: BED file of called peaks, output
    Options:
     -big                   : big input file
     -ss                    : separate the two strands
     --valley-seeking       : find candidate peaks by valley seeking
     --valley-depth [float] : depth of valley if valley seeking (0.9)
     --out-boundary [string]: output cluster boundaries
     --out-half-PH  [string]: output half peak height boundaries
     --dbkey        [string]: species to retrieve the default gene bed file (mm10|hg19)
     --gene         [string]: custom gene bed file for scan statistics (will override --dbkey)
     --use-expr             : use expression levels given in the score column in the custom gene bed file for normalization
     -p             [float] : threshold of p-value to call peak (0.01)
     --multi-test           : do Bonferroni multiple test correction
     -minPH         [int]   : min peak height (2)
     -maxPH         [int]   : max peak height to calculate p-value(-1, no limit if < 0)
     --skip-out-of-range-peaks: skip peaks with PH > maxPH
     -gap           [int]   : merge cluster peaks closer than the gap (-1, no merge if < 0)
     --prefix       [string]: prefix of peak id (Peak)
     -c             [dir]   : cache dir
     --keep-cache           : keep cache when the job done
     -v                     : verbose

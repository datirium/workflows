cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:
  infile:
    type: File
    inputBinding:
      position: 16
      prefix: -i
    doc: |
      The input BAM file

  split:
    type:
    - boolean?
    inputBinding:
      position: 8
      prefix: "-split"
    doc: |
      Report each portion of a “split” BAM (i.e., having an “N” CIGAR operation) alignment as a distinct BED intervals.

  outfile:
    type:
      - string?
    default: ""
    doc: |
      Name for generated output file

outputs:
  output_bed:
    type: stdout

stdout: |
        ${
          return inputs.outfile == "" ? inputs.infile.nameroot + ".bed" : inputs.outfile;
        }

baseCommand: ["bedtools", "bamtobed"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools-metadata.yaml

s:name: "bedtools-bamtobed"
label: "bedtools-bamtobed"
s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/tools/bedtools-bamtobed.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

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
  Tool:    bedtools bamtobed (aka bamToBed)
  Version: v2.26.0
  Summary: Converts BAM alignments to BED6 or BEDPE format.

  Usage:   bedtools bamtobed [OPTIONS] -i <bam>

  Options:
      -bedpe	Write BEDPE format.
          - Requires BAM to be grouped or sorted by query.

      -mate1	When writing BEDPE (-bedpe) format,
          always report mate one as the first BEDPE "block".

      -bed12	Write "blocked" BED format (aka "BED12"). Forces -split.

          http://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1

      -split	Report "split" BAM alignments as separate BED entries.
          Splits only on N CIGAR operations.

      -splitD	Split alignments based on N and D CIGAR operators.
          Forces -split.

      -ed	Use BAM edit distance (NM tag) for BED score.
          - Default for BED is to use mapping quality.
          - Default for BEDPE is to use the minimum of
            the two mapping qualities for the pair.
          - When -ed is used with -bedpe, the total edit
            distance from the two mates is reported.

      -tag	Use other NUMERIC BAM alignment tag for BED score.
          - Default for BED is to use mapping quality.
            Disallowed with BEDPE output.

      -color	An R,G,B string for the color used with BED12 format.
          Default is (255,0,0).

      -cigar	Add the CIGAR string to the BED entry as a 7th column.
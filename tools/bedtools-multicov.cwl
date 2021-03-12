cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.intervals_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.intervals_file.basename+".tsv":root+".tsv";
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  alignment_files:
    type:
    - File
    - File[]
    secondaryFiles:
    - .bai
    inputBinding:
      position: 5
      prefix: "-bams"
    doc: "Sorted and indexed BAM file(s)"

  intervals_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-bed"
    doc: "Intervals file defined in a BED/GFF/VCF format"

  split:
    type: boolean?
    inputBinding:
      position: 7
      prefix: "-split"
    doc: "Treat 'split' BAM or BED12 entries as distinct BED intervals"

  same_strand:
    type: boolean?
    inputBinding:
      position: 8
      prefix: "-s"
    doc: |
      Require same strandedness. That is, only report hits in B
      that overlap A on the _same_ strand. By default, overlaps
      are reported without respect to strand.

  opposite_strand:
    type: boolean?
    inputBinding:
      position: 9
      prefix: "-S"
    doc: |
      Require different strandedness. That is, only report hits in B
      that overlap A on the _opposite_ strand. By default, overlaps
      are reported without respect to strand.

  min_overlap_fraction:
    type: float?
    inputBinding:
      position: 10
      prefix: "-f"
    doc: |
      Minimum overlap required as a fraction of each A.
      Default is 1E-9 (i.e., 1bp).

  reciprocal_overlap:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "-r"
    doc: |
      Require that the fraction overlap be reciprocal for each A and B.
      In other words, if -f is 0.90 and -r is used, this requires that
      B overlap 90% of A and A _also_ overlaps 90% of B.

  min_mapping_quality:
    type: int?
    inputBinding:
      position: 12
      prefix: "-q"
    doc: |
      Minimum mapping quality allowed. Default is 0.

  include_duplicate_reads:
    type: boolean?
    inputBinding:
      position: 13
      prefix: "-D"
    doc: |
      Include duplicate reads.  Default counts non-duplicates only.

  include_failed_qc_reads:
    type: boolean?
    inputBinding:
      position: 14
      prefix: "-F"
    doc: |
      Include failed-QC reads.  Default counts pass-QC reads only.

  only_count_proper_pairs:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "-p"
    doc: |
      Only count proper pairs. Default counts all alignments with
      MAPQ > -q argument, regardless of the BAM FLAG field.

  output_filename:
    type: string?
    default: ""
    doc: "Report file name"


outputs:

  report_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "TSV file with collected counts per inteval"

  stderr_log:
    type: stderr


baseCommand: ["bedtools", "multicov"]
stdout: $(default_output_filename())
stderr: "bedtools_multicov_stderr.log"


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "bedtools-multicov"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bedtools-multicov.cwl
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
  Reports the count of alignments from multiple position-sorted and indexed BAM files that overlap intervals in a BED file


s:about: |
  Usage:   bedtools multicov [OPTIONS] -bams aln.1.bam aln.2.bam ... aln.n.bam -bed <bed/gff/vcf>

  Options: 
    -bams	The bam files.

    -bed	The bed file.

    -split	Treat "split" BAM or BED12 entries as distinct BED intervals.

    -s	Require same strandedness.  That is, only report hits in B
      that overlap A on the _same_ strand.
      - By default, overlaps are reported without respect to strand.

    -S	Require different strandedness.  That is, only report hits in B
      that overlap A on the _opposite_ strand.
      - By default, overlaps are reported without respect to strand.

    -f	Minimum overlap required as a fraction of each A.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)

    -r	Require that the fraction overlap be reciprocal for each A and B.
      - In other words, if -f is 0.90 and -r is used, this requires
        that B overlap 90% of A and A _also_ overlaps 90% of B.

    -q	Minimum mapping quality allowed. Default is 0.

    -D	Include duplicate reads.  Default counts non-duplicates only

    -F	Include failed-QC reads.  Default counts pass-QC reads only

    -p	Only count proper pairs.  Default counts all alignments with
      MAPQ > -q argument, regardless of the BAM FLAG field.
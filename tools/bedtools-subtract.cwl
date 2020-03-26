cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return inputs.reduced_bed_file.basename;
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  reduced_bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-a"
    doc: |
      The input BED file from which the features will be substructed

  subtracted_bed_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-b"
    doc: |
      The input BED file that includes the features to be subtracted from -a

  output_filename:
    type: string?
    default: ""
    doc: |
      Output file name


outputs:

  difference_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Difference BED file


baseCommand: ["bedtools", "subtract"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools-metadata.yaml

s:name: "bedtools-subtract"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bedtools-subtract.cwl
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
  Searches for features in B that overlap A by at least 1 base pair.
  If an overlapping feature is found in B, the overlapping portion is removed from A
  and the remaining portion of A is reported. If a feature in B overlaps all of
  a feature in A, the A feature will not be reported.
  All parameters except -a and -b used by default.

  
s:about: |
  Usage:   bedtools subtract [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

  Options: 
    -A	Remove entire feature if any overlap.  That is, by default,
      only subtract the portion of A that overlaps B. Here, if
      any overlap is found (or -f amount), the entire feature is removed.

    -N	Same as -A except when used with -f, the amount is the sum
      of all features (not any single feature).

    -wb	Write the original entry in B for each overlap.
      - Useful for knowing _what_ A overlaps. Restricted by -f and -r.

    -wo	Write the original A and B entries plus the number of base
      pairs of overlap between the two features.
      - Overlaps restricted by -f and -r.
        Only A features with overlap are reported.

    -s	Require same strandedness.  That is, only report hits in B
      that overlap A on the _same_ strand.
      - By default, overlaps are reported without respect to strand.

    -S	Require different strandedness.  That is, only report hits in B
      that overlap A on the _opposite_ strand.
      - By default, overlaps are reported without respect to strand.

    -f	Minimum overlap required as a fraction of A.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)

    -F	Minimum overlap required as a fraction of B.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)

    -r	Require that the fraction overlap be reciprocal for A AND B.
      - In other words, if -f is 0.90 and -r is used, this requires
        that B overlap 90% of A and A _also_ overlaps 90% of B.

    -e	Require that the minimum fraction be satisfied for A OR B.
      - In other words, if -e is used with -f 0.90 and -F 0.10 this requires
        that either 90% of A is covered OR 10% of  B is covered.
        Without -e, both fractions would have to be satisfied.

    -split	Treat "split" BAM or BED12 entries as distinct BED intervals.

    -g	Provide a genome file to enforce consistent chromosome sort order
      across input files. Only applies when used with -sorted option.

    -nonamecheck	For sorted data, don't throw an error if the file has different naming conventions
        for the same chromosome. ex. "chr1" vs "chr01".

    -sorted	Use the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input.

    -bed	If using BAM input, write output as BED.

    -header	Print the header from the A file prior to results.

    -nobuf	Disable buffered output. Using this option will cause each line
      of output to be printed as it is generated, rather than saved
      in a buffer. This will make printing large output files 
      noticeably slower, but can be useful in conjunction with
      other software tools and scripts that need to process one
      line of bedtools output at a time.

    -iobuf	Specify amount of memory to use for input buffer.
      Takes an integer argument. Optional suffixes K/M/G supported.
      Note: currently has no effect with compressed files.
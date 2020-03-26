cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return inputs.bed_file.basename;
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: |
      The input BED file with features to be shuffled

  chrom_length_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-g"
    doc: |
      Input genome file with chromosome lengths

  excl_bed_file:
    type: File?
    inputBinding:
      position: 7
      prefix: "-excl"
    doc: |
      The BED file with regions where you do not want the permuted features to be placed

  incl_bed_file:
    type: File?
    inputBinding:
      position: 8
      prefix: "-incl"
    doc: |
      The BED file with regions where you want the permuted features to be placed

  seed:
    type: int?
    inputBinding:
      position: 9
      prefix: "-seed"
    doc: |
      Seed for pseudo-random number generation. Default: random

  no_overlapping:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "-noOverlapping"
    doc: |
      Don't allow shuffled intervals to overlap
      
  max_tries:
    type: int?
    inputBinding:
      position: 11
      prefix: "-maxTries"
    doc: |
      Maximum number of attempts to find a home for a shuffled interval. Default: 1000

  output_filename:
    type: string?
    default: ""
    doc: |
      Output file name


outputs:

  shuffled_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Shuffled BED file


baseCommand: ["bedtools", "shuffle"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools-metadata.yaml

s:name: "bedtools-shuffle"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/bedtools-shuffle.cwl
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
 Randomly permutes the genomic locations of a feature file among a genome defined in a genome file.
 One can also provide an “exclusions” BED file that lists regions where you do not want the permuted
 features to be placed. Or instead “inclusions” BED fils that defines coordinates in which features
 in -i should be randomly placed. To make experiment reproducible, set "seed" option.
 NOTE: limited parameters are impelented.


s:about: |
  Usage:   bedtools shuffle [OPTIONS] -i <bed/gff/vcf> -g <genome>

  Options: 
    -excl	A BED/GFF/VCF file of coordinates in which features in -i
      should not be placed (e.g. gaps.bed).

    -incl	Instead of randomly placing features in a genome, the -incl
      options defines a BED/GFF/VCF file of coordinates in which 
      features in -i should be randomly placed (e.g. genes.bed). 
      Larger -incl intervals will contain more shuffled regions. 
      This method DISABLES -chromFirst. 
    -chrom	Keep features in -i on the same chromosome.
      - By default, the chrom and position are randomly chosen.
      - NOTE: Forces use of -chromFirst (see below).

    -seed	Supply an integer seed for the shuffling.
      - By default, the seed is chosen automatically.
      - (INTEGER)

    -f	Maximum overlap (as a fraction of the -i feature) with an -excl
      feature that is tolerated before searching for a new, 
      randomized locus. For example, -f 0.10 allows up to 10%
      of a randomized feature to overlap with a given feature
      in the -excl file. **Cannot be used with -incl file.**
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)

    -chromFirst	
      Instead of choosing a position randomly among the entire
      genome (the default), first choose a chrom randomly, and then
      choose a random start coordinate on that chrom.  This leads
      to features being ~uniformly distributed among the chroms,
      as opposed to features being distribute as a function of chrom size.

    -bedpe	Indicate that the A file is in BEDPE format.

    -maxTries	
      Max. number of attempts to find a home for a shuffled interval
      in the presence of -incl or -excl.
      Default = 1000.
    -noOverlapping	
      Don't allow shuffled intervals to overlap.
    -allowBeyondChromEnd	
      Allow shuffled intervals to be relocated to a position
      in which the entire original interval cannot fit w/o exceeding
      the end of the chromosome.  In this case, the end coordinate of the
      shuffled interval will be set to the chromosome's length.
      By default, an interval's original length must be fully-contained
      within the chromosome.
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - $import: ./metadata/envvar-global.yml
  - class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
  dockerFile: >
    $import: ./dockerfiles/bedtools-Dockerfile

inputs:
  input:
    type: File
    doc: |
      The input file can be in BAM format
          (Note: BAM _must_ be sorted by position)
      or <bed/gff/vcf>
    inputBinding:
      position: 10
      valueFrom: |
          ${
            var prefix = ((/.*\.bam$/i).test(inputs.input.path))?'-ibam':'-i';
            return [prefix,inputs.input.path];
          }
# TODO Need this only to make sure that bam file is sorted.
    secondaryFiles: |
           ${
            if ((/.*\.bam$/i).test(self.location))
               return {"location": self.location+".bai", "class": "File"};
            return [];
           }

# Needed only when -i flag
  genomeFile:
    type: File?
    doc:
      Input genome file.
    inputBinding:
      position: 11
      prefix: "-g"

  dept:
    type:
      name: "JustDepts"
      type: enum
      symbols: ["-bg","-bga","-d", "-dz"]
    inputBinding:
      position: 5

  scale:
    type: float?
    doc: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
    inputBinding:
      position: 5
      prefix: -scale

  mappedreads:
    type: int?
    doc: |
      Optional parameter to calculate scale as 1000000/mappedreads
    inputBinding:
      position: 5
      prefix: -scale
      valueFrom: |
        ${
          if (inputs.scale){
            return null;
          } else {
            return 1000000/inputs.mappedreads;
          }
        }

  split:
    type: ["null",boolean]
    doc: |
      treat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
    inputBinding:
      position: 5
      prefix: "-split"

  strand:
    type: ["null", string]
    doc: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
    inputBinding:
      position: 5
      prefix: "-strand"

  pairchip:
    type: ["null", boolean]
    doc: "pair-end chip seq experiment"
    inputBinding:
      position: 5
      prefix: "-pc"

  du:
    type: ["null", boolean]
    doc: |
      Change strand af the mate read (so both reads from the same strand) useful for strand specific.
      Works for BAM files only
    inputBinding:
      position: 5
      prefix: "-du"

  fragmentsize:
    type: ["null", int]
    doc: "fixed fragment size"
    inputBinding:
      position: 5
      prefix: "-fs"

  max:
    type: ["null",int]
    doc: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
    inputBinding:
      position: 5
      prefix: "-max"

  m5:
    type: ["null",boolean]
    doc: |
      Calculate coverage of 5" positions (instead of entire interval).
    inputBinding:
      position: 5
      prefix: "-5"

  m3:
    type: ["null",boolean]
    doc: |
      Calculate coverage of 3" positions (instead of entire interval).
    inputBinding:
      position: 5
      prefix: "-3"

  genomecoverageout:
    type: string?

outputs:
  genomecoverage:
    type: File
    doc: "The file containing the genome coverage"
    outputBinding:
      glob: |
        ${
          if (inputs.genomecoverageout == null){
            return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".bed";
          } else {
            return inputs.genomecoverageout;
          }
        }

stdout: |
  ${
    if (inputs.genomecoverageout == null){
      return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".bed";
    } else {
      return inputs.genomecoverageout;
    }
  }

baseCommand: ["bedtools", "genomecov"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools-metadata.yaml

s:name: "bedtools-genomecov"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/tools/bedtools-genomecov.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681

doc: |
  Tool is used to calculate statistics on the base of FASTQ file quality scores

s:about: >
  Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

  Options:
  	-ibam		The input file is in BAM format.
  			Note: BAM _must_ be sorted by position

  	-d		Report the depth at each genome position (with one-based coordinates).
  			Default behavior is to report a histogram.

  	-dz		Report the depth at each genome position (with zero-based coordinates).
  			Reports only non-zero positions.
  			Default behavior is to report a histogram.

  	-bg		Report depth in BedGraph format. For details, see:
  			genome.ucsc.edu/goldenPath/help/bedgraph.html

  	-bga		Report depth in BedGraph format, as above (-bg).
  			However with this option, regions with zero
  			coverage are also reported. This allows one to
  			quickly extract all regions of a genome with 0
  			coverage by applying: "grep -w 0$" to the output.

  	-split		Treat "split" BAM or BED12 entries as distinct BED intervals.
  			when computing coverage.
  			For BAM files, this uses the CIGAR "N" and "D" operations
  			to infer the blocks for computing coverage.
  			For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
  			fields (i.e., columns 10,11,12).

  	-strand		Calculate coverage of intervals from a specific strand.
  			With BED files, requires at least 6 columns (strand is column 6).
  			- (STRING): can be + or -

  	-pc		Calculate coverage of pair-end fragments.
  			Works for BAM files only
  	-fs		Force to use provided fragment size instead of read length
  			Works for BAM files only
  	-du		Change strand af the mate read (so both reads from the same strand) useful for strand specific
  			Works for BAM files only
  	-5		Calculate coverage of 5" positions (instead of entire interval).

  	-3		Calculate coverage of 3" positions (instead of entire interval).

  	-max		Combine all positions with a depth >= max into
  			a single bin in the histogram. Irrelevant
  			for -d and -bedGraph
  			- (INTEGER)

  	-scale		Scale the coverage by a constant factor.
  			Each coverage value is multiplied by this factor before being reported.
  			Useful for normalizing coverage by, e.g., reads per million (RPM).
  			- Default is 1.0; i.e., unscaled.
  			- (FLOAT)

  	-trackline	Adds a UCSC/Genome-Browser track line definition in the first line of the output.
  			- See here for more details about track line definition:
  			      http://genome.ucsc.edu/goldenPath/help/bedgraph.html
  			- NOTE: When adding a trackline definition, the output BedGraph can be easily
  			      uploaded to the Genome Browser as a custom track,
  			      BUT CAN NOT be converted into a BigWig file (w/o removing the first line).

  	-trackopts	Writes additional track line definition parameters in the first line.
  			- Example:
  			   -trackopts 'name="My Track" visibility=2 color=255,30,30'
  			   Note the use of single-quotes if you have spaces in your parameters.
  			- (TEXT)

  Notes:
  	(1) The genome file should tab delimited and structured as follows:
  	 <chromName><TAB><chromSize>

  	For example, Human (hg19):
  	chr1	249250621
  	chr2	243199373
  	...
  	chr18_gl000207_random	4262

  	(2) The input BED (-i) file must be grouped by chromosome.
  	 A simple "sort -k 1,1 <BED> > <BED>.sorted" will suffice.

  	(3) The input BAM (-ibam) file must be sorted by position.
  	 A "samtools sort <BAM>" should suffice.

  Tips:
  	One can use the UCSC Genome Browser's MySQL database to extract
  	chromosome sizes. For example, H. sapiens:

  	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
  	"select chrom, size from hg19.chromInfo" > hg19.genome

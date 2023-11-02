cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return "barcodes_" + inputs.selected_chromosome.replace(/ /g,'_') + ".tsv";
          } else {
            return inputs.output_filename;
          }
        };
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.11


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      cp $0 tmp_sorted.bam
      cp $0.bai tmp_sorted.bam.bai
      cp $1 tmp_fb.tar.gz
      tar xzf tmp_fb.tar.gz --strip-components=2
      gzip -d barcodes.tsv.gz -c > tmp_barcodes.tsv
      samtools view -@ $3 tmp_sorted.bam $2 | LC_ALL=C grep "xf:i:25" > tmp_selected_chr.sam
      samtools view -@ $3 -H tmp_sorted.bam > tmp_header.tsv
      cat tmp_header.tsv tmp_selected_chr.sam > tmp_selected_chr_with_header.sam
      samtools view -@ $3 -b tmp_selected_chr_with_header.sam > tmp_selected_chr_with_header.bam
      samtools view -@ $3 tmp_selected_chr_with_header.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c | tr -s " " | tr " " "\t" | cut -f 2-3 | awk '{print $2"-1\t"$1}' | sort -k1,1 > tmp_selected_chr_counts.tsv
      cat tmp_barcodes.tsv | sort -k1,1 > tmp_sorted_barcodes.tsv
      join -1 1 -2 1 -t $'\t' -a1 -e "0" -o "0,2.2" tmp_sorted_barcodes.tsv tmp_selected_chr_counts.tsv
    inputBinding:
      position: 5
    doc: "Bash script to count reads from the selected chromosomes"

  possorted_genome_bam_bai:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 6
    doc: Position-sorted and indexed BAM file

  feature_bc_matrix_folder:
    type: File
    inputBinding:
      position: 7
    doc: Feature barcode matrix, compressed, MEX

  selected_chromosome:
    type: string
    inputBinding:
      position: 8
    doc: |
      Chromosome name(s) to include in the results.
      If space separated list of values provided,
      sum the counts and report in one column

  threads:
    type: int?
    inputBinding:
      position: 9
    default: 1
    doc: "Number of threads to use"

  output_filename:
    type: string?
    default: ""
    doc: "Output filename"

outputs:

  barcodes_data:
    type: File
    outputBinding:
      glob: $(default_output_filename())

baseCommand: [bash, '-c']
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Utilities Chromosome Reads Counts"
s:name: "Cell Ranger Utilities Chromosome Reads Counts"
s:alternateName: "Counts reads from the selected chromosomes of Cell Ranger generated BAM file"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellranger-utils-count-chr.cwl
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
 Cell Ranger Utilities Chromosome Reads Counts

s:about: |
 Counts reads from the selected chromosomes of Cell Ranger generated BAM file
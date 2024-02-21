cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement


inputs:

  genome:
    type: string
    label: "Genome type"
    doc: "Genome type, such as mm10, hg19, hg38, etc"

  genome_label:
    type: string?
    label: "Genome label"
    sd:preview:
      position: 1

  genome_description:
    type: string?
    label: "Genome description"
    sd:preview:
      position: 2

  genome_details:
    type: string?
    label: "Genome details"
    sd:preview:
      position: 3

  genome_file:
    type: File
    format: "http://edamontology.org/format_3009"
    label: "Reference genome file (*.2bit, *.fasta, *.fa, *.fa.gz, *.fasta.gz)"
    doc: "Reference genome file (*.2bit, *.fasta, *.fa, *.fa.gz, *.fasta.gz). All chromosomes are included"

  chromosome_list:
    type:
      - "null"
      - string[]
    label: "Chromosome list to be included into the reference genome FASTA file"
    doc: "Filter chromosomes while extracting FASTA from 2bit"

  chrom_sizes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Genome chromosome length file. All chromosomes are included"
    doc: "Genome chromosome length file. All chromosomes are included"

  annotation_tab:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Compressed tsv.gz annotation file"
    doc: "Compressed tab-separated annotation file. Doesn't include chrM"

  mitochondrial_annotation_tab:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Compressed tsv.gz mitochondrial DNA annotation file"
    doc: "Compressed mitochondrial DNA tab-separated annotation file. Includes only chrM"

  cytoband:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Compressed cytoBand file for IGV browser"
    doc: "Compressed tab-separated cytoBand file for IGV browser"


outputs:

  indices_folder:
    type: Directory
    label: "Bismark indices folder"
    doc: "Bismark generated indices folder"
    outputSource: prepare_indices/indices_folder

  fasta_output:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes only selected chromosomes"
    outputSource: extract_fasta/fasta_file

  fasta_fai_output:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "FAI index for genome FASTA file"
    doc: "Tab-separated FAI index file"
    outputSource: index_fasta/fai_file

  annotation:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "TSV annotation file"
    doc: "Tab-separated annotation file. Includes reference genome and mitochondrial DNA annotations"
    outputSource: prepare_annotation/annotation_tsv_file

  annotation_gtf:
    type: File
    format: "http://edamontology.org/format_2306"
    label: "GTF annotation file"
    doc: "GTF annotation file. Includes reference genome and mitochondrial DNA annotations"
    outputSource: prepare_annotation/annotation_gtf_file

  annotation_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Sorted BED annotation file"
    doc: "Sorted BED annotation file"
    outputSource: sort_annotation_bed/sorted_file

  annotation_bed_tbi:
    type: File
    format: "http://edamontology.org/format_3004"
    label: "Sorted bigBed annotation file"
    doc: "Sorted bigBed annotation file"
    outputSource: annotation_bed_to_bigbed/bigbed_file

  cytoband_output:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "CytoBand file for IGV browser"
    doc: "Tab-separated cytoBand file for IGV browser"
    outputSource: extract_cytoband/output_file

  chrom_length:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: filter_chrom_sizes/output_file
    label: "Genome chromosome length file"
    doc: "Genome chromosome length file"

  stdout_log:
    type: File
    label: "Bismark stdout log"
    doc: "Bismark generated stdout log"
    outputSource: prepare_indices/stdout_log

  stderr_log:
    type: File
    label: "Bismark stderr log"
    doc: "Bismark generated stderr log"
    outputSource: prepare_indices/stderr_log


steps:

  extract_fasta:
    run: ../tools/ucsc-twobit-to-fa.cwl
    in:
      reference_file: genome_file
      chr_list: chromosome_list
    out:
    - fasta_file

  fasta_to_folder:
    in:
      genome_fasta: extract_fasta/fasta_file
      genome_type: genome
    out: [genome_folder]
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      requirements:
      - class: InlineJavascriptRequirement
      inputs:
        genome_fasta: File
        genome_type: string
      outputs:
        genome_folder: Directory
      expression: |
        ${
            return { "genome_folder": {
              "class": "Directory",
              "basename": inputs.genome_type,
              "listing": [inputs.genome_fasta]
            }};
        }

  prepare_indices:
    run: ../tools/bismark-prepare-genome.cwl
    in:
      genome_folder: fasta_to_folder/genome_folder
    out:
    - indices_folder
    - stdout_log
    - stderr_log

  index_fasta:
    run: ../tools/samtools-faidx.cwl
    in:
      fasta_file: extract_fasta/fasta_file
    out:
    - fai_file

  filter_chrom_sizes:
    run: ../tools/custom-bash.cwl
    in:
      input_file: chrom_sizes
      param: chromosome_list
      script:
        default: |
          CHR="$(IFS="|"; echo "$*")"
          cat "$0" | grep -E -w $CHR > chrom_length_file.tsv
    out:
    - output_file

  prepare_annotation:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/ucscuserapps:v358
      inputs:
        script:
          type: string?
          default: |
            #!/bin/bash
            gunzip $0 -c | grep -v "exonCount" > refgene.txt
            gunzip $1 -c | grep -v "exonCount" | awk '{ if ($3=="chrM") print $0 }' >> refgene.txt
            if [ "$#" -ge 2 ]; then
              FILTER=${@:2}
              FILTER=$( IFS=$','; echo "${FILTER[*]}" )
              FILTER=(${FILTER//, / })
              echo "Filtering by" ${FILTER[*]}
              cat refgene.txt | awk -v filter="${FILTER[*]}" 'BEGIN {split(filter, f); for (i in f) d[f[i]]} {if ($3 in d) print $0}' > refgene_filtered.txt  
              mv refgene_filtered.txt refgene.txt
            fi
            cut -f 2- refgene.txt | genePredToGtf file stdin refgene.gtf
            echo -e "bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames" > refgene.tsv
            cat refgene.txt >> refgene.tsv
          inputBinding:
            position: 5
        genome_annotation:
          type: File
          inputBinding:
            position: 6
        mitochondrial_annotation:
          type: File
          inputBinding:
            position: 7
        chromosome_list:
          type:
            - "null"
            - string
            - string[]
          inputBinding:
            position: 8
      outputs:
        annotation_tsv_file:
          type: File
          outputBinding:
            glob: "refgene.tsv"
        annotation_gtf_file:
          type: File
          outputBinding:
            glob: "refgene.gtf"
      baseCommand: ["bash", "-c"]
    in:
      genome_annotation:  annotation_tab
      mitochondrial_annotation: mitochondrial_annotation_tab
      chromosome_list: chromosome_list
    out:
    - annotation_tsv_file
    - annotation_gtf_file

  convert_annotation_to_bed:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: InlineJavascriptRequirement
        expressionLib:
        - var default_output_filename = function() {
                var root = inputs.annotation_tsv_file.basename.split('.').slice(0,-1).join('.');
                return (root == "")?inputs.annotation_tsv_file.basename+".bed":root+".bed";
              };
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        script:
          type: string?
          default: |
            import fileinput
            for line in fileinput.input():
                if "txStart" in line:
                  continue
                cols = line.split("\t")
                refName = cols[1]
                chrom = cols[2]
                txStart = cols[4]
                txEnd = cols[5]
                txStart1 = cols[6]
                txEnd1 = cols[7]
                try: 
                    name = cols[12]
                except Exception:
                    name = cols[11]
                    pass
                strand = cols[3]
                exonCount = cols[8]
                cdsStart = cols[9].split(',')[0:-1]
                cdsEnd = cols[10].split(',')[0:-1]
                startEndPairs = zip(cdsStart, cdsEnd)
                sizes =  ','.join(map(lambda pair: str(int(pair[1])-int(pair[0])), startEndPairs))
                deltas = ','.join(map(lambda offset: str(int(offset)-int(txStart)), cdsStart))
                if 'fix' in chrom or '_' in chrom:
                    continue
                output = [chrom, txStart, txEnd, name, '1000', strand, txStart1, txEnd1, '.', exonCount, sizes, deltas]
                print "\t".join(output)
          inputBinding:
            position: 5
          doc: "Python script to get convert TSV annotation file to BED"
        annotation_tsv_file:
          type: File
          inputBinding:
            position: 6
          doc: "Annotation TSV file"
      outputs:
        annotation_bed_file:
          type: File
          outputBinding:
            glob: "*"
      baseCommand: ["python", "-c"]
      stdout: $(default_output_filename())
    in:
      annotation_tsv_file: prepare_annotation/annotation_tsv_file
    out:
    - annotation_bed_file

  sort_annotation_bed:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: convert_annotation_to_bed/annotation_bed_file
      key:
        default: ["1,1","2,2n"]
    out:
    - sorted_file

  annotation_bed_to_bigbed:
    run: ../tools/ucsc-bedtobigbed.cwl
    in:
      input_bed: sort_annotation_bed/sorted_file
      chrom_length_file: filter_chrom_sizes/output_file
      bed_type:
        default: "bed4+8"
      output_filename:
        source: sort_annotation_bed/sorted_file
        valueFrom: $(self.basename + ".tbi")
    out:
    - bigbed_file

  extract_cytoband:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: InitialWorkDirRequirement
        listing: |
          ${
            return  [
                      {
                        "entry": inputs.input_file,
                        "entryname": inputs.input_file.basename,
                        "writable": true
                      }
                    ]
          }
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        input_file:
          type: File
          inputBinding:
            position: 1
      outputs:
        output_file:
          type: File
          outputBinding:
            glob: "*"
      baseCommand: [gunzip]
    in:
      input_file: cytoband
    out:
    - output_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Build Bismark indices"
label: "Build Bismark indices"
s:alternateName: "Build indices for Bismark Methylation Pipeline. Bowtie2 aligner is used by default"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/bismark-index.cwl
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
        s:name: Michael Kotlia
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

# doc:
#   $include: ../descriptions/bismark-index.md


doc: |
  Copy fasta_file file to the folder and run run bismark_genome_preparation script to prepare indices for Bismark Methylation Analysis.
  Bowtie2 aligner is used by default. The name of the output indices folder is equal to the genome input.
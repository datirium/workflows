cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/altanalyze:v0.0.6


inputs: 

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      
      # Copy altanalyze to the current working directory
      # which is mount with -rw- permissions. Otherwise 
      # we can't download anything because by default
      # container is run by cwltool with --read-only
      
      cp -r /opt/altanalyze .
      
      GENOME_TYPE=$0
      case ${GENOME_TYPE} in
        mm10)
          SPECIES=Mm
          # ENSEMBL_VERSION=EnsMart69
          ENSEMBL_VERSION=EnsMart72
          ;;
        hg19)
          SPECIES=Hs
          # ENSEMBL_VERSION=EnsMart62
          ENSEMBL_VERSION=EnsMart72
          ;;
        hg38)
          SPECIES=Hs
          # ENSEMBL_VERSION=EnsMart78
          ENSEMBL_VERSION=EnsMart72
          ;;
        *)
          echo "Unknown genome"
          exit 1
          ;;
      esac
      echo "Selected ${SPECIES} and ${ENSEMBL_VERSION}"
      python ./altanalyze/AltAnalyze.py --species ${SPECIES} --update Official --version ${ENSEMBL_VERSION}
      mkdir genome_data
      mv ./altanalyze/AltDatabase/EnsMart* ./genome_data/${ENSEMBL_VERSION}__${SPECIES}
      echo "Remove files with ' in the name as we can't mount them to container anyway"
      cd ./genome_data/${ENSEMBL_VERSION}__${SPECIES}
      du -a | grep "'" | cut -f 2 | awk '{print "\""$1"\""}' | xargs rm -f
    inputBinding:
      position: 5
    doc: |
      Bash script to select proper parameters for AltAnalyze based on the genome value

  genome:
    type:
      type: enum
      symbols:
      - "mm10"
      - "hg19"
      - "hg38"
      inputBinding:
        position: 6
    doc: |
      Genome type, such as mm10, hg19 or hg38. Based on the selected value we will set
      --species and --version parameters.


outputs:

  genome_data:
    type: Directory
    outputBinding:
      glob: "genome_data/*"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]


stdout: altanalyze_prepare_genome_stdout.log
stderr: altanalyze_prepare_genome_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "altanalyze-prepare-genome"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/altanalyze-prepare-genome.cwl
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
    - class: s:Organization
      s:legalName: "Salomonis Research Lab"
      s:member:
      - class: s:Person
        s:name: Stuart Hay
        s:email: mailto:haysb91@gmail.com


doc: |
  Downloads data from the Ensembl Database based on the selected genome.
  Output folder name will be always equal to ENSEMBL_VERSION__SPECIES.

  AltAnalyze has problems downloading EnsMart69 and EnsMart78, so only hg19
  genome is currently supported.


s:about: |
  EnsMart72 has GRCh37.p11, GRCm38.p1 and no GRCh38, which is not good for us.

  EnsMart69 has original version of GRCm38. The same we use for our mm10 pipelines.
  EnsMart62 has original version of GRCh37. The same we use for our hg19 pipelines.
  EnsMart78 has original version of GRCh38. The same we use for our hg38 pipelines.

  For reference:
  http://useast.ensembl.org/info/website/archives/assembly.html
  https://m.ensembl.org/Help/Faq?id=286

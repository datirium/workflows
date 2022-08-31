cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for calc_frip function in collect-statistics-frip.cwl tool:\n\n" > "calc_frip.log"
      bam="$0"; bed1="$1"; bed2="$2";
      md="$3"; tsv="$4"; yaml="$5"
      # count of total aligned reads
      tar=$(samtools view -cF0x4 "$bam" 2>> "calc_frip.log")
      # counts of reads in peaks (raw/rip2, norm/rip2)
      rip1=$(samtools view -c "$bam" -L "$bed1" 2>> "calc_frip.log")
      rip2=$(samtools view -c "$bam" -L "$bed2" 2>> "calc_frip.log")
      # frip=rip/tar
      frip=$(printf "$tar" | awk -v rip="$rip1" '{printf("%.3f",rip/$0)}' 2>> "calc_frip.log")
      frip_norm=$(printf "$tar" | awk -v rip="$rip2" '{printf("%.3f",rip/$0)}' 2>> "calc_frip.log")
      printf "$tar, $rip1, $rip2, $frip, $frip_norm\n" >> "calc_frip.log"
      # concatenate formatted frip onto md, tsv, and yaml files
      #   md
      cat "$md" > collected_statistics_report.md
      printf "-" >> collected_statistics_report.md 2>> "calc_frip.log"
      printf "   fraction of (aligned) reads in peaks (raw): $frip\n" >> collected_statistics_report.md 2>> "calc_frip.log"
      printf "-" >> collected_statistics_report.md 2>> "calc_frip.log"
      printf "   fraction of (aligned) reads in peaks (scaled): $frip_norm\n" >> collected_statistics_report.md 2>> "calc_frip.log"
      #   tsv
      headers=$(head -1 "$tsv"); data=$(tail -1 "$tsv")
      printf "%s\t%s\t%s\n%s\t%s\t%s\n" "$headers" "fraction of (aligned) reads in peaks (raw)" "fraction of (aligned) reads in peaks (scaled)" "$data" "$frip" "$frip_norm" > collected_statistics_report.tsv 2>> "calc_frip.log"
      #   yaml
      cat "$yaml" > collected_statistics_report.yaml
      printf "  fraction of (aligned) reads in peaks (raw): $frip\n" >> collected_statistics_report.yaml 2>> "calc_frip.log"
      printf "  fraction of (aligned) reads in peaks (scaled): $frip_norm\n" >> collected_statistics_report.yaml 2>> "calc_frip.log"
    inputBinding:
        position: 4

  bam_file:
    type: File
    label: "Input BAM file"
    inputBinding:
      position: 5
    doc: "Input BAM file, does not have to be coordinates sorted"

  seacr_called_peaks:
    type: File
    inputBinding:
      position: 11

  seacr_called_peaks_norm:
    type: File
    inputBinding:
      position: 12

  collected_statistics_md:
    type: File
    inputBinding:
      position: 20

  collected_statistics_tsv:
    type: File
    inputBinding:
      position: 21

  collected_statistics_yaml:
    type: File
    inputBinding:
      position: 22


outputs:

  modified_file_md:
    type: File
    outputBinding:
      glob: "collected_statistics_report.md"

  modified_file_tsv:
    type: File
    outputBinding:
      glob: "collected_statistics_report.tsv"

  modified_file_yaml:
    type: File
    outputBinding:
      glob: "collected_statistics_report.yaml"

  log_file:
    type: File
    outputBinding:
      glob: "calc_frip.log"


baseCommand: ["bash", "-c"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "collect-statistics-frips"
s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/tools/collect-statistics-frips.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:robert.player@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

doc: |
    Tool processes BAM file and output from SEACR to produce the FRIP (fraction of reads in peaks)
    statistic for ATAC-seq and cut&run type sequencing experiments. This stat is calculated for
    raw and spike-in normalized peak data, then concatentated to the *"_collected_statistics_report"
    files (md, tsv, and yaml).

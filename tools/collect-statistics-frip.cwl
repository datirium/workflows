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
      printf "$(date)\nLog file for collect-statistics-frip.cwl tool:\n\n"
      bam=$0; bed=$1; md=$2; tsv=$3; yaml=$4; spikein=$5
      # count of total aligned reads
      tar=$(samtools view -cF0x4 $bam)
      # order bed coordinates of max bedgraph signal (col6)
      cut -f6 $bed | sed -e 's/:/\t/' -e 's/-/\t/' | awk -F'\t' '{if($3>$2){printf("%s\t%s\t%s\n",$1,$2,$3)}else{printf("%s\t%s\t%s\n",$1,$3,$2)}}' > ordered.bed
      # counts of reads in peaks (split col6 due to start(col2) and end(col3) not always in ascending order - req by samtools)
      rip=$(samtools view -c $bam -L ordered.bed 2> /dev/null)
      # frip=rip/tar
      frip=$(printf $tar | awk -v rip=$rip '{printf("%.3f",rip/$0)}')
      # calculate mean max signal length (mmpl) from end-start sites
      mmpl=$(awk -F'\t' '{if($3>$2){x+=$3-$2}else{$2-$3}}END{printf("%.0f\n",x/NR)}' ordered.bed)
      printf "$tar, $rip, $frip, $mmpl\n"
      # concatenate frip, mmpl, and spikein read count onto md, tsv, and yaml files
      #   md
      cat $md > collected_statistics_report.md
      printf "-" >> collected_statistics_report.md
      printf "   fraction of (aligned) reads in peaks: $frip\n" >> collected_statistics_report.md
      printf "-" >> collected_statistics_report.md
      printf "   mean maximum signal length: $mmpl\n" >> collected_statistics_report.md
      printf "-" >> collected_statistics_report.md
      printf "   spike-in mapped read count (scaling_factor=10,000/x): $spikein\n" >> collected_statistics_report.md
      #   tsv
      headers=$(head -1 $tsv); data=$(tail -1 $tsv)
      printf "%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n" "$headers" "fraction of (aligned) reads in peaks" "mean maximum signal length" "spike-in mapped read count (scaling_factor=10,000/x)" "$data" "$frip" "$mmpl" "$spikein" > collected_statistics_report.tsv
      #   yaml
      cat $yaml > collected_statistics_report.yaml
      printf "  fraction of (aligned) reads in peaks: $frip\n" >> collected_statistics_report.yaml
      printf "  mean maximum signal length: $mmpl\n" >> collected_statistics_report.yaml
      printf "  spike-in mapped read count (scaling_factor=10,000/x): $spikein\n" >> collected_statistics_report.yaml
    inputBinding:
        position: 4

  bam_file:
    type: File
    label: "Input BAM file"
    inputBinding:
      position: 5
    doc: "Input BAM file, does not have to be coordinates sorted"

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

  spikein_reads_mapped:
      type: int
      label: "spike-in mapped reads from get_spikein_bam_statistics"
      inputBinding:
          position: 50


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

  log_file_stdout:
    type: File
    outputBinding:
      glob: "collected_stats_for_vis.log.stdout"
    doc: |
      log for stdout

  log_file_stderr:
    type: File
    outputBinding:
      glob: "collected_stats_for_vis.log.stderr"
    doc: |
      log for stderr


baseCommand: ["bash", "-c"]
stdout: 'collected_stats_for_vis.log.stdout'
stderr: 'collected_stats_for_vis.log.stderr'


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "collect-statistics-frip"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/collect-statistics-frip.cwl
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
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

doc: |
    Tool processes BAM file and output from SEACR to produce the FRIP (fraction of
    reads in peaks) and mean peak length statistics for ATAC-seq and cut&run type
    sequencing experiments. These stats are calculated for spike-in normalized 
    peak data, then concatentated to the *"_collected_statistics_report" files (md,
    tsv, and yaml). Additionally, a re-formatted peakcalled bed file is produced
    with headers per column, and "nearest gene" annotation per peak.

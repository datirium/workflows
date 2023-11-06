cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-bedops:v1.0.0


inputs:

  script_command:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for kallisto-index.cwl tool:\n\n"
      reference_bed=$0; all_other_bed=$1; set_operation=$2;
      printf "INPUTS:\n"
      printf "\$0 - $reference_bed\n"
      printf "\$1 - $all_other_bed\n"
      printf "\$2 - $set_operation\n"
      # commands start

      #   format tsv to bed files and sort bed for input into bedops
      tail -n+2 $reference_bed | cut -f6-8 | sort | uniq | bedtools sort -i - > groupA-list-1.sorted.bed
      x=2; for f in $(echo "$all_other_bed" | sed 's/,/\n/g'); do tail -n+2 $f | cut -f6-8 | sort | uniq | bedtools sort -i - > groupB-list-${x}.sorted.bed; ((x++)); done

      # get peak counts for list A and all unique peaks among list B group
      peak_count_A=$(tail -n+2 $reference_bed | wc -l)
      peak_count_B=$(find ./ -mindepth 1 -maxdepth 1 -name "groupB*.sorted.bed" -exec cat {} + | sort | uniq | wc -l)
      printf "\n\tStarting overview.md file...\n"
      printf "## Inputs:\n" > overview.md
      printf "-" >> overview.md
      printf " peak_list_A, $reference_bed\n" >> overview.md
      printf "-" >> overview.md
      printf " peak_list_B_group, $all_other_bed\n" >> overview.md
      printf "-" >> overview.md
      printf " set_operation, $set_operation\n" >> overview.md
      printf "\n" >> overview.md
      printf "## Results\n" >> overview.md
      printf "-" >> overview.md
      printf " Total peaks in list A: $peak_count_A\n" >> overview.md
      printf "-" >> overview.md
      printf " Total unique peaks in list B group: $peak_count_B\n" >> overview.md

      #   perform set operation with bedops
      # Intersection
      #       The --intersect operation determines genomic regions common to all input sets.
      if [[ "$set_operation" == "Intersection" ]]; then
        bedops --intersect groupA-list-1.sorted.bed $(find ./ -mindepth 1 -maxdepth 1 -name "groupB*.sorted.bed" | sort -V | sed 's/$/ /' | tr -d '\n' | sed 's/ $//') > output.bed
      fi
      # Union
      #       The --merge operation flattens all disjoint, overlapping, and adjoining element regions into contiguous, disjoint regions.
      if [[ "$set_operation" == "Union" ]]; then
        bedops --merge groupA-list-1.sorted.bed $(find ./ -mindepth 1 -maxdepth 1 -name "groupB*.sorted.bed" | sort -V | sed 's/$/ /' | tr -d '\n' | sed 's/ $//') > output.bed
      fi
      # Difference
      #       The --difference operation calculates the genomic regions found within the first (reference, e.g. A) input file, excluding regions in all other (B, C, +) input files.
      if [[ "$set_operation" == "Difference" ]]; then
        bedops --difference groupA-list-1.sorted.bed $(find ./ -mindepth 1 -maxdepth 1 -name "groupB*.sorted.bed" | sort -V | sed 's/$/ /' | tr -d '\n' | sed 's/ $//') > output.bed
      fi
      # Complement
      #       The --complement operation calculates the genomic regions in the gaps between the contiguous per-chromosome ranges defined by one or more inputs. The following example shows the use of two inputs.
      if [[ "$set_operation" == "Complement" ]]; then
        bedops --complement groupA-list-1.sorted.bed $(find ./ -mindepth 1 -maxdepth 1 -name "groupB*.sorted.bed" | sort -V | sed 's/$/ /' | tr -d '\n' | sed 's/ $//') > output.bed
      fi

      # add peak count after running set operation to overview file
      peak_set_count=$(awk 'END{print(NR)}' output.bed)
      printf "\tappending to overview.md file...\n"
      printf "-" >> overview.md
      printf " Peaks after running $set_operation: $peak_set_count\n" >> overview.md

      # format for IGV
      awk -F'\t' '{if($3>$2){printf("%s\t%.0f\t%.0f\t%s\n",$1,$2,$3,"peak_"NR)}}' output.bed > output-for-igv.tsv
      # format for island intersect tool
      awk -F'\t' 'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"};{if($3>$2){printf("%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t0\t0\t0\t%s\n",$1,$2,$3,$3-$2+1,$2+(($3-$2)/2),"0","peak_"NR)}}' output.bed > output-for-iaintersect.tsv
    inputBinding:
      position: 1

  peak_list_A:
    type: File
    inputBinding:
      position: 2
    doc: |
      filtered differential genelist from DESeq or diffbind pipelines

  peak_list_B_group:
    type: File[]
    inputBinding:
      position: 3
      itemSeparator: ","
    doc: |
      filtered differential genelist from DESeq or diffbind pipelines

  set_operator:
    type: string
    inputBinding:
      position: 4
    doc: |
      user selected set operation to perform and return


outputs:

  list_A_bed_file_for_igv:
    type: File
    outputBinding:
      glob: groupA-list-1.sorted.bed
    doc: |
      simple 3 column bed for IGV

  list_B_bed_array_for_igv:
    type: File[]
    outputBinding:
      glob: groupB*.sorted.bed
    doc: |
      simple 3 column bed for IGV

  overview_file:
    type: File
    outputBinding:
      glob: overview.md
    doc: |
      overview file for Overview tab in scidap containing inputs and metrics

  filtered_set_for_igv:
    type: File
    outputBinding:
      glob: output-for-igv.tsv
    doc: |
      peaks based on set operator chosen, formatted as headerless BED file with [chrom start end]

  filtered_set_for_iaintersect:
    type: File
    outputBinding:
      glob: output-for-iaintersect.tsv
    doc: |
      peaks based on set operator chosen, formatted as headered bed file for input into iaintersect.cwl to find nearest gene

  log_file_stdout:
    type: stdout

  log_file_stderr:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: filter-peaks_stdout.log
stderr: filter-peaks_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Set Operations for Called Peaks (ChIP/ATAC/C&R/diffbind) tool"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/filter-peaks-by-overlap.cwl
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
  This tool takes as input multiple peak list TSV files (the `iaintersect_result.tsv` output under the
  "Files" output tab) from the ChIP, ATAC, C&R, or diffbind workflows and performs the user-selected set
  operation on the group. Set operations include intersection, union, difference, and complement. See the
  tooltip for the `set_operator` input for more details.
      
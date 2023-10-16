cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-ubuntu22:v1.0.0


inputs:

  script_command:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for kallisto-index.cwl tool:\n\n"
      list1=$0; list2array=$1; set_operation=$2;
      printf "INPUTS:\n\n"
      printf "\$0 - $list1\n\n"
      printf "\$1 - $list2array\n\n"
      printf "\$2 - $set_operation\n\n"
      # commands start
      # zero out col5, then ensure there is only a single row per gene name, as some genes potentially have slightly different start/end positions (likely due to difference of 0- or 1-based upstream tools)
      #   listA
      awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"0.0",$6)}' $list1 | sort | uniq > list1.tmpx
      awk -F'\t' '{split($4,col4,","); for(i in col4){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,col4[i],$5,$6)}}' list1.tmpx > list1.tmp
      awk -F'\t' '{if(NR==FNR){c1[$4]=$1; c2[$4]=$2; c3[$4]=$3; c6[$4]=$6}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",c1[$4],c2[$4],c3[$4],$4,"0.0",c6[$4])}}' list1.tmp list1.tmp | sort | uniq > list1.tsv
      #   groupB, only keep unique rows (score values will be lost) per input list
      rep=1; echo "$list2array" | sed 's/,/\n/g' | while read filepath; do
        awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"0.0",$6)}' $filepath | sort | uniq > list2.tmpx
        awk -F'\t' '{split($4,col4,","); for(i in col4){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,col4[i],$5,$6)}}' list2.tmpx > list2.tmp
        awk -F'\t' '{if(NR==FNR){c1[$4]=$1; c2[$4]=$2; c3[$4]=$3; c6[$4]=$6}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",c1[$4],c2[$4],c3[$4],$4,"0.0",c6[$4])}}' list2.tmp list2.tmp | sort | uniq > groupB_list${rep}.tsv
        ((rep++))
      done

      # Intersection
      #       list of genes shared between all input lists (have to loop through each list, compare to list A)
      if [[ "$set_operation" == "Intersection" ]]; then
        cp list1.tsv int_list1.tmp
        find ./ -mindepth 1 -maxdepth 1 -name "groupB_list*.tsv" | sort -V | while read rep; do
          cp $rep list2.tmp
          comm -12 <(cut -f4 int_list1.tmp | sort) <(cut -f4 list2.tmp | sort) | while read gene; do grep -P "\t$gene\t" int_list1.tmp; done > int_list1.out
          mv int_list1.out int_list1.tmp
        done
        mv int_list1.tmp intersection.tsv
        #       reclaim score for column 5 from file A (list 1)
        awk -F'\t' '{if(NR==FNR){score[$4]=$5}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,score[$4],$6)}}' $list1 intersection.tsv > genelist-filtered-set.bed
      fi

      # Union
      #       list of genes unique among all genes of all input lists
      if [[ "$set_operation" == "Union" ]]; then
        # concatenate everything from all group B lists
        find ./ -mindepth 1 -maxdepth 1 -name "groupB_list*.tsv" -exec cat {} + | sort | uniq > list2.tsv
        # concatenate all lists
        cp list1.tsv union.tmp
        cat list2.tsv >> union.tmp
        awk -F'\t' '{if(NR==FNR){c1[$4]=$1; c2[$4]=$2; c3[$4]=$3; c6[$4]=$6}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",c1[$4],c2[$4],c3[$4],$4,"0.0",c6[$4])}}' union.tmp union.tmp | sort | uniq > union.tsv
        #       reclaim score for column 5 from list1 genes (for overlaps, scores from list A reported)
        awk -F'\t' '{if(NR==FNR){score[$4]=$5}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,score[$4],$6)}}' $list1 union.tsv > genelist-filtered-set.bed
      fi

      # Symmetric Difference (removing for now, check previous commits if need code)
      #       list of genes that would be left out of an Intersection

      # Relative Complement (the relative completment of set B would be denoted "A / B")
      #       list of genes that are in set A and not in group set B.
      #       so the RC of list 2 are genes in list 1 that are not in list2
      if [[ "$set_operation" == "Relative_Complement" ]]; then
        # concatenate everything from all group B lists
        find ./ -mindepth 1 -maxdepth 1 -name "groupB_list*.tsv" -exec cat {} + | sort | uniq > list2.tsv
        # get rel comp between the 2 lists
        comm -23 <(cut -f4 list1.tsv | sort) <(cut -f4 list2.tsv | sort) | while read gene; do grep -P "\t$gene\t" list1.tsv; done | sort | uniq > relative_complement.tsv
        #       reclaim score for column 5 from list1 genes
        awk -F'\t' '{if(NR==FNR){score[$4]=$5}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,score[$4],$6)}}' $list1 relative_complement.tsv > genelist-filtered-set.bed
      fi

      # make headered bed file for scidap table interface
      printf "Chrom\tStart\tEnd\tName\tScore\tStrand\n" > genelist-filtered-set-w-header.bed
      cat genelist-filtered-set.bed >> genelist-filtered-set-w-header.bed
    inputBinding:
      position: 1

  filtered_list_A:
    type: File
    inputBinding:
      position: 2
    doc: |
      filtered differential genelist from DESeq or diffbind pipelines

  filtered_list_B_group:
    type: File[]
    inputBinding:
      position: 3
      itemSeparator: ","
    doc: |
      filtered differential genelist from DESeq or diffbind pipelines

  set_operation:
    type: string
    inputBinding:
      position: 4
    doc: |
      user selected set operation to perform and return


outputs:

  genelist_filtered_set:
    type: File
    outputBinding:
      glob: genelist-filtered-set.bed
    doc: |
      filtered gene list based on set operator chosen, formatted as headerless BED file with [chrom start end name score strand]

  genelist_filtered_set_with_header:
    type: File
    outputBinding:
      glob: genelist-filtered-set-w-header.bed
    doc: |
      filtered gene list based on set operator chosen, formatted as headered bed file for scidap table interface [chrom start end name score strand]

  log_file_stdout:
    type: stdout

  log_file_stderr:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: genelists-sets_stdout.log
stderr: genelists-sets_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Genelists Set Operations"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/genelists-sets.cwl
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
  This tool takes as input 2 filtered genelists samples and performs the user input set operation on them.
  The output is a single filtered gene list in the same format as the input files (headerless BED file with [chrom start end name score strand]).
  The returned score value (column 5) is always derived from file A.
      
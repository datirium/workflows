cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kallisto:stable


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for kallisto-quant.cwl tool:\n\n"
      INDEX=$0; ANNO=$1; R1=$2; R2=$3; THREADS=$4;
      printf "INPUTS:\n\n"
      printf "\$0 - $INDEX\n\n"
      printf "\$1 - $ANNO\n\n"
      printf "\$2 - $R1\n\n"
      printf "\$3 - $R2\n\n"
      printf "\$4 - $THREADS\n\n"
      # commands start
      kallisto quant -t $THREADS -i $INDEX -o quant_outdir $R1 $R2
      # format output for as deseq upstream (e.g. rpkm_isoforms_cond_1, rpkm_genes_cond_1, rpkm_common_tss_cond_1), only using "genes" in this case
      # using kallisto's "est_counts" output (col4 in abundance.tsv) counts per transcript (as required/expect by deseq tool for diffexp analysis)
      printf "RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm\n" > transcript_counts.tsv
      #   force "est_counts" to integers
      awk -F'\t' '{if(NR==FNR){anno[$3]=$0}else{printf("%s\t%0.f\t%s\n",anno[$1],$4,$5)}}' $ANNO <(tail -n+2 ./quant_outdir/abundance.tsv) >> transcript_counts.tsv

      # print for overview.md
      #   read metrics
      total_aligned=$(tail -n+2 transcript_counts.tsv | awk -F'\t' '{x+=$7}END{printf("%0.f",x)}')
      annotated_aligned=$(tail -n+2 transcript_counts.tsv | grep -v "^na" | awk -F'\t' '{x+=$7}END{printf("%0.f",x)}')
      unannotated_aligned=$(tail -n+2 transcript_counts.tsv | grep "^na" | awk -F'\t' '{x+=$7}END{printf("%0.f",x)}')
      if [[ $(basename $R1 | sed 's/.*\.//') == "gz" ]]; then read_count_r1=$(gunzip -c $R1 | wc -l | awk '{print($0/4)}'); else read_count_r1=$(wc -l $R1 | awk '{print($0/4)}'); fi
      if [[ $(basename $R2 | sed 's/.*\.//') == "gz" ]]; then read_count_r2=$(gunzip -c $R2 | wc -l | awk '{print($0/4)}'); else read_count_r2=$(wc -l $R2 | awk '{print($0/4)}'); fi
      unmapped=$(printf "$read_count_r1" | awk -v x="$total_aligned" '{print($0-x)}')

      #   output stats for pie chart
      printf "\n\n\tgenerating pie_stats.tsv file...\n"
      printf "Total reads\tAnnotated transcript reads\tUnannotated transcript reads\tUnmapped reads\n" > pie_stats.tsv
      printf "$read_count_r1\t$annotated_aligned\t$unannotated_aligned\t$unmapped\n" >> pie_stats.tsv

      #   transcript metrics
      transcriptome_count=$(tail -n+2 transcript_counts.tsv | wc -l)
      transcriptome_gt0=$(tail -n+2 transcript_counts.tsv | awk -F'\t' '{if($7>0){x+=1}}END{printf("%0.f\n",x)}')
      annotated_gt0=$(tail -n+2 transcript_counts.tsv | grep -v "^na" | awk -F'\t' '{if($7>0){x+=1}}END{printf("%0.f\n",x)}')
      annotated_eq0=$(tail -n+2 transcript_counts.tsv | grep -v "^na" | awk -F'\t' '{if($7==0){x+=1}}END{printf("%0.f\n",x)}')
      unannotated_gt0=$(tail -n+2 transcript_counts.tsv | grep "^na" | awk -F'\t' '{if($7>0){x+=1}}END{printf("%0.f\n",x)}')
      unannotated_eq0=$(tail -n+2 transcript_counts.tsv | grep "^na" | awk -F'\t' '{if($7==0){x+=1}}END{printf("%0.f\n",x)}')

      #   format for overview file
      printf "\n\n\tgenerating overview.md file...\n"
      printf "#### INPUTS\n" > overview.md
      printf "-" >> overview.md
      printf " \$INDEX, $INDEX\n" >> overview.md
      printf "-" >> overview.md
      printf " \$ANNO, $ANNO\n" >> overview.md
      printf "-" >> overview.md
      printf " \$R1, $R1\n" >> overview.md
      printf "-" >> overview.md
      printf " \$R2, $R2\n" >> overview.md

      printf "\n" >> overview.md

      printf "#### METRICS\n" >> overview.md
      printf "-" >> overview.md
      printf " Total transcripts in transcriptome: $transcriptome_count\n" >> overview.md
      printf "-" >> overview.md
      printf " Total transcripts with at least 1 read aligned: $transcriptome_gt0\n" >> overview.md
      printf "-" >> overview.md
      printf " Total annotated transcripts with at least 1 read aligned: $annotated_gt0\n" >> overview.md
      printf "-" >> overview.md
      printf " Total unannotated transcripts with at least 1 read aligned: $unannotated_gt0\n" >> overview.md
      printf "-" >> overview.md
      printf " R1 read count: $read_count_r1\n" >> overview.md
      printf "-" >> overview.md
      printf " R2 read count: $read_count_r2\n" >> overview.md
      printf "-" >> overview.md
      printf " Estimated aligned read count, transcriptome: $total_aligned\n" >> overview.md
      printf "-" >> overview.md
      printf " Estimated aligned read count, annotated transcriptome: $annotated_aligned\n" >> overview.md
      printf "-" >> overview.md
      printf " Estimated aligned read count, unannotated transcriptome: $unannotated_aligned\n" >> overview.md

      printf "\n\nWorkflow script complete!\n"
    inputBinding:
      position: 1

  kallisto_index:
    type: File
    inputBinding:
      position: 2
    doc: |
      Kallisto index file

  annotation_tsv:
    type: File
    inputBinding:
      position: 3
    doc: |
      TSV file containing gene annotations for the reference genome. From kallisto index upstream.
      Required columns (include headers as row 1 of TSV): RefseqId, GeneId, Chrom (transcript id/name), TxStart (start of alignment in query), TxEnd (end of alignment in query), Strand (if query start < query end strand +, else -).

  fastq_R1:
    type: File
    label: "R1 fastq"
    inputBinding:
      position: 6
    doc: |
      FASTQ file 1 of paired end read data.

  fastq_R2:
    type: File
    label: "R2 fastq"
    inputBinding:
      position: 7
    doc: |
      FASTQ file 2 of paired end read data.

  threads:
    type: int
    label: "threads"
    inputBinding:
      position: 10
    doc: |
      Number of threads for steps that support multithreading.


outputs:

  overview:
    type: File
    outputBinding:
      glob: overview.md
    doc: |
      Summary of script run and some alignment metrics.

  pie_stats:
    type: File
    outputBinding:
      glob: pie_stats.tsv
    doc: |
      Summary of script run and some alignment metrics.

  kallisto_abundance_file:
    type: File
    outputBinding:
      glob: quant_outdir/abundance.tsv

  kallisto_runinfo_file:
    type: File
    outputBinding:
      glob: quant_outdir/run_info.json

  transcript_counts:
    type: File
    outputBinding:
      glob: transcript_counts.tsv
    doc: |
      Gene expression table formatted for input into DESeq

  log_file_stdout:
    type: File
    outputBinding:
      glob: "log.stdout"
    doc: |
      log for stdout

  log_file_stderr:
    type: File
    outputBinding:
      glob: "log.stderr"
    doc: |
      log for stderr


baseCommand: ["bash", "-c"]
stdout: 'log.stdout'
stderr: 'log.stderr'


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "kallisto-index"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/kallisto-index.cwl
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
    Tool indexes a reference genome fasta using `kallisto index`.

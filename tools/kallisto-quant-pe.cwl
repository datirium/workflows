cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: scrowley1/scidap-kallisto:v1.0.0
inputs:
  script:
    type: string?
    default: "#!/bin/bash\nprintf \"$(date)\\n\\nStdout log file for kallisto-quant.cwl tool:\\n\\n\"\nINDEX=$0; ANNO=$1; R1=$2; R2=$3; THREADS=$4;\nprintf \"INPUTS:\\n\\n\"\nprintf \"\\$0 - $INDEX\\n\\n\"\nprintf \"\\$1 - $ANNO\\n\\n\"\nprintf \"\\$2 - $R1\\n\\n\"\nprintf \"\\$3 - $R2\\n\\n\"\nprintf \"\\$4 - $THREADS\\n\\n\"\n# commands start\nkallisto quant -t $THREADS -i $INDEX -o quant_outdir $R1 $R2\n# format output for use as deseq upstream (e.g. rpkm_isoforms_cond_1, rpkm_genes_cond_1, rpkm_common_tss_cond_1), only using \"genes\" in this case\n#   original file (works with standard deseq) - transcript_counts.tsv\n#   reformatted file (for deseq multi-factor) - transcript_counts_mf.tsv\n# using kallisto's \"est_counts\" output (col4 in abundance.tsv) counts per transcript (as required/expect by deseq tool for diffexp analysis)\nprintf \"RefseqId\\tGeneId\\tChrom\\tTxStart\\tTxEnd\\tStrand\\tTotalReads\\tRpkm\\n\" > transcript_counts.tsv\n#   force \"est_counts\" to integers\nawk -F'\\t' '{if(NR==FNR){anno[$3]=$0}else{printf(\"%s\\t%0.f\\t%s\\n\",anno[$1],$4,$5)}}' $ANNO <(tail -n+2 ./quant_outdir/abundance.tsv) >> transcript_counts.tsv\n\n# making reformatted file for deseq multi-factor (removing unannotated transcripts labeled as \"na\" for col1 [RefseqId] and col2 [GeneId] from the output count table)\nprintf \"RefseqId\\tGeneId\\tChrom\\tTxStart\\tTxEnd\\tStrand\\tTotalReads\\tRpkm\\n\" > transcript_counts_mf.tsv\ntail -n+2 transcript_counts.tsv | grep -vP \"^na\\tna\\t\" >> transcript_counts_mf.tsv\n# convert to csv or 'get_gene_n_tss.R'\nsed 's/\\t/,/g' transcript_counts_mf.tsv > transcript_counts_mf.csv\n\n# REMOVING THIS FOR NOW, REPLACING WITH \n# and if there are duplicate geneIds, only retain the one with the higher read count\n#awk -F'\\t' '{if(NR==FNR){if($7>=tr[$2]){c1[$2]=$1; c3[$2]=$3; c4[$2]=$4; c5[$2]=$5; c6[$2]=$6; tr[$2]=$7; rpkm[$2]=$8}}else{printf(\"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%0.f\\t%s\\n\",c1[$2],$2,c3[$2],c4[$2],c5[$2],c6[$2],tr[$2],rpkm[$2])}}' transcript_counts_mf.tmp transcript_counts_mf.tmp | sort | uniq >> transcript_counts_mf.tsv\n\n# print for overview.md\n#   read metrics\ntotal_aligned=$(tail -n+2 transcript_counts.tsv | awk -F'\\t' '{x+=$7}END{printf(\"%0.f\",x)}')\nannotated_aligned=$(tail -n+2 transcript_counts.tsv | grep -v \"^na\" | awk -F'\\t' '{x+=$7}END{printf(\"%0.f\",x)}')\nunannotated_aligned=$(tail -n+2 transcript_counts.tsv | grep \"^na\" | awk -F'\\t' '{x+=$7}END{printf(\"%0.f\",x)}')\nif [[ $(basename $R1 | sed 's/.*\\.//') == \"gz\" ]]; then read_count_r1=$(gunzip -c $R1 | wc -l | awk '{print($0/4)}'); else read_count_r1=$(wc -l $R1 | awk '{print($0/4)}'); fi\nif [[ $(basename $R2 | sed 's/.*\\.//') == \"gz\" ]]; then read_count_r2=$(gunzip -c $R2 | wc -l | awk '{print($0/4)}'); else read_count_r2=$(wc -l $R2 | awk '{print($0/4)}'); fi\nunmapped=$(printf \"$read_count_r1\" | awk -v x=\"$total_aligned\" '{print($0-x)}')\n\n#   output stats for pie chart\nprintf \"\\n\\n\\tgenerating pie_stats.tsv file...\\n\"\nprintf \"Total reads\\tAnnotated transcript reads\\tUnannotated transcript reads\\tUnmapped reads\\n\" > pie_stats.tsv\nprintf \"$read_count_r1\\t$annotated_aligned\\t$unannotated_aligned\\t$unmapped\\n\" >> pie_stats.tsv\n\n#   transcript metrics\ntranscriptome_count=$(tail -n+2 transcript_counts.tsv | wc -l)\ntranscriptome_gt0=$(tail -n+2 transcript_counts.tsv | awk -F'\\t' '{if($7>0){x+=1}}END{printf(\"%0.f\\n\",x)}')\nannotated_gt0=$(tail -n+2 transcript_counts.tsv | grep -v \"^na\" | awk -F'\\t' '{if($7>0){x+=1}}END{printf(\"%0.f\\n\",x)}')\nannotated_eq0=$(tail -n+2 transcript_counts.tsv | grep -v \"^na\" | awk -F'\\t' '{if($7==0){x+=1}}END{printf(\"%0.f\\n\",x)}')\nunannotated_gt0=$(tail -n+2 transcript_counts.tsv | grep \"^na\" | awk -F'\\t' '{if($7>0){x+=1}}END{printf(\"%0.f\\n\",x)}')\nunannotated_eq0=$(tail -n+2 transcript_counts.tsv | grep \"^na\" | awk -F'\\t' '{if($7==0){x+=1}}END{printf(\"%0.f\\n\",x)}')\n\n#   format for overview file\nprintf \"\\n\\n\\tgenerating overview.md file...\\n\"\n\nprintf \"-\" > overview.md\nprintf \" NOTE: Unannotated transcripts are not shown in the gene expression tab, and will not be used in downstream differential expression analysis.\\n\" >> overview.md\nprintf \"\\n\" >> overview.md\nprintf \"#### INPUTS\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" \\$INDEX, $INDEX\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" \\$ANNO, $ANNO\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" \\$R1, $R1\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" \\$R2, $R2\\n\" >> overview.md\n\nprintf \"\\n\" >> overview.md\n\nprintf \"#### METRICS\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Total transcripts in transcriptome: $transcriptome_count\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Total transcripts with at least 1 read aligned: $transcriptome_gt0\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Total annotated transcripts with at least 1 read aligned: $annotated_gt0\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Total unannotated transcripts with at least 1 read aligned: $unannotated_gt0\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" R1 read count: $read_count_r1\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" R2 read count: $read_count_r2\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Estimated aligned read count, transcriptome: $total_aligned\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Estimated aligned read count, annotated transcriptome: $annotated_aligned\\n\" >> overview.md\nprintf \"-\" >> overview.md\nprintf \" Estimated aligned read count, unannotated transcriptome: $unannotated_aligned\\n\" >> overview.md\n\nprintf \"\\n\\nWorkflow script complete!\\n\"\n"
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
    label: R1 fastq
    inputBinding:
      position: 6
    doc: |
      FASTQ file 1 of paired end read data.
  fastq_R2:
    type: File
    label: R2 fastq
    inputBinding:
      position: 7
    doc: |
      FASTQ file 2 of paired end read data.
  threads:
    type: int
    label: threads
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
      glob: transcript_counts_mf.csv
    doc: |
      Gene expression table formatted for input into DESeq and DESeq multi factor. The na values for unannotated genes have been removed.
  transcript_counts_standard:
    type: File
    outputBinding:
      glob: transcript_counts.tsv
    doc: |
      Gene expression table formatted for input into DESeq (contains na values where annotations are not found). This particular output works with standard deseq, the multi-factor output will be used for both.
  log_file_stdout:
    type: File
    outputBinding:
      glob: log.stdout
    doc: |
      log for stdout
  log_file_stderr:
    type: File
    outputBinding:
      glob: log.stderr
    doc: |
      log for stderr
baseCommand:
- bash
- -c
stdout: log.stdout
stderr: log.stderr
doc: |-
  This tool quantitates RNA-Seq reads using the pseudo aligner Kallisto. Read counts are estimates.

  ### __Data Analysis Steps__
  1. cwl calls dockercontainer robertplayer/scidap-kallisto to pseudo align reads using `kallisto quant`
  2. abundance tsv is merged with annotation data
  3. rows are summed per unique GeneID, and resulting file formatted for use in differential expression analysis
  4. read and alignment metrics are calculated for the sample piechart, and output to the overview.md file
label: kallisto-index

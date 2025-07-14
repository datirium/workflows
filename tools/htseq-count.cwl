cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/htseq:0.13.5--py38h1773678_0
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      shopt -s nocaseglob
      set -- "$0" "$@"
      echo "Run htseq-count with the following parameters"
      for i in "$@";
        do echo $i;
      done;
      htseq-count -f bam -r pos "$@" 1> temp.tsv 2> feature_counts_stderr.log
      cat temp.tsv | grep "^__" > feature_counts_stdout.log
      cat temp.tsv | grep -v "^__" > feature_counts_report.tsv
      rm -f temp.tsv
    inputBinding:
      position: 1
    doc: |
      Script runs htseq-count with the provided parameters splitting
      the report file into summary and the actual counts data.
  alignment_bam_file:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 5
    doc: |
      Path to the coordinate sorted indexed BAM file
  annotation_gtf_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      GTF annotation file
  strand_specific:
    type:
    - 'null'
    - type: enum
      symbols:
      - 'yes'
      - 'no'
      - reverse
    inputBinding:
      position: 7
      prefix: -s
    doc: |
      Whether the data is from a strand-specific assay. For stranded=no, a read is
      considered overlapping with a feature regardless of whether it is mapped to
      the same or the opposite strand as the feature. For stranded=yes and single-end
      reads, the read has to be mapped to the same strand as the feature. For paired-end
      reads, the first read has to be on the same strand and the second read on the
      opposite strand. For stranded=reverse, these rules are reversed.
      Default: "yes"
  feature_type:
    type:
    - 'null'
    - type: enum
      symbols:
      - CDS
      - exon
      - start_codon
      - stop_codon
      - transcript
    inputBinding:
      position: 8
      prefix: -t
    doc: |
      Feature type (3rd column in GFF file) to be used,
      all features of other type are ignored.
      Default: exon
  feature_id:
    type:
    - 'null'
    - type: enum
      symbols:
      - gene_id
      - transcript_id
      - exon_id
      - gene_name
    inputBinding:
      position: 9
      prefix: -i
    doc: |
      GTF attribute to be used as feature ID. Several GTF lines with the same feature
      ID will be considered as parts of the same feature. The feature ID is used to
      identity the counts in the output table.
      FYI: no use of having here "exon_number", so it was excluded from enum
      Default: gene_id
  additional_id:
    type:
    - 'null'
    - type: enum
      symbols:
      - gene_id
      - transcript_id
      - exon_number
      - exon_id
      - gene_name
    inputBinding:
      position: 10
      prefix: --additional-attr
    doc: |
      Additional feature attributes, which will be printed as an additional
      column after the primary attribute column but before the counts column.
      Default: none
  overlapping_mode:
    type:
    - 'null'
    - type: enum
      symbols:
      - union
      - intersection-strict
      - intersection-nonempty
    inputBinding:
      position: 11
      prefix: -m
    doc: |
      Mode to handle reads overlapping more than one feature
      Default: union
  nonunique_mode:
    type:
    - 'null'
    - type: enum
      symbols:
      - none
      - all
      - fraction
      - random
    inputBinding:
      position: 12
      prefix: --nonunique
    doc: |
      Mode to handle reads that align to or are assigned to more than
      one feature in the selected overlapping_mode
      Options:
        "none": the read (or read pair) is counted as ambiguous and not
          counted for any features. Also, if the read (or read pair) aligns
          to more than one location in the reference, it is scored as
          alignment_not_unique.
        "all": the read (or read pair) is counted as ambiguous and is also
          counted in all features to which it was assigned. Also, if the read
          (or read pair) aligns to more than one location in the reference, it
          is scored as alignment_not_unique and also separately for each location.
        "fraction": the read (or read pair) is counted as ambiguous and is also
          counted fractionally in all features to which it was assigned. For
          example, if the read overlaps with 3 features, it will be counted 1/3
          to each of them.
        "random": the read (or read pair) is counted as ambiguous and is also
          counted uniformly at random to one of the features to which it
          was assigned.
      Default: "none"
  secondary_alignments_mode:
    type:
    - 'null'
    - type: enum
      symbols:
      - score
      - ignore
    inputBinding:
      position: 13
      prefix: --secondary-alignments
    doc: |
      Mode to handle secondary alignments (SAM flag 0x100).
      Can be "score" and "ignore"
      Default: "score"
  supplementary_alignments_mode:
    type:
    - 'null'
    - type: enum
      symbols:
      - score
      - ignore
    inputBinding:
      position: 14
      prefix: --supplementary-alignments
    doc: |
      Mode to handle supplementary/chimeric alignments (SAM flag 0x800).
      Can be "score" and "ignore"
      Default: "score"
outputs:
  feature_counts_report_file:
    type: File
    outputBinding:
      glob: feature_counts_report.tsv
  stdout_log:
    type: File
    outputBinding:
      glob: feature_counts_stdout.log
  stderr_log:
    type: File
    outputBinding:
      glob: feature_counts_stderr.log
baseCommand:
- bash
- -c
label: 'HTSeq: Analysing high-throughput sequencing data'
doc: |
  For convenience to use in the workflow that sort and index BAM files by coordinate
  this tools expects coordinate sorted and indexed BAM file as input. For single-read
  dat it won't influence on anything, for paired-end the more memory will be used to
  keep reads while looking for their proper pairs (see --max-reads-in-buffer parameter).

  Current limitations:
    - only one `--additional-attr` is supported
    - skip `--nprocesses` parameter as it's not helpful when we use only one input BAM file

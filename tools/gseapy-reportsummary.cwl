cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for bash script in seacr.cwl tool:\n\n"
      printf "INPUTS:\n"
      echo "\$0 - $0"
      echo "\$1 - $1"
      echo "\$2 - $2"
      # inputs
      deseq_counts="$0"
      deseq_phenotypes="$1"
      report="$2"
      #   DATASET DETAILS
      # total unique genes in rna-seq analysis (remove first 3 rows of header info)
      total_unique_genes=$(cut -f1 $deseq_counts | tail -n+4 | sort | uniq | awk 'END{print(NR)}')
      # total unique genes in marker gene datasets
      feature_genes=$(cut -f8 $report | tail -n+2 | sed 's/;/\n/g' | sort | uniq | awk 'END{print(NR)}')
      # total unique genes among all leading edge genes of all marker datasets
      ledge_genes=$(cut -f9 $report | tail -n+2 | sed 's/;/\n/g' | sort | uniq | awk 'END{print(NR)}')
      # get number of gene sets in file (remove first 1 rows of header info)
      total_gene_sets=$(tail -n+2 $report | awk 'END{print(NR)}')
      #   FORMATTING OUTPUT
      printf "#### DATA DETAILS\n" > reportsummary.md
      printf "-" >> reportsummary.md
      printf " The RNA-Seq sample contained $total_unique_genes total unique genes\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " There are $total_gene_sets marker gene datasets analyzed\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " The marker gene datasets have $feature_genes features (genes)\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " There are $ledge_genes leading edge subset genes among all marker gene datasets analyzed\n" >> reportsummary.md
      printf "\n" >> reportsummary.md
      #   PHENOTYPE AND ANALYSIS DETAILS
      # get phenotype names
      p1=$(tail -1 $deseq_phenotypes | cut -f2 -d$' ')
      p2=$(tail -1 $deseq_phenotypes | cut -f3 -d$' ')
      echo $p1 $p2
      # phenotype 1
      enriched=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1>0){x++}}END{print(x)}')
      sigsets_fdr25=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1>0){if($4<0.25){x++}}}END{print(x)}')
      sigsets_pv5=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1>0){if($3<0.05){x++}}}END{print(x)}')
      sigsets_pv1=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1>0){if($3<0.01){x++}}}END{print(x)}')
      printf "#### Enrichment in phenotype: $p1\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " $enriched / $total_gene_sets marker gene sets are enriched in phenotype $p1\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " %s %s\n" "$sigsets_fdr25" "marker gene sets are significant at FDR <25%" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " %s %s\n" "$sigsets_pv5" "marker gene sets are significantly enriched at nominal pvalue <5%" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " %s %s\n" "$sigsets_pv1" "marker gene sets are significantly enriched at nominal pvalue <1%" >> reportsummary.md
      printf "\n" >> reportsummary.md
      # phenotype 2
      enriched=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1<0){x++}}END{print(x)}')
      sigsets_fdr25=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1<0){if($4<0.25){x++}}}END{print(x)}')
      sigsets_pv5=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1<0){if($3<0.05){x++}}}END{print(x)}')
      sigsets_pv1=$(cut -f2-5 $report | tail -n+2 | awk -F'\t' 'BEGIN{x=0};{if($1<0){if($3<0.01){x++}}}END{print(x)}')
      printf "#### Enrichment in phenotype: $p2\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " $enriched / $total_gene_sets marker gene sets are enriched in phenotype $p2\n" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " %s %s\n" "$sigsets_fdr25" "marker gene sets are significant at FDR <25%" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " %s %s\n" "$sigsets_pv5" "marker gene sets are significantly enriched at nominal pvalue <5%" >> reportsummary.md
      printf "-" >> reportsummary.md
      printf " %s %s\n" "$sigsets_pv1" "marker gene sets are significantly enriched at nominal pvalue <1%" >> reportsummary.md
    inputBinding:
        position: 1

  read_counts_file:
    type: File?
    inputBinding:
      position: 2
    doc: |
      Input gene expression dataset file in txt or gct
      format produced from deseq_experiment upstream.

  phenotypes_file:
    type: File?
    inputBinding:
      position: 3
    doc: |
      Input class vector (phenotype) file in CLS format.

  enrichment_report:
    type: File?
    inputBinding:
      position: 4
    doc: |
      Gene set enrichment report from gseapt.cwl tool.

outputs:

  summary_file:
    type: File
    outputBinding:
      glob: "reportsummary.md"
    doc: |
      Summary report file of GSEAPY results.

  log_file_stderr:
    type: File
    outputBinding:
      glob: "reportsummary.stderr"
    doc: |
      log for BASH stderr

  log_file_stdout:
    type: File
    outputBinding:
      glob: "reportsummary.stdout"
    doc: |
      log for BASH stdout


baseCommand: ["bash", "-c"]
stderr: "reportsummary.stderr"
stdout: "reportsummary.stdout"


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "gseapy-reportsummary"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/gseapy-reportsummary.cwl
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
  Tool runs a custom BASH command to process GSEAPY report
  and associated input files (gene_counts and phenotypes)
  to produce a summary report file of GSEA results similar
  to the one here: http://diverge.hunter.cuny.edu/~weigang/silac-chromatin-gsea/
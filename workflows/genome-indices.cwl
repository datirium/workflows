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

  fasta_ribosomal:
    type: File?
    format: "http://edamontology.org/format_1929"
    label: "Ribosomal DNA file (*.fasta, *.fa)"
    doc: "Ribosomal DNA file (*.fasta, *.fa). Default: hg19"

  chromosome_list:
    type:
      - "null"
      - string[]
    default:
      [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM"
      ]
    label: "Chromosome list to be included into the reference genome FASTA file"
    doc: "Filter chromosomes while extracting FASTA from 2bit"
  
  effective_genome_size:
    type: string
    label: "Effective genome size"
    doc: "MACS2 effective genome sizes: hs, mm, ce, dm or number, for example 2.7e9"

  gtf_annotation:
    type: File?
    format: "http://edamontology.org/format_2306"
    label: "GTF annotation file (gzip compressed, from Gencode)"
    doc: |
      GTF genome annotation file. Primary assembly with
      the reference chromosomes and scaffolds prefixed
      with "chr", excluding the haplotypes and patches.
      This annotation file should include chrM. It can
      be downloaded from the Gencode database (for example,
      gencode.v44.primary_assembly.annotation.gtf.gz).
      If both GTF and TSV annotation files are provided,
      the GTF input will have a higher priority.

  remove_par_y_genes:
    type: boolean?
    default: false
    label: "Exclude PAR locus genes from the GTF annotation file (for hg38 only)"
    doc: |
      Since Ensembl 110, the PAR locus genes
      are included on chrY as copies of chrX.
      Removing the chrY PAR genes is desirable
      so they do not end up as extra entries in
      the output. In accordance to the 10x
      guidelines, it's needed only for hg38.
      Ignored if annotation is loaded not from
      the GTF file.

  annotation_tab:
    type: File?
    format: "http://edamontology.org/format_3475"
    label: "Compressed tsv.gz annotation file"
    doc: |
      Compressed tab-separated annotation file.
      Shouldn't include chrM as it will be loaded
      separately. Ignored if annotation is loaded
      from the GTF file.

  mitochondrial_annotation_tab:
    type: File?
    format: "http://edamontology.org/format_3475"
    label: "Compressed tsv.gz mitochondrial DNA annotation file"
    doc: |
      Compressed mitochondrial DNA tab-separated
      annotation file. Should include chrM as only
      this chromosome will be loaded from it.
      Ignored if annotation is loaded from the
      GTF file.

  cytoband:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Compressed cytoBand file for IGV browser"
    doc: "Compressed tab-separated cytoBand file for IGV browser"

  genome_sa_index_n_bases:
    type: int?
    label: "Length of SA pre-indexing string for reference genome indices"
    doc: |
      Length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory,
      but allow faster searches. For small genomes, the parameter –genomeSAindexNbases must be scaled down to
      min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 9, for 100 kiloBase genome,
      this is equal to 7.
      default: 14
    'sd:layout':
      advanced: true

  genome_sa_index_n_bases_mitochondrial:
    type: int?
    label: "Length of SA pre-indexing string for mitochondrial DNA indices"
    doc: |
      Length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory,
      but allow faster searches. For small genomes, the parameter –genomeSAindexNbases must be scaled down to
      min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 9, for 100 kiloBase genome,
      this is equal to 7.
      default: 14
    'sd:layout':
      advanced: true

  genome_chr_bin_n_bits:
    type: int?
    label: "Number of bins allocated for each chromosome of reference genome"
    doc: |
      If you are using a genome with a large (>5,000) number of references (chrosomes/scaﬀolds), you may need to reduce the
      --genomeChrBinNbits to reduce RAM consumption. For a genome with large number of contigs, it is recommended to scale
      this parameter as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)]).
      default: 18
    'sd:layout':
      advanced: true

  genome_sa_sparse_d:
    type: int?
    label: "Suffix array sparsity for reference genome and mitochondrial DNA indices"
    doc: |
      Suffix array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAMat the cost of mapping speed reduction"
    'sd:layout':
      advanced: true

  limit_genome_generate_ram:
    type: long?
    label: "Limit maximum available RAM (bytes) for reference genome indices generation"
    doc: "Maximum available RAM (bytes) for genome generation. Default 31000000000"
    'sd:layout':
      advanced: true


outputs:

  star_indices:
    type: Directory
    outputSource: star_generate_indices/indices_folder
    label: "STAR genome indices"
    doc: "STAR generated genome indices folder"

  star_indices_stdout_log:
    type: File
    outputSource: star_generate_indices/stdout_log
    label: "STAR stdout log for genome indices"
    doc: "STAR generated stdout log for genome indices"
    
  star_indices_stderr_log:
    type: File
    outputSource: star_generate_indices/stderr_log
    label: "STAR stderr log for genome indices"
    doc: "STAR generated stderr log for genome indices"
    
  bowtie_indices:
    type: Directory
    outputSource: bowtie_generate_indices/indices_folder
    label: "Bowtie genome indices"
    doc: "Bowtie generated genome indices folder"

  bowtie_indices_stdout_log:
    type: File
    outputSource: bowtie_generate_indices/stdout_log
    label: "Bowtie stdout log for genome indices"
    doc: "Bowtie generated stdout log for genome indices"
    
  bowtie_indices_stderr_log:
    type: File
    outputSource: bowtie_generate_indices/stderr_log
    label: "Bowtie stderr log genome indices"
    doc: "Bowtie generated stderr log for genome indices"

  ribosomal_indices:
    type: Directory
    outputSource: ribosomal_generate_indices/indices_folder
    label: "Bowtie ribosomal DNA indices"
    doc: "Bowtie generated ribosomal DNA indices folder"

  ribosomal_indices_stdout_log:
    type: File
    outputSource: ribosomal_generate_indices/stdout_log
    label: "Bowtie stdout log for ribosomal DNA indices"
    doc: "Bowtie generated stdout log for ribosomal DNA indices"
    
  ribosomal_indices_stderr_log:
    type: File
    outputSource: ribosomal_generate_indices/stderr_log
    label: "Bowtie stderr log for ribosomal DNA indices"
    doc: "Bowtie generated stderr log for ribosomal DNA indices"

  mitochondrial_indices:
    type: Directory
    outputSource: mitochondrial_generate_indices/indices_folder
    label: "STAR mitochondrial DNA indices"
    doc: "STAR generated mitochondrial DNA indices folder"
    
  mitochondrial_indices_stdout_log:
    type: File
    outputSource: mitochondrial_generate_indices/stdout_log
    label: "STAR stdout log for mitochondrial DNA indices"
    doc: "STAR generated stdout log for mitochondrial DNA indices"
    
  mitochondrial_indices_stderr_log:
    type: File
    outputSource: mitochondrial_generate_indices/stderr_log
    label: "STAR stderr log for mitochondrial DNA indices"
    doc: "STAR generated stderr log for mitochondrial DNA indices"

  annotation_gtf:
    type: File
    format: "http://edamontology.org/format_2306"
    outputSource: prepare_annotation/annotation_gtf_file
    label: "GTF annotation file"
    doc: "GTF annotation file. Includes reference genome and mitochondrial DNA annotations"
    
  annotation:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: prepare_annotation/annotation_tsv_file
    label: "TSV annotation file"
    doc: "Tab-separated annotation file. Includes reference genome and mitochondrial DNA annotations"
  
  fasta_output:
    type: File
    format: "http://edamontology.org/format_1929"
    outputSource: extract_fasta/fasta_file
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes only selected chromosomes"

  fasta_fai_output:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: index_fasta/fai_file
    label: "FAI index for genome FASTA file"
    doc: "Tab-separated FAI index file"

  cytoband_output:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: extract_cytoband/output_file
    label: "CytoBand file for IGV browser"
    doc: "Tab-separated cytoBand file for IGV browser"
  
  annotation_bed:
    type: File
    format: "http://edamontology.org/format_3003"
    outputSource: sort_annotation_bed/sorted_file
    label: "Sorted BED annotation file"
    doc: "Sorted BED annotation file"
  
  annotation_bed_tbi:
    type: File
    format: "http://edamontology.org/format_3004"
    outputSource: annotation_bed_to_bigbed/bigbed_file
    label: "Sorted bigBed annotation file"
    doc: "Sorted bigBed annotation file"

  genome_size:
    type: string
    outputSource: effective_genome_size
    label: "Effective genome size"
    doc: "MACS2 effective genome sizes: hs, mm, ce, dm or number, for example 2.7e9"
    
  chrom_length:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: star_generate_indices/chrom_length
    label: "Genome chromosome length file"
    doc: "Genome chromosome length file"


steps:

  extract_fasta:
    run: ../tools/ucsc-twobit-to-fa.cwl
    in:
      reference_file: genome_file
      chr_list: chromosome_list
    out:
    - fasta_file

  extract_mitochondrial_fasta:
    run: ../tools/ucsc-twobit-to-fa.cwl
    in:
      reference_file: genome_file
      chr_list:
        default: ["chrM"]
    out:
    - fasta_file

  extract_cytoband:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
        - class: ResourceRequirement
          ramMin: 7620
          coresMin: 1
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

  prepare_annotation:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
        - class: ResourceRequirement
          ramMin: 7620
          coresMin: 1
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/ucscuserapps:v358
      inputs:
        script:
          type: string?
          default: |
            #!/bin/bash
            set -- "$0" "$@"
            CHR=()
            RM_PAR=false
            for arg in "$@"; do
              case $arg in
                --tab=*)
                  TAB="${arg#*=}"
                  ;;
                --mt=*)
                  MT="${arg#*=}"
                  ;;
                --gtf=*)
                  GTF="${arg#*=}"
                  ;;
                --rmpar)
                  RM_PAR=true
                  ;;
                *)
                  CHR+=("$arg")
                  ;;
              esac
            done

            if [[ -n "$GTF" ]]; then

              echo "Extracting GTF annotation"
              echo "Removing version suffix from transcript, gene, and exon IDs"
              echo "Example: gene_id ENSG00000223972.5 -> gene_id ENSG00000223972 gene_version 5"
              ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
              zcat "$GTF" \
                | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
                | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
                | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
                > refgene.gtf

              BIOTYPE_PATTERN="(protein_coding|protein_coding_LoF|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
              GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
              TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
              READTHROUGH_PATTERN="tag \"readthrough_transcript\""

              echo "Filtering GTF annotation to include only specific gene and transcript biotypes"
              cat refgene.gtf \
                | awk '$3 == "transcript"' \
                | grep -E "$GENE_PATTERN" \
                | grep -E "$TX_PATTERN" \
                | grep -Ev "$READTHROUGH_PATTERN" \
                | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
                | sort \
                | uniq \
                > gene_allowlist.txt

              grep -E "^#" refgene.gtf > filtered_refgene.gtf

              if [[ "$RM_PAR" = true ]]; then
                echo "Removing PAR locus genes from the chrY"
                grep -Ff gene_allowlist.txt refgene.gtf \
                  | awk -F "\t" '$1 != "chrY" || $1 == "chrY" && $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
                  >> filtered_refgene.gtf
              else
                grep -Ff gene_allowlist.txt refgene.gtf >> filtered_refgene.gtf
              fi

              rm -f gene_allowlist.txt refgene.gtf
              mv filtered_refgene.gtf refgene.gtf

              echo "Renaming gene_type to gene_biotype, transcript_type to transcript_biotype"
              sed -i 's/gene_type/gene_biotype/g' refgene.gtf
              sed -i 's/transcript_type/transcript_biotype/g' refgene.gtf

              if [[ ${#CHR[@]} -gt 0 ]]; then
                grep -E "^#" refgene.gtf > filtered_refgene.gtf
                echo "Subsetting by" ${CHR[*]}
                cat refgene.gtf | grep -E -v "^#" | awk -v filter="${CHR[*]}" 'BEGIN {split(filter, f); for (i in f) d[f[i]]} {if ($1 in d) print $0}' >> filtered_refgene.gtf
                rm -f refgene.gtf
                mv filtered_refgene.gtf refgene.gtf
              fi

              echo "Converting refgene.gtf to refgene.genepred"
              gtfToGenePred -genePredExt -geneNameAsName2 refgene.gtf refgene.genepred

              echo "Converting refgene.genepred to refgene.tsv (refGene table format)"
              echo -e "bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames" > refgene.tsv
              cat refgene.genepred | awk 'BEGIN{srand()} {print int(rand()*100)+1"\t"$0}' >> refgene.tsv
              rm -f refgene.genepred
            elif [[ -n "$TAB" ]] && [[ -n "$MT" ]]; then
              echo "Extracting TAB annotation"
              gunzip $TAB -c | grep -v "exonCount" > temp_refgene.tsv
              gunzip $MT -c | grep -v "exonCount" | awk '{ if ($3=="chrM") print $0 }' >> temp_refgene.tsv
              if [[ ${#CHR[@]} -gt 0 ]]; then
                echo "Subsetting by" ${CHR[*]}
                cat temp_refgene.tsv | awk -v filter="${CHR[*]}" 'BEGIN {split(filter, f); for (i in f) d[f[i]]} {if ($3 in d) print $0}' > filtered_temp_refgene.tsv
                rm -f temp_refgene.tsv
                mv filtered_temp_refgene.tsv temp_refgene.tsv
              fi
              echo "Converting to GTF format"
              cut -f 2- temp_refgene.tsv | genePredToGtf file stdin refgene.gtf
              echo -e "bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames" > refgene.tsv
              cat temp_refgene.tsv >> refgene.tsv
              rm -f temp_refgene.tsv
            else
              echo "Either a pair of --tab and --mt or --gtf parameter should be provided"
              exit 1
            fi
          inputBinding:
            position: 5
        genome_annotation:
          type: File?
          inputBinding:
            position: 6
            prefix: "--tab="
            separate: false
        mitochondrial_annotation:
          type: File?
          inputBinding:
            position: 7
            prefix: "--mt="
            separate: false
        gtf_annotation:
          type: File?
          inputBinding:
            position: 8
            prefix: "--gtf="
            separate: false
        remove_par_y_genes:
          type: boolean?
          inputBinding:
            position: 9
            prefix: "--rmpar"
        chromosome_list:
          type:
            - "null"
            - string
            - string[]
          inputBinding:
            position: 10
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
      gtf_annotation: gtf_annotation
      remove_par_y_genes: remove_par_y_genes
      chromosome_list: chromosome_list
    out:
    - annotation_tsv_file
    - annotation_gtf_file

  star_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: extract_fasta/fasta_file
      sjdb_gtf_file: prepare_annotation/annotation_gtf_file
      genome_sa_sparse_d: genome_sa_sparse_d
      limit_genome_generate_ram: limit_genome_generate_ram
      genome_sa_index_n_bases: genome_sa_index_n_bases
      genome_chr_bin_n_bits: genome_chr_bin_n_bits
      genome_dir:
        source: genome
        valueFrom: $(self + "_star_genome")
    out:
    - indices_folder
    - chrom_length
    - stdout_log
    - stderr_log

  mitochondrial_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: extract_mitochondrial_fasta/fasta_file
      sjdb_gtf_file: prepare_annotation/annotation_gtf_file
      genome_sa_sparse_d: genome_sa_sparse_d
      genome_sa_index_n_bases: genome_sa_index_n_bases_mitochondrial
      genome_dir:
        source: genome
        valueFrom: $(self + "_star_mitochondrial")
    out:
    - indices_folder
    - chrom_length
    - stdout_log
    - stderr_log

  bowtie_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: extract_fasta/fasta_file
      index_base_name:
        source: genome
        valueFrom: $(self + "_bowtie_genome")
    out:
    - indices_folder
    - stdout_log
    - stderr_log

  ribosomal_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: 
        source: fasta_ribosomal
        valueFrom: |
          ${
            if (self == null){
              var default_ribo = {}
              default_ribo["class"] = "File";
              default_ribo["basename"] = "hg19_default_ribo.fa";
              default_ribo["contents"] = ">gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit\n\
          GCTGACACGCTGTCCTCTGGCGACCTGTCGTCGGAGAGGTTGGGCCTCCGGATGCGCGCGGGGCTCTGGC\n\
          CTCACGGTGACCGGCTAGCCGGCCGCGCTCCTGCCTTGAGCCGCCTGCCGCGGCCCGCGGGCCTGCTGTT\n\
          CTCTCGCGCGTCCGAGCGTCCCGACTCCCGGTGCCGGCCCGGGTCCGGGTCTCTGACCCACCCGGGGGCG\n\
          GCGGGGAAGGCGGCGAGGGCCACCGTGCCCCGTGCGCTCTCCGCTGCGGGCGCCCGGGGCGCCGCACAAC\n\
          CCCACCCGCTGGCTCCGTGCCGTGCGTGTCAGGCGTTCTCGTCTCCGCGGGGTTGTCCGCCGCCCCTTCC\n\
          CCGGAGTGGGGGGTGGCCGGAGCCGATCGGCTCGCTGGCCGGCCGGCCTCCGCTCCCGGGGGGCTCTTCG\n\
          ATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGAGGGGCGGACGTTCGTG\n\
          GCGAACGGGACCGTCCTTCTCGCTCCGCCCGCGCGGTCCCCTCGTCTGCTCCTCTCCCCGCCCGCCGGCC\n\
          GGCGTGTGGGAAGGCGTGGGGTGCGGACCCCGGCCCGACCTCGCCGTCCCGCCCGCCGCCTTCGCTTCGC\n\
          GGGTGCGGGCCGGCGGGGTCCTCTGACGCGGCAGACAGCCCTGCCTGTCGCCTCCAGTGGTTGTCGACTT\n\
          GCGGGCGGCCCCCCTCCGCGGCGGTGGGGGTGCCGTCCCGCCGGCCCGTCGTGCTGCCCTCTCGGGGGGG\n\
          GTTTGCGCGAGCGTCGGCTCCGCCTGGGCCCTTGCGGTGCTCCTGGAGCGCTCCGGGTTGTCCCTCAGGT\n\
          GCCCGAGGCCGAACGGTGGTGTGTCGTTCCCGCCCCCGGCGCCCCCTCCTCCGGTCGCCGCCGCGGTGTC\n\
          CGCGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGAGGCGGTTTGAGTGAGACGAGACGAGACGC\n\
          GCCCCTCCCACGCGGGGAAGGGCGCCCGCCTGCTCTCGGTGAGCGCACGTCCCGTGCTCCCCTCTGGCGG\n\
          GTGCGCGCGGGCCGTGTGAGCGATCGCGGTGGGTTCGGGCCGGTGTGACGCGTGCGCCGGCCGGCCGCCG\n\
          AGGGGCTGCCGTTCTGCCTCCGACCGGTCGTGTGTGGGTTGACTTCGGAGGCGCTCTGCCTCGGAAGGAA\n\
          GGAGGTGGGTGGACGGGGGGGCCTGGTGGGGTTGCGCGCACGCGCGCACCGGCCGGGCCCCCGCCCTGAA\n\
          CGCGAACGCTCGAGGTGGCCGCGCGCAGGTGTTTCCTCGTACCGCAGGGCCCCCTCCCTTCCCCAGGCGT\n\
          CCCTCGGCGCCTCTGCGGGCCCGAGGAGGAGCGGCTGGCGGGTGGGGGGAGTGTGACCCACCCTCGGTGA\n\
          GAAAAGCCTTCTCTAGCGATCTGAGAGGCGTGCCTTGGGGGTACCGGATCCCCCGGGCCGCCGCCTCTGT\n\
          CTCTGCCTCCGTTATGGTAGCGCTGCCGTAGCGACCCGCTCGCAGAGGACCCTCCTCCGCTTCCCCCTCG\n\
          ACGGGGTTGGGGGGGAGAAGCGAGGGTTCCGCCGGCCACCGCGGTGGTGGCCGAGTGCGGCTCGTCGCCT\n\
          ACTGTGGCCCGCGCCTCCCCCTTCCGAGTCGGGGGAGGATCCCGCCGGGCCGGGCCCGGCGCTCCCACCC\n\
          AGCGGGTTGGGACGCGGCGGCCGGCGGGCGGTGGGTGTGCGCGCCCGGCGCTCTGTCCGGCGCGTGACCC\n\
          CCTCCGTCCGCGAGTCGGCTCTCCGCCCGCTCCCGTGCCGAGTCGTGACCGGTGCCGACGACCGCGTTTG\n\
          CGTGGCACGGGGTCGGGCCCGCCTGGCCCTGGGAAAGCGTCCCACGGTGGGGGCGCGCCGGTCTCCCGGA\n\
          GCGGGACCGGGTCGGAGGATGGACGAGAATCACGAGCGACGGTGGTGGTGGCGTGTCGGGTTCGTGGCTG\n\
          CGGTCGCTCCGGGGCCCCCGGTGGCGGGGCCCCGGGGCTCGCGAGGCGGTTCTCGGTGGGGGCCGAGGGC\n\
          CGTCCGGCGTCCCAGGCGGGGCGCCGCGGGACCGCCCTCGTGTCTGTGGCGGTGGGATCCCGCGGCCGTG\n\
          TTTTCCTGGTGGCCCGGCCGTGCCTGAGGTTTCTCCCCGAGCCGCCGCCTCTGCGGGCTCCCGGGTGCCC\n\
          TTGCCCTCGCGGTCCCCGGCCCTCGCCCGTCTGTGCCCTCTTCCCCGCCCGCCGCCCGCCGATCCTCTTC\n\
          TTCCCCCCGAGCGGCTCACCGGCTTCACGTCCGTTGGTGGCCCCGCCTGGGACCGAACCCGGCACCGCCT\n\
          CGTGGGGCGCCGCCGCCGGCCACTGATCGGCCCGGCGTCCGCGTCCCCCGGCGCGCGCCTTGGGGACCGG\n\
          GTCGGTGGCGCGCCGCGTGGGGCCCGGTGGGCTTCCCGGAGGGTTCCGGGGGTCGGCCTGCGGCGCGTGC\n\
          GGGGGAGGAGACGGTTCCGGGGGACCGGCCGCGGCTGCGGCGGCGGCGGTGGTGGGGGGAGCCGCGGGGA\n\
          TCGCCGAGGGCCGGTCGGCCGCCCCGGGTGCCCCGCGGTGCCGCCGGCGGCGGTGAGGCCCCGCGCGTGT\n\
          GTCCCGGCTGCGGTCGGCCGCGCTCGAGGGGTCCCCGTGGCGTCCCCTTCCCCGCCGGCCGCCTTTCTCG\n\
          CGCCTTCCCCGTCGCCCCGGCCTCGCCCGTGGTCTCTCGTCTTCTCCCGGCCCGCTCTTCCGAACCGGGT\n\
          CGGCGCGTCCCCCGGGTGCGCCTCGCTTCCCGGGCCTGCCGCGGCCCTTCCCCGAGGCGTCCGTCCCGGG\n\
          CGTCGGCGTCGGGGAGAGCCCGTCCTCCCCGCGTGGCGTCGCCCCGTTCGGCGCGCGCGTGCGCCCGAGC\n\
          GCGGCCCGGTGGTCCCTCCCGGACAGGCGTTCGTGCGACGTGTGGCGTGGGTCGACCTCCGCCTTGCCGG\n\
          TCGCTCGCCCTCTCCCCGGGTCGGGGGGTGGGGCCCGGGCCGGGGCCTCGGCCCCGGTCGCTGCCTCCCG\n\
          TCCCGGGCGGGGGCGGGCGCGCCGGCCGGCCTCGGTCGCCCTCCCTTGGCCGTCGTGTGGCGTGTGCCAC\n\
          CCCTGCGCCGGCGCCCGCCGGCGGGGCTCGGAGCCGGGCTTCGGCCGGGCCCCGGGCCCTCGACCGGACC\n\
          GGCTGCGCGGGCGCTGCGGCCGCACGGCGCGACTGTCCCCGGGCCGGGCACCGCGGTCCGCCTCTCGCTC\n\
          GCCGCCCGGACGTCGGGGCCGCCCCGCGGGGCGGGCGGAGCGCCGTCCCCGCCTCGCCGCCGCCCGCGGG\n\
          CGCCGGCCGCGCGCGCGCGCGCGTGGCCGCCGGTCCCTCCCGGCCGCCGGGCGCGGGTCGGGCCGTCCGC\n\
          CTCCTCGCGGGCGGGCGCGACGAAGAAGCGTCGCGGGTCTGTGGCGCGGGGCCCCCGGTGGTCGTGTCGC\n\
          GTGGGGGGCGGGTGGTTGGGGCGTCCGGTTCGCCGCGCCCCGCCCCGGCCCCACCGGTCCCGGCCGCCGC\n\
          CCCCGCGCCCGCTCGCTCCCTCCCGTCCGCCCGTCCGCGGCCCGTCCGTCCGTCCGTCCGTCGTCCTCCT\n\
          CGCTTGCGGGGCGCCGGGCCCGTCCTCGCGAGGCCCCCCGGCCGGCCGTCCGGCCGCGTCGGGGGCTCGC\n\
          CGCGCTCTACCTTACCTACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCA\n\
          TGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGT\n\
          CGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCC\n\
          CCTTCGCGGGGGGGATGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCTCTCCGGCCCCGGCC\n\
          GGGGGGCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCCCCGTGGCGGC\n\
          GACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGG\n\
          GTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAG\n\
          CAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCG\n\
          AGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGC\n\
          CAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGG\n\
          ATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGCCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCG\n\
          CCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTC\n\
          AAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTG\n\
          GTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAA\n\
          TTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGA\n\
          AAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACCGGCGATGCG\n\
          GCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTAT\n\
          GGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTT\n\
          GACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGAT\n\
          TCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAAC\n\
          GAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACA\n\
          AGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACG\n\
          CGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCA\n\
          TTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGGAATTCCCGAGTAAGTGCGGGTCAT\n\
          AAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAG\n\
          TGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAAC\n\
          TTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAACG\n\
          GAGCCCGGAGGGCGAGGCCCGCGGCGGCGCCGCCGCCGCCGCGCGCTTCCCTCCGCACACCCACCCCCCC\n\
          ACCGCGACGCGGCGCGTGCGCGGGCGGGGCCCGCGTGCCCGTTCGTTCGCTCGCTCGTTCGTTCGCCGCC\n\
          CGGCCCCGCCGCCGCGAGAGCCGAGAACTCGGGAGGGAGACGGGGGGGAGAGAGAGAGAGAGAGAGAGAG\n\
          AGAGAGAGAGAGAGAGAGAAAGAAGGGCGTGTCGTTGGTGTGCGCGTGTCGTGGGGCCGGCGGGCGGCGG\n\
          GGAGCGGTCCCCGGCCGCGGCCCCGACGACGTGGGTGTCGGCGGGCGCGGGGGCGGTTCTCGGCGGCGTC\n\
          GCGGCGGGTCTGGGGGGGTCTCGGTGCCCTCCTCCCCGCCGGGGCCCGTCGTCCGGCCCCGCCGCGCCGG\n\
          CTCCCCGTCTTCGGGGCCGGCCGGATTCCCGTCGCCTCCGCCGCGCCGCTCCGCGCCGCCGGGCACGGCC\n\
          CCGCTCGCTCTCCCCGGCCTTCCCGCTAGGGCGTCTCGAGGGTCGGGGGCCGGACGCCGGTCCCCTCCCC\n\
          CGCCTCCTCGTCCGCCCCCCCGCCGTCCAGGTACCTAGCGCGTTCCGGCGCGGAGGTTTAAAGACCCCTT\n\
          GGGGGGATCGCCCGTCCGCCCGTGGGTCGGGGGCGGTGGTGGGCCCGCGGGGGAGTCCCGTCGGGAGGGG\n\
          CCCGGCCCCTCCCGCGCCTCCACCGCGGACTCCGCTCCCCGGCCGGGGCCGCGCCGCCGCCGCCGCCGCG\n\
          GCGGCCGTCGGGTGGGGGCTTTACCCGGCGGCCGTCGCGCGCCTGCCGCGCGTGTGGCGTGCGCCCCGCG\n\
          CCGTGGGGGCGGGAACCCCCGGGCGCCTGTGGGGTGGTGTCCGCGCTCGCCCCCGCGTGGGCGGCGCGCG\n\
          CCTCCCCGTGGTGTGAAACCTTCCGACCCCTCTCCGGAGTCCGGTCCCGTTTGCTGTCTCGTCTGGCCGG\n\
          CCTGAGGCAACCCCCTCTCCTCTTGGGCGGGGGGGGCGGGGGGACGTGCCGCGCCAGGAAGGGCCTCCTC\n\
          CCGGTGCGTCGTCGGGAGCGCCCTCGCCAAATCGACCTCGTACGACTCTTAGCGGTGGATCACTCGGCTC\n\
          GTGCGTCGATGAAGAACGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATCGACAC\n\
          TTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGCTTGCCGATCAATC\n\
          GCCCCGGGGGTGCCTCCGGGCTCCTCGGGGTGCGCGGCTGGGGGTTCCCTCGCAGGGCCCGCCGGGGGCC\n\
          CTCCGTCCCCCTAAGCGCAGACCCGGCGGCGTCCGCCCTCCTCTTGCCGCCGCGCCCGCCCCTTCCCCCT\n\
          CCCCCCGCGGGCCCTGCGTGGTCACGCGTCGGGTGGCGGGGGGGAGAGGGGGGCGCGCCCGGCTGAGAGA\n\
          GACGGGGAGGGCGGCGCCGCCGCCGGAAGACGGAGAGGGAAAGAGAGAGCCGGCTCGGGCCGAGTTCCCG\n\
          TGGCCGCCGCCTGCGGTCCGGGTTCCTCCCTCGGGGGGCTCCCTCGCGCCGCGCGCGGCTCGGGGTTCGG\n\
          GGTTCGTCGGCCCCGGCCGGGTGGAAGGTCCCGTGCCCGTCGTCGTCGTCGTCGCGCGTCGTCGGCGGTG\n\
          GGGGCGTGTTGCGTGCGGTGTGGTGGTGGGGGAGGAGGAAGGCGGGTCCGGAAGGGGAAGGGTGCCGGCG\n\
          GGGAGAGAGGGTCGGGGGAGCGCGTCCCGGTCGCCGCGGTTCCGCCGCCCGCCCCCGGTGGCGGCCCGGC\n\
          GTCCGGCCGACCGGCCGCTCCCCGCGCCCCTCCTCCTCCCCGCCGCCCCTCCTCCGAGGCCCCGCCCGTC\n\
          CTCCTCGCCCTCCCCGCGCGTACGCGCGCGCGCCCGCCCGCCCGGCTCGCCTCGCGGCGCGTCGGCCGGG\n\
          GCCGGGAGCCCGCCCCGCCGCCCGCCCGTGGCCGCGGCGCCGGGGTTCGCGTGTCCCCGGCGGCGACCCG\n\
          CGGGACGCCGCGGTGTCGTCCGCCGTCGCGCGCCCGCCTCCGGCTCGCGGCCGCGCCGCGCCGCGCCGGG\n\
          GCCCCGTCCCGAGCTTCCGCGTCGGGGCGGCGCGGCTCCGCCGCCGCGTCCTCGGACCCGTCCCCCCGAC\n\
          CTCCGCGGGGGAGACGCGCCGGGGCGTGCGGCGCCCGTCCCGCCCCCGGCCCGTGCCCCTCCCTCCGGTC\n\
          GTCCCGCTCCGGCGGGGCGGCGCGGGGGCGCCGTCGGCCGCGCGCTCTCTCTCCCGTCGCCTCTCCCCCT\n\
          CGCCGGGCCCGTCTCCCGACGGAGCGTCGGGCGGGCGGTCGGGCCGGCGCGATTCCGTCCGTCCGTCCGC\n\
          CGAGCGGCCCGTCCCCCTCCGAGACGCGACCTCAGATCAGACGTGGCGACCCGCTGAATTTAAGCATATT\n\
          AGTCAGCGGAGGAAAAGAAACTAACCAGGATTCCCTCAGTAACGGCGAGTGAACAGGGAAGAGCCCAGCG\n\
          CCGAATCCCCGCCCCGCGGGGCGCGGGACATGTGGCGTACGGAAGACCCGCTCCCCGGCGCCGCTCGTGG\n\
          GGGGCCCAAGTCCTTCTGATCGAGGCCCAGCCCGTGGACGGTGTGAGGCCGGTAGCGGCCGGCGCGCGCC\n\
          CGGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAA\n\
          ATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTT\n\
          CAAGAGGGCGTGAAACCGTTAAGAGGTAAACGGGTGGGGTCCGCGCAGTCCGCCCGGAGGATTCAACCCG\n\
          GCGGCGGGTCCGGCCGTGTCGGCGGCCCGGCGGATCTTTCCCGCCCCCCGTTCCTCCCGACCCCTCCACC\n\
          CGCCCTCCCTTCCCCCGCCGCCCCTCCTCCTCCTCCCCGGAGGGGGCGGGCTCCGGCGGGTGCGGGGGTG\n\
          GGCGGGCGGGGCCGGGGGTGGGGTCGGCGGGGGACCGTCCCCCGACCGGCGACCGGCCGCCGCCGGGCGC\n\
          ATTTCCACCGCGGCGGTGCGCCGCGACCGGCTCCGGGACGGCTGGGAAGGCCCGGCGGGGAAGGTGGCTC\n\
          GGGGGGCCCCGTCCGTCCGTCCGTCCTCCTCCTCCCCCGTCTCCGCCCCCCGGCCCCGCGTCCTCCCTCG\n\
          GGAGGGCGCGCGGGTCGGGGCGGCGGCGGCGGCGGCGGTGGCGGCGGCGGCGGGGGCGGCGGGACCGAAA\n\
          CCCCCCCCGAGTGTTACAGCCCCCCCGGCAGCAGCACTCGCCGAATCCCGGGGCCGAGGGAGCGAGACCC\n\
          GTCGCCGCGCTCTCCCCCCTCCCGGCGCCCACCCCCGCGGGGAATCCCCCGCGAGGGGGGTCTCCCCCGC\n\
          GGGGGCGCGCCGGCGTCTCCTCGTGGGGGGGCCGGGCCACCCCTCCCACGGCGCGACCGCTCTCCCACCC\n\
          CTCCTCCCCGCGCCCCCGCCCCGGCGACGGGGGGGGTGCCGCGCGCGGGTCGGGGGGCGGGGCGGACTGT\n\
          CCCCAGTGCGCCCCGGGCGGGTCGCGCCGTCGGGCCCGGGGGAGGTTCTCTCGGGGCCACGCGCGCGTCC\n\
          CCCGAAGAGGGGGACGGCGGAGCGAGCGCACGGGGTCGGCGGCGACGTCGGCTACCCACCCGACCCGTCT\n\
          TGAAACACGGACCAAGGAGTCTAACACGTGCGCGAGTCGGGGGCTCGCACGAAAGCCGCCGTGGCGCAAT\n\
          GAAGGTGAAGGCCGGCGCGCTCGCCGGCCGAGGTGGGATCCCGAGGCCTCTCCAGTCCGCCGAGGGCGCA\n\
          CCACCGGCCCGTCTCGCCCGCCGCGCCGGGGAGGTGGAGCACGAGCGCACGTGTTAGGACCCGAAAGATG\n\
          GTGAACTATGCCTGGGCAGGGCGAAGCCAGAGGAAACTCTGGTGGAGGTCCGTAGCGGTCCTGACGTGCA\n\
          AATCGGTCGTCCGACCTGGGTATAGGGGCGAAAGACTAATCGAACCATCTAGTAGCTGGTTCCCTCCGAA\n\
          GTTTCCCTCAGGATAGCTGGCGCTCTCGCAGACCCGACGCACCCCCGCCACGCAGTTTTATCCGGTAAAG\n\
          CGAATGATTAGAGGTCTTGGGGCCGAAACGATCTCAACCTATTCTCAAACTTTAAATGGGTAAGAAGCCC\n\
          GGCTCGCTGGCGTGGAGCCGGGCGTGGAATGCGAGTGCCTAGTGGGCCACTTTTGGTAAGCAGAACTGGC\n\
          GCTGCGGGATGAACCGAACGCCGGGTTAAGGCGCCCGATGCCGACGCTCATCAGACCCCAGAAAAGGTGT\n\
          TGGTTGATATAGACAGCAGGACGGTGGCCATGGAAGTCGGAATCCGCTAAGGAGTGTGTAACAACTCACC\n\
          TGCCGAATCAACTAGCCCTGAAAATGGATGGCGCTGGAGCGTCGGGCCCATACCCGGCCGTCGCCGGCAG\n\
          TCGAGAGTGGACGGGAGCGGCGGGGGCGGCGCGCGCGCGCGCGCGTGTGGTGTGCGTCGGAGGGCGGCGG\n\
          CGGCGGCGGCGGCGGGGGTGTGGGGTCCTTCCCCCGCCCCCCCCCCCACGCCTCCTCCCCTCCTCCCGCC\n\
          CACGCCCCGCTCCCCGCCCCCGGAGCCCCGCGGACGCTACGCCGCGACGAGTAGGAGGGCCGCTGCGGTG\n\
          AGCCTTGAAGCCTAGGGCGCGGGCCCGGGTGGAGCCGCCGCAGGTGCAGATCTTGGTGGTAGTAGCAAAT\n\
          ATTCAAACGAGAACTTTGAAGGCCGAAGTGGAGAAGGGTTCCATGTGAACAGCAGTTGAACATGGGTCAG\n\
          TCGGTCCTGAGAGATGGGCGAGCGCCGTTCCGAAGGGACGGGCGATGGCCTCCGTTGCCCTCGGCCGATC\n\
          GAAAGGGAGTCGGGTTCAGATCCCCGAATCCGGAGTGGCGGAGATGGGCGCCGCGAGGCGTCCAGTGCGG\n\
          TAACGCGACCGATCCCGGAGAAGCCGGCGGGAGCCCCGGGGAGAGTTCTCTTTTCTTTGTGAAGGGCAGG\n\
          GCGCCCTGGAATGGGTTCGCCCCGAGAGAGGGGCCCGTGCCTTGGAAAGCGTCGCGGTTCCGGCGGCGTC\n\
          CGGTGAGCTCTCGCTGGCCCTTGAAAATCCGGGGGAGAGGGTGTAAATCTCGCGCCGGGCCGTACCCATA\n\
          TCCGCAGCAGGTCTCCAAGGTGAACAGCCTCTGGCATGTTGGAACAATGTAGGTAAGGGAAGTCGGCAAG\n\
          CCGGATCCGTAACTTCGGGATAAGGATTGGCTCTAAGGGCTGGGTCGGTCGGGCTGGGGCGCGAAGCGGG\n\
          GCTGGGCGCGCGCCGCGGCTGGACGAGGCGCGCGCCCCCCCCACGCCCGGGGCACCCCCCTCGCGGCCCT\n\
          CCCCCGCCCCACCCGCGCGCGCCGCTCGCTCCCTCCCCACCCCGCGCCCTCTCTCTCTCTCTCTCCCCCG\n\
          CTCCCCGTCCTCCCCCCTCCCCGGGGGAGCGCCGCGTGGGGGCGCGGCGGGGGGAGAAGGGTCGGGGCGG\n\
          CAGGGGCCGCGCGGCGGCCGCCGGGGCGGCCGGCGGGGGCAGGTCCCCGCGAGGGGGGCCCCGGGGACCC\n\
          GGGGGGCCGGCGGCGGCGCGGACTCTGGACGCGAGCCGGGCCCTTCCCGTGGATCGCCCCAGCTGCGGCG\n\
          GGCGTCGCGGCCGCCCCCGGGGAGCCCGGCGGCGGCGCGGCGCGCCCCCCACCCCCACCCCACGTCTCGG\n\
          TCGCGCGCGCGTCCGCTGGGGGCGGGAGCGGTCGGGCGGCGGCGGTCGGCGGGCGGCGGGGCGGGGCGGT\n\
          TCGTCCCCCCGCCCTACCCCCCCGGCCCCGTCCGCCCCCCGTTCCCCCCTCCTCCTCGGCGCGCGGCGGC\n\
          GGCGGCGGCAGGCGGCGGAGGGGCCGCGGGCCGGTCCCCCCCGCCGGGTCCGCCCCCGGGGCCGCGGTTC\n\
          CGCGCGCGCCTCGCCTCGGCCGGCGCCTAGCAGCCGACTTAGAACTGGTGCGGACCAGGGGAATCCGACT\n\
          GTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCGGGTGTTGACGCGATGTGATTTCTGCCCAGTGCT\n\
          CTGAATGTCAAAGTGAAGAAATTCAATGAAGCGCGGGTAAACGGCGGGAGTAACTATGACTCTCTTAAGG\n\
          TAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAATGGATGAACGAGATTCCCACTGTCCCTACC\n\
          TACTATCCAGCGAAACCACAGCCAAGGGAACGGGCTTGGCGGAATCAGCGGGGAAAGAAGACCCTGTTGA\n\
          GCTTGACTCTAGTCTGGCACGGTGAAGAGACATGAGAGGTGTAGAATAAGTGGGAGGCCCCCGGCGCCCC\n\
          CCCGGTGTCCCCGCGAGGGGCCCGGGGCGGGGTCCGCGGCCCTGCGGGCCGCCGGTGAAATACCACTACT\n\
          CTGATCGTTTTTTCACTGACCCGGTGAGGCGGGGGGGCGAGCCCGAGGGGCTCTCGCTTCTGGCGCCAAG\n\
          CGCCCGCCCGGCCGGGCGCGACCCGCTCCGGGGACAGTGCCAGGTGGGGAGTTTGACTGGGGCGGTACAC\n\
          CTGTCAAACGGTAACGCAGGTGTCCTAAGGCGAGCTCAGGGAGGACAGAAACCTCCCGTGGAGCAGAAGG\n\
          GCAAAAGCTCGCTTGATCTTGATTTTCAGTACGAATACAGACCGTGAAAGCGGGGCCTCACGATCCTTCT\n\
          GACCTTTTGGGTTTTAAGCAGGAGGTGTCAGAAAAGTTACCACAGGGATAACTGGCTTGTGGCGGCCAAG\n\
          CGTTCATAGCGACGTCGCTTTTTGATCCTTCGATGTCGGCTCTTCCTATCATTGTGAAGCAGAATTCGCC\n\
          AAGCGTTGGATTGTTCACCCACTAATAGGGAACGTGAGCTGGGTTTAGACCGTCGTGAGACAGGTTAGTT\n\
          TTACCCTACTGATGATGTGTTGTTGCCATGGTAATCCTGCTCAGTACGAGAGGAACCGCAGGTTCAGACA\n\
          TTTGGTGTATGTGCTTGGCTGAGGAGCCAATGGGGCGAAGCTACCATCTGTGGGATTATGACTGAACGCC\n\
          TCTAAGTCAGAATCCCGCCCAGGCGAACGATACGGCAGCGCCGCGGAGCCTCGGTTGGCCTCGGATAGCC\n\
          GGTCCCCCGCCTGTCCCCGCCGGCGGGCCGCCCCCCCCTCCACGCGCCCCGCCGCGGGAGGGCGCGTGCC\n\
          CCGCCGCGCGCCGGGACCGGGGTCCGGTGCGGAGTGCCCTTCGTCCTGGGAAACGGGGCGCGGCCGGAAA\n\
          GGCGGCCGCCCCCTCGCCCGTCACGCACCGCACGTTCGTGGGGAACCTGGCGCTAAACCATTCGTAGACG\n\
          ACCTGCTTCTGGGTCGGGGTTTCGTACGTAGCAGAGCAGCTCCCTCGCTGCGATCTATTGAAAGTCAGCC\n\
          CTCGACACAAGGGTTTGTCCGCGCGCGCGTGCGTGCGGGGGGCCCGGCGGGCGTGCGCGTTCGGCGCCGT\n\
          CCGTCCTTCCGTTCGTCTTCCTCCCTCCCGGCCTCTCCCGCCGACCGCGGCGTGGTGGTGGGGTGGGGGG\n\
          GAGGGCGCGCGACCCCGGTCGGCCGCCCCGCTTCTTCGGTTCCCGCCTCCTCCCCGTTCACGCCGGGGCG\n\
          GCTCGTCCGCTCCGGGCCGGGACGGGGTCCGGGGAGCGTGGTTTGGGAGCCGCGGAGGCGCCGCGCCGAG\n\
          CCGGGCCCCGTGGCCCGCCGGTCCCCGTCCCGGGGGTTGGCCGCGCGGCGCGGTGGGGGGCCACCCGGGG\n\
          TCCCGGCCCTCGCGCGTCCTTCCTCCTCGCTCCTCCGCACGGGTCGACCGACGAACCGCGGGTGGCGGGC\n\
          GGCGGGCGGCGAGCCCCACGGGCGTCCCCGCACCCGGCCGACCTCCGCTCGCGACCTCTCCTCGGTCGGG\n\
          CCTCCGGGGTCGACCGCCTGCGCCCGCGGGCGTGAGACTCAGCGGCGTCTCGCCGTGTCCCGGGTCGACC\n\
          GCGGCCTTCTCCACCGAGCGGCGGTGTAGGAGTGCCCGTCGGGACGAACCGCAACCGGAGCGTCCCCGTC\n\
          TCGGTCGGCACCTCCGGGGTCGACCAGCTGCCGCCCGCGAGCTCCGGACTTAGCCGGCGTCTGCACGTGT\n\
          CCCGGGTCGACCAGCAGGCGGCCGCCGGACGCAGCGGCGCACGCACGCGAGGGCGTCGATTCCCCTTCGC\n\
          GCGCCCGCGCCTCCACCGGCCTCGGCCCGCGGTGGAGCTGGGACCACGCGGAACTCCCTCTCCCACATTT\n\
          TTTTCAGCCCCACCGCGAGTTTGCGTCCGCGGGACCTTTAAGAGGGAGTCACTGCTGCCGTCAGCCAGTA\n\
          CTGCCTCCTCCTTTTTCGCTTTTAGGTTTTGCTTGCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTT\n\
          CTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCGCTTGTCTTCTTCTTGTGTTCTCTTCTTG\n\
          CTCTTCCTCTGTCTGTCTCTCTCTCTCTCTCTCTCTCTGTCTCTCGCTCTCGCCCTCTCTCTCTTCTCTC\n\
          TCTCTCTCTCTCTCTCTCTGTCTCTCGCTCTCGCCCTCTCTCTCTCTCTTCTCTCTGTCTCTCTCTCTCT\n\
          CTCTCTCTCTCTCTCTCTCTGTCGCTCTCGCCCTCTCGCTCTCTCTCTGTCTCTGTCTGTGTCTCTCTCT\n\
          CTCCCTCCCTCCCTCCCTCCCTCCCTCCCTCCCTCCCCTTCCTTGGCGCCTTCTCGGCTCTTGAGACTTA\n\
          GCCGCTGTCTCGCCGTACCCCGGGTCGACCGGCGGGCCTTCTCCACCGAGCGGCGTGCCACAGTGCCCGT\n\
          CGGGACGAGCCGGACCCGCCGCGTCCCCGTCTCGGTCGGCACCTCCGGGGTCGACCAGCTGCCGCCCGCG\n\
          AGCTCCGGACTTAGCCGGCGTCTGCACGTGTCCCGGGTCGACCAGCAGGCGGCCGCCGGACGCAGCGGCG\n\
          CACCGACGGAGGGCGCTGATTCCCGTTCACGCGCCCGCGCCTCCACCGGCCTCGGCCCGCCGTGGAGCTG\n\
          GGACCACGCGGAACTCCCTCTCCTACATTTTTTTCAGCCCCACCGCGAGTTTGCGTCCGCGGGACCTTTA\n\
          AGAGGGAGTCACTGCTGCCGTCAGCCAGTACTGCCTCCTCCTTTTTCGCTTTTAGGTTTTGCTTGCCTTT\n\
          TTTTTTTTTTTTTTTTTTTTTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT\n\
          CTTTCGCTCTCGCTCTCTCGCTCTCTCCCTCGCTCGTTTCTTTCTTTCTCTTTCTCTCTCTCTCTCTCTC\n\
          TCTCTCTCTCTCTGTCTCTCGCTCTCGCCCTCTCTCTCTCTTTCTCTCTCTCTCTGTCTCTCTCTCTCTC\n\
          TCTCTCTCTCTCTCTCTCTCCCTCCCTCCCTCCCCCTCCCTCCCTCTCTCCCCTTCCTTGGCGCCTTCTC\n\
          GGCTCTTGAGACTTAGCCGCTGTCTCGCCGTGTCCCGGGTCGACCGGCGGGCCTTCTCCACCGAGCGGCG\n\
          TGCCACAGTGCCCGTCGGGACGAGCCGGACCCGCCGCGTCCCCGTCTCGGTCGGCACCTCCGGGGTCGAC\n\
          CAGCTGCCGCCCGCGAGCTCCGGACTTAGCCGGCGTCTGCACGTGTCCCGGGTCGACCAGCAGGCGGCCG\n\
          CCGGACGCTGCGGCGCACCGACGCGAGGGCGTCGATTCCGGTTCACGCGCCGGCGACCTCCACCGGCCTC\n\
          GGCCCGCGGTGGAGCTGGGACCACGCGGAACTCCCTCTCCCACATTTTTTTCAGCCCCACCGCGAGTTTG\n\
          CGTCCGCGGGACTTTTAAGAGGGAGTCACTGCTGCCGTCAGCCAGTAATGCTTCCTCCTTTTTTGCTTTT\n\
          TGGTTTTGCCTTGCGTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTCTCTCTCTC\n\
          TCTCTCTCTCTCTCTGTCTCTCTCTCTCTGTCTCTCTCCCCTCCCTCCCTCCTTGGTGCCTTCTCGGCTC\n\
          GCTGCTGCTGCTGCCTCTGCCTCCACGGTTCAAGCAAACAGCAAGTTTTCTATTTCGAGTAAAGACGTAA\n\
          TTTCACCATTTTGGCCGGGCTGGTCTCGAACTCCCGACCTAGTGATCCGCCCGCCTCGGCCTCCCAAAGA\n\
          CTGCTGGGAGTACAGATGTGAGCCACCATGCCCGGCCGATTCCTTCCTTTTTTCAATCTTATTTTCTGAA\n\
          CGCTGCCGTGTATGAACATACATCTACACACACACACACACACACACACACACACACACACACACACACA\n\
          CACACACCCCGTAGTGATAAAACTATGTAAATGATATTTCCATAATTAATACGTTTATATTATGTTACTT\n\
          TTAATGGATGAATATGTATCGAAGCCCCATTTCATTTACATACACGTGTATGTATATCCTTCCTCCCTTC\n\
          CTTCATTCATTATTTATTAATAATTTTCGTTTATTTATTTTCTTTTCTTTTGGGGCCGGCCCGCCTGGTC\n\
          TTCTGTCTCTGCGCTCTGGTGACCTCAGCCTCCCAAATAGCTGGGACTACAGGGATCTCTTAAGCCCGGG\n\
          AGGAGAGGTTAACGTGGGCTGTGATCGCACACTTCCACTCCAGCTTACGTGGGCTGCGGTGCGGTGGGGT\n\
          GGGGTGGGGTGGGGTGGGGTGCAGAGAAAACGATTGATTGCGATCTCAATTGCCTTTTAGCTTCATTCAT\n\
          ACCCTGTTATTTGCTCGTTTATTCTCATGGGTTCTTCTGTGTCATTGTCACGTTCATCGTTTGCTTGCCT\n\
          GCTTGCCTGTTTATTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCCTTA\n\
          CTGGCAGGGTCTTCCTCTGTCTCTGCCGCCCAGGATCACCCCAACCTCAACGCTTTGGACCGACCAAACG\n\
          GTCGTTCTGCCTCTGATCCCTCCCATCCCCATTACCTGAGACTACAGGCGCGCACCACCACACCGGCTGA\n\
          CTTTTATGTTGTTTCTCATGTTTTCCGTAGGTAGGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT\n\
          GTGTGTGTGTGTGTGTGTGTGTGTGTATCTATGTATGTACGTATGTATGTATGTATGTGAGTGAGATGGG\n\
          TTTCGGGGTTCTATCATGTTGCCCACGCTGGTCTCGAACTCCTGTCCTCAAGCAATCCGCCTGCCTGCCT\n\
          CGGCCGCCCACACTGCTGCTATTACAGGCGTGAGACGCTGCGCCTGGCTCCTTCTACATTTGCCTGCCTG\n\
          CCTGCCTGCCTGCCTGCCTATCAATCGTCTTCTTTTTAGTACGGATGTCGTCTCGCTTTATTGTCCATGC\n\
          TCTGGGCACACGTGGTCTCTTTTCAAACTTCTATGATTATTATTATTGTAGGCGTCATCTCACGTGTCGA\n\
          GGTGATCTCGAACTTTTAGGCTCCAGAGATCCTCCCGCATCGGCCTCCCGGAGTGCTGTGATGACACGCG\n\
          TGGGCACGGTACGCTCTGGTCGTGTTTGTCGTGGGTCGGTTCTTTCCGTTTTTAATACGGGGACTGCGAA\n\
          CGAAGAAAATTTTCAGACGCATCTCACCGATCCGCCTTTTCGTTCTTTCTTTTTATTCTCTTTAGACGGA\n\
          GTTTCACTCTTGTCGCCCAGGGTGGAGTACGATGGCGGCTCTCGGCTCACCGCACCCTCCGCCTCCCAGG\n\
          TTCAAGTGATTCTCCTGCCTCAGCCTTCCCGAGTAGCTGGAATGACAGAGATGAGCCATCGTGCCCGGCT\n\
          AATTTTTCTATTTTTAGTACAGATGGGGTTTCTCCATCTTGGTCAGGCTGGTCTTCAACTTCCGACCGTT\n\
          GGAGAATCTTAACTTTCTTGGTGGTGGTTGTTTTCCTTTTTCTTTTTTTTTCTTTTCTTTTCTTTCCTTC\n\
          TCCTCCCCCCCCCACCCCCCTTGTCGTCGTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTC\n\
          CTCCTCCTCCTCTTTCATTTCTTTCAGCTGGGCTCTCCTACTTGTGTTGCTCTGTTGCTCACGCTGGTCT\n\
          CAAACTCCTGGCCTTGACTCTTCTCCCGTCACATCCGCCGTCTGGTTGTTGAAATGAGCATCTCTCGTAA\n\
          AATGGAAAAGATGAAAGAAATAAACACGAAGACGGAAAGCACGGTGTGAACGTTTCTCTTGCCGTCTCCC\n\
          GGGGTGTACCTTGGACCCGGAAACACGGAGGGAGCTTGGCTGAGTGGGTTTTCGGTGCCGAAACCTCCCG\n\
          AGGGCCTCCTTCCCTCTCCCCCTTGTCCCCGCTTCTCCGCCAGCCGAGGCTCCCACCGCCGCCCCTGGCA\n\
          TTTTCCATAGGAGAGGTATGGGAGAGGACTGACACGCCTTCCAGATCTATATCCTGCCGGACGTCTCTGG\n\
          CTCGGCGTGCCCCACCGGCTACCTGCCACCTTCCAGGGAGCTCTGAGGCGGATGCGACCCCCACCCCCCC\n\
          GTCACGTCCCGCTACCCTCCCCCGGCTGGCCTTTGCCGGGCGACCCCAGGGGAACCGCGTTGATGCTGCT\n\
          TCGGATCCTCCGGCGAAGACTTCCACCGGATGCCCCGGGTGGGCCGGTTGGGATCAGACTGGACCACCCC\n\
          GGACCGTGCTGTTCTTGGGGGTGGGTTGACGTACAGGGTGGACTGGCAGCCCCAGCATTGTAAAGGGTGC\n\
          GTGGGTATGGAAATGTCACCTAGGATGCCCTCCTTCCCTTCGGTCTGCCTTCAGCTGCCTCAGGCGTGAA\n\
          GACAACTTCCCATCGGAACCTCTTCTCTTCCCTTTCTCCAGCACACAGATGAGACGCACGAGAGGGAGAA\n\
          ACAGCTCAATAGATACCGCTGACCTTCATTTGTGGAATCCTCAGTCATCGACACACAAGACAGGTGACTA\n\
          GGCAGGGACACAGATCAAACACTATTTCCGGGTCCTCGTGGTGGGATTGGTCTCTCTCTCTCTCTCTCTC\n\
          TCTCTCTCTCTCTCTCTCTCTCTCGCACGCGCACGCGCGCACACACACACACAATTTCCATATCTAGTTC\n\
          ACAGAGCACACTCACTTCCCCTTTTCACAGTACGCAGGCTGAGTAAAACGCGCCCCACCCTCCACCCGTT\n\
          GGCTGACGAAACCCCTTCTCTACAATTGATGAAAAAGATGATCTGGGCCGGGCACGCTAGCTCACGCCTG\n\
          TCACTCCGGCACTTTGGGAGGCCGAGGCGGGTGGATCGCTTGGGGCCGGGAGTTCGAGACCAGGCTGGCC\n\
          GACGTGGCGAAACCCCGTCTCTCTGAAAAATAGAACGATTAGCCGGGCCTGGTGGCGTGGGCTTGGAATC\n\
          ACGACCGCTCGGGAGACTGGGGCGGGCGACTTGTTCCAACCGGGGAGGCCGAGGCCGCGATGAGCTGAGA\n\
          TCGTGCCGTGGCGATGCGGCCTGGATGACGGAGCGAGACCCCGTCTCGAGAGAATCATGATGTTATTATA\n\
          AGATGAGTTGTGCGCGGTGATGGCCGCCTGTAGTCGCGGCTACTCGGGAGGCTGAGACGAGGAGAAGATC\n\
          ACTTGAGGCCCCACAGGTCGAGGCTTCGGTCGGCCGTGACCCACTGTATCCTGGGCAGTCACCGGTCAAG\n\
          GAGATATGCCCCTTCCCCGTTTGCTTTTCTTTTCTTCCCTTCTCTTTTCTTCTTTTTGCTTCTCTTTTCT\n\
          TTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTCTTTTTCTCTCTTCCCCTCTTTCTTT\n\
          CCTGCCTTCCTGCCTTTCTTCTTTTCTTCTTTCCTCCCTTCCTCCCTTCCTTCTTTCCTCCCGCCTCAGC\n\
          CTCCCAAAGTGCTGGGATGACTGGCGGGAGGCACCATGCCTGCTTGGCCCAAAGAGACCCTCTTGGAAAG\n\
          TGAGACGCAGAGAGCGCCTTCCAGTGATCTCATTGACTGATTTAGAGACGGCATCTCGCTCCGTCACCCC\n\
          GGCAGTGGTGCCGTCGTAACTCACTCCCTGCAGCGTGGACGCTCCTGGACTCGAGCGATCCTTCCACCTC\n\
          AGCCTCCAGAGTACAGAGCCTGGGACCGCGGGCACGCGCCACTGTGCCCACACCGTTTTTAATTGTTTTT\n\
          TTTTCCCCCGAGACAGAGTTTCACTCTCGTGGCCTAGACTGCAGTGCGGTGGCGCGATCTTGGCTCACCG\n\
          CAACCTCTGCCTCCCGGTTTCAAGCGATTCTCCTGCATCGGCCTCCTGAGTAGCCGGGATTGCGGGCATG\n\
          CGCTGCCACGTCTGGCTGATTTCGTATTTTTAGTGGAGACGGGGCTTCTCCATGTCGATCGGGCTGGTTT\n\
          CGAACTCCCGACCTCAGGTGATCCGCCCTCCCCGGCCTCCGGAAGTGCTGGGATGACAGGCGTGAGCCAC\n\
          CGCGCCCGGCCTTCATTTTTAAATGTTTTCCCACAGACGGGGTCTCATCATTTCTTTGCAACCCTCCTGC\n\
          CCGGCGTCTCAAAGTGCTGGCGTGACGGGCGTGAGCCACTGCGCCTGGACTCCGGGGAATGACTCACGAC\n\
          CACCATCGCTCTACTGATCCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTT\n\
          TCTTTCTTGATGAATTATCTTATGATTTATTTGTGTACTTATTTTCAGACGGAGTCTCGCTCTGGGCGGG\n\
          GCGAGGCGAGGCGAGGCACAGCGCATCGCTTTGGAAGCCGCGGCAACGCCTTTCAAAGCCCCATTCGTAT\n\
          GCACAGAGCCTTATTCCCTTCCTGGAGTTGGAGCTGATGCCTTCCGTAGCCTTGGGCTTCTCTCCATTCG\n\
          GAAGCTTGACAGGCGCAGGGCCACCCAGAGGCTGGCTGCGGCTGAGGATTAGGGGGTGTGTTGGGGCTGA\n\
          AAACTGGGTCCCCTATTTTTGATACCTCAGCCGACACATCCCCCGACCGCCATCGCTTGCTCGCCCTCTG\n\
          AGATCCCCCGCCTCCACCGCCTTGCAGGCTCACCTCTTACTTTCATTTCTTCCTTTCTTGCGTTTGAGGA\n\
          GGGGGTGCGGGAATGAGGGTGTGTGTGGGGAGGGGGTGCGGGGTGGGGACGGAGGGGAGCGTCCTAAGGG\n\
          TCGATTTAGTGTCATGCCTCTTTCACCACCACCACCACCACCGAAGATGACAGCAAGGATCGGCTAAATA\n\
          CCGCGTGTTCTCATCTAGAAGTGGGAACTTACAGATGACAGTTCTTGCATGGGCAGAACGAGGGGGACCG\n\
          GGGACGCGGAAGTCTGCTTGAGGGAGGAGGGGTGGAAGGAGAGACAGCTTCAGGAAGAAAACAAAACACG\n\
          AATACTGTCGGACACAGCACTGACTACCCGGGTGATGAAATCATCTGCACACTGAACACCCCCGTCACAA\n\
          GTTTACCTATGTCACAATCTTGCACATGTATCGCTTGAACGACAAATAAAAGTTAGGGGGGAGAAGAGAG\n\
          GAGAGAGAGAGAGAGAGAGAGACAGAGAGAGACAGAGAGAGAGAGAGAGGAGGGAGAGAGGAAAACGAAA\n\
          CACCACCTCCTTGACCTGAGTCAGGGGGTTTCTGGCCTTTTGGGAGAACGTTCAGCGACAATGCAGTATT\n\
          TGGGCCCGTTCTTTTTTTTTCTTCTTCTTTTCTTTCTTTTTTTTTGGACTGAGTCTCTCTCGCTCTGTCA\n\
          CCCAGGCTGCGGTCGCGGTGGCGCTCTCTCGGCTCACTGAAACCTCTGCTTCCCGGGTTCCAGTGATTCT\n\
          TCTTCGGTAGCTGGGATTACAGGCGCACACCATGACGGCGGGCTCATATTCCTATTTTCAGTAGAGACGG\n\
          GGTTTCTCCACGTTGGCCACGCTGGTCTCGAACTCCTGACCTCAAATGATCCGCCTTCCTGGGCCTCCCA\n\
          AAGTGCTGGAAACGACAGGCCTGAGCCGCCGGGATTTCAGCCTTTAAAAGCGCGGCCCTGCCACCTTTCG\n\
          CTGTGGCCCTTACGCTCAGAATGACGTGTCCTCTCTGCCGTAGGTTGACTCCTTGAGTCCCCTAGGCCAT\n\
          TGCACTGTAGCCTGGGCAGCAAGAGCCAAACTCCGNNCCCCCACCTCCTCGCGCACATAATAACTAACTA\n\
          ACAAACTAACTAACTAACTAAACTAACTAACTAACTAAAATCTCTACACGTCACCCATAAGTGTGTGTTC\n\
          CCGTGAGAGTGATTTCTAAGAAATGGTACTGTACACTGAACGCAGTGGCTCACGTCTGTCATCCCGAGGT\n\
          CAGGAGTTCGAGACCAGCCCGGCCAACGTGGTGAAACCCCGTCTCTACTGAAAATACGAAATGGAGTCAG\n\
          GCGCCGTGGGGCAGGCACCTGTAACCCCAGCTACTCGGGAGGCTGGGGTGGAAGAATTGCTTGAACCTGG\n\
          CAGGCGGAGGCTGCAGTGACCCAAGATCGCACCACTGCACTACAGCCTGGGCGACAGAGTGAGACCCGGT\n\
          CTCCAGATAAATACGTACATAAATAAATACACACATACATACATACATACATACATACATACATACATAC\n\
          ATCCATGCATACAGATATACAAGAAAGAAAAAAAGAAAAGAAAAGAAAGAGAAAATGAAAGAAAAGGCAC\n\
          TGTATTGCTACTGGGCTAGGGCCTTCTCTCTGTCTGTTTCTCTCTGTTCGTCTCTGTCTTTCTCTCTGTG\n\
          TCTCTTTCTCTGTCTGTCTGTCTCTTTCTTTCTCTCTGTCTCTGTCTCTGTCTTTGTCTCTCTCTCTCCC\n\
          TCTCTGCCTGTCTCACTGTGTCTGTCTTCTGTCTTACTCTCTTTCTCTCCCCGTCTGTCTCTCTCTCTCT\n\
          CTCTCCCTCCCTGTTTGTTTCTCTCTCTCCCTCCCTGTCTGTTTCTCTCTCTCTCTTTCTGTCTGTTTCT\n\
          GTCTCTCTCTGTCTGTCTATGTCTTTCTCTGTCTGTCTCTTTCTCTGTCTGTCTGCCTCTCTCTTTCTTT\n\
          TTCTGTGTCTCTCTGTCGGTCTCTCTCTCTCTGTCTGTCTGTCTGTCTCTCTCTCTCTCTCTCTGTGCCT\n\
          ATCTTCTGTCTTACTCTCTTTCTCTGCCTGTCTGTCTGTCTCTCCCTCCCTTTCTGTTTCTCTCTCTCTC\n\
          TCTCTCTCTCTCCCCCTCTCCCTGTCTGTTTCTCTCCGTCTCTCTCTCTTTCTGTCTGTTTCTCACTGTC\n\
          TCTCTCTGTCCATCTCTCTCTCTCTCTGTCTGTCTCTTTCGTTCTCTCTGTCTGTCTGTCTCTCTCTCTC\n\
          TCTCTCTCTCTCTCTCTCTCTCCCTGTCTGTCTGTTTCTCTCTATCTCTCGCTGTCCATCTCTGTCTTTC\n\
          TATGTCTGTCTCTTTCTCTGTCAGTCTGTCAGACACCCCCGTGCCGGGTAGGGCCCTGCCCCTTCCACGA\n\
          AAGTGAGAAGCGCGTGCTTCGGTGCTTAGAGAGGCCGAGAGGAATCTAGACAGGCGGGCCTTGCTGGGCT\n\
          TCCCCACTCGGTGTATGATTTCGGGAGGTCGAGGCCGGGTCCCCGCTTGGATGCGAGGGGCATTTTCAGA\n\
          CTTTTCTCTCGGTCACGTGTGGCGTCCGTACTTCTCCTATTTCCCCGATAAGCTCCTCGACTTCAACATA\n\
          AACGGCGTCCTAAGGGTCGATTTAGTGTCATGCCTCTTTCACCGCCACCACCGAAGATGAAAGCAAAGAT\n\
          CGGCTAAATACCGCGTGTTCTCATCTAGAAGTGGGAACTTACAGATGACAGTTCTTGCATGGGCAGAACG\n\
          AGGGGGACCGGGNACGCGGAAGCCTGCTTGAGGGRGGAGGGGYGGAAGGAGAGACAGCTTCAGGAAGAAA\n\
          ACAAAACACGAATACTGTCGGACACAGCACTGACTACCCGGGTGATGAAATCATCTGCACACTGAACACC\n\
          CCCGTCACAAGTTTACCTATGTCACAGTCTTGCTCATGTATGCTTGAACGACAAATAAAAGTTCGGGGGG\n\
          GAGAAGAGAGGAGAGAGAGAGAGAGACGGGGAGAGAGGGGGGAGAGGGGGGGGGAGAGAGAGAGAGAGAG\n\
          AGAGAGAGAGAGAGAGAGAGAGAAAGAGAAGTAAAACCAACCACCACCTCCTTGACCTGAGTCAGGGGGT\n\
          TTCTGGCCTTTTGGGAGAACGTTCAGCGACAATGCAGTATTTGGGCCCGTTCTTTTTTTCTTCTTCTTCT\n\
          TTTCTTTCTTTTTTTTTGGACTGAGTCTCTCTCGCTCTGTCACCCAGGCTGCGGTGCGGTGGCGCTCTCT\n\
          CGGCTCACTGAAACCTCTGCTTCCCGGGTTCCAGTGATTCTTCTTCGGTAGCTGGGATTACAGGTGCGCA\n\
          CCATGACGGCCGGCTCATCGTTCTATTTTTAGTAGAGACGGGGTTTCTCCACGTTGGCCACGCTGGTCTC\n\
          GAACTCCTGACCACAAATGATCCACCTTCCTGGGCCTCCCAAAGTGCTGGAAACGACAGGCCTGAGCCGC\n\
          CGGGATTTCAGCCTTTAAAAGCGCGCGGCCCTGCCACCTTTCGCTGCGGCCCTTACGCTCAGAATGACGT\n\
          GTCCTCTCTGCCATAGGTTGACTCCTTGAGTCCCCTAGGCCATTGCACTGTAGCCTGGGCAGCAAGAGCC\n\
          AAACTCCGTCCCCCCACCTCCCCGCGCACATAATAACTAACTAACTAACTAACTAACTAAAATCTCTACA\n\
          CGTCACCCATAAGTGTGTGTTCCCGTGAGGAGTGATTTCTAAGAAATGGTACTGTACACTGAACGCAGGC\n\
          TTCACGTCTGTCATCCCGAGGTCAGGAGTTCGAGACCAGCCCGGCCCACGTGGTGAAACCCCCGTCTCTA\n\
          CTGAAAATACGAAATGGAGTCAGGCGCCGTGGGGCAGGCACCTGTAACCCCAGCTACTCGGGAGGCTGGG\n\
          GTGGAAGAATTGCTTGAACCTGGCAGGCGGAGGCTGCAGTGACCCAAGATCGCACCACTGCACTACAGCC\n\
          TGGGCGACAGAGTGAGACCCGGTCTCCAGATAAATACGTACATAAATAAATACACACATACATACATACA\n\
          TACATACAACATACATACATACAGATATACAAGAAAGAAAAAAAGAAAAGAAAAGAAAGAGAAAATGAAA\n\
          GAAAAGGCACTGTATTGCTACTGGGCTAGGGCCTTCTCTCTGTCTGTTTCTCTCTGTTCGTCTCTGTCTT\n\
          TCTCTCTGTGTCTCTTTCTCTGTCTGTCTGTCTGTCTGTCTGTCTGTCTCTTTCTTTCTTTCTGTCTCTG\n\
          TCTTTGTCCCTCTCTCTCCCTCTCTGCCCTGTCTCACTGTGTCTGTCTTCTATCTTACTCTCTTTCTCTC\n\
          CCCGTCTGTCTCTCTCTCACTCCCTCCCTGTCTGTTTCTCTCTCTCTCTCTTTCTGTCTGTTTCTGTCTC\n\
          TCTCTGTCTGCCTCTCTCTTTCTCTATCTGTCTCTTTCTCTGTCTGTCTGCCCCTCTCTTTCTTTTTCTG\n\
          TGTCTCTCTGTCTGTCTCTCTCTCTCTCTGTGCCTATCTTCTGTCTTACTCTCTTTCTCTGCCTGTCTGT\n\
          CTGTCTCTCTCTGTCTCTCCCTCCCTTTCTGCTTCTCTCTCTCTCTCTCTCTCTNNNCCCTCCCTGTCTG\n\
          TTTCTCTCTGTCTCCCTCTCTTTCTGTCTGTTTCTCACTGTCTCTCTCTGTCTGTCTGTTTCATTCTCTC\n\
          TGTCTCTGTCTCTGTCTCTCTCTCTCTCTGTCTCTCCCTCTCTGTGTGTATCTTTTGTCTTACTCTCCTT\n\
          CTCTGCCTGTCCGTCTGTCTGTCTGTCTCTCTCTCTCCCTGTCCCTCTCTCTTTCTGTCTGTTTCTCTCT\n\
          CTCTCTCTCTCTCTCTCTCTCTGTCTCTGTCTTTCTCTGTCTGTCCCTTTCTCTGTCTGTCTGCCTCTCT\n\
          CTTTCTCTTTCTGTGTCTCTCTGTCTCTCTCTCTGTGCCTATCTTCTGTCTTACTCTCTTTCTCTGCCTG\n\
          TCTATCTGTCTGTCTCTCTCTGTCTCTCTCCCTGCCTTTCTGTTTCTCTCTCTCTCCCTCTCTCGCTCTC\n\
          TCTGTCTTTCTCTCTTTCTCTCTGTTTCTCTGTCTCTCTCTGTCCGTCTCTGTCTTTTTCTGTCTGTCTG\n\
          TCTCTCTCTTTCTTTCTGTCGTCTGTCTCTGTCTCTGTCTCTGTCTCTCTCTCTCTCTCTCTCCTTGTCT\n\
          CTCTCACTGTGTCTGTCTTCTGTCTTACTCTCCTTCTCTGCCTGTCCATCTGTCTGTCTGTCTCTCTCTC\n\
          TCTCTCCCTACCTTTCTGTTTCTCTCTCGCTAGCTCTCTCTCTCTCTGCCTGTTTCTCTCTTTCTCTCTC\n\
          TGTCTTTCTCTGTCTGTCTCTTTCTCTGTCTGTCTGTCTCTTTCTCTCTGTCTCTGTCTCTGTCTCTCTC\n\
          TCTCTCTCTCTCTCTCTCTCTGCCTCTCTCACTGTGTCTGTCTTCTGTCTTATTCTCTTTCTCTCTCTGT\n\
          CTCTCTCTCTCTCTCCTTTACTGTCTGTTTCTCTCTCTCTCTCTCTCTTTCTGCCTGTTTCTCTCTGTCT\n\
          GTCTCTGTCTTTCTCTGTCTGTCTGCCTCTCTCTTTCTTTTTCTGCGTCTCTCTGTCTCTCTCTCTCTCT\n\
          CTCTGTTCCTATCTTCTGTCTTACTCTGTTTCCTTGCCTGCCTGCCTGTCTGTGTGTCTGTCTCTCTCTC\n\
          TCTCTCTCTCTCTCTCTCCCTCCCTTTCTCTTTCTCTGTCTCTCTCTCTCTTTCTGGGTGTTTCTCTCTG\n\
          TCTCTCTGTCCATCTCTGTCTTTCTATGTCTGTCTCTCTCTTTCTCTCTGTCTCTGTCTCTGCCTCTCTC\n\
          TCTCTCTCTCTCTCTCTCTCTCTGTCTGTCTCTCTCACTGTGTGTGTCTGTCTTCTGTCTTACTCTCCTT\n\
          CTCTGCCTGTCCGTCTGTCTGTCTGTCTCTCCCTCTCTCTCCCTCCCTTTCTGTTTCTCTCTCTCTCTCT\n\
          TTCTGTCTGTTTCTCTCTTTCTCTCTCTGTCTGTCTCTTTCTCTGTCTGTCTGTCTCTCTCTTTCTTTTT\n\
          CTCTGTCTCTCTGTCTCTCTCTGTGTCTGTCTCTCTGTCTGTGCCTATCTTCTGTCTTACTCTCTTTCTC\n\
          TGGCTGTCTGCCTGTCTCTCTCTCTCTCTCTGTCTGTCTCCGTCCCTCTCTCCCTGTCTGTCTGTTTCTC\n\
          TCTCTGCCTCTCTCTCTCTCTGTCTGTCTCTTTCTCTGTCTGTCTGTCTCTCTCTTTCTTTTTCTCTGTC\n\
          TCTCTGTCTCTCTCTGTGTCTGTCTCTCTTTCTGTGCCTATCTTCTGTCTTACTCTCTTTCTCTGGCTGT\n\
          CTGCCTGTCTCTCTCTCTCTGCCTGTCTCCGTCCCTCCCTCCCTGTCTGTCTGTTTCTCTCTCTGTCTCT\n\
          GTCTCTCTGTCCATCTCTGTCTGTCTCTTTCTCTTTCTCTCTCTCTGTCTCTGTCTCTCTCTCTCTCTGC\n\
          CTGTCTCTCTCACTGTGTCTGTCTTCTGTCTTACTCTCTTTCTCTTGCCTGCCTCTCTGTCTGTCTGTCT\n\
          CTCTCCCTCCATGTCTCTCTCTCTCTCTCACTCACTCTCTCTCCGTCTCTCTCTCTTTCTGTCTGTTTCT\n\
          CTCTCTGTCTGTCTCTCTCCCTCCATGTCTCTCTCTCTCTCTCTCACTCACTCTCTCTCCGTCTCTCTCT\n\
          CTCTTTCTGTCTGTTTCTCTCTCTGTCTGTCTCTCTCCCTCCATGTCTCTCTCTCTCCCTCTCACTCACT\n\
          CTCTCTCCGTCTCTCTCTCTCTTTCTGTCTGTTTCTTTGTCTGTCTGTCTGTCTGTCTGTCTGTCTCTCT\n\
          CTCTCTCTCTCTCTCTCTCTCTCTCTGTTTGTCTTTCTCCCTCCCTGTCTGTCTGTCTGTCTCTCTCTCT\n\
          CTGTCTCTGTCTCTGTCTCTCTCTCTTTCTCTTTCTGTCTGTTTCTCTCTATCTCTCGCTGTCCATCTCT\n\
          GTCTTTCTATGTCTGTCTCTTTCTCTGTCAGTCTGTCAGACACACCCGTGCCGGTAGGGCCCTGCCCTTC\n\
          CACGAGAGTGAGAAGCGCGTGCTTCGGTGCTTAGAGAGGCCGAGAGGAATCTAGACAGGCGGGCCTTGCT\n\
          GGGCTTCCCCACTCGGTGTACGATTTCGGGAGGTCGAGGCCGGGTCCCCGCTTGGATGCGAGGGGCATTT\n\
          TCAGACTTTTCTCTCGGTCACGTGTGGCGTCCGTACTTCTCCTATTTCCCCGATAAGTCTCCTCGACTTC\n\
          AACATAAACTGTTAAGGCCGGACGCCAACACGGCGAAACCCCGTCTCTACTAAAAATACAAAGCTGAGTC\n\
          GGGAGCGGTGGGGCAGGCCCTGTAATGCCAGCTCCTCGGGAGGCTGAGGCGGGAGAATCGCTTGAACCAG\n\
          GGAAGCGGAGGCTGCAGGGAGCCGAGATCGCGCCACTGCACTACGGCCCAGGCTGTAGAGTGAGTGAGAC\n\
          TCGGTCTCTAAATAAATACGGAAATTAATTAATTCATTAATTCTTTTCCCTGCTGACGGACATTTGCAGG\n\
          CAGGCATCGGTTGTCTTCGGGCATCACCTAGCGGCCACTGTTATTGAAAGTCGACGTTGACACGGAGGGA\n\
          GGTCTCGCCGACTTCACCGAGCCTGGGGCAACGGGTTTCTCTCTCTCCCTTCTGGAGGCCCCTCCCTCTC\n\
          TCCCTCGTTGCCTAGGGAACCTCGCCTAGGGAACCTCCGCCCTGGGGGCCCTATTGTTCTTTGATCGGCG\n\
          CTTTACTTTTCTTTGTGTTTTGGCGCCTAGACTCTTCTACTTGGGCTTTGGGAAGGGTCAGTTTAATTTT\n\
          CAAGTTGCCCCCCGGCTCCCCCCACTACCCACGTCCCTTCACCTTAATTTAGTGAGNCGGTTAGGTGGGT\n\
          TTCCCCCAAACCGCCCCCCCCCCCCCGCCTCCCAACACCCTGCTTGGAAACCTTCCAGAGCCACCCCGGT\n\
          GTGCCTCCGTCTTCTCTCCCCTTCCCCCACCCCTTGCCGGCGATCTCATTCTTGCCAGGCTGACATTTGC\n\
          ATCGGTGGGCGTCAGGCCTCACTCGGGGGCCACCGTTTTTGAAGATGGGGGCGGCACGGTCCCACTTCCC\n\
          CGGAGGCAGCTTGGGCCGATGGCATAGCCCCTTGACCCGCGTGGGCAAGCGGGCGGGTCTGCAGTTGTGA\n\
          GGCTTTTCCCCCCGCTGCTTCCCGCTCAGGCCTCCCTCCCTAGGAAAGCTTCACCCTGGCTGGGTCTCGG\n\
          TCACCTTTTATCACGATGTTTTAGTTTCTCCGCCCTCCGGCCAGCAGAGTTTCACAATGCGAAGGGCGCC\n\
          ACGGCTCTAGTCTGGGCCTTCTCAGTACTTGCCCAAAATAGAAACGCTTTCTGAAAACTAATAACTTTNC\n\
          TCACTTAAGATTTCCAGGGACGGCGCCTTGGCCCGTGTTTGTTGGCTTGTTTTGTTTCGTTCTGTTTTGT\n\
          TTTGTTCGTGTTTTTCCTTTCTCGTATGTCTTTCTTTTCAGGTGAAGTAGAAATCCCCAGTTTTCAGGAA\n\
          GACGTCTATTTTCCCCAAGACACGTTAGCTGCCGTTTTTTCCTGTTGTGAACTAGCGCTTTTGTGACTCT\n\
          CTCAACGCTGCAGTGAGAGCCGGTTGATGTTTACNATCCTTCATCATGACATCTTATTTTCTAGAAATCC\n\
          GTAGGCGAATGCTGCTGCTGCTCTTGTTGCTGTTGTTGTTGTTGTTGTTGTCGTCGTTGCTGTTGTCGTT\n\
          GTCGTTGTTGTTGTCGTTGTCGTTGTTTTCAAAGTATACCCCGGCCACCGTTTATGGGATCAAAAGCATT\n\
          ATAAAATATGTGTGATTATTTCTTGAGCACGCCCTTCCTCCCCCTCTCTCTGTCTCTCTGTCTGTCTCTG\n\
          TCTCTCTCTTTCTCTGTCTGTCTTCTCTCTCTCTCTCTCTCTGTGTCTCTCTCTCTCTGCCTGTCTGTTT\n\
          CTCTCTCTCTGCCTCTCTCTCTCTCTCTCTCTCTGCCTGTCTCTCTCACTGTGTCTGTCTTCTGTCTTAC\n\
          TCCCTTTCTCTGTCTGTCTGTCGGTCTCTCTCTCTCTCTCTCCCTGTCTGTATGTTTCTCTCTGTCTCTG\n\
          TCTCTCTCTCTCTTTCTGTTTCTCTCTCTCCGTCTCTGTCTTTCTCTGACTGTCTCTCTCTTTCCTTCTC\n\
          TCTGTCTCTCTCTGCCTGTCTCTCTCACTCTGTCTTCTGTCTTATCTCTCTCTCTGCCTGCCTGTCTCTC\n\
          TCACTCTCTCTCTCTGTGTGTCTCTCTCTCTCTTTCTGTTTCTCTCTGTCTCTCTGTCCGTCTCTGTCTT\n\
          TCTCTGTCTGTCTCTTTGTCTGTCTGTCTTTGTCTTTCCTTCTCTCTGTCTCTGTCTCTCTCACTGTGTC\n\
          TGTCTTCTGTCTTAGTCTCTCTCTCTCTCTCTCCCTGTCTGTCTGTCTCTCTCTCTCTCTCCCCCTGTCT\n\
          GTTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTGTCTTTGTCTTTCTTTCTGTCTCTGTCTCTCTCTCT\n\
          CTCTCTGTGTGTCTGTCTTCTGTCTTACTGTCTTTCTCTGCCTGTCTGTCTGTCTGTCTCTCTCTGTCTG\n\
          TCTCTCTCTCTCTCTCCCCCTGTCGGCTGTTTCTCTGTCTCTGTCTGTGTCTCTCTTTCTGTCTGTTTCT\n\
          CTCTGTCTGTCTTTCTCTCTCTGTCTCTTTCTCTCTGTCTCTCTGTCTGTCTCTGTCTCTCTCTCTGTCT\n\
          CTCTCTCTCTGTGGGGGTGTGTGTGTGTGTGTGTATGTGTGTGTGTGTGTGTGTGTGTGTCTGCCTTCTG\n\
          TCTTACTCTCTTTCTCTGCCTGTCTGTCTGCCTGTCTGTTTGTCTCTCTCTCTCTGCCTGTCTCTCTCCC\n\
          TTCCTGTCTGTTTCTCTCTCTTTCTGTTTCTCTCTGTCTCTGTCCATCTCTGTCTTTCTCCGTCTGTCTC\n\
          TTTATCTGTCTCTCTCCGTCTGTCTCTTTATCTGTCTCTCTCTCTCTTTCTGTCTTTCTCTCTCTGTGTA\n\
          TCGTTGTCTCTCTCTGTCTGTCTCTGTCTCTGTCTCTCTGTCTCTCTCTCTCTCTCTCTCTCTCTGTCTG\n\
          TCTGTCCGTCTGTCTGTCTCGGTCTCTGCGTCTCGCTATCTCCCGCCCTCTCTTTTTTTGCAAAAGAAGC\n\
          TCAAGTACATCTAATCTAATCCCTTACCAAGGCCTGAATTCTTCACTTCTGACATCCCAGATTTGATCTC\n\
          CCTACAGAATGCTGTACAGAACTGGCGAGTTGATTTCTGGACTTGGATACCTCATAGAAACTACATATGA\n\
          ATAAAGATCCAATCCTAAAATCTGGGGTGGCTTCTCCCTCGACTGTCTCGAAAAATCGTACCTCTGTTCC\n\
          CCTAGGATGCCGGAAGAGTTTTCTCAATGTGCATCTGCCCGTGTCCTAAGTGATCTGTGACCGAGCCCTG\n\
          TCCGTCCTGTCTCAAATATGTACGTGCAAACACTTCTCTCCATTTCCACAACTACCCACGGCCCCTTGTG\n\
          GAACCACTGGCTCTTTGAAAAAAATCCCAGAAGTGGTTTTGGCTTTTTGGCTAGGAGGCCTAAGCCTGCT\n\
          GAGAACTTTCCTGCCCAGGATCCTCGGGACCATGCTTGCTAGCGCTGGATGAGTCTCTGGAAGGACGCAC\n\
          GGGACTCCGCAAAGCTGACCTGTCCCACCGAGGTCAAATGGATACCTCTGCATTGGCCCGAGGCCTCCGA\n\
          AGTACATCACCGTCACCAACCGTCACCGTCAGCATCCTTGTGAGCCTGCCCAAGGCCCCGCCTCCGGGGA\n\
          GACTCTTGGGAGCCCGGCCTTCGTCGGCTAAAGTCCAAAGGGATGGTGACTTCCACCCACAAGGTCCCAC\n\
          TGAACGGCGAAGATGTGGAGCGTAGGTCAGAGAGGGGACCAGGAGGGGAGACGTCCCGACAGGCGACGAG\n\
          TTCCCAAGGCTCTGGCCACCCCACCCACGCCCCACGCCCCACGTCCCGGGCACCCGCGGGACACCGCCGC\n\
          TTTATCCCCTCCTCTGTCCACAGCCGGCCCCACCCCACCACGCAACCCACGCACACACGCTGGAGGTTCC\n\
          AAAACCACACGGTGTGACTAGAGCCTGACGGAGCGAGAGCCCATTTCACGAGGTGGGAGGGGTGGGGGTG\n\
          GGGTGGGTTGGGGGTTGTGGGGTCTGTGGCGAGCCCGATTCTCCCTCTTGGGTGGCTACAGGCTAGAAAT\n\
          GAATATCGCTTCTTGGGGGGAGGGGCTTCCTTAGGCCATCACCGCTTGCGGGACTACCTCTCAAACCCTC\n\
          CCTTGAGGCCACAAAATAGATTCCACCCCACCCATCGACGTTTCCCCCGGGTGCTGGATGTATCCTGTCA\n\
          AGAGACCTGAGCCTGACACCGTCGAATTAAACACCTTGACTGGCTTTGTGTGTTTGTTTGTTTCTGAGAT\n\
          GGAGTCTTGCTCTGTCCCCCAGGCTGGAGTGCAGTGGCGTGATCTCAGCTCACTGGAACCTCTGCCTCCT\n\
          GGGTTCAAGTGATTCTCCTGTCTCAGCGCCACCATGGCCGGCTCATTTTTTTTTTTTTTTTTTTTGGTAG\n\
          ACACGGGGTTTCACCCTCTTTCATTGGTTTTCACTGGAGATTCTAGATTCGAGCCACACCTCATTCCGTG\n\
          CCACAGAGAGACTTCTTTTTTTTTTTTTTTTTTTTAAGCGCAACGCAACATGTCTGCCTTATTTGAGTGG\n\
          CTTCCTATATCATTATAATTGTGTTATAGATGAAGAAACGGTATTAAACACTGTGCTAATGATAGTGAAA\n\
          GTGAAGACAAAAGAAAGGCTATCTATTTTGTGGTTAGAATAAAGTTGCTCAGTATTTAGAAGCTACCTAA\n\
          ATACGTCAGCATTTACACTCTTCCTAGTAAAAGCTGGCCGATCTGAATAATCCTCCTTTAAACAAACACA\n\
          ATTTTTGATAGGGTTAAGATTTTTTTAAGAATGCGACTCCTGCAAAATAGCTGAACAGACGATACACATT\n\
          TAAAAAAATAACAACACAAGGATCAACCAGACTTGGGAAAAAATCGAAAACCACACAAGTCTTATGAAGA\n\
          ACTGAGTTCTTAAAATAGGACGGAGAACGTAGCTATCGGAAGAGAAGGCAGTATTGGCAAGTTGATTGTT\n\
          ACGTTGGTCAGCAGTAGCTGGCACTATCTTTTTGGCCATCTTTCGGGCAATGTAACTACTACAGCAAAAT\n\
          GAGATATGATCCATTAAACAACATATTCGCAAATCAAAAAGTGTTTCAGTAATATAATGCTTCAGATTTA\n\
          GAAGCAAATCAAATGATAGAACTCCACTGCTGTAATAAGTCACCCCAAAGATCACCGTATCTGACAAAAT\n\
          AACTACCACAGGGTTATGACTTCAGAATCATACTTTCTTCTTGATATTTACTTATGTATTTATTTTTTTT\n\
          AATTTATTTCTCTTGAGACGCGTCTCGCTCTGTCGCCCAGGCTGGAGTGCGATGGTGTGATCTCGGCTCA\n\
          CTGCAACCGCCACCTCCCTGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGG\n\
          TGCCCGCCACCACGCCCAGCTAATCTTTATACTTTTAATAGAGACGGGGTTTCACCGTGTCGGCCCGGAT\n\
          GGTCTCGATCTCTTGACCTCGTGACCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATGACAGGCGTGAGC\n\
          CACTGAGCCCGGCCTTCTCTTGACGTTTAAACTATGAAGTCAGTCCAGAGAAACGCAATAAATGTCAACG\n\
          GTGAGGATGGTGTTGAGGCAGAAGTAGGACCACACTTTTTCCTATCTTATTCAGTTGATAACAATATGAC\n\
          CTAGGTAGTAATTTCCTATGTGCCTACTTATACACGAGTACAAAAGAGTAAAACAGAGAGACTGCTAAAT\n\
          TAAAGGGTACGTGAAGTTCTTCATAGTAACTCCGTAAACTGGAACACTGTCAAAAAGCAGCAGCTAGTGA\n\
          ATTGTTTCCATGTATTTTTCTATTATCCAATAAGTGAACTATGCTATTCCTTTCCAGTCTCCCAAGCACT\n\
          TCTTGTCCCCATCACCACTTCGGTGCTCGAAGAAAAAGTAAGCAAATCAAGGAACACAAGCTAAAGAAAC\n\
          ACACACACAAACCAAAGACAACTACAGCGTCTGCAAAAGTTTGCTAGAAGACTGAAACTGTTGAGTATAA\n\
          GGATCTGGTATTCTACGATCATGAGTTCACTTCAGAGTTTGTTCAAGACATACGTTTCGTAAGGAAACAT\n\
          CTTAGTTAGAAGTTATTCAGCAGTAGGTACCATCCCTAAGTATTTTTCACCAAATCCGTGACAATAAAGA\n\
          GCTATCTAACCAGAAAAATTAGCGAGTACGGGCACCATCCATAGGGCTTTGTCTTTACGCTTCATTAGCA\n\
          CTTACCATGCCTTACAATGTCTAGGATTGACCCTGATAGCATTTCGAAAACAAGCTAATGCTTTGTCCAG\n\
          TTCTTCAGTGAAGACAACTCACGCCCTAATGCGCTATAGGCATAAGCATCATTTGGATCCACTTCGAGAG\n\
          TTCTCTGGAAGAATTGAATCGCAATATCGTGTTCCCGTTTGCAGACCGAAACAGTTTCCCTGCAGCACAC\n\
          CAGGCCTCTGGCTGGCGAATTTTTATCCATGTCTGTGAAGTCTTTGGACAGAACTGAAAGAGCAACCTCT\n\
          TTCGGAGGATGCCAAAGTGTTGTAGAGTAGATCTCCATGCCTTCGACTCTGTAATTCTCAATCCTCCTAA\n\
          CCTCTGAGAATTGTCTTTCAGCTTGCGTGGACTCTGAAAGTTTACAATAGGCCNTTTCCGATTTGGCACA\n\
          GTACCCAACCGGTATTGCAGTGGTGAGAAGCTAGATGGCTCAAGATGCTGATAGCTTCTTTGCCGTGGTA\n\
          AGAACACAAAGCTAAATAACCTTTCCCCCTTTCACGAAGAAGGCTCATCAAGCCTTCCGCTGCTGCTTTT\n\
          TGTAGATTAAAAGCCTGAATCTGAGGCGCGATTGCGGCTATTTTCCCTTCTGAAATGACGGAAGAGTCCA\n\
          ATTTTGTCACTTCCAGGCTATCACTTATGTTCGGTGGAGTTATTGCTCCTTTATTAGTTTTACTTTTGGT\n\
          TCTTCTGTTTGGGATTTTAGGTGGAAACTTCATTTTTAATTTTCTCCTAATTCTCCTCGGTTGTGGAGCT\n\
          GTCACTAGTCAAGAGTCGTGAATTTCTTCGAGGNCGGTGCATTTGGGGGAGATGCCATAGTGGGGCTCAA\n\
          TACCTGAGGTGTTGCCCTTGTCGGCGGACCAGAACTTTGTGTTTTTGCAAGGACTGGAGTTACCTTTCGG\n\
          CTCTTTCCCCTCTGCGAGAAGACAGACGGTGTTCCGGTTTGGCCGATTCTGGCAACAGGCTTTTCTGAAG\n\
          GGGCTCCGGTGGATGGCACGTCAGTGACAGACGGTGTCTCATACCAGTGCAGTTTTGTCAATAGGGTCCG\n\
          TCTCCGGGACTTGGGGTTTCTAATGGCAAAATGCCAACACTTGGGGTTAATGGACTAACAGCTGCTGGTC\n\
          CTCCTAATAAACTTCGACCAGTTTTTGGTTTATGTTGAACCTGTTTAGATCATATGGAAGTTCCTGTTCC\n\
          CAGTGGGACAGTATCAGGTGAAAGGACAGCTGAATCGATAGAAGACACTGGGGAGTCTGTATTCAAGGAG\n\
          TACTTTGAATTGGAAGATTCTAAATTCCATCCGTTTCATTCGACGGTGTCCTGGGGTGTTTCCGTAAGAA\n\
          CGGTCTCGGGCTGTCTGTGACATAAACTAGGACGAGGTCCAAGTGTTGTGGCGCAACACTTGGACAGGCA\n\
          GTTGCTAAAGCTCTCTAGAGAGGTGAATCAAAATGTTTGGTCAGGATCTGGCTTTTCCCCCCTATTTCAC\n\
          ATCATGATTCAAAGGGACACCAGAGGAAAGGATTTCAACGAAGGCTCTTTTGGTCACATTCTGATCCTTT\n\
          GGTAAGCCGATCTGTCTTGCAATATACATGTCCCGACGATGGAAGGGGAAAGCGAGCTGAATCACCAAAC\n\
          TCAGGAACGATAATATCATCGTGGCTTTTCTGCTTATGAAACACTCCACCCGATAAGATTTGATCCCCTT\n\
          CTGCAAGCTTGCTGAGATCAACACAACATTTCGCAAGCAGGCATTTGCATTGCGGGGTAGTACAACTGTG\n\
          TCCTTTCAAGAGTCTATATGTTTTATAGGCCTTTCCTGAGCGGTAAGAACAGGTCGCCAGTAAGAACAAG\n\
          GCTTCTTCTGAGTGTACTTCTGCATAAAGGCGTTCTGCGGGGGAAACCGCATCTCGGTAGGCATAGTGGT\n\
          TTAGTGCTTGCCATATAGCAGCCTGGACGGGTCCCTGCAGCACCGCCATCCTCGAGGCTCAGGCCCACTT\n\
          TCTGCAGTGCCACAGGCACCCCCCCCCCCCCATAGCGGCTCCGGCCCGGCCAGCCCCGGCTCATTTAAAG\n\
          GCACCAGCCGCCGTTACCGGGGGATGGGGGAGTCCGAGACAGAATGACTTCTTTATCCTGCTGACTCTGG\n\
          AAAGCCCGGCGCCTTGTGATCCATTGCAAACCGAGAGTCACCTCGTGTTTAGAACACGGATCCACTCCCA\n\
          AGTTCAGTGGGGGGATGTGAGGGGTGTGGCAGGTAGGACGAAGGACTCTCTTCCTTCTGATTCGGTCTGC\n\
          ACAGTGGGGCCTAGGGCTGGAGCTCTCTCCGTGCGGACCGCTGACTCCCTCTACCTTGGGTTCCCTCGGC\n\
          CCCACCCTGGAACGCCGGGCCTTGGCAGATTCTGGCCCTTTCTGGCCCTTCAGTCGCTGTCAGAAACCCC\n\
          ATCTCATGCTCGGATGCCCCGAGTGACTGTGGCTCGCACCTCTCCGGAAACATTGGAAATCTCTCCTCTA\n\
          CGCGCGGCCACCTGAAACCACAGGAGCTCGGGACACACGTGCTTTCGGGAGAGAATGCTGAGAGTCTCTC\n\
          GCCGACTCTCTCTTGACTTGAGTTCTTCGTGGGTGCGTGGTTAAGACGTAGTGAGACCAGATGTATTAAC\n\
          TCAGGCCGGGTGCTGGTGGCTCACGCCTGTAACCCCAACACTTTGGGAGGCCGAGGCCGTAGGATCCCTC\n\
          GAGGAATCGCCTAACCCTGGGGAGGTTGAGGTTGCAGTGAGTGAGCCATAGTTGTGTCACTGTGCTCCAG\n\
          TCTGGGCGAAAGACAGAATGAGGCCCTGCCACAGGCAGGCAGGCAGGCAGGCAGGCAGAAAGACAACAGC\n\
          TGTATTATGTTCTTCTCAGGGTAGGAAGCAAAAATAACAGAATACAGCACTTAATTAATTTTTTTTTTTT\n\
          CCTTCGGACGGAGTTTCACTCTTGGTGCCCACGCTGGAGTGCAGTGGCACCATCTCGGCTCACCGCAACC\n\
          TCCACCTCCCGCGTTCAAGCGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGGAGGAGCCA\n\
          CCACACCCAGCTGATTTTGTATTGTTAGTAGAGACGGCATTTCTCCATGTGGGTCAGGCTGGTCTCGAAC\n\
          TGGCGACCCCAGTGGATCTGCCCGCCCCGGCCTCCCAAAGTGCTGGGGTGACAGGCGTGAGCCATCGTGA\n\
          CTGGCCGGCTACGTTTATTTATTTATTTTTTTAATTATTTTACTTTTTTTTAGTTTTCCATTTTAATCTA\n\
          TTTATTTATTTACATTTATTTATTTATTTATTTATTTACTTATTTATTTATTTTCGAGACAGACTCTCGC\n\
          TCTGCTGCCCAGGCTGGAGTGCAGCGGCGTGATCTCGGCTCACTGCAACGTCCGCCTCCCGGGTTCACGC\n\
          CATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACCGTGCCCGGCTAACTTTTT\n\
          GTATTTTGAGTAGAGATGGGGTTTCACTGTGGTAGCCAGGATGGTCTCGATCTCCTGACCCCGTGATCCG\n\
          TCCACCTCGGCCTCCCAAAGTGCTGGGATGACAGGCGTGAGCCACCGGCCCCGGCCTATTTATCTATTTA\n\
          TTAACTTTGAGTCCAGGTTATGAAACCAGTTAGTTTTTGTAATTTTTTTTTTTTTTTTTTTTTTTTGAGA\n\
          CGAGGTTTCACCGTGTTGCCAAGGCTTGGACCGAGGGATCCACCGGCCCTCGGCCTCCCAAAAGTGCGGG\n\
          GATGACAGGCGCGAGCCTACCGCGCCCGGACCCCCCCTTTCCCCTTCCCCCGCTTGTCTTCCCGACAGAC\n\
          AGTTTCACGGCAGAGCGTTTGGCTGGCGTGCTTAAACTCATTCTAAATAGAAATTTGGGACGTCAGCTTC\n\
          TGGCCTCACGGACTCTGAGCCGAGGAGTCCCCTGGTCTGTCTATCACAGGACCGTACACGTAAGGAGGAG\n\
          AAAAATCGTAACGTTCAAAGTCAGTCATTTTGTGATACAGAAATACACGGATTCACCCAAAACACAGAAA\n\
          CCAGTCTTTTAGAAATGGCCTTAGCCCTGGTGTCCGTGCCAGTGATTCTTTTCGGTTTGGACCTTGACTG\n\
          AGAGGATTCCCAGTCGGTCTCTCGTCTCTGGACGGAAGTTCCAGATGATCCGATGGGTGGGGGACTTAGG\n\
          CTGCGTCCCCCCAGGAGCCCTGGTCGATTAGTTGTGGGGATCGCCTTGGAGGGCGCGGTGACCCACTGTG\n\
          CTGTGGGAGCCTCCATCCTTCCCCCCACCCCCTCCCCAGGGGGATCCCAATTCATTCCGGGCTGACACGC\n\
          TCACTGGCAGGCGTCGGGCATCACCTAGCGGTCACTGTTACTCTGAAAACGGAGGCCTCACAGAGGAAGG\n\
          GAGCACCAGGCCGCCTGCGCACAGCCTGGGGCAACTGTGTCTTCTCCACCGCCCCCGCCCCCACCTCCAA\n\
          GTTCCTCCCTCCCTTGTTGCCTAGGAAATCGCCACTTTGACGACCGGGTCTGATTGACCTTTGATCAGGC\n\
          AAAAACGAACAAACAGATAAATAAATAAAATAACACAAAAGTAACTAACTAAATAAAATAAGTCAATACA\n\
          ACCCATTACAATACAATAAGATACGATACGATAGGATGCGATAGGATACGATAGGATACAATACAATAGG\n\
          ATACGATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACG\n\
          CCGGGCGCGGTGGCTCATGCCTGTCATCCCGTCACTTTGGGATGCCGAGGTGGACGCATCACCTGAAGTC\n\
          GGGAGTTGGAGACAAGCCCGACCAACATGGAGAAATCCCGTCTCAATTGAAAATACAAAACTAGCCGGGC\n\
          GCGGTGGCACATGCCTATAATCCCAGCTGCTAGGAAGGCTGAGGCAGGAGAATCGCTTGAACCTGGGAAG\n\
          CGGAGGTTGCAGTGAGCCGAGATTGCGCCATCGCACTCCAGTCTGAGCAACAAGAGCGAAACTCCGTCTC\n\
          AAAAATAAATACATAAATAAATACATACATACATACATACATACATACATACATACATACATAAATTAAA\n\
          ATAAATAAATAAAATAAAATAAATAAATGGGCCCTGCGCGGTGGCTCAAGCCTGTCATCCCCTCACTTTG\n\
          GGAGGCCAAGGCCGGTGGATCAAGAGGCGGTCAGACCAACAGGGCCAGTATGGTGAAACCCCGTCTCTAC\n\
          TCACAATACACAACATTAGCCGGGCGCTGTGCTGTGCTGTACTGTCTGTAATCCCAGCTACTCGGGAGGC\n\
          CGAGCTGAGGCAGGAGAATCGCTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGC\n\
          AACCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCCAAAAAATGAAAATGAAAATGAAACGCAACAAAAT\n\
          AATTAAAAAGTGAGTTTCTGGGGAAAAAGAAGAAAAGAAAAAAGAAAAAAACAACAAAACAGAACAACCC\n\
          CACCGTGACATACACGTACGCTTCTCGCCTTTCGAGGCCTCAAACACGTTAGGAATTATGCGTGATTTCT\n\
          TTTTTTAACTTCATTTTATGTTATTATCATGATTGATGTTTCGAGACGGAGTCTCGGAGGCCCGCCCTCC\n\
          CTGGTTGCCCAGACAACCCCGGGAGACAGACCCTGGCTGGGCCCGATTGTTCTTCTCCTTGGTCAGGGGT\n\
          TTCCTTGTCTTTCTTCGTGTCTTTAACCCGCGTGGACTCTTCCGCCTCGGGTTTGACAGATGGCAGCTCC\n\
          ACTTTAGGCCTTGTTGTTGTTGGGGACTTTCCTGATTCTCCCCAGATGTAGTGAAAGCAGGTAGATTGCC\n\
          TTGCCTGGCCTTGCCTGGCCTTGCCTTTTCTTTCTTTCTTTCTTTCTTTATTACTTTCTCTTTTTCTTCT\n\
          TCTTCTTCTTCTTTTTTTTGAGACAGAGTTTCACTCTTGTTGCCCAGGCTAGAGGGCAATGGCGCGATCT\n\
          CGGCTCACCGCACCCTCCGCCTCCCAGGTTCAAGCGATTCTCCTGCCTCAGCCTCCTGATTAGCTGGGAT\n\
          TACAGGCATGGGCCACCGTGCTGGCTGATGTTTGTACTTTTAGTAGAGACGGTGTTTTTCCATGTTGGTC\n\
          AGGCTGGTCTCCCACTCCCAACCTCAGGTGGTCCGCCTGCCTTAGCCTCCCAAAGTGCTGGGATGACAGG\n\
          CGTGCAACCGCGCCCAGCCTCTCTCTCTCTCTCTCTCTCTCTCGCTCGCTTGCTTGCTTGCTTTCGTGCT\n\
          TTCTTGCTTTCCCGTTTTCTTGCTTTCTTTCTTTCTTTCGTTTCTTTCATGCTTGCTTTCTTGCTTGCTT\n\
          GCTTGCTTTCGTGCTTTCTTGCTTTCCTGTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTGTTTCTTT\n\
          CTTGCTTGCTTTCTTGCTTGCTTGCTTGCTTTCGTGCTTTCTTGCTTTCCTGTTTTCTTTCTTTCTTTCT\n\
          TTCTTTTCTTTCTTTCTTGCTTGCTTTCCTGCTTGCTTGCTTTCGTGCTTTCTTGTTTTCTCGATTTCTT\n\
          TCTTTCTTTTGTTTCTTTCCTGCTTGCTTTCTTGCTTGCTTGCTTTCGTGCTTCTTGCTTTCCTGTTTTC\n\
          TTTCTTTCTTTCTTTCTTTTGTTTCTTTCTTGCTTGCTTTCTTGCTTGCTTGCTTTCGTGCTGTCTTGTT\n\
          TCTCGATTTCTTTCTTTCTTTTGTTTCTTTCCTGCTTGCTTTCTTGCTTGATTGCTTTCGTGCTTTCTTG\n\
          CTTTCTTGTTTTCTTTCTTTCTTTTGTTTCTTTCTTTCTTGCTTCCTTGTTTTCTTGCTTTCTTGCTTGC\n\
          TTGCTTTCGTGCTTTCTTGTTTTCTTGCTTTCTTTCTTTTGTTTCTTTCTTGCTTGCTTTCTTGCTTCCT\n\
          TGTTTTCTTGCTTTCTTGCTTGCTTGCTTTCGTGCTTTCTTTCTTGCTTTCTTTTCTTTCTTTCTTTTCT\n\
          TTTTCTTTCTTTCTTGCTTTCTTTTCTTTCATCATCATCTTTCTTTCTTTCCTTTCTTTCTTTCTTTCTT\n\
          TCTATCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTGTTTCGTCCTTTTGAGACAGAGT\n\
          TTCACTCTTGTTTCCACGGCTAGAGTGCAATGGCGCGATCTTGGCTCACCGCACCTTCCGCCTCCCGGGT\n\
          TCGAGCGCTTCTCCTGCCTCCAGCCTCCCGATTAGCGGGGATTGACAGGGAGGCACCCCCACGCCTGGCT\n\
          TGGCTGATGTTTGTGTTTTTAGTAGGCACGCCGTGTCTCTCCATGTTGCTCAGGCTGGTCTCCAACTCCC\n\
          GACCTCCTGTGATGCGCCCACCTCGGCCTCTCGAAGTGCTGGGATGACGGGCGTGACGACCGTGCCCGGC\n\
          CTGTTGACTCATTTCGCTTTTTTATTTCTTTCGTTTCCACGCGTTTACTTATATGTATTAATGTAAACGT\n\
          TTCTGTACGCTTATATGCAAACAACGACAACGTGTATCTCTGCATTGAATACTCTTGCGTATGGTAAATA\n\
          CGTATCGGTTGTATGGAAATAGACTTCTGTATGATAGATGTAGGTGTCTGTGTTATACAAATAAATACAC\n\
          ATCGCTCTATAAAGAAGGGATCGTCGATAAAGACGTTTATTTTACGTATGAAAAGCGTCGTATTTATGTG\n\
          TGTAAATGAACCGAGCGTACGTAGTTATCTCTGTTTTCTTTCTTCCTCTCCTTCGTGTTTTTCTTCCTTC\n\
          CTTTCTTCCTTTCTCTCCTTCTTTAGGTTTTTCTTCCTCTCTTCCTTTCCTTCTTTCTCTCTTTCTGTCC\n\
          TTTTTTCCTTCGTGCTTTATTTCTCTTTCGTTCCCTGTGTTTCCTTCTTTTTTCTTTCCTCTCTGTTTCT\n\
          TTTTCCCTTCTTTCCTTCGTTTCTTTCCTCATTCTTTCTCTCTTTTTCGTTGTTTCTTTCCTTCCCGTCT\n\
          GTCTTTTAAAAAATTGGAGTGTTTCAGAAGTTTACTTTGTGTATCTACGTTTTCTAAATTGTCTCTCTTT\n\
          TCTCCATTTTCTTCCTCCCTCCCTCCCTCCCTCCCTGCTCCCTTCCCTCCCTCCTTCCCTTTCGCCATCT\n\
          GTCTCTTTTCCCCACTCCCCTCCCCCCGTCTGTCTCTGCGTGGATTCCGGAAGAGCCTACCGATTCTGCC\n\
          TCTCCGTGTGTCTGCAGCGACCCCGCGACCGAGTCCTTGTGTGTTCTTTCTCCCTCCCTCCCTCCCTCCC\n\
          TCCCTCCCTCCCTCCCTGCTTCCGAGAGGCATCTCCAGAGACCGCGCCGTGGGTTGTCTTCTGACTCTGT\n\
          CGCGGTCGAGGCAGAGACGCGTTTTGGGCACCGTTTGTGTGGGGTTGGGGCAGAGGGGCTGCGTTTTCGG\n\
          CCTCGGGAAGAGCTTCTCGACTCACGGTTTCGCTTTCGCGGTCCACGGGCCGCCCTGCCAGCCGGATCTG\n\
          TCTCGCTGACGTCCGCGGCGGTTGTCGGGCTCCATCTGGCGGCCGCTTTGAGATCGTGCTCTCGGCTTCC\n\
          GGAGCTGCGGTGGCAGCTGCCGAGGGAGGGGACCGTCCCCGCTGTGAGCTAGGCAGAGCTCCGGAAAGCC\n\
          CGCGGTCGTCAGCCCGGCTGGCCCGGTGGCGCCAGAGCTGTGGCCGGTCGCTTGTGAGTCACAGCTCTGG\n\
          CGTGCAGGTTTATGTGGGGGAGAGGCTGTCGCTGCGCTTCTGGGCCCGCGGCGGGCGTGGGGCTGCCCGG\n\
          GCCGGTCGACCAGCGCGCCGTAGCTCCCGAGGCCCGAGCCGCGACCCGGCGGACCCGCCGCGCGTGGCGG\n\
          AGGCTGGGGACGCCCTTCCCGGCCCGGTCGCGGTCCGCTCATCCTGGCCGTCTGAGGCGGCGGCCGAATT\n\
          CGTTTCCGAGATCCCCGTGGGGAGCCGGGGACCGTCCCGCCCCCGTCCCCCGGGTGCCGGGGAGCGGTCC\n\
          CCGGGCCGGGCCGCGGTCCCTCTGCCGCGATCCTTTCTGGCGAGTCCCCGTGGCCAGTCGGAGAGCGCTC\n\
          CCTGAGCCGGTGCGGCCCGAGAGGTCGCGCTGGCCGGCCTTCGGTCCCTCGTGTGTCCCGGTCGTAGGAG\n\
          GGGCCGGCCGAAAATGCTTCCGGCTCCCGCTCTGGAGACACGGGCCGGCCCCTGCGTGTGGCCAGGGCGG\n\
          CCGGGAGGGCTCCCCGGCCCGGCGCTGTCCCCGCGTGTGTCCTTGGGTTGACCAGAGGGACCCCGGGCGC\n\
          TCCGTGTGTGGCTGCGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCCGTGTCGCGCGTCGCCTGGGC\n\
          CGGCGGCGTGGTCGGTGACGCGACCTCCCGGCCCCGGGGGAGGTATATCTTTCGCTCCGAGTCGGCAATT\n\
          TTGGGCCGCCGGGTTATAT"
              return default_ribo;
            } else {
              return self;
            }
          }
      index_base_name:
        source: genome
        valueFrom: $(self + "_bowtie_ribosomal")
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

  convert_annotation_to_bed:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
        - class: ResourceRequirement
          ramMin: 7620
          coresMin: 1
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
    out: [sorted_file]

  annotation_bed_to_bigbed:
    run: ../tools/ucsc-bedtobigbed.cwl
    in:
      input_bed: sort_annotation_bed/sorted_file
      chrom_length_file: star_generate_indices/chrom_length
      bed_type:
        default: "bed4+8"
      output_filename:
        source: sort_annotation_bed/sorted_file
        valueFrom: $(self.basename + ".tbi")
    out: [bigbed_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Generate genome indices for STAR & bowtie"
label: "Generate genome indices for STAR & bowtie"
s:alternateName: "Generates genome indices for STAR v2.5.3a (03/17/2017) & bowtie v1.2.0 (12/30/2016)."

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/genome-indices.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium, LLC"
  s:email: mailto:support@datirium.com
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45226"
    s:streetAddress: "3559 Kroger Ave"


# doc:
#   $include: ../descriptions/genome-indices.md


doc: |
  Creates indices for:
   * [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
   * [bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016)

  It performs the following steps:

  1. `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
  2. Outputs indices as [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
  3. Separates *chrNameLength.txt* file from Directory output
  4. `bowtie-build` to generate indices requires genome [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file as input, returns results as a group of main and secondary files
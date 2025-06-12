cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-mirdeep2:v1.0.1
inputs:
  threads:
    type: int
    inputBinding:
      prefix: -t
    doc: |
      Number of threads for parallel processing.
  bt2dir:
    type: Directory
    inputBinding:
      prefix: -b
    doc: |
      BWA index directory containing reference genome FASTA with associated indices.
  refgenome:
    type: File
    inputBinding:
      prefix: -r
    doc: |
      Reference genome FASTA file to be used for alignment.
  genome:
    type: string
    inputBinding:
      prefix: -g
    doc: |
      Name used for setting organism name, genus, species, and tax ID.
  adapter:
    type: string
    inputBinding:
      prefix: -a
    doc: |
      Adapter sequence to be trimmed from miRNA sequence reads.
  fastq:
    type: File
    inputBinding:
      prefix: -f
    doc: |
      Single-end read data in FASTQ format, received after sequencing.
outputs:
  known_novel_mir_pdfs:
    type: File
    outputBinding:
      glob: known_novel_mir_pdfs.tar.gz
    doc: |
      output directory gzip tarball for result html references
  pdfs_directory:
    type: Directory
    outputBinding:
      glob: pdfs_*
    doc: |
      output directory for column 1 hyperlinks in mirdeep2_result html
  mirdeep2_result:
    type: File
    outputBinding:
      glob: mirdeep2_result.html
    doc: |
      summary of mirdeep2 results, "miRDeep2 Results" tab
  mirs_known:
    type: File
    outputBinding:
      glob: mirs_known.tsv
    doc: |
      known mature miRNAs detected by mirdeep2
  mirs_known_bed:
    type: File
    outputBinding:
      glob: mirs_known.bed
    doc: |
      known mature miRNAs detected by mirdeep2 bed file for igv annotation
  deseq_input_isoforms:
    type: File
    outputBinding:
      glob: deseq_input.tsv
    doc: |
      novel mature miRNAs detected by mirdeep2 formatted for input into DESeq
  deseq_input_genes:
    type: File
    outputBinding:
      glob: deseq_input.tsv
    doc: |
      novel mature miRNAs detected by mirdeep2 formatted for input into DESeq
  deseq_input_common_tss:
    type: File
    outputBinding:
      glob: deseq_input.tsv
    doc: |
      novel mature miRNAs detected by mirdeep2 formatted for input into DESeq
  mirs_novel:
    type: File
    outputBinding:
      glob: mirs_novel.tsv
    doc: |
      novel mature miRNAs detected by mirdeep2
  mirs_known_exocarta_deepmirs:
    type: File
    outputBinding:
      glob: mirs_known_exocarta_deepmirs.tsv
    doc: |
      list of detected miRNA also in ExoCarta's exosome database, "Detected Exosome miRNAs" tab
  mirs_known_gene_targets:
    type: File
    outputBinding:
      glob: mirs_known_gene_targets.tsv
    doc: |
      pre-computed gene targets of known mature mirs, downloadable
  known_mirs_mature:
    type: File
    outputBinding:
      glob: known_mirs_mature.fa
    doc: |
      known mature mir sequences, downloadable
  known_mirs_precursor:
    type: File
    outputBinding:
      glob: known_mirs_precursor.fa
    doc: |
      known precursor mir sequences, downloadable
  novel_mirs_mature:
    type: File
    outputBinding:
      glob: novel_mirs_mature.fa
    doc: |
      novel precursor mir sequences, downloadable
  novel_mirs_precursor:
    type: File
    outputBinding:
      glob: novel_mirs_precursor.fa
    doc: |
      novel precursor mir sequences, downloadable
  overview:
    type: File
    outputBinding:
      glob: overview.md
    doc: |
      input list, alignment & mir metrics, "Overview" tab
  log_file_stdout:
    type: stdout
  log_file_stderr:
    type: stderr
baseCommand:
- run_mirdeep2.sh
stdout: mirna-mirdeep2_stdout.log
stderr: mirna-mirdeep2_stderr.log
doc: "A CWL tool for discovering known or novel miRNAs from deep sequencing data using the miRDeep2 tool.\nCalls the shell wrapper for the SciDAP miRDeep2 pipeline. The ExoCarta exosome database is also used\nfor identifying exosome-related miRNAs, and TargetScan's organism-specific databases are used for\nidentifying miRNA gene targets.\n\nPrimary Output files:\n- mirs_known.tsv, known mature miRNAs detected by mirdeep2, \"Known miRNAs\" tab\n- mirs_novel.tsv, known novel miRNAs detected by mirdeep2, \"Novel miRNAs\" tab\nSecondary Output files:\n- mirs_known_exocarta_deepmirs.tsv, list of detected miRNA also in ExoCarta's exosome database, \"Detected Exosome miRNAs\" tab\n- mirs_known_gene_targets.tsv, pre-computed gene targets of known mature mirs, downloadable\n- known_mirs_mature.fa, known mature mir sequences, downloadable\n- known_mirs_precursor.fa, known precursor mir sequences, downloadable\n- novel_mirs_mature.fa, novel mature mir sequences, downloadable\n- novel_mirs_precursor.fa, novel precursor mir sequences, downloadable\nReports:\n- overview.md (input list, alignment & mir metrics), \"Overview\" tab\n- mirdeep2_result.html, summary of mirdeep2 results, \"miRDeep2 Results\" tab\n\nPARAMS:\n  -h\thelp\tshow this message\n  -t  INT\t\tnumber of threads\n  -b  DIR\t\tpath to bowtie2 indices folder of genome reference\n  -r  FILE\tpath to genome reference FASTA file\n  -g  STRING\tgenome short name (hg19, hg38, mm10, rn7, dm3)\n  -a  STRING\tsequencing adapter for clipping reads (default: TCGTAT)\n  -f  FILE\tpath to sequence read FASTQ file\n\nGENOME:\nFor filtering mirbase by organism.\n  genome\t#organism   #division   #name   \t\t\t\t  #tree   \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t                                    #NCBI-taxid\n  hg19\thsa     \tHSA     \tHomo sapiens    \t\t    Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae; 9606\n  hg38\thsa     \tHSA     \tHomo sapiens    \t\t    Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae; 9606\n  mm10\tmmu     \tMMU     \tMus musculus    \t\t    Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia;   \t\t    10090\n  rn7\t\trno     \tRNO     \tRattus norvegicus   \t  Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia;   \t\t    10116\n  dm3\t\tdme     \tDME     \tDrosophila melanogaster Metazoa;Bilateria;Ecdysozoa;Arthropoda;Hexapoda;        \t\t\t\t\t\t            7227\n\nNOTES:\nFor the identification of novel miRNA candidates, the following may be used as a filtering guideline:\n1. miRDeep score > 4 (some authors use 1)\n2. not present a match with rfam\n3. should present a significant RNAfold (\"yes\")\n4. a number of mature reads > 10\n5. if applicable, novel mir must be expressed in multiple samples\n\n____________________________________________________________________________________________________\nReferences:\n- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245920\n- https://github.com/rajewsky-lab/mirdeep2\n- https://biocontainers.pro/tools/mirdeep2\n- https://www.mirbase.org/\n- http://exocarta.org/index.html\n- https://www.targetscan.org/vert_80/\n"
label: mirna-mirdeep2

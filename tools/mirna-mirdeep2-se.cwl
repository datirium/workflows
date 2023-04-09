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
      prefix: "-t"
    doc: |
      Number of threads for parallel processing.

  bt2dir:
    type: Directory
    inputBinding:
      prefix: "-b"
    doc: |
      BWA index directory containing reference genome FASTA with associated indices.

  refgenome:
    type: File
    inputBinding:
      prefix: "-r"
    doc: |
      Reference genome FASTA file to be used for alignment.

  genome:
    type: string
    inputBinding:
      prefix: "-g"
    doc: |
      Name used for setting organism name, genus, species, and tax ID.

  adapter:
    type: string
    inputBinding:
      prefix: "-a"
    doc: |
      Adapter sequence to be trimmed from miRNA sequence reads.

  fastq:
    type: File
    inputBinding:
      prefix: "-f"
    doc: |
      Single-end read data in FASTQ format, received after sequencing.


outputs:

  known_novel_mir_pdfs:
    type: File
    outputBinding:
      glob: known_novel_mir_pdfs.tar.gz
    doc: |
      output directory gzip tarball for result html references

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


baseCommand: ["run_mirdeep2.sh"]
stdout: mirna-mirdeep2_stdout.log
stderr: mirna-mirdeep2_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "mirna-mirdeep2"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/mirna-mirdeep2-se.cwl
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
  A CWL tool for discovering known or novel miRNAs from deep sequencing data using the miRDeep2 tool.
  Calls the shell wrapper for the SciDAP miRDeep2 pipeline. The ExoCarta exosome database is also used
  for identifying exosome-related miRNAs, and TargetScan's organism-specific databases are used for
  identifying miRNA gene targets.

  Primary Output files:
  - mirs_known.tsv, known mature miRNAs detected by mirdeep2, "Known miRNAs" tab
  - mirs_novel.tsv, known novel miRNAs detected by mirdeep2, "Novel miRNAs" tab
  Secondary Output files:
  - mirs_known_exocarta_deepmirs.tsv, list of detected miRNA also in ExoCarta's exosome database, "Detected Exosome miRNAs" tab
  - mirs_known_gene_targets.tsv, pre-computed gene targets of known mature mirs, downloadable
  - known_mirs_mature.fa, known mature mir sequences, downloadable
  - known_mirs_precursor.fa, known precursor mir sequences, downloadable
  - novel_mirs_mature.fa, novel mature mir sequences, downloadable
  - novel_mirs_precursor.fa, novel precursor mir sequences, downloadable
  Reports:
  - overview.md (input list, alignment & mir metrics), "Overview" tab
  - mirdeep2_result.html, summary of mirdeep2 results, "miRDeep2 Results" tab

  PARAMS:
    -h	help	show this message
    -t  INT		number of threads
    -b  DIR		path to bowtie2 indices folder of genome reference
    -r  FILE	path to genome reference FASTA file
    -g  STRING	genome short name (hg19, hg38, mm10, rn7, dm3)
    -a  STRING	sequencing adapter for clipping reads (default: TCGTAT)
    -f  FILE	path to sequence read FASTQ file

  GENOME:
  For filtering mirbase by organism.
    genome	#organism   #division   #name   				  #tree   																		                                    #NCBI-taxid
    hg19	hsa     	HSA     	Homo sapiens    		    Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae; 9606
    hg38	hsa     	HSA     	Homo sapiens    		    Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae; 9606
    mm10	mmu     	MMU     	Mus musculus    		    Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia;   		    10090
    rn7		rno     	RNO     	Rattus norvegicus   	  Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia;   		    10116
    dm3		dme     	DME     	Drosophila melanogaster Metazoa;Bilateria;Ecdysozoa;Arthropoda;Hexapoda;        						            7227

  NOTES:
  For the identification of novel miRNA candidates, the following may be used as a filtering guideline:
  1. miRDeep score > 4 (some authors use 1)
  2. not present a match with rfam
  3. should present a significant RNAfold ("yes")
  4. a number of mature reads > 10
  5. if applicable, novel mir must be expressed in multiple samples

  ____________________________________________________________________________________________________
  References:
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245920
  - https://github.com/rajewsky-lab/mirdeep2
  - https://biocontainers.pro/tools/mirdeep2
  - https://www.mirbase.org/
  - http://exocarta.org/index.html
  - https://www.targetscan.org/vert_80/

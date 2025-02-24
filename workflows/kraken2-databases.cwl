cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


inputs:

  alias:
    type: string
    label: "Sample short name/Alias"
    sd:preview:
      position: 1

  database_name:
    type:
    - "null"
    - type: enum
      name: "name of Kraken2 database to download"
      symbols:
      - Viral
      - Standard
      - Standard-16
      - MinusB
      - PlusPFP-16
      - EuPathDB46
      - 16S_Greengenes
      - 16S_Silva_138
    label: "Select Kraken2 database for download:"
    'sd:localLabel': true
    doc: "Database details:\n
       - Viral (500 MB), all refseq viral genomes\n
       - Standard (80 GB), Refeq archaea, bacteria, viral, plasmid, human, UniVec_Core
       - Standard-16 (15 GB), Refeq archaea, bacteria, viral, plasmid, human, UniVec_Core, capped at 16 GB (shrunk via random kmer downselect)
       - MinusB (9 GB), Standard minus bacteria\n
       - PlusPFP-16 (15 GB), Standard + protozoa, fungi & plant, capped at 16 GB (shrunk via random kmer downselect)\n
       - EuPathDB46 (11 GB), eukaryotic pathogen genomes with contaminants removed (https://veupathdb.org/veupathdb/app)\n
       - 16S_Greengenes (73 MB), Greengenes 16S rRNA database (release 13.5, 20200326)\n
       - 16S_Silva_138 (112 MB), SILVA 16S rRNA database (release 138, 20200326)"
    sd:preview:
      position: 2


outputs:

  k2db:
    type: Directory
    label: "decompressed and untarred kraken2 database directory used as input for kraken2 classification"
    outputSource: download_db/k2db

  compressed_k2db_tar:
    type: File
    label: "compressed and tarred kraken2 database directory file for download and use outside of scidap"
    outputSource: download_db/compressed_k2db_tarball


steps:

  download_db:
    run: ../tools/k2-download-db.cwl
    in:
      user_selection:
        source: database_name
        valueFrom: $(self)
    out:
      - k2db
      - compressed_k2db_tarball


label: "Kraken2 Database installation pipeline"
doc: |
  This workflow downloads the user-selected pre-built kraken2 database from: https://benlangmead.github.io/aws-indexes/k2

  ### __Inputs__
  Select a pre-built Kraken2 database to download and use for metagenomic classification:
   - Available options comprised of various combinations of RefSeq reference genome sets:
     - [Viral (0.5 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz), all refseq viral genomes
     - [MinusB (8.7 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20221209.tar.gz), standard minus bacteria (archaea, viral, plasmid, human1, UniVec_Core)
     - [PlusPFP-16 (15.0 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz), standard (archaea, bacteria, viral, plasmid, human1, UniVec_Core) + (protozoa, fungi & plant) capped at 16 GB (shrunk via random kmer downselect)
     - [EuPathDB46 (34.1 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20201113.tar.gz), eukaryotic pathogen genomes with contaminants removed (https://veupathdb.org/veupathdb/app)
     - [16S_gg_13_5 (73 MB)](https://genome-idx.s3.amazonaws.com/kraken/16S_Greengenes13.5_20200326.tgz), Greengenes 16S rRNA database ([release 13.5](https://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/gg_13_5/), 20200326)\n
     - [16S_silva_138 (112 MB)](https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz), SILVA 16S rRNA database ([release 138.1](https://www.arb-silva.de/documentation/release-1381/), 20200827)

  ### __Outputs__
   - k2db, an upstream database used by kraken2 classification tool
   - compressed_k2db_tar, compressed and tarred kraken2 database directory file for download and use outside of scidap

  ### __Data Analysis Steps__
  1. download selected pre-built kraken2 database.
  2. make available as upstream source for kraken2 metagenomic taxonomic classification.

  ### __References__
    - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0

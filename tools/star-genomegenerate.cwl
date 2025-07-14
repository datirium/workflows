cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 30510
  coresMin: 8
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {
          "class": "Directory",
          "basename": inputs.genome_dir,
          "listing": [],
          "writable": true}
      ]
    }
hints:
- class: DockerRequirement
  dockerPull: scidap/star:v2.7.5c
inputs:
  parameters_files:
    type: string?
    inputBinding:
      position: 1
      prefix: --parametersFiles
    doc: |
      string: name of a user-defined parameters file, "-": none. Can only be
      defined on the command line.
  run_dir_perm:
    type: string?
    inputBinding:
      position: 1
      prefix: --runDirPerm
    doc: |
      User_RWX
      string: permissions for the directories created at the run-time.
      User_RWX ... user-read/write/execute
      All_RWX  ... all-read/write/execute (same as chmod 777)
  sjdb_gtf_feature_exon:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFfeatureExon
    doc: |
      exon
      string: feature type in GTF file to be used as exons for building
      transcripts
  genome_dir:
    type: string
    default: star_indices
    inputBinding:
      position: 1
      prefix: --genomeDir
    doc: Directory where genome index files are stored
  genome_fasta_files:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 1
      itemSeparator: ' '
      prefix: --genomeFastaFiles
    doc: |
      string(s): path(s) to the fasta files with genomic sequences for genome
      generation, separated by spaces. Only used if runMode==genomeGenerate.
      These files should be plain text FASTA files, they *cannot* be zipped.
  limit_io_buffer_size:
    type: long?
    inputBinding:
      position: 1
      prefix: --limitIObufferSize
    doc: |
      150000000
      int>0: max available buffers size (bytes) for input/output, per thread
  sjdb_gtf_chr_prefix:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFchrPrefix
    doc: |
      string: prefix for chromosome names in a GTF file (e.g. 'chr' for using ENSMEBL annotations with UCSC genomes)
  genome_chr_bin_n_bits:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeChrBinNbits
    doc: |
      int: =log2(chrBin), where chrBin is the size of the bins for genome
      storage: each chromosome will occupy an integer number of bins
  genome_sa_sparse_d:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeSAsparseD
    doc: |
      int>0: suffux array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAM at the cost of mapping speed reduction
  genome_suffix_length_max:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeSuffixLengthMax
    doc: |
      int: maximum length of the suffixes, has to be longer than read length. -1 = infinite.
  genome_chain_files:
    type:
    - 'null'
    - File
    - type: array
      items: File
    inputBinding:
      position: 1
      prefix: --genomeChainFiles
    doc: |
      string: chain files for genomic liftover
  limit_genome_generate_ram:
    type: long?
    inputBinding:
      position: 1
      prefix: --limitGenomeGenerateRAM
    doc: |
      31000000000
      int>0: maximum available RAM (bytes) for genome generation
  sjdb_gtf_tag_exon_parent_transcript:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFtagExonParentTranscript
    doc: |
      transcript_id
      string: tag name to be used as exons' transcript-parents (default
      "transcript_id" works for GTF files)
  sjdb_gtf_tag_exon_parent_gene:
    type: string?
    inputBinding:
      position: 1
      prefix: --sjdbGTFtagExonParentGene
    doc: |
      gene_id
      string: tag name to be used as exons'' gene-parents (default "gene_id" works
      for GTF files)
  sjdb_gtf_file:
    type: File?
    inputBinding:
      position: 1
      prefix: --sjdbGTFfile
    doc: |
      string: path to the GTF file with annotations
  sys_shell:
    type: string?
    inputBinding:
      position: 1
      prefix: --sysShell
    doc: |
      string: path to the shell binary, preferrably bash, e.g. /bin/bash.
      - ... the default shell is executed, typically /bin/sh. This was reported to fail on some Ubuntu systems - then you need to specify path to bash.
  sjdb_overhang:
    type: int?
    inputBinding:
      position: 1
      prefix: --sjdbOverhang
    doc: |
      default: 100
      int>0: length of the donor/acceptor sequence on each side of the junctions,
      ideally = (mate_length - 1)
  genome_sa_index_n_bases:
    type: int?
    inputBinding:
      position: 1
      prefix: --genomeSAindexNbases
    doc: |
      int: length (bases) of the SA pre-indexing string. Typically between 10 and
      15. Longer strings will use much more memory, but allow faster searches.
outputs:
  indices_folder:
    type: Directory
    outputBinding:
      glob: $(inputs.genome_dir)
  chrom_length:
    type: File
    outputBinding:
      glob: $(inputs.genome_dir + "/chrNameLength.txt")
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- STAR
- --runMode
- genomeGenerate
- --runThreadN 8
stdout: star_build_stdout.log
stderr: star_build_stderr.log
doc: |
  Tool returns directory with indices generated by STAR. If genome_dir input is not provided,
  use default output directory name star_indices.
  Output chr_name_length should not be moved outside the indices folder.
label: star-genomegenerate

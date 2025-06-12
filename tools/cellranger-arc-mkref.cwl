cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_folder_name = function() { if (inputs.output_folder_name == ""){ var root = inputs.genome_fasta_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.genome_fasta_file.basename:root; } else { return inputs.output_folder_name; } };
- class: InitialWorkDirRequirement
  listing: |
    ${
      var exclude_chr = "[]";
      if (inputs.exclude_chr && inputs.exclude_chr.length > 0){
        exclude_chr = '["' + inputs.exclude_chr.join('", "') + '"]'
      }
      var entry = `
      {
          genome: ["${get_output_folder_name()}"]
          input_fasta: ["${inputs.genome_fasta_file.path}"]
          input_gtf: ["${inputs.annotation_gtf_file.path}"]
          non_nuclear_contigs: ${exclude_chr}
      }`
      return [{
        "entry": entry,
        "entryname": "config.txt"
      }];
    }
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/cellranger-arc:v0.0.1
inputs:
  genome_fasta_file:
    type: File
    doc: |
      Genome FASTA file
  annotation_gtf_file:
    type: File
    doc: |
      GTF annotation file.
  exclude_chr:
    type:
    - 'null'
    - string[]
    doc: |
      Contigs that do not have any chromatin structure,
      for example, mitochondria or plastids. These
      contigs are excluded from peak calling since the
      entire contig will be "open" due to a lack of
      chromatin structure.
  output_folder_name:
    type: string?
    default: ''
    doc: |
      Unique genome name, used to name output folder
  threads:
    type: int?
    inputBinding:
      position: 5
      prefix: --nthreads
    doc: |
      Number of threads used during STAR genome indexing
      Default: 1
  memory_limit:
    type: int?
    inputBinding:
      position: 6
      prefix: --memgb
    doc: |
      Maximum memory (GB) used when aligning reads with STAR
      Defaults: 16
outputs:
  indices_folder:
    type: Directory
    outputBinding:
      glob: $(get_output_folder_name())
    doc: |
      Compatible with Cell Ranger ARC reference
      folder that includes STAR and BWA indices.
  chrom_length_file:
    type: File
    outputBinding:
      glob: $(get_output_folder_name() + "/star/chrNameLength.txt")
    doc: |
      Chromosome length file
      in TSV format.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- cellranger-arc
- mkref
- --config
- config.txt
stdout: cellranger_arc_mkref_stdout.log
stderr: cellranger_arc_mkref_stderr.log
label: Cell Ranger Reference (RNA+ATAC)
doc: |
  Cell Ranger Reference (RNA, ATAC, RNA+ATAC)

  Builds a reference genome of a selected species for quantifying
  gene expression and chromatin accessibility

  Notes:
  - `input_motifs` parameter in the `config.txt` file is not
    implemented.
  - if GTF file provided in `annotation_gtf_file` has records with
    duplicate gene_id, they should be grouped together. Applicable to
    USCS RefGene annotations.

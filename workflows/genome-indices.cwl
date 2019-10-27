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

  fasta:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes all chromosomes"

  fasta_ribosomal:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Ribosomal DNA FASTA file"
    doc: "Ribosomal DNA FASTA file"

  fasta_mitochondrial:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Mitochondrial DNA FASTA file"
    doc: "Mitochondrial DNA FASTA file"

  effective_genome_size:
    type: string
    label: "Effective genome size"
    doc: "MACS2 effective genome sizes: hs, mm, ce, dm or number, for example 2.7e9"

  input_annotation_gtf:
    type: File
    format: "http://edamontology.org/format_2306"
    label: "GTF annotation file"
    doc: "GTF annotation file. Includes reference genome and mitochondrial DNA annotations"

  annotation_tab:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "TSV annotation file"
    doc: "Tab-separated annotation file. Includes reference genome and mitochondrial DNA annotations"

  cytoband:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "CytoBand file for IGV browser"
    doc: "Tab-separated cytoBand file for IGV browser"
  
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

  threads:
    type: int?
    label: "Number of threads to run tools"
    doc: "Number of threads for those steps that support multithreading"
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
    outputSource: input_annotation_gtf
    label: "GTF annotation file"
    doc: "GTF annotation file. Includes reference genome and mitochondrial DNA annotations"
    
  annotation:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: annotation_tab
    label: "TSV annotation file"
    doc: "Tab-separated annotation file. Includes reference genome and mitochondrial DNA annotations"
  
  fasta_output:
    type: File
    format: "http://edamontology.org/format_1929"
    outputSource: fasta
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes all chromosomes"

  fasta_fai_output:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: index_fasta/fai_file
    label: "FAI index for genome FASTA file"
    doc: "Tab-separated FAI index file"

  cytoband_output:
    type: File
    format: "http://edamontology.org/format_3475"
    outputSource: cytoband
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

  star_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: fasta
      sjdb_gtf_file: input_annotation_gtf
      genome_sa_sparse_d: genome_sa_sparse_d
      limit_genome_generate_ram: limit_genome_generate_ram
      genome_sa_index_n_bases: genome_sa_index_n_bases
      genome_chr_bin_n_bits: genome_chr_bin_n_bits
      genome_dir:
        source: genome
        valueFrom: $(self + "_star_genome")
      threads: threads
    out:
    - indices_folder
    - chrom_length
    - stdout_log
    - stderr_log

  mitochondrial_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: fasta_mitochondrial
      sjdb_gtf_file: input_annotation_gtf
      genome_sa_sparse_d: genome_sa_sparse_d
      genome_sa_index_n_bases: genome_sa_index_n_bases_mitochondrial
      genome_dir:
        source: genome
        valueFrom: $(self + "_star_mitochondrial")
      threads: threads
    out:
    - indices_folder
    - chrom_length
    - stdout_log
    - stderr_log

  bowtie_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta
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
      fasta_file: fasta_ribosomal
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
      fasta_file: fasta
    out:
    - fai_file

  convert_annotation_to_bed:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
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
      annotation_tsv_file: annotation_tab
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
- http://schema.org/docs/schema_org_rdfa.html

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

doc: |
  Creates indices for:
   * [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
   * [bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016)

  It performs the following steps:

  1. `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
  2. Outputs indices as [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
  3. Separates *chrNameLength.txt* file from Directory output
  4. `bowtie-build` to generate indices requires genome [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file as input, returns results as a group of main and secondary files

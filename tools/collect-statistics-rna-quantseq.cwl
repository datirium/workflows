cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function() {
        if (inputs.output_prefix) {
          return inputs.output_prefix;
        }
        var root = inputs.star_alignment_report.basename.split('.').slice(0,-1).join('.');
        var suffix = "_collected_statistics_report";
        return (root == "")?inputs.star_alignment_report.basename+suffix:root+suffix;
    };


hints:
- class: DockerRequirement
  dockerPull: rackspacedot/python37


inputs:

  script:
    type: string?
    default: |
        #!/usr/bin/env python
        import os
        import sys
        import argparse
        import yaml
        import math
        import re


        def cut_int(s):
            return int(s.strip().split()[0])

        def cut_float(s):
            return float(s.strip().split()[0])

        def cut_percent(s):
            return float(s.strip().replace("%", "").split()[0])


        ADAPTER_TRIMMING = {
            "Trimming mode": {
                "alias": "trimming mode",
                "function": str,
                "pair_end_specific": False
            },
            "Quality Phred score cutoff": {
                "alias": "quality phred score cutoff",
                "function": int,
                "pair_end_specific": False
            },
            "Quality encoding type selected": {
                "alias": "quality encoding type",
                "function": str,
                "pair_end_specific": False
            },
            "Adapter sequence": {
                "alias": "adapter sequence",
                "function": str,
                "pair_end_specific": False
            },
            "Minimum required adapter overlap": {
                "alias": "minimum required adapter overlap",
                "function": cut_int,
                "pair_end_specific": False
            },
            "Maximum trimming error rate": {
                "alias": "maximum trimming error rate",
                "function": cut_float,
                "pair_end_specific": False
            },
            "Minimum required sequence length": {
                "alias": "minimum required read length",
                "function": cut_int,
                "pair_end_specific": False
            },
            "Total number of sequences analysed": {
                "alias": "number of reads/pairs analysed for length validation",
                "function": int,
                "pair_end_specific": False
            },
            "Number of sequence pairs removed": {
                "alias": "reads/pairs removed because of length cutoff",
                "function": cut_int,
                "pair_end_specific": False
            },
            "Sequences removed because they became shorter": {
                "alias": "reads/pairs removed because of length cutoff",
                "function": cut_int,
                "pair_end_specific": False
            },
            "order": ["fastq",
                    "trimming mode",
                    "adapter sequence",
                    "number of reads/pairs analysed for length validation",
                    "reads/pairs removed because of length cutoff",
                    "minimum required read length",
                    "quality phred score cutoff",
                    "quality encoding type",
                    "minimum required adapter overlap",
                    "maximum trimming error rate"]
        }


        ALIGNMENT = {
            "Number of input reads": {
                "alias": "total reads/pairs processed",
                "function": int,
                "pair_end_specific": False
            },
            "Uniquely mapped reads number": {
                "alias": "uniquely mapped reads/pairs number",
                "function": int,
                "pair_end_specific": False
            },
            "Mismatch rate per base, %": {
                "alias": "mismatch rate per base, %",
                "function": cut_percent,
                "pair_end_specific": False
            },
            "Deletion rate per base": {
                "alias": "deletion rate per base, %",
                "function": cut_percent,
                "pair_end_specific": False
            },
            "Number of reads mapped to multiple loci": {
                "alias": "reads/pairs mapped to multiple loci",
                "function": int,
                "pair_end_specific": False
            },
            "Number of reads mapped to too many loci": {
                "alias": "reads/pairs suppressed due to mapping to too many loci",
                "function": int,
                "pair_end_specific": False
            },
            "reads unmapped: too many mismatches": {
                "alias": "reads/pairs unmapped due to too many mismatches, %",
                "function": cut_percent,
                "pair_end_specific": False
            },
            "reads unmapped: too short": {
                "alias": "reads/pairs unmapped due to too short, %",
                "function": cut_percent,
                "pair_end_specific": False
            },
            "reads unmapped: other": {
                "alias": "reads/pairs unmapped due to other reasons, %",
                "function": cut_percent,
                "pair_end_specific": False
            },
            "reads with at least one reported alignment": {
                "alias": "number of reads/pairs in ribosomal dna",
                "function": cut_int,
                "pair_end_specific": False
            },
            "order": ["total reads/pairs processed",
                    "uniquely mapped reads/pairs number",
                    "number of reads/pairs in transcriptome",
                    "number of reads/pairs in ribosomal dna",
                    "reads/pairs mapped to multiple loci",
                    "reads/pairs suppressed due to mapping to too many loci",
                    "reads/pairs unmapped due to too many mismatches, %",
                    "reads/pairs unmapped due to too short, %",
                    "reads/pairs unmapped due to other reasons, %",
                    "mismatch rate per base, %",
                    "deletion rate per base, %"]
        }


        BAMSTATS = {
            "raw total sequences": {
                "alias": "total reads/pairs",
                "function": int,
                "pair_end_specific": True
            },
            "reads mapped": {
                "alias": "reads/pairs mapped",
                "function": int,
                "pair_end_specific": True
            },
            "reads unmapped": {
                "alias": "reads/pairs unmapped",
                "function": int,
                "pair_end_specific": True
            },
            "average length": {
                "alias": "reads average length",
                "function": float,
                "pair_end_specific": False
            },
            "maximum length": {
                "alias": "reads maximum length",
                "function": int,
                "pair_end_specific": False
            },
            "average quality": {
                "alias": "reads average quality",
                "function": float,
                "pair_end_specific": False
            },
            "insert size average": {
                "alias": "insert size average",
                "function": float,
                "pair_end_specific": False
            },
            "insert size standard deviation": {
                "alias": "insert size standard deviation",
                "function": float,
                "pair_end_specific": False
            },
            "order": ["total reads/pairs",
                    "reads/pairs mapped",
                    "reads/pairs unmapped",
                    "reads average length",
                    "reads maximum length",
                    "reads average quality",
                    "insert size average",
                    "insert size standard deviation"]
        }


        def arg_parser():
            general_parser = argparse.ArgumentParser()
            general_parser.add_argument("--trim1",           help="Path to Trimgalore report file for FASTQ 1",           required=False)
            general_parser.add_argument("--trim2",           help="Path to Trimgalore report file for FASTQ 2",           required=False)
            general_parser.add_argument("--star",            help="Path to STAR report file",                             required=False)
            general_parser.add_argument("--bowtie",          help="Path to Bowtie report file",                           required=False)
            general_parser.add_argument("--bamstats",        help="Path to bam statistics report file",                   required=False)
            general_parser.add_argument("--isoforms",        help="Path to isoforms file",                                required=False)
            general_parser.add_argument("--paired",          help="Process as paired-end",                           action="store_true")
            general_parser.add_argument("--output",          help="Output filename prefix",               default="collected_statistics")
            return general_parser


        def normalize_args(args, skip_list=[]):
            normalized_args = {}
            for key,value in args.__dict__.items():
                if key not in skip_list:
                    normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
                else:
                    normalized_args[key]=value
            return argparse.Namespace (**normalized_args)


        def open_file(filename):
            lines = []
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.strip():
                        lines.append(line.strip())
            return lines


        def split_line(l, delimiter):
            return [i.strip() for i in l.split(delimiter)]


        def get_correspondent_key(data, long_k):
            for short_k, v in data.items():
                if short_k in long_k:
                    return v["alias"], v["function"], v["pair_end_specific"]
            raise Exception


        def process_trimgalore_report(filepath, collected_results, header, key_dict, pair_end, delimiter):
            if not collected_results.get(header, None):
                collected_results[header] = {"fastq": []}
            fastq = {}
            for line in open_file(filepath):
                try:
                    key, value = split_line(line, delimiter)
                    if "Total reads processed" in line:
                        fastq["total reads processed"] = int(value.strip().replace(",",""))
                    elif "Reads with adapters" in line:
                        fastq["reads with adapters"] = int(value.strip().replace(",","").split()[0])
                    else:
                        res_key, res_function, pair_end_specific = get_correspondent_key(key_dict, key)
                        if not collected_results[header].get(res_key, None):
                            collected_results[header][res_key] = res_function(value)
                except Exception:
                    pass
            if not pair_end:
                collected_results[header]["number of reads/pairs analysed for length validation"] = fastq["total reads processed"]
            collected_results[header]["fastq"].append(fastq)
            if key_dict.get("order", None):
                collected_results[header] = {k: collected_results[header][k] for k in key_dict["order"] if k in collected_results[header]}


        def process_custom_report(filepath, collected_results, header, key_dict, pair_end, delimiter):
            if not collected_results.get(header, None):
                collected_results[header] = {}
            for line in open_file(filepath):
                try:
                    key, value = split_line(line, delimiter)
                    res_key, res_function, pair_end_specific = get_correspondent_key(key_dict, key)
                    if not collected_results[header].get(res_key, None):
                        if pair_end_specific and pair_end:
                            collected_results[header][res_key] = res_function(res_function(value)/2)
                        else:
                            collected_results[header][res_key] = res_function(value)
                except Exception:
                    pass
            if key_dict.get("order", None):
                collected_results[header] = {k: collected_results[header][k] for k in key_dict["order"] if k in collected_results[header]}


        def process_isoforms_report(filepath, collected_results, header, key_dict, pair_end):
            if not collected_results.get(header, None):
                collected_results[header] = {}
            total_reads_index = None
            used_reads = 0
            for line in open_file(filepath):
                if re.match ('.*RefseqId.*|.*GeneId.*|.*Chrom.*|.*TotalReads.*', line) and total_reads_index is None:
                    total_reads_index = line.split(',').index('TotalReads')
                    continue
                line_splitted = line.split(',')
                used_reads += int(line_splitted[total_reads_index])
            if pair_end:
                used_reads = int(used_reads/2)
            collected_results[header]["number of reads/pairs in transcriptome"] = used_reads
            if key_dict.get("order", None):
                collected_results[header] = {k: collected_results[header][k] for k in key_dict["order"] if k in collected_results[header]}


        def collect_stats(args):
            collected_results = {}

            if args.trim1:
                process_trimgalore_report(filepath=args.trim1,
                                        collected_results=collected_results,
                                        header="adapter trimming statistics",
                                        key_dict=ADAPTER_TRIMMING,
                                        pair_end=args.paired,
                                        delimiter=":")

            if args.trim2:
                process_trimgalore_report(filepath=args.trim2,
                                        collected_results=collected_results,
                                        header="adapter trimming statistics",
                                        key_dict=ADAPTER_TRIMMING,
                                        pair_end=args.paired,
                                        delimiter=":")

            if args.star:
                process_custom_report(filepath=args.star,
                                    collected_results=collected_results,
                                    header="alignment statistics",
                                    key_dict=ALIGNMENT,
                                    pair_end=args.paired,
                                    delimiter="|")

            if args.isoforms:
                process_isoforms_report(filepath=args.isoforms,
                                        collected_results=collected_results,
                                        header="alignment statistics",
                                        key_dict=ALIGNMENT,
                                        pair_end=args.paired)

            if args.bowtie:
                process_custom_report(filepath=args.bowtie,
                                    collected_results=collected_results,
                                    header="alignment statistics",
                                    key_dict=ALIGNMENT,
                                    pair_end=args.paired,
                                    delimiter=":")

            if args.bamstats:
                process_custom_report(filepath=args.bamstats,
                                    collected_results=collected_results,
                                    header="BAM statistics",
                                    key_dict=BAMSTATS,
                                    pair_end=args.paired,
                                    delimiter=":")

            return (collected_results)


        def export_results_yaml(collected_data, filepath):
            with open(filepath+".yaml", 'w') as output_stream:
                output_stream.write(yaml.dump(collected_data, width=1000, sort_keys=False))


        def export_results_markdown(collected_data, filepath):
            with open(filepath+".md", 'w') as output_stream:
                for line in yaml.dump(collected_data, width=1000, sort_keys=False).split("\n"):
                    if not line.strip():
                        continue
                    if line.startswith("  - "):
                        output_stream.write(line+"\n")
                    elif line.startswith("    "):
                        output_stream.write("<br>"+line+"\n")
                    elif line.startswith("  "):
                        output_stream.write("- "+line+"\n")
                    else:
                        output_stream.write("### "+line+"\n")


        def export_results_table(collected_data, filepath):
            with open(filepath+".tsv", 'w') as output_stream:
                total_reads= collected_data["alignment statistics"]["total reads/pairs processed"]
                uniquely_mapped_reads = collected_data["alignment statistics"]["uniquely mapped reads/pairs number"]
                uniquely_mapped_reads_dedup = collected_data["BAM statistics"]["reads/pairs mapped"]
                transcriptome_reads = collected_data["alignment statistics"]["number of reads/pairs in transcriptome"]
                multimapped_reads = collected_data["alignment statistics"]["reads/pairs suppressed due to mapping to too many loci"]
                not_transcriptome_reads = uniquely_mapped_reads_dedup - transcriptome_reads
                unmapped_reads = total_reads - uniquely_mapped_reads - multimapped_reads
                ribosomal_reads = collected_data["alignment statistics"]["number of reads/pairs in ribosomal dna"]
                reads_dedup = uniquely_mapped_reads - uniquely_mapped_reads_dedup

                header = [
                            "Tags total",
                            "Transcriptome",
                            "Multi-mapped",
                            "Outside annotation",
                            "Unmapped",
                            "Deduplicated",
                            "Ribosomal contamination",

                            "alignment statistics",
                            "total reads/pairs processed",
                            "uniquely mapped reads/pairs number",
                            "reads/pairs mapped to multiple loci",
                            "reads/pairs suppressed due to mapping to too many loci",
                            "mismatch rate per base, %",
                            "deletion rate per base, %",
                            "reads/pairs unmapped due to too many mismatches, %",
                            "reads/pairs unmapped due to too short, %",
                            "reads/pairs unmapped due to other reasons, %",
                            "number of reads/pairs in transcriptome",
                            "number of reads/pairs in ribosomal dna",

                            "BAM statistics",
                            "total reads/pairs",
                            "reads/pairs mapped",
                            "reads/pairs unmapped",
                            "insert size average",
                            "insert size standard deviation",
                            "reads average length",
                            "reads average quality",
                            "reads maximum length"]

                if collected_data.get("adapter trimming statistics", None):
                    header.extend(["adapter trimming statistics",
                                "trimming mode",
                                "adapter sequence",
                                "quality phred score cutoff",
                                "minimum required adapter overlap",
                                "maximum trimming error rate",
                                "minimum required read length",
                                "number of reads/pairs analysed for length validation",
                                "reads/pairs removed because of length cutoff",
                                "fastq1: total reads processed",
                                "fastq1: reads with adapters"])
                    if len(collected_data["adapter trimming statistics"]["fastq"]) > 1:
                        header.extend(["fastq2: total reads processed", "fastq2: reads with adapters"])

                output_stream.write("\t".join(header)+"\n")

                data = [
                        total_reads,
                        transcriptome_reads,
                        multimapped_reads,
                        not_transcriptome_reads,
                        unmapped_reads,
                        reads_dedup,
                        ribosomal_reads,

                        "",
                        collected_data["alignment statistics"]["total reads/pairs processed"],
                        collected_data["alignment statistics"]["uniquely mapped reads/pairs number"],
                        collected_data["alignment statistics"]["reads/pairs mapped to multiple loci"],
                        collected_data["alignment statistics"]["reads/pairs suppressed due to mapping to too many loci"],
                        collected_data["alignment statistics"]["mismatch rate per base, %"],
                        collected_data["alignment statistics"]["deletion rate per base, %"],
                        collected_data["alignment statistics"]["reads/pairs unmapped due to too many mismatches, %"],
                        collected_data["alignment statistics"]["reads/pairs unmapped due to too short, %"],
                        collected_data["alignment statistics"]["reads/pairs unmapped due to other reasons, %"],
                        collected_data["alignment statistics"]["number of reads/pairs in transcriptome"],
                        collected_data["alignment statistics"]["number of reads/pairs in ribosomal dna"],

                        "",
                        collected_data["BAM statistics"]["total reads/pairs"],
                        collected_data["BAM statistics"]["reads/pairs mapped"],
                        collected_data["BAM statistics"]["reads/pairs unmapped"],
                        collected_data["BAM statistics"]["insert size average"],
                        collected_data["BAM statistics"]["insert size standard deviation"],
                        collected_data["BAM statistics"]["reads average length"],
                        collected_data["BAM statistics"]["reads average quality"],
                        collected_data["BAM statistics"]["reads maximum length"]]

                if collected_data.get("adapter trimming statistics", None):
                    data.extend(["",
                                collected_data["adapter trimming statistics"]["trimming mode"],
                                collected_data["adapter trimming statistics"]["adapter sequence"],
                                collected_data["adapter trimming statistics"]["quality phred score cutoff"],
                                collected_data["adapter trimming statistics"]["minimum required adapter overlap"],
                                collected_data["adapter trimming statistics"]["maximum trimming error rate"],
                                collected_data["adapter trimming statistics"]["minimum required read length"],
                                collected_data["adapter trimming statistics"]["number of reads/pairs analysed for length validation"],
                                collected_data["adapter trimming statistics"]["reads/pairs removed because of length cutoff"],
                                collected_data["adapter trimming statistics"]["fastq"][0]["total reads processed"],
                                collected_data["adapter trimming statistics"]["fastq"][0]["reads with adapters"]])
                    if len(collected_data["adapter trimming statistics"]["fastq"]) > 1:
                        data.extend([collected_data["adapter trimming statistics"]["fastq"][1]["total reads processed"],
                                    collected_data["adapter trimming statistics"]["fastq"][1]["reads with adapters"]])

                output_stream.write("\t".join(str(l) for l in data))


        def main(argsl=None):
            if argsl is None:
                argsl = sys.argv[1:]
            args,_ = arg_parser().parse_known_args(argsl)
            args = normalize_args(args, skip_list=["paired"])
            collected_data = collect_stats(args)
            export_results_yaml(collected_data, args.output)
            export_results_table(collected_data, args.output)
            export_results_markdown(collected_data, args.output)


        if __name__ == "__main__":
            sys.exit(main(sys.argv[1:]))
    inputBinding:
      position: 5

  trimgalore_report_fastq_1:
    type: File?
    inputBinding:
      position: 6
      prefix: "--trim1"

  trimgalore_report_fastq_2:
    type: File?
    inputBinding:
      position: 7
      prefix: "--trim2"

  star_alignment_report:
    type: File?
    inputBinding:
      position: 8
      prefix: "--star"

  bowtie_alignment_report:
    type: File?
    inputBinding:
      position: 9
      prefix: "--bowtie"

  bam_statistics_report:
    type: File?
    inputBinding:
      position: 10
      prefix: "--bamstats"

  isoforms_file:
    type: File?
    inputBinding:
      position: 11
      prefix: "--isoforms"

  paired_end:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "--paired"

  output_prefix:
    type: string?
    inputBinding:
      position: 13
      prefix: "--output"
      valueFrom: $(get_output_prefix())
    default: ""


outputs:

  collected_statistics_yaml:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".yaml")

  collected_statistics_tsv:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".tsv")

  collected_statistics_md:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".md")


baseCommand: [python, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "collect-statistics-rna-seq"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/collect-statistics-rna-seq.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool processes and combines log files generated by Trimgalore, Bowtie, Samtools and MACS2.

  `get_output_prefix` function returns output file prefix equal to `output_prefix`+`_collected_statistics_report` (if this input is provided) or
  generated on the base of bowtie log basename with `_collected_statistics_report` extension.


s:about: |
  Runs python code from the script input

cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement


inputs:

    script:
        type: string?
        default: |
            #!/bin/bash
            if [[ "$0" -eq 0 ]]; then
                printf "1" > scaling_factor.txt
            else
                printf "$0" | awk '{printf("%.3f",10000/$0)}' > scaling_factor.txt
            fi
        inputBinding:
            position: 4

    reads_mapped:
        type: int
        label: "spike-in mapped reads from get_spikein_bam_statistics"
        inputBinding:
            position: 5


outputs:

    scaling_factor:
        type: float?
        outputBinding:
            loadContents: true
            glob: scaling_factor.txt
            outputEval: |
                ${
                    var s = self[0].contents;
                    return parseInt(s);
                }
        doc: |
            Scaling factor (sf) for seq library normalization.
            sf=C/mapped reads where C is a constant (10000 used here)

    scaling_factor_file:
        type: File
        outputBinding:
            glob: "scaling_factor.txt"
        doc: |
            Text file output of scaling factor. If contents is null, this
            means zero reads were mapped to the spike-in reference genome.


baseCommand: ["bash", "-c"]


doc: |
  Takes as input the *_bam_statistics*.report.txt from `samtools_stat.cwl` tool.
  Returns a float scale factor for normalizing a bam-to-bedgraph, using the
  `bam-bedgraph-bigwig.cwl` tool.

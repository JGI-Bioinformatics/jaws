workflow jgi_dap_leo {

    ### from inputs.json:

    Array[File] expt_raw_fastqs
    Array[File] ctl_raw_fastqs 
    File adapters
    File genome_fasta
    Array[File] bt2index_list
    Int effgsize
    File genes_gff
    File bgmodel
    Map[File, String] library_names_map
    Map[File, String] sample_names_map
    File outdir
    String amplified
    Int maxfrags
    Boolean find_motifs


### Process negative controls in parallel, then merge bams

    if(length(ctl_raw_fastqs) > 0) {
        scatter (raw_fastq in ctl_raw_fastqs) {

            String ctl_basename = sample_names_map[raw_fastq] + "_" + library_names_map[raw_fastq]

            call trimAlign as trimAlign_ctl {
                input:  raw_fastq=raw_fastq,
                        basename=ctl_basename,
                        adapters=adapters,
                        effgsize=effgsize,
                        bt2index_list=bt2index_list,
                        maxfrags=maxfrags
            }
            
            call copyOutput as copyOutput_ctl {
                input:  basename=ctl_basename,
                        genes_gff=genes_gff,
                        effgsize=effgsize,
                        genome_fasta=genome_fasta,
                        bam=trimAlign_ctl.bam,
                        bai=trimAlign_ctl.bai,
                        bigwig=trimAlign_ctl.bw,
                        trim_stats=trimAlign_ctl.trim_stats,
                        align_stats=trimAlign_ctl.align_stats,
                        destination=outdir+"/"+ctl_basename
            }
        }

        call mergeBams {
            input:  bams=trimAlign_ctl.bam
        }
    }

### Process experimental samples in parallel
### using merged negative control bam as background in findPeaks

    scatter (raw_fastq in expt_raw_fastqs) {
    
    String basename = sample_names_map[raw_fastq] + "_" + library_names_map[raw_fastq]

        call trimAlign as trimAlign_expt {
            input:  raw_fastq=raw_fastq,
                    basename=basename,
                    adapters=adapters,
                    effgsize=effgsize,
                    bt2index_list=bt2index_list,
                    maxfrags=maxfrags
        }

        call findPeaks {
            input:  expt_bam=trimAlign_expt.bam,
                    expt_bai=trimAlign_expt.bai,
                    ctl_bam=mergeBams.merged_bam, # optional
                    ctl_bai=mergeBams.merged_bai, # optional
                    basename=basename,
                    effgsize=effgsize
        }

        if(findPeaks.numpeaks_filt > 0) {

            if(find_motifs) {

                call motifInputs {
                    input:  peaks_narrow=findPeaks.peaks_narrow_filt,
                            basename=basename,
                            genome_fasta=genome_fasta
                }

                call findMotifs {
                    input:  summit_seqs=motifInputs.summit_seqs,
                            peak_seqs=motifInputs.peak_seqs,
                            bgmodel=bgmodel,
                            genome_fasta=genome_fasta,
                            basename=basename
                }
            }

            call assignGenes {
                input:  peaks_narrow=findPeaks.peaks_narrow_filt,
                        trim=round(trimAlign_expt.avg_fragsize/3),
                        basename=basename,
                        genes_gff=genes_gff
            }
        }
            
        call dapStats {
            input:  peaks_narrow=findPeaks.peaks_narrow,
                    expt_bam=trimAlign_expt.bam,
                    expt_bai=trimAlign_expt.bai,
                    basename=basename,
                    avg_fragsize=trimAlign_expt.avg_fragsize,
                    align_rate=trimAlign_expt.align_rate,
                    amplified=amplified
        }
        
        call copyOutput as copyOutput_expt {
            input:  basename=basename,
                    genes_gff=genes_gff,
                    effgsize=effgsize,
                    genome_fasta=genome_fasta,
                    bam=trimAlign_expt.bam,
                    bai=trimAlign_expt.bai,
                    bigwig=trimAlign_expt.bw,
                    peaks_narrow=findPeaks.peaks_narrow,
                    peaks_narrow_filt=findPeaks.peaks_narrow_filt,
                    dap_stats=dapStats.dap_stats,
                    peak_plot=dapStats.peak_plot,
                    motif1=findMotifs.motifs_summits,
                    motif2=findMotifs.motifs_summits_pal,
                    motif3=findMotifs.motifs_peaks,
                    motif4=findMotifs.motifs_peaks_pal,
                    trim_stats=trimAlign_expt.trim_stats,
                    align_stats=trimAlign_expt.align_stats,
                    macs2_stats=findPeaks.macs2_stats,
                    overlap_genes=assignGenes.overlap_genes,
                    assigned_genes=assignGenes.assigned_genes,
                    destination=outdir+"/"+basename
        }
    }
}

### Task definitions

# bbduk to trim and quality filter raw .fastq.gz
# bowtie2 to align trimmed fastq to genome
# samtools to convert sam to bam, sort, and index
# deeptools bamCoverage to make bigwig
task trimAlign {
    File raw_fastq
    String basename
    File adapters
    Int effgsize
    Array[File] bt2index_list
    Int maxfrags
    Int threads = 4
    Int memory_gb = 7
    
    runtime {
		docker: "jfroula/dap_py:3.2"
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
		bt2index_list_path=${bt2index_list[0]}
        echo "bt2index_list_path: $bt2index_list_path"

        bt2index_dir=$(dirname $bt2index_list_path)
        echo "bt2index_dir: $bt2index_dir"

        bt2index_name=$(basename $bt2index_list_path .1.bt2)
        echo "bt2index_name: $bt2index_name"

        bt2index=$bt2index_dir/$bt2index_name
        echo "bt2index: $bt2index"

        trim_align.sh \
        -i ${raw_fastq} \
        -n ${basename} \
        -a ${adapters} \
        -e ${effgsize} \
        -r $bt2index \
        -s ${maxfrags} \
        -p ${threads} \
        -m ${memory_gb}
    }
    output {
        File trim_fastq = "${basename}_sub_trim.fastq"
        File bam = "${basename}.bam"
        File bai = "${basename}.bam.bai"
        File bw = "${basename}_dedup_norm.bw"
        File trim_stats = "${basename}_trim_stats.txt"
        File align_stats = "${basename}_align_stats.txt"
        Float align_rate = read_float("${basename}_alignRate.txt")
        Float avg_fragsize = read_float("${basename}_avgFragSize.txt")
    }
}

# samtools to merge all control bams, and index merged.bam
task mergeBams {
    Array[File] bams

    runtime {
		docker: "jfroula/dap_py:3.2"
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        /bin/bash -c \
        " \
        samtools merge \
        merged.bam \
        ${sep=" " bams} \
        && \
        samtools index \
        merged.bam \
        "
    }
    output {
        File merged_bam = "merged.bam"
        File merged_bai = "merged.bam.bai"
    }
}

# macs2 to call peaks, using optional control.bam for background
task findPeaks {
    File expt_bam
    File expt_bai
    File? ctl_bam
    File? ctl_bai
    String basename
    Int effgsize
    Int min_foldch = 5

    runtime {
		docker: "jfroula/dap_py:2"
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        find_peaks.sh \
        -i ${expt_bam} \
        ${"-c " + ctl_bam} \
        -n ${basename} \
        -e ${effgsize} \
        -f ${min_foldch}
    }
    output {
        File peaks_narrow = "${basename}_peaks.narrowPeak"
        File peaks_narrow_filt = "${basename}_peaks_filt.narrowPeak"
        File summits_bed = "${basename}_summits.bed"
        Int numpeaks = read_int("${basename}_numpeaks.txt")
        Int numpeaks_filt = read_int("${basename}_numpeaks_filt.txt")
        File macs2_stats = "${basename}_macs2_stats.txt"
    }
}

# Write fasta files corresponding to the sequences in peak regions
# for inputs in the motif caller (meme).
# Call the same script twice, once to get entire peak region seqs,
# and once to get only seqs the summits +/- 30bp
task motifInputs {
    File peaks_narrow
    String basename
    File genome_fasta

    runtime {
		docker: "jfroula/dap_py:3.2"
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        /bin/bash -c \
        " \
        narrowPeak_to_fasta.py \
        -narrowPeak ${peaks_narrow} \
        -out ${basename}_summits.fasta \
        -ref ${genome_fasta} \
        -maxpeaks 100 \
        -extend 30 \
        -fimocoords \
        && \
        narrowPeak_to_fasta.py \
        -narrowPeak ${peaks_narrow} \
        -out ${basename}_peaks.fasta \
        -ref ${genome_fasta} \
        -maxpeaks 100 \
        -extend all \
        -fimocoords \
        "
    }
    output {
        File summit_seqs="${basename}_summits.fasta"
        File peak_seqs="${basename}_peaks.fasta"
    }
}

# meme (from meme-suite) to find motifs in peak regions
# fimo to map the found motifs back to the genome
task findMotifs {
    File summit_seqs
    File peak_seqs
    File bgmodel
    File genome_fasta
    String basename

    runtime {
		docker: "jfroula/dap_py:2"
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        /global/projectb/scratch/leo/dap_analysis/scripts/find_motifs.sh \
        -s ${summit_seqs} \
        -p ${peak_seqs} \
        -b ${bgmodel} \
        -r ${genome_fasta} \
        -n ${basename}
    }
    output {
        File motifs_summits="${basename}_meme_summits"
        File motifs_summits_pal="${basename}_meme_summits_pal"
        File motifs_peaks="${basename}_meme_peaks"
        File motifs_peaks_pal="${basename}_meme_peaks_pal"
    }
}

# Filter for intergenic peaks, then assign each to up to 2 genes
# Also report peaks that overlap genes
task assignGenes {
    File peaks_narrow
    Int trim
    String basename
    File genes_gff

    runtime {
		docker: "jfroula/dap_py:3.2"
        shared: 0
        time: "05:00:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        assign_genes.sh \
        -p ${peaks_narrow} \
        -t ${trim} \
        -n ${basename} \
        -r ${genes_gff}
    }
    output {
        File assigned_genes = "${basename}_assigned_genes.bed"
        File overlap_genes = "${basename}_overlap_genes.bed"
    }
}

# calculate fraction of reads in peaks
task dapStats {
    File expt_bam
    File expt_bai
    File peaks_narrow
    String basename
    Float avg_fragsize
    Float align_rate
    String amplified

    runtime {
		docker: "jfroula/dap_py:3.2"
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        dap_stats.sh \
        -i ${expt_bam} \
        -j ${peaks_narrow} \
        -n ${basename} \
        -s ${avg_fragsize} \
        -a ${align_rate} \
        -c ${amplified}
    }
    output {
        File dap_stats = "${basename}_dap_stats.tsv"
        File? peak_plot = "${basename}_peak_distribution.png"
    }
}

task copyOutput {
    String basename
    String genes_gff
    String genome_fasta
    Int effgsize
    File trim_stats
    File align_stats
    File? macs2_stats
    File bam
    File bai
    File bigwig
    File? peaks_narrow
    File? peaks_narrow_filt
    File? dap_stats
    File? peak_plot
    File? motif1
    File? motif2
    File? motif3
    File? motif4
    File? overlap_genes
    File? assigned_genes
    
    String destination

    runtime {
        shared: 0
        time: "00:15:00"
        memory: "115G"
        poolname: "dapseq_leo"
        node: 2
        nwpn: 12
    }
    command {
        mkdir -p ${destination} \
        && \
        rsync -avzh \
        ${trim_stats} \
        ${align_stats} \
        ${macs2_stats} \
        ${bam} \
        ${bai} \
        ${bigwig} \
        ${peaks_narrow} \
        ${peaks_narrow_filt} \
        ${dap_stats} \
        ${peak_plot} \
        ${motif1} \
        ${motif2} \
        ${motif3} \
        ${motif4} \
        ${overlap_genes} \
        ${assigned_genes} \
        ${destination} \
        && \
        echo "${genome_fasta}" > "${destination}/${basename}_reference.txt" \
        && \
        echo "${genes_gff}" >> "${destination}/${basename}_reference.txt" \
        && \
        echo "effectiveGenomeSize: ${effgsize}" >> "${destination}/${basename}_reference.txt" \
        && \
        chmod -R 0777 ${destination}
    }
    output {
    }
}

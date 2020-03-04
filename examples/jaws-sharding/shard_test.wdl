workflow shard_test_wf {
    meta {
        version: '0.1.0'
        author: 'Jeff Froula <jlfroula@lbl.gov>'
    }

    File reads
    File reference
    Int chunk_size=100000000

    call shard {
        input: reads = reads,
        chunk_size=chunk_size
    }
    call bbmap_indexing {
        input: reference = reference
    }

    Array[String] shard_array = shard.shards

    scatter(coords in shard_array) {
        call alignment {
        input: reads = reads,
               ref = bbmap_indexing.ref,
               coords = coords
        }
    }

	call merge_bams {
		input: bams = alignment.bam
	}
}

task shard {
    File reads
    String bname = basename(reads)
    Int chunk_size

    command {
        set -e -o pipefail
        shifter --image=jfroula/jaws-sharding:1.0.10 fastq_indexer.py --input ${reads} --output ${bname}.index

        shifter --image=jfroula/jaws-sharding:1.0.10 create_blocks.py -if ${bname}.index -ff ${reads} -of ${bname}.sharded -bs ${chunk_size}
    }

    output {
        Array[String] shards = read_lines("${bname}.sharded")
    }
}

task bbmap_indexing {
    File reference

    command {
        shifter --image=jfroula/aligner-bbmap:1.0.0 bbmap.sh ref=${reference} 
    }

    output {
        File ref = "ref"
    }
}

task alignment {
    File reads
    File ref
    String coords
    Int threads=4
    String bname = basename(reads)
    String path_to_ref = sub(ref, "/ref", "")

    command<<<
        start=$(echo ${coords} | awk '{print $1}')
        end=$(echo ${coords} | awk '{print $2}')

		# we are piping a block of the fastq sequence to the aligner 
        shifter --image=jfroula/jaws-sharding:1.0.10 \
		  shard_reader.py -i ${reads} -s $start -e $end | \
        shifter --image=jfroula/aligner-bbmap:1.0.0 \
		  bbmap.sh int in=stdin.fq path=${path_to_ref} out=${bname}.sam overwrite keepnames mappedonly threads=${threads}

		# create a sorted bam file from the sam file
		shifter --image=jfroula/aligner-bbmap:1.0.0 \
		  samtools view -uS ${bname}.sam | \
		shifter --image=jfroula/aligner-bbmap:1.0.0 \
		  samtools sort - -o ${bname}.sorted.bam
    >>>

    output {
        File bam = "${bname}.sorted.bam"
    }
}
task merge_bams {
    Array[File] bams

    command {
	  shifter --image=jfroula/aligner-bbmap:1.0.0 \
	    picard MergeSamFiles I=${sep=' I=' bams} OUTPUT=merged.sorted.bam SORT_ORDER=coordinate ASSUME_SORTED=true USE_THREADING=true
    }

    output {
       File merged = "merged.sorted.bam"
    }
}


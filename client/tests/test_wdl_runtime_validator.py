#!/usr/bin/env python
"""
These functions are to test the "runtime{}" validation tests are working.

This library uses "fixtures" from conftest.py which should be located in the same directory.
There is no need to import conftest.py as it is done automatically.

Author: Jeff Froula jlfroula@lbl.gov>
Updated: 04/18/21
"""

import pytest
from jaws_client.wdl_runtime_validator import (
    validate_wdl_runtime,
    WdlRuntimeError,
)  # noqa


def test_good_wdl():
    validate_wdl_runtime(GOOD_WDL, 500)


def test_allRequiredParams():
    """
    Make sure there is an error message when a WDL has a runtime missing
    time or memory, which are the minimum settings.
    """
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(NO_TIME_AND_MEMORY, 500)

    validate_wdl_runtime(TIME_AND_MEMORY_ARE_VARIABLES, 500)


def test_timeParam():
    """
    1. if constraint: knl, then limit is 48hrs
    2. You are limited to 48hrs where constraint=knl"
    3. You are limited to 72hrs when constraint is the default value "haswell"
    4. You are limited to 72hrs where constraint=haswell"
    """
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(BAD_TIME_LIMITS, 500)


def test_memoryParam():
    """Does the test for over-memory request give valid user error."""
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(BAD_MEMORY_VALUE, 100)

    validate_wdl_runtime(BAD_MEMORY_VALUE, 500)


def test_runtimeCombinations():
    """testing the test that validates the user included the correct combination of params for various resource allocations"""  # noqa
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(BAD_COMBINATIONS, 500)


BAD_COMBINATIONS = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
}
task bbcms {
     File infile
     String container

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        qos: "regular"
        account: "wrong"
    }

}
task read_mapping_pairs{

    runtime {
        docker: container
        time: "03:00:00"
        memory: "115G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
        qos: "regular"

    }
}
"""


BAD_MEMORY_VALUE = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
}
task bbcms {
    runtime {
        docker: container
        time: "03:00:00"
        memory: "200G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task assy {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "300G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }
}

"""


BAD_TIME_LIMITS = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
    call assy {
         input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    call read_mapping_pairs {
     input: reads=input_file, ref=create_agp.outcontigs, container=bbtools_container
    }

}
task bbcms {
    runtime {
        docker: container
        time: "200:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "skylake"
    }

}

task assy {
    runtime {
        docker: container
        time: "80:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "knl"
    }

}

task create_agp {

    runtime {
        docker: container
        time: "80:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}


task read_mapping_pairs{
    runtime {
        docker: container
        time: "80:00:00"
        memory: "115G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
    }

}
"""


NO_TIME_AND_MEMORY = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
    call assy {
         input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    call read_mapping_pairs {
     input: reads=input_file, ref=create_agp.outcontigs, container=bbtools_container
    }

}
task bbcms {
    runtime {
        docker: container
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task assy {
    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task create_agp {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task read_mapping_pairs{

    runtime {
        docker: container
        time: "03:00:00"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
    }

}
"""


GOOD_WDL = """
workflow jgi_meta {
    File input_file
    Float uniquekmer=1000
    String bbtools_container="bryce911/bbtools:38.44"
    String spades_container="bryce911/spades:3.13.0"
    String basic_container="bryce911/bbtools:38.44"
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
    call assy {
         input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    call read_mapping_pairs {
     input: reads=input_file, ref=create_agp.outcontigs, container=bbtools_container
    }

}
task bbcms {
     File infile
     String container

     String filename_resources="resources.log"
     String filename_outfile="input.corr.fastq.gz"
     String filename_outfile1="input.corr.left.fastq.gz"
     String filename_outfile2="input.corr.right.fastq.gz"
     String filename_readlen="readlen.txt"
     String filename_outlog="stdout.log"
     String filename_errlog="stderr.log"
     String filename_kmerfile="unique31mer.txt"
     String filename_counts="counts.metadata.json"
     String dollar="$"

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

     command <<<
        touch ${filename_resources} && \
        bbcms.sh -Xmx105g metadatafile=${filename_counts} mincount=2 highcountfraction=0.6 \
        in=${infile} out=${filename_outfile} > >(tee -a ${filename_outlog}) \
        2> >(tee -a ${filename_errlog} >&2) && grep Unique ${filename_errlog} | \
        rev |  cut -f 1 | rev  > ${filename_kmerfile} && \
        reformat.sh -Xmx105g in=${filename_outfile} out1=${filename_outfile1} out2=${filename_outfile2} && \
        readlength.sh -Xmx105g in=${filename_outfile} out=${filename_readlen}
     >>>

     output {
            File out = filename_outfile
            File out1 = filename_outfile1
            File out2 = filename_outfile2
            File outreadlen = filename_readlen
            File stdout = filename_outlog
            File stderr = filename_errlog
            File outcounts = filename_counts
            File outkmer = filename_kmerfile
            File outresources = filename_resources
     }
}

task assy {
     File infile1
     File infile2
     String container

     String filename_resources="resources.log"
     String outprefix="spades3"
     String filename_outfile="${outprefix}/scaffolds.fasta"
     String filename_spadeslog ="${outprefix}/spades.log"
     String dollar="$"

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

     command <<<
        touch ${filename_resources} && \
        spades.py -m 2000 --tmp-dir /tmp -o ${outprefix} --only-assembler -k 33,55,77,99,127  \
        --meta -t ${dollar}(grep "model name" /proc/cpuinfo | wc -l) -1 ${infile1} -2 ${infile2}
     >>>

     output {
            File out = filename_outfile
            File outlog = filename_spadeslog
            File outresources = filename_resources
     }
}


task create_agp {
    File scaffolds_in
    String container

    String filename_resources="resources.log"
    String prefix="assembly"
    String filename_contigs="${prefix}.contigs.fasta"
    String filename_scaffolds="${prefix}.scaffolds.fasta"
    String filename_agp="${prefix}.agp"
    String filename_legend="${prefix}.scaffolds.legend"

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

    command{
        touch ${filename_resources} && \
        fungalrelease.sh -Xmx105g in=${scaffolds_in} out=${filename_scaffolds} \
        outc=${filename_contigs} agp=${filename_agp} legend=${filename_legend} \
        mincontig=200 minscaf=200 sortscaffolds=t sortcontigs=t overwrite=t
  }

    output{
        File outcontigs = filename_contigs
        File outscaffolds = filename_scaffolds
        File outagp = filename_agp
        File outlegend = filename_legend
        File outresources = filename_resources
    }
}


task read_mapping_pairs{
    File reads
    File ref
    String container

    String filename_resources="resources.log"
    String filename_unsorted="pairedMapped.bam"
    String filename_outsam="pairedMapped.sam.gz"
    String filename_sorted="pairedMapped_sorted.bam"
    String filename_sorted_idx="pairedMapped_sorted.bam.bai"
    String filename_bamscript="to_bam.sh"
    String filename_cov="covstats.txt"
    String dollar="$"

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

    command{
        touch ${filename_resources} && \
        bbmap.sh -Xmx105g threads=${dollar}(grep "model name" /proc/cpuinfo | wc -l) \
        nodisk=true interleaved=true ambiguous=random in=${reads} ref=${ref} \
        out=${filename_unsorted} covstats=${filename_cov} bamscript=${filename_bamscript} && \
        samtools sort -m100M -@ ${dollar}(grep "model name" /proc/cpuinfo | wc -l) \
        ${filename_unsorted} -o ${filename_sorted} && \
        samtools index ${filename_sorted} && \
        reformat.sh -Xmx105g in=${filename_unsorted} out=${filename_outsam} overwrite=true
    }


  output {
      File outbamfile = filename_sorted
      File outbamfileidx = filename_sorted_idx
      File outcovfile = filename_cov
      File outsamfile = filename_outsam
      File outresources = filename_resources
  }
}
"""
TIME_AND_MEMORY_ARE_VARIABLES =  """
workflow jgi_meta {
    String time = "00:30:00"
    String memory = "5G"
    String cpu = 5

    call task1 {
	  input: time=time, memory=memory, cpu=cpu
    }
}
task task1 {
	String time
	String memory
	String cpu

    runtime {
		docker: "doejgi/jaws-debian:latest"
        node: 1
        nwpn: 1
        poolname: "dashboard_test"
        shared: 0
		time: time 
		memory: memory
		cpu: cpu
    }

     command {
        echo "task one gobble-di-gook" > output.txt
        sleep 50
     }

     output {
            File out = "output.txt"
     }
}
"""

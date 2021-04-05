workflow jgi_rqcfilter {
    Array[File] input_files
    String bbtools_container="microbiomedata/bbtools:38.44"
    String database="/global/cscratch1/sd/jaws/refdata/nmdc/RQCFilterData"

    scatter(file in input_files) {
        call rqcfilter{
             input:  input_file=file, container=bbtools_container, database=database
        }
    }

    parameter_meta {
        input_files: "illumina paired-end interleaved fastq files"
  	    outdir: "The final output directory path"
        database : "database path to RQCFilterData directory"
        clean_fastq_files: "after QC fastq files"
    }
    meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.0"
    }
}

task rqcfilter {
     File input_file
     String container
     String database
     String filename_outlog="stdout.log"
     String filename_errlog="stderr.log"
     String filename_stat="filtered/filterStats.txt"
     String filename_stat2="filtered/filterStats2.txt"
     String dollar="$"
     String num_threads=8

     runtime {
        docker: container
        time: "01:00:00"
        memory: "40G"
        poolname: "nmdc_readqc_test"
        shared: 0
        node: 1
        nwpn: 1
     }

     command {
        #sleep 30
        export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
        set -eo pipefail
        rqcfilter2.sh -Xmx35g threads=${num_threads} jni=t in=${input_file} path=filtered rna=f trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t tmpdir= barcodefilter=f trimpolyg=5 usejni=f rqcfilterdata=${database} > >(tee -a ${filename_outlog}) 2> >(tee -a ${filename_errlog} >&2)
     }
     output {
            File stdout = filename_outlog
            File stderr = filename_errlog
            File stat = filename_stat
            File stat2 = filename_stat2
     }
}

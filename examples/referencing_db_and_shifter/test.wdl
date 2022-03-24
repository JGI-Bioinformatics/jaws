workflow runblastplus_sub {
    File reads
    String ncbi_nt
    
    call task1 { input: ncbi_nt = ncbi_nt }
    call task2 { input: tmpfile = task1.outfile }
}

### ------------------------------------ ###
task task1 {
	String ncbi_nt

    command {
      # how to access reference data
      ls ${ncbi_nt} > tmp.txt
    }

    runtime {   
        docker: "ubuntu:16.04"
        poolname: "referece_db_pool"
        node: 1
        nwpn: 1
        memory: "5G"
        time: "00:10:00"
    }

    output { File outfile = "tmp.txt" }
}

task task2 {
    File tmpfile

    command {
		cat ${tmpfile}
    }

    output { String outfile = stdout() }
}


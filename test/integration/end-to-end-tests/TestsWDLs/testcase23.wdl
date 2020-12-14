workflow runblastplus_sub {
    File reads
    String ncbi_nt="/refdata"
    
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
        poolname: "useforalltests"
        shared: 0
        node: 1
        nwpn: 1
        mem: "5G"
        time: "99:99:99"
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


workflow runblastplus_sub {
    File reads
<<<<<<< HEAD
    String ncbi_nt
    
    call task1 { input: ncbi_nt = ncbi_nt }
=======
    call task1 { }
>>>>>>> 8fdd6c85f29e50f51ba78169b488b7c9a681c11d
    call task2 { input: tmpfile = task1.outfile }
}

### ------------------------------------ ###
task task1 {
<<<<<<< HEAD
	String ncbi_nt

    command {
      # how to access reference data
      ls ${ncbi_nt} > tmp.txt
    }

    runtime {   
        docker: "ubuntu:16.04"
        poolname: "mysmall"
=======
    command {
      # how to access reference data
      ls /refdata > tmp.txt
    }

    runtime {
        docker: "jfroula/jaws-blastplus:1.0.18"
        poolname: "muysmall"
        shared: 1
>>>>>>> 8fdd6c85f29e50f51ba78169b488b7c9a681c11d
        node: 1
        nwpn: 1
        mem: "5G"
        time: "00:10:00"
    }

    output { File outfile = "tmp.txt" }
}

task task2 {
    File tmpfile

    command {
<<<<<<< HEAD
      cat ${tmpfile}
=======
        cat ${tmpfile}
>>>>>>> 8fdd6c85f29e50f51ba78169b488b7c9a681c11d
    }

    output { String outfile = stdout() }
}


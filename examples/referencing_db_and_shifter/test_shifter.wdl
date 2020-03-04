workflow runblastplus_sub {
    File reads
    call task1 { }
    call task2 { input: tmpfile = task1.outfile }
}

### ------------------------------------ ###
task task1 {
    command {
      # how to access reference data
      ls /refdata > tmp.txt
    }

    runtime {
        docker: "jfroula/jaws-blastplus:1.0.18"
        poolname: "muysmall"
        shared: 1
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
        cat ${tmpfile}
    }

    output { String outfile = stdout() }
}


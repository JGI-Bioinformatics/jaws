workflow runblastplus_sub {
    call task1 { }
    call task2 { input: tmpfile = task1.outfile }
}

### ------------------------------------ ###
task task1 {

    command {
      # how to access reference data
      ls /refdata/i_dont_exist > tmp.txt
    }

    runtime {   
        poolname: "referece_db_pool"
        shared: 0
        node: 1
        nwpn: 1
        memory: "5G"
        time: "00:20:00"
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


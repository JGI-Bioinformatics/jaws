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
        memory: "5G"
        time: "00:20:00"
        cpu: 1
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


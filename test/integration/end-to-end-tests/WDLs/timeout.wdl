workflow timeout {
    String label
    call task_timeout { 
        input: label = label 
    }
    output {
        File outfile = task_timeout.outfile
    }
}

task task_timeout {
    String label
    command <<<
        echo $label > timeout.txt
        date >> timeout.txt
        sleep 10m
        date >> timeout.txt
    >>>
    output {
        File outfile = "timeout.txt"
    }
    runtime {
        poolname: "test_small"
        node: 1
        nwpn: 1
        memory: "5G"
        time: "00:02:00"
        shared: 0
    }
}


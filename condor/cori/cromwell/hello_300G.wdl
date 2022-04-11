workflow myWorkflow {
    call myTaskOut
    output { File out = myTaskOut.out }
}

task myTask {
    command {
        echo "hello world"
    }
    output {
        String out = read_string(stdout())
    }

    runtime {
        #poolname: "sulsjhello"
        #node: 1
        #nwpn: 1
        cpu: 4
        memory: "300G"
        #time: "00:10:00"
        #shared: 0
    }
}

task myTaskOut {
    command {
        echo "hello workld!"
        echo "hello world" > hello.txt && wc -l hello.txt > wc.out
    }
    output {
        File out = "wc.out"
    }
    runtime {
        cpu: 4
        memory: "300G"
    }
}

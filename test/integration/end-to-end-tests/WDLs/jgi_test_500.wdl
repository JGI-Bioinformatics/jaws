workflow myWorkflow {
    call myTaskOut
    output { File out = myTaskOut.out }
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
        poolname: "helloworldtest"
        node: 1
        nwpn: 1
        memory: "500G"
        time: "00:10:00"        
        constraint: "lr3_c32,jgi_m512"
        account: "lr_jgicloud"
        qos: "condo_jgicloud"
        partition: "lr3"
    }
}

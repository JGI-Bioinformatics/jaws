workflow ConfigurableScatter {
    Int num_tasks
    Int task_time

    call makeArray { input: num = num_tasks }
    
    scatter (one in makeArray.out) {
        call stepA { input: in=one, sleep_time=task_time }
    }
    
    call stepB { input: files=stepA.out }
  
}

task makeArray {
    Int num

    command {
        for ((i=0;i<${num};i++))
        do
            echo "Line $i"
        done
    }

    output { Array [String] out = read_lines(stdout()) }
}

task stepA {

  String in
  Int sleep_time

  command <<< 
       echo "$(date) - ${in} job started" 
       echo "sleeping for ${sleep_time} seconds ..."
       sleep ${sleep_time}
       echo "$(date) - job ended"
       echo "-----------------"
  >>> 

  runtime {
    time: "02:00:00"
    memory: "5G"
    poolname: "configurable_scatter"
    shared: 0
    node: 1
    nwpn: 1           
  }
  
  output { File out = stdout() }
}

task stepB {

  Array[File] files

  command { cat ${sep=' ' files} > final.txt }  # simalar to python's ( ' '.join(array) )
  output { File out = "final.txt" }
}
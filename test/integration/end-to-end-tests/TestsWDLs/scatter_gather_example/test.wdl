workflow ScatterGather {

  Array[String] test_input

  scatter (one in test_input) {
    call stepA { input: in=one }
  }
  call stepB { input: files=stepA.out }
}

task stepA {

  String in

  command <<< 
       echo "${in} is sleeping for 5 seconds on $(date)" > file;
       sleep 5
       echo "job ended $(date)" >> file;
       echo "-----------------" >> file;
  >>> 

  runtime {
    time: "00:20:00"
    mem: "5G"
    poolname: "smallworkerpoolforscattergather"
    shared: 0
    node: 1
    nwpn: 1           
  }

  output { File out = "file" }
}

task stepB {

  Array[File] files

  command { cat ${sep=' ' files} > final }  # simalar to python's ( ' '.join(array) )
  output { String out = "final" }
}

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
    cluster: "jaws_lbl_gov"
    time: "00:10:00"
    mem: "5G"
    poolname: "smallworkerpoolforscattergather"
    shared: 1
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

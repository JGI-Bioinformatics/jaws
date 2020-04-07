workflow ScatterGather {

  Array[String] test_input

  scatter (one in test_input) {
    call stepA { input: in=one }
  }
  call stepB { input: files=stepA.out }
}

task stepA {

  String in

  command { 
       echo ${in} is sleeping for 5 seconds > file;
       sleep 5
  }
  output { File out = "file" }
}

task stepB {

  Array[File] files

  command { cat ${sep=' ' files} }  # simalar to python's ( ' '.join(array) )
  output { String out = read_string(stdout()) }
}

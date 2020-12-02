workflow ScatterGather {

  Array[String] test_input

  # we'll run stepA in parallel
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

  # we'll create a worker pool of two, so each task can run in parallel
  runtime {
	poolname="anything_unique"
  	node=1
    nwpn=2
  }

  # an out file gets created for each shard so in the next step we'll have
  # an array of files as input
  output { File out = "file" }
}

task stepB {

  # we have an array of files created from all the previous sharding jobs.
  Array[File] files

  # this ${sep=' ' files} is a special wdl function that makes a string out
  # out of an array simalar to python's ( ' '.join(array) )
  command { cat ${sep=' ' files} }  

  output { String out = read_string(stdout()) }
}

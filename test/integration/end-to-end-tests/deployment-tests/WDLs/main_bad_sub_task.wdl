import "bad_sub_task.wdl" as sub

workflow main {
    String in1
    String out1
    String in2
    String out2
    String out3
    
    call sub.sub_workflow as echo1 {
      input: in=in1,
          out=out1
    }
    
    call sub.sub_workflow as echo2 {
      input: in=in2,
          out=out2
    }
    
    call cat {
        input: file1=echo1.sub_out,
            file2=echo2.sub_out,
            outFile=out3
    }
}

task cat {
    File file1
    File file2
    String outFile
    command {
      echo  "cat task msg to STDOUT"
      1>&2 echo "cat task msg to STDERR"  
      cat ${file1} ${file2} >> ${outFile}.txt
    }
    runtime {
      poolname: "alkmain"
      shared: 0
      time: "0:20:00"
      memory: "5G"
      node: 1
      nwpn: 1
    }
    output {
        File out = "${outFile}.txt"
    }
}


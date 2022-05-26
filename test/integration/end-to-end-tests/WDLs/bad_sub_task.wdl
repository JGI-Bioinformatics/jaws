workflow sub_workflow {
    String in
    String out

    call echo {
      input: text=in,
          outFile=out
    }

    output {
      File sub_out=echo.out
    }
}

task echo {
    String text
    String outFile
    command {
        echo "echo task in sub, msg to STDOUT" 
        1>&2 echo "echo task in sub, msg to STDERR"
        echoooo ${text} > ${outFile}.txt
    }

    runtime {
      time: "0:10:00"
      memory: "5G"
      cpu: 1
    }

    output {
        File out = "${outFile}.txt"
    }
}


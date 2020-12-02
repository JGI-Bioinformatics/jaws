import "sub.wdl" as sub

workflow main_wdl { 
    call sub.sub_workflow { 
           input: sub_input = "we are running sub-workflow"
     }
    call task2 {
           input: in = sub_workflow.run_task_output
    }
}

task task2 {
    String in
    command {
        # does something with "in"
    }
}


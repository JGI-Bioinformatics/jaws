workflow sub1_workflow {
    String sub1_input

    call run_task1 { input: sub1_input = "we are running sub-workflow 1" }

    output { String run_task1_output = run_task1.status  }
}

task run_task1 {
    String sub1_input

    command
    {    
        echo ${sub1_input}
	}
    output {
        String status = "task1 was run successfully"
    }
}




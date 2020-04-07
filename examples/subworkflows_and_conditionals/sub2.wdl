workflow sub2_workflow {
    String sub2_input

    call run_task2 { input: sub2_input = "we are running sub-workflow 2" }

    output { String run_task2_output = run_task2.status  }
}

task run_task2 {
    String sub2_input

    command
    {    
        echo ${sub2_input}
	}
    output {
        String status = "task2 was run successfully"
    }
}




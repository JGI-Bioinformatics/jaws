workflow sub_workflow { 
    String sub_input 
    
    call run_task { 
        input: sub_input = sub_input
    } 

    output { 
        String run_task_output = run_task.status 
     }
} 

task run_task { 
    String sub_input 
    
    command { 
        echo ${sub_input}
    } 
    output { 
        String status = read_string(stdout())

    }
} 

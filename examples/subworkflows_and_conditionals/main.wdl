import "sub1.wdl" as sub1
import "sub2.wdl" as sub2

workflow main_wdl {
    Boolean flag 

    call run_preprocess { input: preprocess_input = "we are running preprocessing steps" }

    if( flag ){
        call sub1.sub1_workflow { input: sub1_input = "we are running sub-workflow 1" }
    }
    if( ! flag ){
        call sub2.sub2_workflow { input: sub2_input = "we are running sub-workflow 2"}
    }

    output { 
        String? main_output_wdl1 = sub1_workflow.run_task1_output
        String? main_output_wdl2 = sub2_workflow.run_task2_output
    }
}

task run_preprocess {
    String preprocess_input

    command
    {    
        echo ${preprocess_input}
	}
}




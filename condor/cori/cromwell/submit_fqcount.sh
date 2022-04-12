curl -X POST --header "Accept: application/json"    -v "localhost:50011/api/workflows/v1" -F workflowSource=@fq_count.wdl -F workflowInputs=@fq_count.json -F workflowOptions=@options.json


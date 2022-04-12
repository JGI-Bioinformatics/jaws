curl -X POST --header "Accept: application/json"    -v "localhost:50011/api/workflows/v1" -F workflowSource=@hello.wdl -F workflowInputs=@inputs.json -F workflowOptions=@options.json

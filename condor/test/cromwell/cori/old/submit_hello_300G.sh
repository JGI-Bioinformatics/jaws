curl -X POST --header "Accept: application/json"    -v "localhost:50011/api/workflows/v1" -F workflowSource=@hello_300G.wdl -F workflowInputs=@inputs.json

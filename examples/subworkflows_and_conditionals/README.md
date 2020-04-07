# Running sub-workflows and using conditional statements

## this example includes two wdl sub-workflows, and which one is run depends on some flag (using conditional statements).

### log onto cori or denovo
```
ssh denovo.nersc.gov or cori.nersc.gov
module load java # this may only be required on cori to load java v1.8
```  

### validate wdl main.wdl
```
java -jar /global/cfs/projectdirs/jaws/cromwell/womtool.jar validate main.wdl
```  

### run it
```
java -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run main.wdl --inputs inputs.json
```  


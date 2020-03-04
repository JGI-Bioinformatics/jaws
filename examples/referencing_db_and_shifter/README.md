# set dockerRoot variable
You need to make sure the config file's `dockerRoot` is set to your current working directory and you've added `cromwell-executions` at the end.  For example:
`dockerRoot = /path/to/current/working/dir/cromwell-executions`

# How to Access Reference Databases
in your commands section in the WDL, you use `/refdata` as a root, and then add whatever specific database directory you want, for example to use  the blast nt database from NCBI you would have a command
`blastn -db /refdata/nt/nt`
where the first `nt` is the directory with all the index files and the second `nt` is the prefix to the index files (i.e. `nt.nih`).

# Run Example 
This workflow just lists the contents of the nt database.

```
java -Dconfig.file=shifter.conf -jar /global/dna/projectdirs/DSI/workflows/cromwell/java/cromwell.jar run test_shifter.wdl
```

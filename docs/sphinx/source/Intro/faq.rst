====
FAQs
====

##################
JAWS command line
##################
    
Does Cromwell offer checkpointing?
    sort of; Cromwell has call caching instead which accomplishes the same thing. When a task completes successfully, it's results are capable of being reused if the same task and inputs are run again. Use `jaws submit --no-cache` to turn caching off.

|

Why didn't call caching work for me?
    Changes to the WDL, the name contents of the inputs.json, or the name of the inputs.json will prevent call-caching.

    For example, if you set your task's runtime attributes using input variables, changes to the values of these variables count as changes to the inputs, resulting in a different hash for the task (the wdl and inputs.json are hashed).

    Call caching may have failed if your files are being fed in as String rather than File inputs. The hashes of two identical Files stored in different locations would be the same. The hashes of the String values for the different locations would be different, even though the contents of the file are the same.

    Call caching also requires consistency in the outputs of the task, both the count (number of outputs) and the output expressions. If you publish a new version of your WDL that has one extra or one fewer output, it will not be able to benefit from a previously successful run of the same task, even if the inputs are the same.
    
    The only things you can change are the (1) filepaths inside the inputs.json if they are declared as "File" and not "String" and (2) hard-coded values inside the runtime{} section (except for docker).

|

Can I copy only a select set of output files so I'm not copying uneccessary files.
    When JAWS is finished (i.e. "download-complete") you should be able to copy your output files using the 'jaws get' command. When you use the command without the --complete flag, you will only get the files that were listed in the output{} section of the WDL.  The --complete flag will get you all the files in cromwell's 'execution' directory.

|

What flavor of linux do the compute nodes run?
    JAWS makes multiple computing resources available, using various linux distros.  Thus, we recommend that a docker container be specified for every task; if not, the default container is Debian.

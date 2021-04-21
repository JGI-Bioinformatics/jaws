====
FAQs
====

##################
JAWS command line
##################
    
Does Cromwell offer checkpointing?
    sort of; Cromwell has call caching instead which accomplishes the same thing. When a task completes successfully, it's results are capable of being reused if the same task and inputs are run again. Use `jaws submit --no-cache` to turn caching off.
    
|


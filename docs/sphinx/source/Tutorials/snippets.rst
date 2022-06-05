==================================================
Code Snippets to Answer Common WDL Design Problems
==================================================

.. role:: bash(code)
    :language: bash

OpenWDL provides the WDL functions in `specs for version 1.0 <https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md>`_

#############
Building WDLs
#############

.. raw:: html

  <details>
  <summary><a>How can I use bash commands that require curly braces?</a></summary>

  <br>
  If you ever need to use curly braces in bash to strip a suffix txt or set a default, there are two ways: 1) make your WDL use "version 1.0", or 2) write a hack as shown below:
  
  <br><br>This is the prefered way; have <a href=https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#versioning>version 1.0</a> as the first line in your WDL. You'll specify the command section like "command <<< >>>" instead of using curly braces. You'll have to make some other formatting changes too, see the link to the version 1.0 spec.

  <code>
    <pre>
        command <<<
            # setting a default value in bash
            VAR=${VAR:=25}

            # strip a suffix
            myvar=${somefile%.txt}
        >>>
    </pre>
  </code>

  <br><br>This is the hack if you want to keep your WDL the default version which is "draft"

  <code>
    <pre>

        task somthing {
            String dollar='$'
            command { 
                ${dollar}{parameter:=default} 
            }
        }
    </pre>
  </code>
  </details>


  <details>
  <summary><a> How can I output a file that has been named dynamically as a bash variable </a></summary>
  <br>
  Bash variables created in the command{} block cannot be seen outside the block, for example, in the output {} section. Therefore, you can write the name(s) of any output files to another file which will be read inside the output {} block.
  <br>
  <br>
  This is the official WDL way, using glob
  <code>
    <pre>
        output {
            Array[File] output_bams = glob("*.bam")
        }
    </pre>
  </code>

  This is another method

  <code>
    <pre>
        command{
           echo $lib.bam > list_of_files
        }
        output {
           Array[File] = read_lines("list_of_files")
        }
        

    To see more about read_lines() and other WDL functions, see `openwdl/wdl <https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md>`_
    </pre>
  </code>
  </details>

  <details>
  <summary><a> Using conditionals</a></summary>

  <code><pre>

        workflow conditional_example {
          File infile

          call wc as wc_before { input: infile = infile }

          Int num_lines = wc_before.num_lines

          if (num_lines > 10) {
            call truncate { input: infile = infile }
          }

          # This function will return false if the defined() argument is an 
          # unset optional value. It will return true in all other cases.
          Boolean has_head_file = defined(truncate.outfile)

          if (has_head_file) {
            call wc as wc_after { input: infile = truncate.outfile }
          }

          # notice the '?' after File. These are required since these files may not exist.
          output {
            File wc_before_file = wc_before.outfile
            File? head_file = truncate.outfile
            File? wc_after_file = wc_after.outfile
          }
        }

        task wc {
          File infile
          command { wc -l < ${infile} | tee wc.txt }
          output {
            Int num_lines = read_int(stdout())
            File outfile = "wc.txt"
          }
        }
  </pre></code></details>

  <details>
  <summary><a> How to scatter over arrays and maps </a></summary>
  <br>
    Although you can scatter over arrays and maps, there is different syntax for each.
    You can only scatter over an array with this syntax
  <br>
    
  <code><pre>
        Array[String] some_array
        scatter (e in some_array) {
            String value = some_array[e]
            call some_task {input: value = value}
        }
  </pre></code>

  But you can iterate over a map by using the 'pair' keyword and then '.left' and '.right' as such

  <code><pre>
        Map[String,String] some_map
        scatter (pair in some_map) {
            String key= pair.left
            String value = pair.right # or String val = some_map[key]
            call some_task {input: value = value}
        }

    You can see working examples for <a href=https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/tree/master/scatter_gather_example>scattering an array and scattering a map</a> 
  </pre></code></details>

  <details>
  <summary><a> Custom data structures </a></summary>
  <br>
    Besides Map, Array, Pair you can create a custom data structure using "struct". This will be similar to a hash but can contain any combination of data types. 
    <br>
    <ul>
      <li>Documentation for <a href=https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#custom-types-structs>Custom Type "Struct"</a></li>
      <li>Example <a href=https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/blob/main/custom_datastructure/main.wdl>main.wdl</a> && <a href=https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/blob/main/custom_datastructure/inputs.json>inputs.json</a></li>
    </ul>
    <br>
  </details>

  <details>
  <summary><a> How to copy a whole directory that is listed in my inputs.json</a></summary>
  <br>
	Sometimes you may want to copy all the contents of a directory. Unfortunately Cromwell doesn't allow for this (it's a limitation of the variable declaration "File"). 
    <br><br>1. One solution would be to include a tar file in the inputs.json and then untar it inside the task.
	<br><br>2. Another solution is to list all the files in the inputs.json.  The files will be put into a cromwell generated folder inside the "inputs" directory of that task (i.e. inputs/-697750178/). See the example WDL and inputs json that shows you how you would access that folder.
  <br><a href=https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/blob/main/copy-refdata-as-inputs/refdata.wdl>refdata.wdl</a> && <a href=https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/blob/main/copy-refdata-as-inputs/refdata.json>refdata.json</a>
  </details>
|
|

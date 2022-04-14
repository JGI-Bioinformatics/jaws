==================================================
Code Snippets to Answer Common WDL Design Problems
==================================================

.. role:: bash(code)
    :language: bash

OpenWDL provides the WDL functions in `specs for version 1.0 <https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md>`_

#############
Building WDLs
#############

How can I use bash commands that require curly braces?
    If you ever need to use curly braces in bash to strip a suffix txt or set a default:

    .. code-block:: text

        # strip a suffix
        myvar=${somefile%.txt}
        or 
        # set defaults
        ${parameter:=default}


    Then you need to do

    .. code-block:: text

        task somthing {
            String dollar='$'
            command { 
                ${dollar}{parameter:=default} 
            }
        }


How can I output a file that has been named dynamically as a bash variable
    Bash variables created in the command{} block cannot be seen outside the block, for example, in the output {} section. Therefore, you can write the name(s) of any output files to another file which will be read inside the output {} block.

    This is the official WDL way, using glob

    .. code-block:: text

        output {
            Array[File] output_bams = glob("*.bam")
        }

    This is another method

    .. code-block:: text

        command{
           echo $lib.bam > list_of_files
        }
        output {
           Array[File] = read_lines("list_of_files")
        }
        

    To see more about read_lines() and other WDL functions, see `openwdl/wdl <https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md>`_


Using Conditionals

    .. code-block:: text

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


How to scatter over arrays and maps
    Although you can scatter over arrays and maps, there is different syntax for each.
    You can only scatter over an array with this syntax
    
    .. code-block:: text

        Array[String] some_array
        scatter (e in some_array) {
            String value = some_array[e]
            call some_task {input: value = value}
        }

    But you can iterate over a map by using the :bash:`pair` keyword and then :bash:`.left` and :bash:`.right` as such

    .. code-block:: text

        Map[String,String] some_map
        scatter (pair in some_map) {
            String key= pair.left
            String value = pair.right # or String val = some_map[key]
            call some_task {input: value = value}
        }

    You can see working examples for `scattering an array and scattering a map <https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/tree/master/scatter_gather_example>`_.


Custom data structures
    Besides Map, Array, Pair you can create a custom data structure using "struct". This will be similar to a hash but can contain any combination of data types.  See `WDL Spec for v1.1: Custom Type "Struct" <https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#custom-types-structs>`_.

    Example 
    `main.wdl <https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/blob/main/custom_datastructure/main.wdl>`_ && `inputs.json <https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/blob/main/custom_datastructure/inputs.json>`_

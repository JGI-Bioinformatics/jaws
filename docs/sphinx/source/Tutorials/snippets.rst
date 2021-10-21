==================================================
Code Snippets to Answer Common WDL Design Problems
==================================================

#############
Building WDLs
#############

How do I use Arrays and Maps in my WDL. 
    Specifically, how do I dereference the contents of the array or map so I can use them in my commands?
    This example was copied from github:gist `scottfrazer/style_guide.md <https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211>`_.  Also, you can see more about the functions used here on the official Broad Institute's `WDL spec.md <https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md>`_.
    
    WDL allows compound types like Array[String] or Map[String, Int] or Array[Array[String]]. There are two ways to get these data types into a form that the command can use:
    
    1. Serialization by concatenation (only for Array)
    2. Serialization by write-to-file

    Use WDL functions for common transformations
    
    
    .. code-block:: text

        task example {
          Array[String] array
          Map[String, File] map
          Array[Array[Int]] matrix
          
          command {
            echo ${sep=',' array}
            cat ${write_lines(array)}
            python script.py --map=${write_map(map)}
            python process.py ${write_tsv(matrix)}
          }
        }
        
        workflow test {
          call example
        }
        {
          "test.example.array": ["a", "b", "c"],
          "test.example.map": {
            "key0": "/path/to/file0",
            "key1": "/path/to/file1",
            "key2": "/path/to/file2",
          },
          "test.example.matrix": [
            [0, 1, 2],
            [3, 4, 5],
            [6, 7, 8]
          ]
        }

        Produces this command:
        
        echo a,b,c
        cat /tmp/array.txt
        python script.py --map=<cromwell-execution/path>/map_<hash_id>.txt
        python process.py <cromwell-execution/path>/matrix_<hash_id>.txt

        array.txt would contain
        a
        b
        c

        map.txt would contain
        key0  /path/to/file0
        key1  /path/to/file1
        key2  /path/to/file2

        matrix.txt would contain
        0 1 2
        3 4 5
        6 7 8

        use read_* functions go to from files output by your command into WDL values

        task example {
          command {
            echo 'first' > file
            echo 'second' >> file
            echo 'third' >> file
          }
          output {
            Array[String] out = read_lines("file")
          }
        }
    
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



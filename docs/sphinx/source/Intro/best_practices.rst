================================
Best Practices for Creating WDLs
================================

.. role:: listsize
.. role:: textborder
.. role:: bash(code)

There are opportunities to participate in code reviews with other WDL developers. `ContactUs <contact_us.html>`_ 

----------------------

.. raw:: html

   <font class="listsize";>1) organize workflows into tasks</font>

   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    Workflows are easily converted to WDLs when they are composed of self contained tasks. Each task should be runnable on its own with explicitly defined inputs and outputs. Creating WDLs are trivial when the original workflows are organized thus.
    </p>

   </details>

   <br><font class="listsize";>2) set -euo pipefail</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    This command should be used at the begining of the command{} section in your WDL. This command will help capture errors at the point where they occur in your unix code, rather than having the commands run beyond where the error happened, since this makes debugging more difficult. However, the 'set -euo pipefail' command can cause strange crashes with no error. (e.g. this command returns a non zero error code even if it is successful, java -Xmx1536m -jar /opt/omics/bin/CRT-CLI.jar -version), so it is not always appropriate to use.
   </p>

   </details>

   <br><font class="listsize">3) Use Docker containers</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    The running environment and required scripts should be encapsulated in a docker image which is pused to hub.docker.com and have a versioned Dockerfile. Once testing is done, avoid commands that use shifter or singularity, but instead refer to images by adding them to the runtime{} section of the WDL. 
    JAWS makes multiple computing resources available, using various linux distros.  Thus, we recommend that a docker container be specified for every task; if not, the default container is Debian.
    <br><br>
    Save all versions of your docker images at hub.docker.com. Jaws will pull images from there by default.
   </p>
   
   Example

.. code-block:: text

    # You should be able to do this inside your WDL (cori example).
    "docker: bryce911/bbtools:38.86" 

    # instead of
    shifter --image=bryce911/bbtools:38.02 filterFastq.py in=${fastq}

.. raw:: html

   </details>

.. raw:: html

   <br><font class="listsize">4) Provide everything required for a test run</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    The developer should provide a git repository with testing data and instructions so the JAWS team can validate that the pipeline works in JAWS before it is “ready for production use.”  The README.md should include information on the testing environment (i.e Cromwell version), testing data set, and execution time and output size.
   </p>
   
   </details>


   <br><font class="listsize">5) Avoid hard-coding paths in the WDL</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    File/directory paths should be put into the inputs.json file, not the WDL. The exeption to this rule are docker images which should be hard-coded so the WDL contains information about the version of the docker container.
   </p>
   
   </details>


   <br><font class="listsize">6) WDL tasks should be self-sufficient</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    Imagine the WDL task as a wrapper script, it should be able to run independently of the pipeline. This means that a script should explicitly list all required input files as arguments and not assume some input files already exist in the current working directory. Scripts should also specify output files as arguments and shouldn't write them somewhere other than the current working directory if they will be needed for the next task. These rules make writing the WDL trivial.
    <br><br>
    The WDL should be expected to handle minimal logic.  Use wrapper scripts to deal with logic if need be.

    Also, scripts should return a code of 0 if it was successfull. And don't write anything but errors to stderr. Cromwell depends on seeing a return code of 0 on success and JAWS depends on seeing errors written to stderr. Sometimes, scripts write errors to stdout and these will be missed if you try and see the errors via running the JAWS command (jaws errors).
   </p>
   
   Example

.. code-block:: text

    # This explicitly lists all input files, and output file.
    filterFastq.py in=${fastq} ref=${refdata} huseq=${hu_fasta} out=myout

    # This script expects the files to exist implicitly
    filterFastq.py ref=${refdata} 

.. raw:: html

    </details>


.. raw:: html

   <br><font class="listsize">7) Use subworkflows</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
   Consider using sub-workflows if grouping tasks makes the complete workflow more understandable, reusable, and maintainable. Even a task can be its own workflow.
   <br>
    You use sub-workflows tasks as if they were regular tasks <br>(example copied from https://cromwell.readthedocs.io/en/stable/SubWorkflows/).
   </p>
   
   Example


.. code-block:: text
   
    # main.wdl
    
    import "sub_wdl.wdl" as sub

    workflow main_workflow {

        call sub.hello_and_goodbye { input: hello_and_goodbye_input = "sub world" }

        # call myTask { input: hello_and_goodbye.hello_output }

        output {
            String main_output = hello_and_goodbye.hello_output
        }
    }
    

.. code-block:: text
    
    # sub_wdl.wdl

    workflow hello_and_goodbye {
    String hello_and_goodbye_input

    call hello {input: addressee = hello_and_goodbye_input }
    call goodbye {input: addressee = hello_and_goodbye_input }

    output {
        String hello_output = hello.salutation
        String goodbye_output = goodbye.salutation
      }
    }
  
    task hello {
        String addressee
        command {
            echo "Hello ${addressee}!"
        }
        output {
            String salutation = read_string(stdout())
        }
    }

    task goodbye {
        String addressee
        command {
            echo "Goodbye ${addressee}!"
        }
        output {
            String salutation = read_string(stdout())
        }
    }

.. raw:: html

   </details>

.. raw:: html

   <br><font class="listsize">8) Documenting your WDLs</font>
   <details>
   <summary style="color: #448ecf";>expand</summary>

   <p class="textborder">
    The best way to document your WDLs is with a README.md that is in the same repository as the WDL. However, adding "metadata" sections in the WDL is also best practice since you will hard-code some relevant information this way, like author, contact info, etc.  See the WDL template as an example.
   </p>
   
.. raw:: html

   </details>
|

|

=========
Templates
=========


.. raw:: html

    <font class="listsize">WDL Best Practices Template</font>
    <details>
    <summary style="color: #448ecf";>example</summary>

.. code-block:: text

    # By versioning your WDL, you specify which specification cromwell uses to decifer the WDL.
    # New features come with new versions.
    version 1.0 
    
    # import any subworkflows
    import "subworkflow.wdl" as firstStep
    
    workflow bbtools {
        meta {
		    developer: "Jackson Brown jbrown@my-inst"
			institution: "JGI"
			version: "2222.2.0"
			notes: "this is the official release version"
        }
    
        # you must have this input section within the "workflow" stanza if you are using version 1
        input {
            File reads
            File ref
            String bbtools_docker = "jfroula/bbtools:1.0.4"
        }
    
        call firstStep {
          input: fastq=reads,
                 container=bbtools_docker
        }
        
        call alignment {
           input: fastq=reads,
                  fasta=ref,
                  container=bbtools_docker
        }
    
        call samtools {
           input: sam=alignment.sam
       }
    }
    
    #
    # below are task definitions
    #
    task alignment {
        # Metadata is good for helping the next guy understand your code. 
        # This meta section can also be used for documentation generated by wdl-aid.
        # You can run "wdl-aid <workflow.wdl>" if it is installed, see https://wdl-aid.readthedocs.io/en/latest/usage.html)
        meta {
            metaParameter1: "Some meta Data I"
            metaParameter2: "Some meta Data II"
            description: "Add a brief description of what this task does in this optional block. One can add as much text as one wants in this section to inform an outsider to understand the mechanics of this task."
        }
    
        input {
            File fastq
            File fasta
        }
    
        command {
            # Use this command to help debug your bash code (i.e. prevents hidden bugs).
            # For a description, see https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
            set -euo pipefail
    
            # Note that ~{} is prefered over the old ${} syntax
            bbmap.sh in=~{fastq} ref=~{fasta} out=test.sam
        }
        
        runtime {
            docker: "jfroula/bbtools:1.0.4"
            time: "12:00:00"      
            poolname: "medium"    
            shared: 0         
            constraint: "haswell"
            nodes: 1
            nwpn: 1
        }
    
        output {
           File sam = "test.sam"
        }
    
        # This section is optional and used to create documentation using the wdl-aid tool. 
        # see https://wdl-aid.readthedocs.io/en/latest/usage.html
        # You can run "wdl-aid <workflow.wdl>" if it is installed.
        parameter_meta {
            WDL_AID: {
              exclude: ["input_name", "call.input_name"]
            }
            fastq: {description: "henryInputFile Description", category: "advanced"}
            fasta: {description: "henryInputFile Description", category: "advanced"}
            dockerImage:    {description: "dockerImage Description", category: "advanced"}
        }
        
    }

.. raw:: html

    </details>

|

.. raw:: html

    <font class="listsize">Dockerfile template</font>
    <details>
    <summary style="color: #448ecf";>example</summary>

.. code-block:: text

    FROM ubuntu:16.04

    # install stuff with apt-get
    RUN apt-get update && apt-get install -y wget bzip2
    
    # install miniconda
    # There is a good reason to install miniconda in a path other than its default.  
    # The default intallation directory is /root/miniconda3 but this path will not be 
    # accessible by shifter or singularity so we'll install under /usr/local/bin/miniconda3.
    RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh \
    && bash ./Miniconda3-4.5.11-Linux-x86_64.sh -b -p /usr/local/bin/miniconda3 \
    && rm Miniconda3-4.5.11-Linux-x86_64.sh
    
    # point to all the future conda installations you are going to do
    ENV PATH=/usr/local/bin/miniconda3/bin:$PATH
    
    # Install stuff with conda
    # Remember to use versions of everything you install with conda as shown in example.
    RUN conda install -y -c bioconda bowtie2=2.3.4.3
    RUN conda install -y -c anaconda biopython=1.72
    
    # copy bash/python scripts specific to your pipeline
    COPY scripts/* /usr/local/bin/

.. raw:: html

    </details>

|
|

Additional helpful notes when building Docker images:
-----------------------------------------------------

* The dockerfile template uses the strategy of installing miniconda so you can use :bash:`conda install` for probably, most of your tools.  However, :bash:`pip install` and :bash:`apt-get install` work in addition to, or instead of miniconda.

* Also, remember to use versions of everything you install with conda as shown in example.

* There is a good reason to install miniconda in a path other than its default.  The default installation directory is :bash:`/root/miniconda3` but this path will not be accessible by shifter or singularity.

* When you build your docker (i.e. :bash:`docker build --tag <somename> -f ./Dockerfile3 .`) run this in a CLEAN directory with only the essential files in there because everything in your local dir will become part of the image.

* One helpful thing you can do when developing docker images is to create a bare essentials image with your favorite editor installed (i.e. vim). Then you can go into the container interactively :bash:`docker run --it <image>` and see if you can install stuff manually, then just copy those same commands into the final dockerfile.


For more see the docker official docs on `best practices <https://docs.docker.com/develop/develop-images/dockerfile_best-practices/>`_


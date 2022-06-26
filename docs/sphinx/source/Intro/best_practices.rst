================================
Best Practices for Creating WDLs
================================

.. role:: listsize
.. role:: textborder
.. role:: bash(code)

There are opportunities to participate in code reviews with other WDL developers. `ContactUs <contact_us.html>`_ 

----------------------

.. raw:: html

   <details>
   <summary style="color: #448ecf";>set -euo pipefail</summary>

   <p class="textborder">
    The <b>set -euo pipefail</b> command is actually a composition of three tests.
    <br>
    <br>For example:
    <br>use <b>set -e</b> to make your script exit when a command fails.
    <br>use <b>set -u</b> to exit when your script tries to use undeclared variables.
    <br>use <b>set -o pipefail</b> in scripts to catch failures in "cat myfile" in e.g. "cat myfile | grep id". Instead of the successful error code from grep id getting returned, we get a non-zero exit code from cat myfile
    <br>use <b>set -x</b> to trace what gets executed. Useful for debugging.

    <br><br>
    This command can be useful when used at the begining of the command{} section in your WDL. This command will help capture errors at the point where they occur in your unix code, rather than having the commands run beyond where the error happened, since this makes debugging more difficult.  Another way of saying it is that, without set -e, the wdl-task will use the error code from the last command even if an ealier command failed.  However, the <b>set -euo pipefail</b> command can cause the task to exit without any error printed stderr, so it is not always appropriate to use.
   </p>

   </details>

   <details>
   <summary style="color: #448ecf";>Use Docker containers with SHA256 instead of tags</summary>

   <p class="textborder">
    <br>1. The running environment and required scripts should be encapsulated in a docker image. 
    <br>2. The image should be pushed to hub.docker.com and have a versioned Dockerfile. JAWS will pull images from there by default. 
    <br>3. We recommend that a docker container be specified for every task; if not, the default container is ubuntu.
    <br>4. It is important to reference containers by their SHA256 instead of tag (e.g. doejgi/bbtools@sha256:64088.. instead of doejgi/bbtools:latest) for both reproducability (a container can change and have the same tag) and because call-caching only works when the container is referenced by SHA256 version.
   </p>
   
   SHA Example

.. code-block:: text

    # call-caching will not work
    runtime { "docker: ubuntu:20.04" }

    # call-caching will work
    runtime { "docker: ubuntu@sha256:47f14534bda344d9fe6ffd6effb95eefe579f4be0d508b7445cf77f61a0e5724" }

    # find the sha
    docker pull ubuntu:20.04
    Digest: sha256:47f14534bda344d9fe6ffd6effb95eefe579f4be0d508b7445cf77f61a0e5724

    # or 
    docker inspect --format='{{.RepoDigests}}' ubuntu:20.04
    ubuntu@sha256:47f14534bda344d9fe6ffd6effb95eefe579f4be0d508b7445cf77f61a0e5724

.. raw:: html

   </details>

   <details>
   <summary style="color: #448ecf";>Avoid hard-coding paths in the WDL</summary>

   <p class="textborder">
    Paths to files or directories should be put into the inputs.json file, not the WDL. The exeption to this rule are docker images which <i>should</i> be hard-coded so the WDL contains information about the version of the docker container.
   </p>
   
   </details>

   <details>
   <summary style="color: #448ecf";>WDL tasks should be self-sufficient</summary>

   <p class="textborder">
    <br>1. Imagine the WDL task as a wrapper script, it should be able to run independently of the pipeline. This means that a script should explicitly list all required input files as arguments and not assume some input files already exist in the current working directory. 
    <br>2. Scripts should also specify output files as arguments and shouldn't write them somewhere other than the current working directory if they will be needed for the next task. These rules make writing the WDL trivial.
    <br>3. The WDL should be expected to handle minimal logic.  Use wrapper scripts to deal with logic if need be.
    <br>4. Also, scripts should return a code of 0 if it was successfull. And don't write anything but errors to stderr. Cromwell depends on seeing a return code of 0 on success and JAWS depends on seeing errors written to stderr. Sometimes, scripts write errors to stdout and these will be missed if you try and see the errors via running the JAWS command (jaws errors).
   </p>
   
   Example

.. code-block:: text

    # This explicitly lists all input files, and output file.
    filterFastq.py in=${fastq} ref=${refdata} huseq=${hu_fasta} out=myout

    # This script expects the files to exist implicitly
    filterFastq.py ref=${refdata} 

.. raw:: html

    </details>

   <details>
   <summary style="color: #448ecf";>Use subworkflows</summary>

   <p class="textborder">
   Consider using subworkflows if organizing tasks that way makes the main workflow more understandable, reusable, and maintainable. Even a single task can be its own workflow.
   <br>
    Subworkflows are imported and used as if they were normal tasks, see the example below that was copied from https://cromwell.readthedocs.io/en/stable/SubWorkflows/.
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

   <details>
   <summary style="color: #448ecf";>Documenting your WDLs</summary>

   <p class="textborder">
    The best way to document your WDLs is with a README.md that is in the same repository as the WDL. However, adding "metadata" sections in the WDL is also best practice since you will hard-code some relevant information this way, like author, contact info, etc.  See the WDL template as an example.
   </p>
   
.. raw:: html

   </details>
|

Templates
-----------------------------------------------------

.. raw:: html

    <details>
    <summary style="color: #448ecf";>WDL Best Practices Template</summary>

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

    <details>
    <summary style="color: #448ecf";>Dockerfile template</summary>

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

* Also, remember to use versions of everything you install with conda as shown in above docker template example.

* There is a good reason to install miniconda in a path other than its default.  The default installation directory is :bash:`/root/miniconda3` but this path will not be accessible by shifter or singularity.

* When you build your docker (i.e. :bash:`docker build --tag <somename> -f ./Dockerfile3 .`) run this in a CLEAN directory with only the essential files in there because everything in your local dir will become part of the image.

* One helpful thing you can do when developing docker images is to create a bare essentials image with your favorite editor installed (i.e. vim). Then you can go into the container interactively :bash:`docker run --it <image>` and see if you can install stuff manually, then just copy those same commands into the final dockerfile.


For more see the docker official docs on `best practices <https://docs.docker.com/develop/develop-images/dockerfile_best-practices/>`_


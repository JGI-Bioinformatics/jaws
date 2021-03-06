################
Backend Examples
################


.. _singularity_backend:

Singularity Backend Example
---------------------------

From file on lawrencium
Cromwell_NERSC_dev.conf

.. code-block:: text

    JTM
    {
      actor-factory = "Cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

      config {
        runtime-attributes = """
          String? docker
          String time = "00:00:00"
          Int cpu = 1
          Float? memory_gb = 5
          Int node = 1
          Int nwpn = 1
        """

        submit = "jtm submit \
          -cr '/bin/bash ${script}' \
          -cl ${cluster} \
          -t ${time} \
          -c ${cpu} \
          -m ${memory_gb} \
          -p ${poolname} \
          -C ${constraint} \
          -N ${node} \
          -nwpn ${nwpn} \
          -jid ${job_name} \
          --qos ${qos} \
          -A ${account}"
        kill = "jtm kill -tid ${job_id}"
        check-alive = "jtm isalive -tid ${job_id}"

        job-id-regex = "JTM task ID (\\d+)"

        # Submit string when there is a "docker" runtime attribute.
        submit-docker = """
            LOOKUP=$(shifterimg lookup ${docker})
            if [[ ! $LOOKUP ]]; then
                shifterimg pull ${docker}
            fi

		    jtm submit
              -cr 'shifter_exec.sh ${docker} ${job_shell} ${script}' \
              -cl ${cluster} \
              -t ${time} \
              -c ${cpu} \
              -m ${memory_gb} \
              -p ${poolname} \
              -C ${constraint} \
              -N ${node} \
              -nwpn ${nwpn} \
              -jid ${job_name} \
              --qos ${qos} \
              -A ${account}"
        """

        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above AND in the generation of the "script" file.
        dockerRoot = /<working_dir_for_Cromwell>/Cromwell-executions
    }


If your wdl task is being run in a docker image, then the "submit-docker" command is run, otherwise the "submit" command is run.  Note that "submit-docker" runs "singularity_exec.sh"

**singularity_exec.sh**

.. code-block:: text

	export SINGULARITY_CACHEDIR=/<somepath>/sif_files
	export SINGULARITY_PULLFOLDER=/<somepath>/sif_files
	export SINGULARITY_TMPDIR=/<somepath>/sif_files
	export SINGULARITY_LOCALCACHEDIR=/<somepath>/sif_files
	singularity exec --bind $1:$2 --bind /<somepath_to_db>:/refdata docker://$3 $4 $5


.. _shifter_backend:

Shifter Backend Example
---------------------------

From file on cori
Cromwell.conf

.. code-block:: text

	# this is required for shifter to find image from its registry.
	docker {
		hash-lookup {
			enabled = false
		}
	}

    JTM
    {
      actor-factory = "Cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

      config {
        runtime-attributes = """
          String? docker
          String time = "00:00:00"
          Int cpu = 1
          Float? memory_gb = 5
          String cluster = "cori"
          String poolname = "small"
          String constraint = "haswell"
          String qos = "genepool_special"
          String account = "fungalp"
          Int node = 1
          Int nwpn = 1
        """

        submit = "jtm submit \
          -cr '/bin/bash ${script}' \
          -cl ${cluster} \
          -t ${time} \
          -c ${cpu} \
          -m ${memory_gb} \
          -p ${poolname} \
          -C ${constraint} \
          -N ${node} \
          -nwpn ${nwpn} \
          -jid ${job_name} \
          --qos ${qos} \
          -A ${account}"
        kill = "jtm kill -tid ${job_id}"
        check-alive = "jtm isalive -tid ${job_id}"
        job-id-regex = "JTM task ID (\\d+)"

        # Submit string when there is a "docker" runtime attribute.
        submit-docker = """
            LOOKUP=$(shifterimg lookup ${docker})
            if [[ ! $LOOKUP ]]; then
                shifterimg pull ${docker}
            fi

            jtm submit \
              -cr 'shifter_exec.sh ${docker} ${job_shell} ${script}' \
              -cl ${cluster} \
              -t ${time} \
              -c ${cpu} \
              -m ${memory_gb} \
              -p ${poolname} \
              -C ${constraint} \
              -N ${node} \
              -nwpn ${nwpn} \
              -jid ${job_name} \
              --qos ${qos} \
              -A ${account}"
        """

        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above AND in the generation of the "script" file.
        dockerRoot = <working_dir_for_Cromwell>/Cromwell-executions
      }
    }

If your wdl task is being run in a docker image, then the "submit-docker" command is run, otherwise the "submit" command is run.  Note that "submit-docker" runs "shifter_exec.sh"

**shifter_exec.sh**

.. code-block:: text

	#!/bin/bash
	shifter --image=$1 -V <full_path_to_db>:/refdata $2 $3

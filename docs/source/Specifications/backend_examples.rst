################
Backend Examples
################


.. _singularity_backend:

Singularity Backend Example
---------------------------

From file on lawrencium
/global/home/groups-sw/lr_jgicloud/dev/jtm/opt/cromwell/cromwell_nersc_dev.conf

.. code-block:: bash

    JTM
    {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      # this is required for shifter to find image from its registry.
      docker {
        hash-lookup {
          enabled = false
        }
      }

      config {
        runtime-attributes = """
        String? docker
        String time = "00:00:00"
        Int cpu = 1
        String mem = "0G"
            String cluster = "cori"
            String poolname = "small"
            #Int poolsize = 1
        String constraint = "haswell"
        Int node = 1
        Int nwpn = 1
        Int shared = 1
        """

        submit = "jtm-submit -cr '/bin/bash ${script}' -cl ${cluster} -t ${time} -c ${cpu} -m ${mem} -p ${poolname} -C ${constraint} -N ${node} -nwpn ${nwpn} -jid ${job_name} --shared ${shared}"
        kill = "jtm-kill ${job_id}"
        check-alive = "jtm-isalive ${job_id}"

        job-id-regex = "JTM task ID (\\d+)"

        # Submit string when there is a "docker" runtime attribute.
        submit-docker = """
            LOOKUP=$(shifterimg lookup ${docker})
            if [[ ! $LOOKUP ]]; then
                shifterimg pull ${docker}
            fi

		jtm-submit -cr 'shifter_exec.sh ${docker} ${job_shell} ${script}' \
                -cl ${cluster} -t ${time} -c ${cpu} -m ${mem} -p ${poolname} -C ${constraint} \
            -N ${node} -nwpn ${nwpn} -jid ${job_name} --shared ${shared}
        """

        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above AND in the generation of the "script" file.
        dockerRoot = /global/cscratch1/sd/jaws_jtm/dev/cromwell-executions
    }
 }





.. _shifter_backend:

Shifter Backend Example
---------------------------

From file on cori
/global/project/projectdirs/jaws_jtm/dev/etc/cromwell.conf

.. code-block:: bash

    JTM
    {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      # this is required for shifter to find image from its registry.
      docker {
        hash-lookup {
          enabled = false
        }
      }

      config {
        runtime-attributes = """
        String? docker
        String time = "00:00:00"
        Int cpu = 1
        String mem = "0G"
        String cluster = "cori"
        String poolname = "small"
        #Int poolsize = 1
        String constraint = "haswell"
        String qos = "genepool_special"
        String account = "fungalp"
        Int node = 1
        Int nwpn = 1
        Int shared = 1
        """

        submit = "source /global/project/projectdirs/jaws_jtm/jtm/venv/bin/activate && jtm-submit -cr '/bin/bash ${script}' -cl ${cluster} -t ${time} -c ${cpu} -m ${mem} -p ${poolname} -C ${constraint} -N ${node} -nwpn ${nwpn} -jid ${job_name} --shared ${shared} --qos ${qos} --account ${account}"
        kill = "jtm-kill ${job_id}"
        check-alive = "jtm-isalive ${job_id}"
        job-id-regex = "JTM task ID (\\d+)"

        # Submit string when there is a "docker" runtime attribute.
        submit-docker = """
            LOOKUP=$(shifterimg lookup ${docker})
            if [[ ! $LOOKUP ]]; then
                shifterimg pull ${docker}
            fi

            jtm-submit -cr 'shifter_exec.sh ${docker} ${job_shell} ${script}' \
            -cl ${cluster} -t ${time} -c ${cpu} -m ${mem} -p ${poolname} -C ${constraint} \
            -N ${node} -nwpn ${nwpn} -jid ${job_name} --shared ${shared}
        """

        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above AND in the generation of the "script" file.
        dockerRoot = /global/cscratch1/sd/jaws_jtm/dev/cromwell-executions
      }
    }

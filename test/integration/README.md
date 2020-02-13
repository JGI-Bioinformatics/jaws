# JAWS Integration Tests

The purpose of the integration tests is to validate the functionality of JAWS, from cradle to grave.

To do so, the pipeline spins up a Slurm cluster, using the [ElastiCluster](https://github.com/elasticluster/elasticluster) tool.
It then goes on to install Docker on the head-node of the cluster. After this the deploy step can take necessary steps to
deploy all the necessary components and run the integration test suite against JAWS. Once the pipeline is done, it destroys the AWS instances.

Important: Should the destroy step not run, eg because it was interrupted or skipped, the cluster will _not_ be destroyed and will
keep on incurring cost. If need arises it can be destroyed manually by setting $CI_PROJECT_NAME and $CI_PIPELINE_ID to the corresponding
values of the failed pipeline and executing the destroy_cluster script.

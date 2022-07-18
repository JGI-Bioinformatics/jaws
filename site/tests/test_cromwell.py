import json
import os
import glob
from deepdiff import DeepDiff
from jaws_site import cromwell

tests_dir = os.path.dirname(os.path.abspath(__file__))

example_cromwell_url = "http://localhost:8000"
example_workflows_url = f"{example_cromwell_url}/api/workflows/v1"
crom = cromwell.Cromwell(example_cromwell_url)

# simple successful run
example_cromwell_run_id_1 = "ee30d68f-39d4-4fde-85c2-afdecce2bad3"

# successful run with subs
example_cromwell_run_id_2 = "c720836c-0931-4ddc-8366-774160e05531"

# simple failed run
example_cromwell_run_id_3 = "dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15"

# failed, missing infile
example_cromwell_run_id_4 = "469bcdd6-d67f-455f-9475-438349f41631"

# successful run with subs
example_cromwell_run_id_5 = "5c7ca4aa-6bc2-4f90-a80c-6232ee64e1f0"

# failed simple run
example_cromwell_run_id_6 = "cbbbc75f-8920-495c-a290-0a1a5f0d1c20"

# failed run, with subs (sub's docker container not found)
example_cromwell_run_id_7 = "dcb85c55-bf74-4f63-bce3-fe61f7b84ebb"

# successful run, with outputs list
example_cromwell_run_id_8 = "5a2cbafe-56ed-42aa-955d-fef8cb5014bb"

# successful AWS run
example_cromwell_run_id_9 = "f4f5afd1-79f5-497a-9612-baed76dc365d"

# running Run
example_cromwell_run_id_10 = "dcc24ca7-c303-4e8e-ad26-7b2644308fab"


def __load_example_output_from_file(cromwell_run_id, output_type):
    with open(f"{tests_dir}/{cromwell_run_id}.{output_type}.json", "r") as fh:
        output = json.load(fh)
    return output


def test_get_metadata(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expectedWorkflowRoot = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # noqa
    metadata1 = crom.get_metadata(example_cromwell_run_id_1)
    assert expectedWorkflowRoot == metadata1.get("workflowRoot")


def test_workflow_root(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_9}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_9, "metadata"),
    )
    expected = "s3://jaws-site/cromwell-execution/jgi_meta/f4f5afd1-79f5-497a-9612-baed76dc365d"
    metadata = crom.get_metadata(example_cromwell_run_id_9)
    actual = metadata.workflow_root()
    assert expected == actual


def test_task(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_1)
    task = metadata.tasks["fq_count.count_seqs"]
    shard_index = -1  # not a sharded task
    attempt = 1  # generally only 1 attempt
    assert task.calls[shard_index][attempt].execution_status == "Done"


def test_task_subworkflow(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata2 = crom.get_metadata(example_cromwell_run_id_2)
    task = metadata2.tasks["main_workflow.hello_and_goodbye_1"]
    shard_index = -1  # not a sharded task
    attempt = 1  # generally only 1 attempt
    taskMetadata = task.subworkflows[shard_index][attempt]
    assert taskMetadata.get("rootWorkflowId") == "c720836c-0931-4ddc-8366-774160e05531"


def test_task_summary(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    actual = metadata.task_summary()
    expected = __load_example_output_from_file(
        example_cromwell_run_id_2, "task-summary"
    )
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_5}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_5, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_5)
    actual = metadata.task_summary()
    expected = __load_example_output_from_file(
        example_cromwell_run_id_5, "task-summary"
    )
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


def test_task_log(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    actual = metadata.task_log()
    expected = __load_example_output_from_file(example_cromwell_run_id_2, "task-log")
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


def test_task_stdout(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_1)
    task = metadata.tasks["fq_count.count_seqs"]
    shard_index = -1  # not a sharded task
    attempt = 1  # generally only 1 attempt
    call = task.calls[shard_index][attempt]
    expected_stderr = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr"  # noqa
    expected_stdout = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout"  # noqa
    assert call.stderr == expected_stderr
    assert call.stdout == expected_stdout


def test_metadata_tasks(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )

    metadata = crom.get_metadata(example_cromwell_run_id_2)
    for task_name, task in metadata.tasks.items():
        assert len(task.data) > 0
        call = task.data[0]
        shard_index = call["shardIndex"]
        attempt = call["attempt"]
        if "subWorkflowMetadata" in call:
            subworkflow = task.subworkflows[shard_index][attempt]
            assert isinstance(subworkflow, cromwell.Metadata)
        else:
            assert "jobId" in call


def test_errors(requests_mock, monkeypatch):
    def mock_read_file(path):
        return None

    monkeypatch.setattr(cromwell, "_read_file", mock_read_file)

    # simple workflow with no errors
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_errors_report_1 = {}
    metadata_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_errors_report_1 = metadata_1.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_1, expected_errors_report_1, ignore_order=True
            )
        )
        is False
    )

    # simple workflow with errors
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_3}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_3, "metadata"),
    )
    expected_errors_report_3 = __load_example_output_from_file(
        example_cromwell_run_id_3, "errors"
    )
    metadata_3 = crom.get_metadata(example_cromwell_run_id_3)
    actual_errors_report_3 = metadata_3.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_3, expected_errors_report_3, ignore_order=True
            )
        )
        is False
    )

    # failed, missing infile
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_4}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_4, "metadata"),
    )
    expected_errors_report_4 = __load_example_output_from_file(
        example_cromwell_run_id_4, "errors"
    )
    metadata_4 = crom.get_metadata(example_cromwell_run_id_4)
    actual_errors_report_4 = metadata_4.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_4, expected_errors_report_4, ignore_order=True
            )
        )
        is False
    )

    # failed simple run
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_6}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_6, "metadata"),
    )
    expected_errors_report_6 = __load_example_output_from_file(
        example_cromwell_run_id_6, "errors"
    )
    metadata_6 = crom.get_metadata(example_cromwell_run_id_6)
    actual_errors_report_6 = metadata_6.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_6, expected_errors_report_6, ignore_order=True
            )
        )
        is False
    )

    # failed run with subworkflows
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_7}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_7, "metadata"),
    )
    expected_errors_report_7 = __load_example_output_from_file(
        example_cromwell_run_id_7, "errors"
    )
    metadata_7 = crom.get_metadata(example_cromwell_run_id_7)
    actual_errors_report_7 = metadata_7.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_7, expected_errors_report_7, ignore_order=True
            )
        )
        is False
    )


def test_running(requests_mock, monkeypatch):
    def mock_read_file(path):
        return None

    monkeypatch.setattr(cromwell, "_read_file", mock_read_file)

    def mock_glob(path):
        return []

    monkeypatch.setattr(glob, "glob", mock_glob)

    # completed workflow
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_running_report_1 = {}
    metadata_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_running_report_1 = metadata_1.running()
    assert (
        bool(
            DeepDiff(
                actual_running_report_1, expected_running_report_1, ignore_order=True
            )
        )
        is False
    )

    # running workflow
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_10}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_10, "metadata"),
    )
    expected_running_report_10 = __load_example_output_from_file(
        example_cromwell_run_id_10, "running"
    )
    metadata_10 = crom.get_metadata(example_cromwell_run_id_10)
    actual_running_report_10 = metadata_10.running()
    print(actual_running_report_10)
    assert (
        bool(
            DeepDiff(
                actual_running_report_10, expected_running_report_10, ignore_order=True
            )
        )
        is False
    )


def test_get_outputs(requests_mock):
    # test 1 : outputs scalar
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_outputs_1 = {
        "fq_count.outfile": "./call-count_seqs/execution/num_seqs.txt"
    }
    ex_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_outputs_1 = ex_1.outputs(relpath=True)
    assert (
        bool(DeepDiff(actual_outputs_1, expected_outputs_1, ignore_order=True)) is False
    )

    # test 2 : outputs list
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_8}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_8, "metadata"),
    )
    expected_outputs_8 = {
        "create_lastdb.dbFiles": [
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.bck",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.des",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.prj",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.sds",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.ssp",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.suf",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.tis",
        ]
    }
    ex_8 = crom.get_metadata(example_cromwell_run_id_8)
    actual_outputs_8 = ex_8.outputs(relpath=True)
    assert (
        bool(DeepDiff(actual_outputs_8, expected_outputs_8, ignore_order=True)) is False
    )


def test_outfiles(requests_mock):
    # test 1 : outputs scalar
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_outfiles_1 = ["./call-count_seqs/execution/num_seqs.txt"]
    ex_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_outfiles_1 = ex_1.outfiles(relpath=True)
    assert (
        bool(DeepDiff(actual_outfiles_1, expected_outfiles_1, ignore_order=True))
        is False
    )

    # test 2 : outputs list
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_8}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_8, "metadata"),
    )
    expected_outfiles_8 = [
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.bck",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.des",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.prj",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.sds",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.ssp",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.suf",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.tis",
    ]
    ex_8 = crom.get_metadata(example_cromwell_run_id_8)
    actual_outfiles_8 = ex_8.outfiles(relpath=True)
    assert (
        bool(DeepDiff(actual_outfiles_8, expected_outfiles_8, ignore_order=True))
        is False
    )

    # test 9 : aws outputs
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_9}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_9, "metadata"),
    )
    expected_outfiles_9 = [
        "./call-create_agp/cacheCopy/assembly.contigs.fasta",
        "./call-read_mapping_pairs/covstats.txt",
        "./call-read_mapping_pairs/pairedMapped.sam.gz",
        "./call-bbcms/cacheCopy/resources.log",
        "./call-create_agp/cacheCopy/assembly.scaffolds.fasta",
        "./call-bbcms/cacheCopy/input.corr.fastq.gz",
        "./call-assy/resources.log",
        "./call-read_mapping_pairs/resources.log",
        "./call-assy/spades3/spades.log",
        "./call-create_agp/cacheCopy/resources.log",
        "./call-read_mapping_pairs/pairedMapped_sorted.bam.bai",
        "./call-bbcms/cacheCopy/counts.metadata.json",
        "./call-bbcms/cacheCopy/stdout.log",
        "./call-assy/spades3/scaffolds.fasta",
        "./call-bbcms/cacheCopy/input.corr.left.fastq.gz",
        "./call-bbcms/cacheCopy/stderr.log",
        "./call-bbcms/cacheCopy/readlen.txt",
        "./call-bbcms/cacheCopy/unique31mer.txt",
        "./call-create_agp/cacheCopy/assembly.scaffolds.legend",
        "./call-bbcms/cacheCopy/input.corr.right.fastq.gz",
        "./call-read_mapping_pairs/pairedMapped_sorted.bam",
        "./call-create_agp/cacheCopy/assembly.agp",
    ]
    ex_9 = crom.get_metadata(example_cromwell_run_id_9)
    actual_outfiles_9 = ex_9.outfiles(relpath=True)
    assert (
        bool(DeepDiff(actual_outfiles_9, expected_outfiles_9, ignore_order=True))
        is False
    )


def test_job_summary(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    actual = metadata.job_summary()
    expected = {
        "12129": "main_workflow.goodbye",
        "12130": "main_workflow.hello",
        "12131": "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
        "12132": "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
        "12133": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
        "12134": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
    }
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


def test_parse_cromwell_task_dir():

    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            {
                "call_root": "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9",  # noqa
                "cached": False,
                "shard": 9,
                "name": "jgi_dap_leo.trimAlign_expt[9]",
                "call_root_rel_path": "jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
                "wdl_name": "jgi_dap_leo",
                "cromwell_run_id": "cda3cb3f-535c-400d-ab61-2e41aeb35a80",
                "task_name": "trimAlign_expt",
            },
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            {
                "call_root": "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello",  # noqa
                "cached": False,
                "shard": -1,
                "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
                "call_root_rel_path": "main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
                "wdl_name": "main_workflow",
                "cromwell_run_id": "e7f02164-2d3d-4cfb-828a-f3da23c43280",
                "task_name": "hello_and_goodbye_1",
                "subworkflow_name": "hello_and_goodbye",
                "subworkflow_cromwell_run_id": "3327f701-769a-49fe-b407-eb4be3a4a373",
                "sub_task_name": "hello",
            },
        ],
        [
            "s3://jaws-site-prod/cromwell-execution/jgi_meta/bbfad8f4-f5de-43c2-94ef-4bd43f1de4d3/call-bbcms",  # noqa
            {
                "call_root": "s3://jaws-site-prod/cromwell-execution/jgi_meta/bbfad8f4-f5de-43c2-94ef-4bd43f1de4d3/call-bbcms",  # noqa
                "cached": False,
                "shard": -1,
                "name": "jgi_meta.bbcms",
                "call_root_rel_path": "jgi_meta/bbfad8f4-f5de-43c2-94ef-4bd43f1de4d3/call-bbcms/execution",  # noqa
                "wdl_name": "jgi_meta",
                "cromwell_run_id": "bbfad8f4-f5de-43c2-94ef-4bd43f1de4d3",
                "task_name": "bbcms",
            },
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/inhomo/8cc6f043-13b1-4e18-b685-b9533e6704cf/call-processTaxon/shard-5",  # noqa
            {
                "call_root": "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/inhomo/8cc6f043-13b1-4e18-b685-b9533e6704cf/call-processTaxon/shard-5",  # noqa
                "cached": False,
                "shard": 5,
                "name": "inhomo.processTaxon[5]",
                "call_root_rel_path": "inhomo/8cc6f043-13b1-4e18-b685-b9533e6704cf/call-processTaxon/shard-5/execution",  # noqa
                "wdl_name": "inhomo",
                "cromwell_run_id": "8cc6f043-13b1-4e18-b685-b9533e6704cf",
                "task_name": "processTaxon",
            },
        ],
        [
            "/tahoma/mscjgi/scratch/jaws-prod/cromwell-executions/nmdc_metag/c16846ff-3a7b-444e-a26b-ce484eb205b5/call-annotation/awf.annotation/e0910a3c-6ba1-43e3-8b4b-d275fb0601fb/call-s_annotate/shard-0/sa.s_annotate/1671df94-89d9-4418-a949-737038f458a0/call-fasta_merge",  # noqa
            {
                "call_root": "/tahoma/mscjgi/scratch/jaws-prod/cromwell-executions/nmdc_metag/c16846ff-3a7b-444e-a26b-ce484eb205b5/call-annotation/awf.annotation/e0910a3c-6ba1-43e3-8b4b-d275fb0601fb/call-s_annotate/shard-0/sa.s_annotate/1671df94-89d9-4418-a949-737038f458a0/call-fasta_merge",  # noqa
                "cached": False,
                "shard": -1,
                "name": "nmdc_metag.annotation:annotation.s_annotate:s_annotate.fasta_merge",
                "call_root_rel_path": "nmdc_metag/c16846ff-3a7b-444e-a26b-ce484eb205b5/call-annotation/awf.annotation/e0910a3c-6ba1-43e3-8b4b-d275fb0601fb/call-s_annotate/shard-0/sa.s_annotate/1671df94-89d9-4418-a949-737038f458a0/call-fasta_merge/execution",  # noqa
                "wdl_name": "nmdc_metag",
                "cromwell_run_id": "c16846ff-3a7b-444e-a26b-ce484eb205b5",
                "task_name": "annotation",
                "subworkflow_name": "s_annotate",
                "subworkflow_cromwell_run_id": "1671df94-89d9-4418-a949-737038f458a0",
                "sub_task_name": "fasta_merge",
                "sub_shard": 0,
            },
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/bbmap_shard_wf/be518d2a-6232-4f50-b6cf-7e1a3a995ad3/call-alignment/shard-0",  # noqa
            {
                "call_root": "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/bbmap_shard_wf/be518d2a-6232-4f50-b6cf-7e1a3a995ad3/call-alignment/shard-0",  # noqa
                "cached": False,
                "shard": 0,
                "name": "bbmap_shard_wf.alignment[0]",
                "call_root_rel_path": "bbmap_shard_wf/be518d2a-6232-4f50-b6cf-7e1a3a995ad3/call-alignment/shard-0/execution",  # noqa
                "wdl_name": "bbmap_shard_wf",
                "cromwell_run_id": "be518d2a-6232-4f50-b6cf-7e1a3a995ad3",
                "task_name": "alignment",
            },
        ],
    ]

    for task_dir, expected in test_data:
        actual = cromwell.parse_cromwell_task_dir(task_dir)
        print(actual, flush=True)
        assert actual == expected


def test_sort_table():
    test_table = [
        ["bbtools.samtools", False, "Queued", "2022-07-11 17:52:37", "", "", "", ""],
        [
            "bbtools.alignment",
            False,
            "Done",
            "2022-07-11 17:51:13",
            "2022-07-11 17:51:14",
            "2022-07-11 17:52:35",
            "0:00:01",
            "0:01:21",
        ],
    ]
    expected = [
        [
            "bbtools.alignment",
            False,
            "Done",
            "2022-07-11 17:51:13",
            "2022-07-11 17:51:14",
            "2022-07-11 17:52:35",
            "0:00:01",
            "0:01:21",
        ],
        ["bbtools.samtools", False, "Queued", "2022-07-11 17:52:37", "", "", "", ""],
    ]
    actual = cromwell.sort_table(test_table, 3)
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False


def test_sort_table_dict():
    test_table = [{"id": 41, "key": "B"}, {"id": 51, "key": "A"}]
    expected = [{"id": 51, "key": "A"}, {"id": 41, "key": "B"}]
    actual = cromwell.sort_table_dict(test_table, "key")
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False

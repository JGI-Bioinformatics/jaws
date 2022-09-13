import os
from jaws_site import perf_metrics

tests_dir = os.path.dirname(os.path.abspath(__file__))


def test_extract_jaws_info():
    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            "cda3cb3f-535c-400d-ab61-2e41aeb35a80",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            "e7f02164-2d3d-4cfb-828a-f3da23c43280",
        ],
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.extract_jaws_info(task_dir)
        assert result == expected


def test_remove_beginning_path():
    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            "cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            "cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
        ],
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.remove_beginning_path(task_dir)
        assert result == expected


def test_parse_perf_metrics_task_dir():
    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            "jgi_dap_leo.trimAlign_expt[9]",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
        ],
        [
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-bbmap_indexing/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.bbmap_indexing",
        ],
        [
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-shard/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.shard",
        ],
        [
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-alignment/shard-0/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.alignment[0]",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/nmdc_metag/9969e560-7f3e-4305-b66e-324d199d3b33/call-annotation/awf.annotation/8ea94de6-e56f-4b4a-9015-25eb998e68cc/call-f_annotate/shard-0/fa.f_annotate/0e16887c-11ec-48d3-a343-27a9876a7c70/call-smart",  # noqa
            "nmdc_metag.annotation:annotation.f_annotate:f_annotate.smart",
        ],
        [
            "s3://jaws-site-prod/cromwell-execution/jgi_meta/bbfad8f4-f5de-43c2-94ef-4bd43f1de4d3/call-bbcms",  # noqa
            "jgi_meta.bbcms",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/inhomo/8cc6f043-13b1-4e18-b685-b9533e6704cf/call-processTaxon/shard-5",  # noqa
            "inhomo.processTaxon[5]",
        ],
        [
            "/tahoma/mscjgi/scratch/jaws-prod/cromwell-executions/nmdc_metag/c16846ff-3a7b-444e-a26b-ce484eb205b5/call-annotation/awf.annotation/e0910a3c-6ba1-43e3-8b4b-d275fb0601fb/call-s_annotate/shard-0/sa.s_annotate/1671df94-89d9-4418-a949-737038f458a0/call-fasta_merge",  # noqa
            "nmdc_metag.annotation:annotation.s_annotate:s_annotate.fasta_merge",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/bbmap_shard_wf/be518d2a-6232-4f50-b6cf-7e1a3a995ad3/call-alignment/shard-0",  # noqa
            "bbmap_shard_wf.alignment[0]",
        ],
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.parse_cromwell_task_dir_name(task_dir)
        assert result == expected

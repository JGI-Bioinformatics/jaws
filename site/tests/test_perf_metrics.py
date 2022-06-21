
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
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-alignment/shard-0/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.alignment[0]"
        ]
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.parse_cromwell_task_dir(task_dir)
        assert result == expected

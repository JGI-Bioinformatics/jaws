import pytest
from jaws_site.perf_metrics_es import extract_jaws_info


@pytest.mark.parametrize(
    "working_dir, expect",
    [
        (
            "/var/udiMount/global/cscratch1/sd/jaws_jtm/\
jaws-prod/cromwell-executions/rna_count_wrapper/a304bf3c-f042-4af4-8cc5-20bfba245dc9\
/call-rna_count_DGE/DGE.rna_count_DGE/\
d04e8fb6-4a37-4f87-8e4e-8a7107523df9/call-diffGeneExp/execution/DGE_files",
            "a304bf3c-f042-4af4-8cc5-20bfba245dc9"
        ),
        (
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/\
cromwell-executions/jgi_dap_leo/f92649f2-31f9-4fc4-b1fd-18ed4de94c12\
/call-trimAlign_expt/shard-37/execution",
            "f92649f2-31f9-4fc4-b1fd-18ed4de94c12"
        ),
        (
            "cromwell-executions/\
                fq_count/8b1669c3-9bca-4047-99b4-e3556b71c4f4/\
                    call-count_seqs/execution",
            "8b1669c3-9bca-4047-99b4-e3556b71c4f4"
        ),
        (
            "/global/cfs/cdirs/jaws/jaws-install/jaws-supervisord-prod",
            "naw"
        ),
        (
            "None",
            "naw"
        )
    ],
)
def test_extract_cromwell_es(working_dir, expect):
    out = extract_jaws_info(working_dir=working_dir)
    assert out == expect

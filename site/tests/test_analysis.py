import json

from jaws_site import analysis


def test_tail_files(log_file):
    logs, is_truncated = analysis.__tail(log_file)
    print(logs)
    assert is_truncated
    assert len(logs) == 1000
    assert logs[0] == "this is line number 2000\n"
    assert logs[-1] == "this is line number 2999\n"


def test_find_rc_failed_files(cromwell_run_dir):
    run_dir = cromwell_run_dir
    out_json = analysis._find_outfiles(run_dir, failed_only=True)
    obs_output = json.dumps(out_json, sort_keys=True, indent=4)
    exp_output = """
{
    "asm_1": {
        "stderr": "This is standard error. This call-asm_1 had an error",
        "stderr.submit": "This is submit stderr from call-asm_1",
        "stdout": "This is standard output from call-asm_1"
    },
    "circularizeAssembly": {
        "stderr": "This is standard error. This call-circularizeAssembly had an error",
        "stdout": "This is standard output from call-circularizeAssembly"
    },
    "filterHighGc": {
        "stderr": "This is standard error. This call-filterHighGc had an error",
        "stderr.submit": "This is submit stderr from call-filterHighGc",
        "stdout": "This is standard output from call-filterHighGc"
    }
}
    """.strip()
    assert obs_output == exp_output

import json

from jaws_site import central_rpc_operations


def testtail_files(log_file):
    logs, is_truncated = central_rpc_operations._tail(log_file)
    print(logs)
    assert is_truncated
    assert len(logs) == 1000
    assert logs[0] == "this is line number 2000\n"
    assert logs[-1] == "this is line number 2999\n"


def test_find_rc_failed_files(cromwell_run_dir):
    run_dir = cromwell_run_dir
    print(f"run_dir={run_dir}")
    out_json = central_rpc_operations._find_outfiles(run_dir, failed_only=True)
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


def test_find_rc_all_files(cromwell_run_dir):
    run_dir = cromwell_run_dir
    print(f"run_dir={run_dir}")
    out_json = central_rpc_operations._find_outfiles(run_dir, failed_only=False)
    obs_output = json.dumps(out_json, sort_keys=True, indent=4)
    exp_output = """
{
    "asm_1": {
        "stderr": "This is standard error. This call-asm_1 had an error",
        "stderr.submit": "This is submit stderr from call-asm_1",
        "stdout": "This is standard output from call-asm_1"
    },
    "asm_2": {
        "stderr": "This is standard error. This call-asm_2 had no errors",
        "stderr.submit": "",
        "stdout": "This is standard output from call-asm_2"
    },
    "asm_3": {
        "stderr": "This is standard error. This call-asm_3 had no errors",
        "stderr.submit": "",
        "stdout": "This is standard output from call-asm_3"
    },
    "callGenes": {
        "stderr": "This is standard error. This call-callGenes had no errors",
        "stderr.submit": "",
        "stdout": "This is standard output from call-callGenes"
    },
    "circularizeAssembly": {
        "stderr": "This is standard error. This call-circularizeAssembly had an error",
        "stdout": "This is standard output from call-circularizeAssembly"
    },
    "createNonMitoReads": {
        "stderr": "This is standard error. This call-createNonMitoReads had no errors",
        "stderr.submit": "",
        "stdout": "This is standard output from call-createNonMitoReads"
    },
    "doHmmSearch": {
        "stderr": "This is standard error. This call-doHmmSearch had no errors",
        "stderr.submit": "",
        "stdout": "This is standard output from call-doHmmSearch"
    },
    "filterHighGc": {
        "stderr": "This is standard error. This call-filterHighGc had an error",
        "stderr.submit": "This is submit stderr from call-filterHighGc",
        "stdout": "This is standard output from call-filterHighGc"
    },
    "getContigsMatchingHmmSearch": {
        "stderr": "This is standard error. This call-getContigsMatchingHmmSearch had no errors",
        "stderr.submit": "",
        "stdout": "This is standard output from call-getContigsMatchingHmmSearch"
    }
}
    """.strip()
    assert obs_output == exp_output

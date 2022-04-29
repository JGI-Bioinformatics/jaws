import jaws_condor.utils


def test_extract_cromwell_id():
    """
    /bin/bash /..../cromwell-executions/x/74a668ad-958e-437e-a941-6e5e23f8716d/call-y/shard-01/execution/script
    ==> extract 74a668ad-958e-437e-a941-6e5e23f8716d
    """
    msg = '/x/y/z/cromwell-executions/1-234-5/74a668ad-958e-437e-a941-6e5e23f8716d/call-y/shard-01/execution/script'
    expected_uuid = '74a668ad-958e-437e-a941-6e5e23f8716d'
    actual_uuid = jaws_condor.utils.extract_cromwell_id(msg)
    assert actual_uuid == expected_uuid

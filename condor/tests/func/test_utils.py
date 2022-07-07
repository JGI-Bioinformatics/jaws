from jaws_condor import utils


def test_extract_cromwell_id():
    """
    /bin/bash /..../cromwell-executions/x/74a668ad-958e-437e-a941-6e5e23f8716d/call-y/shard-01/execution/script
    ==> extract 74a668ad-958e-437e-a941-6e5e23f8716d
    """
    msg = "/x/y/z/cromwell-executions/1-234-5/74a668ad-958e-437e-a941-6e5e23f8716d/call-y/shard-01/execution/script"
    expected_uuid = "74a668ad-958e-437e-a941-6e5e23f8716d"
    actual_uuid = utils.extract_cromwell_id(msg)
    assert actual_uuid == expected_uuid


def test_mem_unit_to_g():
    tests = [["G", 1, 1], ["T", 1, 1024], ["M", 1024, 1], ["K", 1024 * 1024, 1]]
    for unit, mem, expected in tests:
        actual = utils.mem_unit_to_g(unit, mem)
        assert actual == expected


def test_pad_string_path():
    tests = [["11", 8, "00/00/00/11/"], ["11", 2, "11/"], [12345, 6, "01/23/45/"]]
    for my_string, pad_length, expected in tests:
        actual = utils.pad_string_path(my_string, pad_length)
        assert actual == expected

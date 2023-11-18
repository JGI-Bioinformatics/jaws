import pytest
from jaws_rpc import jsonrpc_utils


def test_validate_request():
    valid_requests = [
        {"jsonrpc": "2.0", "method": "hello", "params": {"name": "world"}},
        {"jsonrpc": "2.0", "method": "hello"},  # params is optional
    ]
    invalid_requests = [
        {"method": "hello", "params": {"name": "world"}},  # missing "jsonrpc"
        {"jsonrpc": "2.0", "params": {"name": "world"}},  # missing "method"
        {"jsonrpc": "2.0", "method": "hello", "name": "world"},  # missing "params"
        {},  # empty
        [],  # not a dict
        "hello",  # not a dict
        None,  # undefined
        {"jsonrpc": "2.1", "method": "hello"},  # not jsonrpc 2.0
    ]
    for request in valid_requests:
        assert jsonrpc_utils.is_valid_request(request)
    for request in invalid_requests:
        with pytest.raises(Exception):
            assert jsonrpc_utils.is_valid_request(request)


def test_validate_response():
    valid_responses = [
        {
            "jsonrpc": "2.0",
            "error": {"code": 500, "message": "Server error"},
        },  # correct failure response
        {"jsonrpc": "2.0", "result": {"foo": "bar"}},  # correct success response
        {"jsonrpc": "2.0", "result": 42},  # result may be scalar
        {"jsonrpc": "2.0", "result": None},  # result may be None
    ]
    invalid_responses = [
        {
            "jsonrpc": "2.0",
            "result": 42,
            "error": None,
        },  # cannot have both "error" and "result" keys
        {"jsonrpc": "2.0", "error": {"message": "Invalid input"}},  # missing "code"
        {},  # empty
        [],  # not a dict
        42,  # not a dict
        None,  # undefined
        {
            "jsonrpc": "2.0",
        },  # must have either "error" or "result"
        {"result": {"foo": "bar"}},  # missing "jsonrpc"
        {"jsonrpc": "2.1", "result": 42},  # not jsonrpc 2.0
    ]
    for response in valid_responses:
        assert jsonrpc_utils.is_valid_response(response)
    for response in invalid_responses:
        with pytest.raises(Exception):
            jsonrpc_utils.is_valid_response(response)

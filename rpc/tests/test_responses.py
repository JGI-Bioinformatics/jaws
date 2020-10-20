from jaws_rpc import responses


def test_success():
    response = responses.success({"apple": "red"})
    assert response["result"]["apple"] == "red"


def test_failure_with_custom_exception():
    error = responses.ErrorWithCode(message="test error message", code=999)
    print(f"{error}")
    response = responses.failure(error)
    assert response["error"]["message"] == "test error message"
    assert response["error"]["code"] == 999


def test_failure_with_standard_exception():
    response = responses.failure(Exception("test error message"))
    assert response["error"]["message"] == "test error message"
    assert response["error"]["code"] == 500

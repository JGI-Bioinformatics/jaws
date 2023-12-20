"""Utilities related to JSON-RPC2 specification per: https://www.jsonrpc.org/specification,
except we neither use nor require the "id" field because AMQP provides the correlation_id."""


request_valid_keys = ["jsonrpc", "id", "method", "params"]
response_valid_keys = ["jsonrpc", "id", "error", "result"]
response_error_valid_keys = ["code", "message", "data"]


def is_valid_request(request):
    """Verifies the request dict conforms to JSON-RPC2 spec.  Raises error if invalid.

    :param request: The request object to validate
    :type request: dict
    :return: True
    :rtype: bool
    """
    if request is None:
        raise InvalidRequest("request is undefined")
    elif type(request) is not dict:
        raise InvalidRequest("request is not a dict")
    elif "method" not in request:
        raise InvalidRequest("method is undefined")
    elif "jsonrpc" not in request:
        raise InvalidRequest("jsonrpc is undefined")
    elif request["jsonrpc"] != "2.0":
        raise InvalidRequest("jsonrpc is not 2.0")
    else:
        for key in request.keys():
            if key in request_valid_keys:
                pass
            else:
                raise InvalidRequest(f"invalid key, {key}")
    return True


def is_valid_response(response):
    """Verifies the response dict conforms to JSON-RPC2 spec.  Raises error if invalid.

    :param response: The response object to validate
    :type response: dict
    :return: True
    :rtype: bool
    """
    if "jsonrpc" not in response:
        raise InvalidResponse("jsonrpc is undefined")
    elif response["jsonrpc"] != "2.0":
        raise InvalidResponse("jsonrpc is not 2.0")
    elif "error" in response:
        if type(response["error"]) is not dict:
            raise InvalidResponse("error is not a dict")
        elif "code" not in response["error"]:
            raise InvalidResponse("error code is undefined")
        elif "message" not in response["error"]:
            raise InvalidResponse("error message is undefined")
        elif "result" in response:
            raise InvalidResponse("cannot have error and result")
        else:
            for key in response["error"].keys():
                if key not in response_error_valid_keys:
                    raise InvalidResponse(f"invalid error key, {key}")
    elif "result" not in response:
        raise InvalidResponse("neither result or code is defined")
    else:
        for key in response.keys():
            if key not in response_valid_keys:
                raise InvalidResponse(f"invalid key, {key}")
    return True


class InvalidRequest(Exception):
    pass


class InvalidResponse(Exception):
    pass

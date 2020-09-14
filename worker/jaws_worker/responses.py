"""
Convenience functions to prepare JSON-RPC2 compliant responses.
"""

import http.client.responses as http_responses


def success(result={}, task=None):
    """Return a JSON-RPC2 successful result message.

    :param result: The result returned by a successful RPC call.
    :type result: str|int|dict|list
    :param task: Optional task (Popen obj and rc_file path)
    :type task: dict
    :return: JSON-RPC2 formatted response
    :rtype: str
    """
    return {"jsonrpc": "2.0", "result": result}, task


def failure(code, message=None):
    """Return a JSON-RPC2 error message.

    :param code: The error code returned by the RPC call
    :type code: int
    :param message: The error message returned by the RPC call
    :type message: str, optional
    :return: JSON-RPC2 formatted response
    :rtype: dict
    """
    if message is None:
        message = (
            http_responses["status_code"]
            if "status_code" in http_responses
            else "Unknown error"
        )
    return {"jsonrpc": "2.0", "error": {"code": code, "message": message}}, None

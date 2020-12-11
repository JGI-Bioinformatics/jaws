"""
Convenience functions to prepare JSON-RPC2 compliant responses.
"""


class ErrorWithCode(Exception):
    def __init__(self, message, code=500):
        self.message = message
        self.code = code
        super().__init__(self.message)


def success(result={}):
    """Return a JSON-RPC2 successful result message.

    :param result: The result returned by a successful RPC call.
    :type result: str|int|dict|list
    :param task: Optional task (Popen obj and rc_file path)
    :type task: dict
    :return: JSON-RPC2 formatted response
    :rtype: str
    """
    return {"jsonrpc": "2.0", "result": result}


def failure(error: Exception):
    """Return a JSON-RPC2 error message.

    :param error: Exception object with code and message
    :type error: ErrorWithCode
    :return: JSON-RPC2 formatted response
    :rtype: dict
    """
    code = 500
    if hasattr(error, "code"):
        code = error.code
    message = f"{error}"
    return {"jsonrpc": "2.0", "error": {"code": code, "message": message}}

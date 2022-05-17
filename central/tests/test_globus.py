import pytest
from jaws_central import globus
import globus_sdk


class MockClient:
    def __init__(self):
        pass


class MockAuthorizer:
    def __init__(self):
        pass


class MockNativeAppAuthClient:
    def __init__(self, client_id):
        pass


class MockConfidentialAppAuthClient:
    def __init__(self, client_id, client_secret):
        pass


class MockRefreshTokenAuthorizer:
    def __init__(self, token, client):
        pass


class MockClientCredentialsAuthorizer:
    def __init__(self, client, scopes):
        pass


class MockGlobusTransferClient:
    def __init__(self, authorizer=None):
        self.authorizer = authorizer

    def cancel_task(self, task_id):
        if task_id == "error":
            raise globus_sdk.GlobusAPIError("error")


@pytest.fixture()
def mock_globus(monkeypatch):
    monkeypatch.setattr(globus_sdk, "NativeAppAuthClient", MockNativeAppAuthClient)
    monkeypatch.setattr(
        globus_sdk, "ConfidentialAppAuthClient", MockConfidentialAppAuthClient
    )
    monkeypatch.setattr(
        globus_sdk, "RefreshTokenAuthorizer", MockRefreshTokenAuthorizer
    )
    monkeypatch.setattr(
        globus_sdk, "ClientCredentialsAuthorizer", MockClientCredentialsAuthorizer
    )
    monkeypatch.setattr(globus_sdk, "TransferClient", MockGlobusTransferClient)


def test_virtual_transfer_path(configuration):
    globus_client = globus.GlobusService()
    test_data = [
        (
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
            "/",
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
        ),  # noqa: E501
        (
            "/global/scratch/jaws/jaws-dev/inputs/mmelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
            "/global/scratch/jaws",
            "/jaws-dev/inputs/mmelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
        ),  # noqa: E501
        (
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",  # noqa
            "/",
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",  # noqa
        ),  # noqa: E501
        (
            "/global/scratch/jaws/jaws-dev/inputs/mmelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",  # noqa
            "/global/scratch/jaws",
            "/jaws-dev/inputs/mmelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",
        ),  # noqa: E501
        (
            "/global/cfs/cdirs/jaws/data-repository-dev/mmelara/CORI/5247",
            "/",
            "/global/cfs/cdirs/jaws/data-repository-dev/mmelara/CORI/5247",
        ),  # noqa: E501
        (
            "/global/scratch/jaws/data-repository-dev/mmelara/JGI/3593",
            "/global/scratch/jaws",
            "/data-repository-dev/mmelara/JGI/3593",
        ),  # noqa: E501
    ]
    for (full_path, host_path, expected_path) in test_data:
        result = globus_client.virtual_transfer_path(full_path, host_path)
        assert result == expected_path

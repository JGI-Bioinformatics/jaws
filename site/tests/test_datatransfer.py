import pytest
from jaws_site.datatransfer_protocol import DataTransferFactory


@pytest.mark.parametrize(
    'transfer_type, expect',
    [
        ('globus_transfer', "<class 'jaws_site.datatransfer_plugins.globus_transfer.DataTransfer'>"),
        ('aws_s3_transfer', "<class 'jaws_site.datatransfer_plugins.aws_s3_transfer.DataTransfer'>"),
    ]
)
def test_datatransfer_factory(transfer_type, expect):

    data_transfer = DataTransferFactory(transfer_type)
    assert str(data_transfer.__class__) == expect

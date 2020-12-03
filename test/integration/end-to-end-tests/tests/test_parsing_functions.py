import pytest

from JawsTestWrappers.parsing_functions import source_environment

def test_source_environment(mocker):
    mocker.patch('JawsTestWrappers.parsing_functions.submit_cmd',
                 return_value='jaws-dev activated; see "jaws --help" for more.')
    assert source_environment('dev') == 'source ~/jaws-dev.sh'
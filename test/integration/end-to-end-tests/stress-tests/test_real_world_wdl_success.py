#!/usr/bin/env python

import pytest
import submission_utils as util


@pytest.mark.parametrize(
    "wdl, input_json",
    (
        #  ('../../../examples/marcels_meta_annot/annotation.wdl',
        #   '../../../examples/marcels_meta_annot/JsonFiles/marcel.2.lrcs.clean.json'),
        (
            "../../../examples/bfoster_meta_assem/jgi_meta.jaws.wdl",
            "../../../examples/bfoster_meta_assem/inputs.json",
        ),
#        (
#            "../../../examples/leo_dapseq/Azospirillum_brasilense.wdl",
#            "../../../examples/leo_dapseq/shortened-100.json",
#        ),
    ),
)
def test_real_world_success(site, wdl, input_json):
    # run the test against all the wdl/json in the @pytest.mark.parametrize
    util.run_success(site, wdl, input_json)

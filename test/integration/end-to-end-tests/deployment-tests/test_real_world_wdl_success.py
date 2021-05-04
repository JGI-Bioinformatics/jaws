#!/usr/bin/env python

import pytest
import test_tutorial_wdl_success as tts


class TestRealWorldWdlSuccess(tts.TestRunSuccess):

    @pytest.mark.parametrize(
        'wdl, input_json',
        (
              #  ('../../../examples/marcels_meta_annot/annotation.wdl',
              #   '../../../examples/marcels_meta_annot/JsonFiles/marcel.2.lrcs.clean.json'),
                ('../../../examples/bfoster_meta_assem/jgi_meta.jaws.wdl',
                 '../../../examples/bfoster_meta_assem/inputs.json'),
                ('../../../examples/leo_dapseq/Azospirillum_brasilense.wdl',
                 '../../../examples/leo_dapseq/shortened-100.json'),

        )
    )
    def test_tutorial_success(self, env, site, wdl, input_json):
        # run the test against all the wdl/json in the @pytest.mark.parametrize
        self.run_success(env, site, wdl, input_json)

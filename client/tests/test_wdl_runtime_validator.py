#!/usr/bin/env python
"""
These functions are to test the "runtime{}" validation tests are working.

This library uses "fixtures" from conftest.py which should be located in the same directory.
There is no need to import conftest.py as it is done automatically.

Author: Jeff Froula jlfroula@lbl.gov>
Updated: 04/18/21
"""

import pytest
from jaws_client.wdl_runtime_validator import ( validate_wdl_runtime, WdlRuntimeError, WdlRuntimeMemoryError, WdlRuntimeTimeError,) # noqa


def test_good_wdl():
    validate_wdl_runtime(GOOD_WDL)


def test_allRequiredParams():
    """
    Make sure there is an error message when a WDL has a runtime missing
    time or memory, which are the minimum settings.
    """
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(NO_TIME_AND_MEMORY)


def test_timeParam():
    """
    1. if constraint: knl, then limit is 48hrs
    2. You are limited to 48hrs where constraint=knl"
    3. You are limited to 72hrs when constraint is the default value "haswell"
    4. You are limited to 72hrs where constraint=haswell"
    """
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(BAD_TIME_LIMITS)


def test_memoryParam():
    """
    1. If you are using skylake, you must have account: set to fungalp."
    2. if constraint: skylake, then qos must be jgi_exvivo(250G) or jgi_shared(758G)
    3. if constraint: knl, then limit is 96G memory
    """
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(BAD_MEMORY_VALUE)


def test_runtimeCombinations():
    """testing the test that validates the user included the correct combination of params for various resource allocations""" # noqa
    with pytest.raises(WdlRuntimeError):
        validate_wdl_runtime(BAD_COMBINATIONS)


# TODO no_qa
#
BAD_COMBINATIONS = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
}
task bbcms {
     File infile
     String container

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        qos: "regular"
        account: "wrong"
    }

}
task read_mapping_pairs{

    runtime {
        docker: container
        time: "03:00:00"
        memory: "115G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
        qos: "regular"

    }
}
"""


BAD_MEMORY_VALUE = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
}
task bbcms {
    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "skylake"
    }

}

task assy {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "knl"
    }
}


task create_agp {
    runtime {
        docker: container
        time: "03:00:00"
        memory: "800G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "skylake"
        qos: "jgi_exvivo"
    }

}

task read_mapping_pairs{
    runtime {
        docker: container
        time: "03:00:00"
        memory: "300G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "skylake"
        qos: "jgi_shared"

    }

}
"""


BAD_TIME_LIMITS = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
    call assy {
         input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    call read_mapping_pairs {
     input: reads=input_file, ref=create_agp.outcontigs, container=bbtools_container
    }

}
task bbcms {
    runtime {
        docker: container
        time: "200:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "skylake"
    }

}

task assy {
    runtime {
        docker: container
        time: "80:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "knl"
    }

}

task create_agp {

    runtime {
        docker: container
        time: "80:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}


task read_mapping_pairs{
    runtime {
        docker: container
        time: "80:00:00"
        memory: "115G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
    }

}
"""


NO_TIME_AND_MEMORY = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
    call assy {
         input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    call read_mapping_pairs {
     input: reads=input_file, ref=create_agp.outcontigs, container=bbtools_container
    }

}
task bbcms {
    runtime {
        docker: container
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task assy {
    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task create_agp {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task read_mapping_pairs{

    runtime {
        docker: container
        time: "03:00:00"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
    }

}
"""


GOOD_WDL = """
workflow jgi_meta {
    call bbcms {
          input: infile=input_file, container=bbtools_container
    }
    call assy {
         input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    call read_mapping_pairs {
     input: reads=input_file, ref=create_agp.outcontigs, container=bbtools_container
    }

}
task bbcms {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task assy {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task create_agp {

    runtime {
        docker: container
        time: "03:00:00"
        memory: "118G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
    }

}

task read_mapping_pairs{

    runtime {
        docker: container
        time: "03:00:00"
        memory: "115G"
        node: 1
        nwpn: 1
        poolname: "bfostersmall"
        shared: 0
        constraint: "haswell"
    }
}
"""

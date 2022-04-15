import pytest
import os
import shutil
import uuid
import json

# flake8: noqa

def join_path(*args):
    return os.path.join(*args)

def recur_dict_comparison(expected_val, actual_val):
    if isinstance(actual_val, str) or isinstance(actual_val, int):
        assert actual_val == expected_val
    elif isinstance(actual_val, list):
        for i, elem in enumerate(actual_val):
            recur_dict_comparison(expected_val[i], elem)
    elif isinstance(actual_val, dict):
        for k in actual_val:
            assert k in expected_val
            recur_dict_comparison(expected_val[k], actual_val[k])


def dict_comparison(expected_dict, actual_dict):
    for k in expected_dict:
        recur_dict_comparison(expected_dict[k], actual_dict[k])

def test_create_destination_json(configuration, input_file):
    root_dir = input_file
    staging_dir = os.path.join(root_dir, "test_user", "CORI")
    expected = {
        "test.file1": f"{staging_dir}{root_dir}/test.fasta"
    }

    uuid = "12345"
    json_filepath = os.path.join(root_dir, "test.json")
    wdl_filepath = os.path.join(root_dir, "test.wdl")
    import jaws_client.workflow
    wf_inputs = jaws_client.workflow.WorkflowInputs(json_filepath, uuid, \
        wdl_loc=wdl_filepath)
    actual = wf_inputs.prepend_paths_to_json(staging_dir)
    assert actual
    dict_comparison(expected, actual.inputs_json)

def test_files_inputs(configuration, files_inputs):
    uuid = "1234"
    import jaws_client.workflow
    inputs_json = os.path.join(files_inputs, "files.json")
    wdl = os.path.join(files_inputs, "files.wdl")
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, uuid, wdl_loc = wdl)
    expected = ["/path/1", "2", "/path/3", "https://path4", "5", \
        "6", "/path/7", "/path/8", "9", "/path/10", "http://abf", "ftp://xyz"]
    assert 12 == len(inputs.src_file_inputs)
    for expect in expected:
        if not expect.startswith(("/", "http://", "ftp://", "https://")):
            expect = os.path.abspath(join_path(inputs.basedir, expect))
        assert expect in inputs.src_file_inputs

@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_validation(configuration, simple_wdl_example):
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(
        os.path.join(simple_wdl_example, "align.wdl"), "1234"
    )
    wdl.validate(1000)


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_validation_with_no_subworkflows(configuration, no_subworkflows_present):
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(no_subworkflows_present, "1234")
    with pytest.raises(jaws_client.workflow.WdlError):
        wdl.validate(1000)


@pytest.mark.skipif(
    shutil.which('womtool') is None, reason="WOMTool needs to be installed"
)
def test_bad_syntax_wdl(configuration, incorrect_wdl):
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(incorrect_wdl, "1234")
    with pytest.raises(jaws_client.workflow.WdlError):
        wdl.validate(1000)


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_subworkflows(configuration, subworkflows_example):
    basedir = subworkflows_example
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "main.wdl"), "1234")
    sub1 = os.path.join(basedir, "sub1.wdl")
    sub2 = os.path.join(basedir, "sub2.wdl")

    subworkflows = [
        jaws_client.workflow.WdlFile(sub1, "1234"),
        jaws_client.workflow.WdlFile(sub2, "1234"),
    ]

    assert len(wdl.subworkflows) == 2

    for expected in subworkflows:
        assert expected in wdl.subworkflows
    assert wdl.max_ram_gb == 6


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed"
)
def test_calculate_wdl_max_ram_gb(configuration, dap_seq_example):
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(
        os.path.join(dap_seq_example, "test.wdl"), "1234"
    )
    assert 5 == wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed"
)
def test_calculate_wdl_max_ram_gb_warn_on_mem_keyword(configuration, deprecated_mem_example):
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(
        os.path.join(deprecated_mem_example, "deprecated_mem.wdl"), "1234"
    )
    with pytest.raises(jaws_client.workflow.WdlError):
        gb = wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_calculate_wdl_max_ram_gb_with_subworkflows(
    configuration, subworkflows_example
):
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(
        os.path.join(subworkflows_example, "main.wdl"), "1234"
    )
    assert 6 == wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_appropriate_staging_dir_for_all_wdls(configuration, subworkflows_example):
    basedir = subworkflows_example
    staging = os.path.join(basedir, "staging")
    os.mkdir(staging)

    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "main.wdl"), "1234")

    new_wdl_path = os.path.join(staging, wdl.name)
    wdl.copy_to(new_wdl_path)

    assert staging not in wdl.file_location
    assert os.path.exists(os.path.join(staging, wdl.name))

    zip_path = os.path.join(staging, wdl.submission_id)
    os.mkdir(zip_path)

    for sub in wdl.subworkflows:
        new_path = os.path.join(zip_path, sub.name)
        sub.copy_to(new_path)
        assert os.path.exists(new_path)


#  We are temporarily deactivating this check
#def test_fail_invalid_backend(wdl_with_invalid_backend):
#    import jaws_client.workflow
#    wdl = jaws_client.workflow.WdlFile(wdl_with_invalid_backend, "1234")
#    with pytest.raises(jaws_client.workflow.WdlError) as e_info:
#        wdl.verify_wdl_has_no_backend_tags()


def test_move_input_files_to_destination(configuration, sample_workflow):
    inputs = os.path.join(sample_workflow, "workflow", "sample.json")
    wdl = os.path.join(sample_workflow, "workflow", "sample.wdl")
    staging_dir = os.path.join(sample_workflow, "staging")
    import jaws_client.workflow
    inputs_json = jaws_client.workflow.WorkflowInputs(inputs, uuid.uuid4(), \
        wdl_loc=wdl)
    inputs_json.move_input_files(os.path.join(staging_dir, "NERSC"))


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_zipping_up_of_subworkflow_files(configuration, subworkflows_example):
    basedir = os.path.basename(subworkflows_example)
    shutil.copytree(subworkflows_example, basedir)
    staging_dir = os.path.join(basedir, "staging")
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "main.wdl"), "1234")
    staged_wdl, zip_file = wdl.prepare_wdls(staging_dir)
    #assert os.path.exists(staged_wdl)
    #assert os.path.exists(zip_file)


def test_no_zip_file_in_manifest_if_no_subworkflows(simple_wdl_example):
    basedir = simple_wdl_example
    staging_dir = os.path.join(basedir, "staging")
    compute_dir = os.path.join(basedir, "compute")

    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "align.wdl"), "1234")
    staged_wdl, zip_file = wdl.prepare_wdls(basedir)
    manifest_file = jaws_client.workflow.Manifest(staging_dir, compute_dir)
    manifest_file.add(staged_wdl, zip_file)

    for infiles in manifest_file.manifest:
        assert ".zip" not in infiles[0]  # tuple where first entry is the filename


def test_manifest_file(staged_files):
    stage_dir, dest_dir = staged_files
    import jaws_client.workflow
    manifest_file = jaws_client.workflow.Manifest(stage_dir, dest_dir)

    wdl = os.path.join(stage_dir, "example.wdl")
    json_file = os.path.join(stage_dir, "example.json")

    moved_wdl = os.path.join(dest_dir, "example.wdl")
    moved_json = os.path.join(dest_dir, "example.json")

    expected_items = [(wdl, moved_wdl, "F"), (json_file, moved_json, "F")]
    manifest_file.add(wdl, json_file)

    for expected, actual in zip(expected_items, manifest_file.manifest):

        expected_src_file, expected_moved_file, expected_inode_type = expected
        actual_src_file, actual_moved_file, actual_inode_type = actual

        assert expected_src_file == actual_src_file
        assert expected_moved_file == actual_moved_file
        assert expected_inode_type == actual_inode_type


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_same_submission_id_in_workflow_files(subworkflows_example):
    submission_id = "1234567890"
    basedir = os.path.basename(subworkflows_example)
    shutil.copytree(subworkflows_example, basedir)
    wdl_file = os.path.join(basedir, "main.wdl")
    import jaws_client.workflow
    wdl = jaws_client.workflow.WdlFile(wdl_file, submission_id)
    zip_path = os.path.join(subworkflows_example, "zip_directory")
    zip_wdl, _ = wdl.prepare_wdls(zip_path)
    assert os.path.basename(zip_wdl).strip(".wdl") == submission_id


def test_refdata_not_translated(refdata_inputs, refdata_wdl):
    inputs_json = os.path.join(refdata_inputs, "inputs.json")
    text_file = os.path.join(refdata_inputs, "file1.txt")
    wdl = os.path.join(refdata_wdl, "refdata.wdl")
    import jaws_client.workflow
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, "1231231", wdl_loc=wdl)
    modified_json = inputs.prepend_paths_to_json("/remote/uploads/NERSC/staging")
    expected = {"refdata.file1": "/remote/uploads/NERSC/staging" + text_file,
                "refdata.runblastplus_sub.ncbi_nt": "/refdata/nt"}
    dict_comparison(modified_json.inputs_json, expected)


def test_refdata_in_different_form(refdata_inputs_missing_slash, refdata_wdl):
    inputs_json = os.path.join(refdata_inputs_missing_slash, "inputs.json")
    text_file = os.path.join(refdata_inputs_missing_slash, "file1.txt")
    wdl = os.path.join(refdata_wdl, "refdata.wdl")
    import jaws_client.workflow
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, "1231231", wdl_loc=wdl)
    modified_json = inputs.prepend_paths_to_json("/remote/uploads/NERSC/staging")
    expected = {"refdata.file1": "/remote/uploads/NERSC/staging" + text_file,
                "refdata.runblastplus_sub.ncbi_nt": "/refdata"}
    dict_comparison(modified_json.inputs_json, expected)


def test_refdata_in_inputs_json(refdata_inputs, monkeypatch, refdata_wdl):
    inputs_json = os.path.join(refdata_inputs, "inputs.json")
    wdl = os.path.join(refdata_wdl, "refdata.wdl")
    import jaws_client.workflow
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, "12312", wdl_loc=wdl)
    assert "/refdata/nt" in inputs.src_file_inputs


def test_rel_path_in_input_files(relative_inputs):
    """test_json = {
        "fileX": "./fileX.txt",
        "fileY": "../fileY.txt",
        "fileZ": "/home/profx/fileZ.txt"
    }"""
    root_dir = relative_inputs
    inputs_json = os.path.join(root_dir, "test", "relative.json")
    wdl = os.path.join(root_dir, "test", "relative.wdl")
    import jaws_client.workflow
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, 'ABCDEF', wdl_loc=wdl)
    for path in inputs.src_file_inputs:
        assert path.startswith(root_dir)


"""
# commented out because nested file feature is being checked in test_src_json_inputs
def test_nested_files_are_in_src_file_inputs():
    test_json = {
        "jgi_dap_leo.adapters": "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/adapters.fa",
        "jgi_dap_leo.genome_fasta": "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245.fasta",
        "jgi_dap_leo.bt2index_file1": "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index",
        "jgi_dap_leo.bt2index_list": [
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.1.bt2",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.2.bt2",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.3.bt2",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.4.bt2",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.rev.1.bt2",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.rev.2.bt2"
        ],
        "jgi_dap_leo.effgsize": 7682393,
        "jgi_dap_leo.genes_gff": "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_filt.gff",
        "jgi_dap_leo.bgmodel": "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245.bgmodel",
        "jgi_dap_leo.outdir": "nersc_out",
        "jgi_dap_leo.amplified": "10cyc",
        "jgi_dap_leo.maxfrags": 1000011,
        "jgi_dap_leo.find_motifs": "false",
        "jgi_dap_leo.ctl_raw_fastqs": [],
        "jgi_dap_leo.expt_raw_fastqs": [
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-CATTGGAC+TCAAGCAC.fastq.gz",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TACGTGAC+TCAAGCAC.fastq.gz",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TGAGGATG+TCAAGCAC.fastq.gz"
        ],
        "jgi_dap_leo.library_names_map": {
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-CATTGGAC+TCAAGCAC.fastq.gz": "A-CATTGGAC+TCAAGCAC",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TACGTGAC+TCAAGCAC.fastq.gz": "A-TACGTGAC+TCAAGCAC",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TGAGGATG+TCAAGCAC.fastq.gz": "A-TGAGGATG+TCAAGCAC"
        },
        "jgi_dap_leo.sample_names_map": {
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-CATTGGAC+TCAAGCAC.fastq.gz": "accB",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TACGTGAC+TCAAGCAC.fastq.gz": "araC",
            "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TGAGGATG+TCAAGCAC.fastq.gz": "birA"
        }
    }
    expected_file_paths = {
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/adapters.fa",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245.fasta",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.1.bt2",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.2.bt2",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.3.bt2",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.4.bt2",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.rev.1.bt2",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_bt2index/Azospirillum_brasilense_Sp245.rev.2.bt2",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245_filt.gff",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/Azospirillum_brasilense_Sp245.bgmodel",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-CATTGGAC+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TACGTGAC+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TGAGGATG+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-CATTGGAC+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TACGTGAC+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TGAGGATG+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-CATTGGAC+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TACGTGAC+TCAAGCAC.fastq.gz",
        "/global/cfs/cdirs/jaws/test/tutorial_test_data/LeosData/48sp01/A-TGAGGATG+TCAAGCAC.fastq.gz"
    }

    import jaws_client.workflow
    wf_inputs = jaws_client.workflow.WorkflowInputs('test.json', "ABCDEF", test_json)

    assert len(wf_inputs.src_file_inputs) == len(expected_file_paths)

    for src_file in wf_inputs.src_file_inputs:
        assert src_file in expected_file_paths

# commented out because looks_like_file_path function no longer exists
def test_looks_like_file_path():
    test_inputs = [
        ("./fileX", True),
        ("../fileY", True),
        ("/opt/fileZ", True),
        ("file0", False),
        ("http://some-service.lbl.gov", False)
    ]
    import jaws_client.workflow
    for (input, expected) in test_inputs:
        assert jaws_client.workflow.looks_like_file_path(input) is expected
"""


def test_rsync_excludes(configuration, output_example):
    base_dir = output_example
    src = f"{base_dir}/run1"
    dest = f"{base_dir}/test_copy"
    import jaws_client.workflow
    result = jaws_client.workflow.rsync(
        src,
        dest,
        ["-rLtq", "--chmod=Du=rwx,Dg=rwx,Do=,Fu=rw,Fg=rw,Fo=", "--exclude=inputs"],
    )
    assert result.returncode == 0
    assert os.path.exists(f"{dest}/run1/task1/execution/stdout") is True
    assert os.path.exists(f"{dest}/run1/task1/inputs/infile") is False

def test_struct(relative_inputs):
    root_dir = relative_inputs
    inputs_json = os.path.join(root_dir, "test", "struct.json")
    wdl = os.path.join(root_dir, "test", "struct.wdl")
    import jaws_client.workflow
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, 'ABCDEF', wdl_loc=wdl)
    assert len(inputs.inputs_json["test_struct.product_list"]) == 3
    assert inputs.inputs_json["test_struct.product_list"][0]["name"] == "Apple"
import pytest
import os
import shutil


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
    expected = {"test.file1": f"{staging_dir}{root_dir}/test.fasta"}

    json_filepath = os.path.join(root_dir, "test.json")
    wdl_filepath = os.path.join(root_dir, "test.wdl")
    import jaws_client.workflow

    wf_inputs = jaws_client.workflow.WorkflowInputs(json_filepath, wdl_filepath)
    actual = wf_inputs.prepend_paths_to_json(staging_dir)
    assert actual
    dict_comparison(expected, actual)
    dict_comparison(expected, wf_inputs.inputs_json)


def test_files_inputs(configuration, files_inputs):
    import jaws_client.workflow

    inputs_json = os.path.join(files_inputs, "files.json")
    wdl_file = os.path.join(files_inputs, "files.wdl")
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, wdl_file)
    expected = [
        "/path/1",
        "2",
        "/path/3",
        "https://path4",
        "5",
        "6",
        "/path/7",
        "/path/8",
        "9",
        "/path/10",
        "http://abf",
        "ftp://xyz",
    ]
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

    wdl = jaws_client.workflow.WdlFile(os.path.join(simple_wdl_example, "align.wdl"))
    wdl.validate()


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_validation_with_no_subworkflows(configuration, no_subworkflows_present):
    import jaws_client.workflow

    wdl = jaws_client.workflow.WdlFile(no_subworkflows_present)
    with pytest.raises(jaws_client.workflow.WdlError):
        wdl.validate()


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed"
)
def test_bad_syntax_wdl(configuration, incorrect_wdl):
    import jaws_client.workflow

    wdl = jaws_client.workflow.WdlFile(incorrect_wdl)
    with pytest.raises(jaws_client.workflow.WdlError):
        wdl.validate()


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_subworkflows(configuration, subworkflows_example):
    basedir = subworkflows_example
    import jaws_client.workflow

    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "main.wdl"))
    wdl.validate()

    sub1 = os.path.join(basedir, "sub1.wdl")
    sub2 = os.path.join(basedir, "sub2.wdl")
    subworkflows = [
        jaws_client.workflow.WdlFile(sub1),
        jaws_client.workflow.WdlFile(sub2),
    ]

    assert len(wdl.subworkflows) == 2

    for expected in subworkflows:
        assert expected in wdl.subworkflows
    assert wdl.max_ram_gb == 6


def test_calculate_wdl_max_ram_gb(configuration, dap_seq_example):
    import jaws_client.workflow

    wdl_file = os.path.join(dap_seq_example, "test.wdl")
    wdl = jaws_client.workflow.WdlFile(wdl_file)
    assert 5 == wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_calculate_wdl_max_ram_gb_with_subworkflows(
    configuration, subworkflows_example
):
    import jaws_client.workflow

    wdl_file = os.path.join(subworkflows_example, "main.wdl")
    wdl = jaws_client.workflow.WdlFile(wdl_file)
    wdl.validate()
    assert 6 == wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_copy_to(configuration, subworkflows_example):
    basedir = subworkflows_example
    staging = os.path.join(basedir, "staging")
    os.mkdir(staging)

    import jaws_client.workflow

    wdl_file = os.path.join(basedir, "main.wdl")
    wdl = jaws_client.workflow.WdlFile(wdl_file)

    new_wdl_path = os.path.join(staging, wdl.name)
    wdl.copy_to(new_wdl_path)
    assert staging not in wdl.file_location
    assert os.path.exists(new_wdl_path)


def test_fail_invalid_backend(wdl_with_invalid_backend):
    import jaws_client.workflow

    with pytest.raises(jaws_client.workflow.WdlError):
        jaws_client.workflow.WdlFile(wdl_with_invalid_backend)


def test_copy_input_files_to_destination(configuration, sample_workflow):
    inputs_file = os.path.join(sample_workflow, "workflow", "sample.json")
    wdl_file = os.path.join(sample_workflow, "workflow", "sample.wdl")
    staging_dir = os.path.join(sample_workflow, "staging")
    import jaws_client.workflow

    inputs = jaws_client.workflow.WorkflowInputs(inputs_file, wdl_file)
    copied_files = inputs.copy_input_files(os.path.join(staging_dir, "NERSC"))
    assert len(copied_files) == 7


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_zipping_up_of_subworkflow_files(configuration, subworkflows_example):
    staging_dir = os.path.join(subworkflows_example, "staging")
    import jaws_client.workflow

    wdl_file = os.path.join(subworkflows_example, "main.wdl")
    wdl = jaws_client.workflow.WdlFile(wdl_file)
    submission_id = "SUBMISSIONID"
    staged_wdl, zip_file = wdl.prepare_wdls(staging_dir, submission_id)

    expected_staged_wdl = os.path.join(staging_dir, f"{submission_id}.wdl")
    expected_staged_zip = os.path.join(staging_dir, f"{submission_id}.zip")

    assert staged_wdl
    assert zip_file
    assert staged_wdl == expected_staged_wdl
    assert zip_file == expected_staged_zip
    assert os.path.exists(staged_wdl)
    assert os.path.exists(zip_file)


def test_manifest_add_nonexistant_file(simple_wdl_example):
    basedir = simple_wdl_example
    nonexistant_file = f"{basedir}/foo.bar"

    import jaws_client.workflow

    manifest = jaws_client.workflow.Manifest(basedir)
    manifest.add(nonexistant_file)

    for file_path in manifest.files:
        assert "foo.bar" not in file_path


def test_manifest_file(staged_files):
    basedir, file1, file2 = staged_files
    import jaws_client.workflow

    manifest = jaws_client.workflow.Manifest(basedir)
    expected = ["file1.txt", "file2.txt"]
    file1 = os.path.abspath(file1)
    file2 = os.path.abspath(file2)
    manifest.add(file1, file2)

    for actual in manifest.files:
        assert actual in expected


def test_refdata_not_translated(refdata_inputs, refdata_wdl):
    inputs_json = os.path.join(refdata_inputs, "inputs.json")
    text_file = os.path.join(refdata_inputs, "file1.txt")
    wdl_file = os.path.join(refdata_wdl, "refdata.wdl")
    import jaws_client.workflow

    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, wdl_file)
    staging_dir = "/STAGING"
    modified_json = inputs.prepend_paths_to_json(staging_dir)

    expected_translated = f"{staging_dir}{text_file}"
    expected_not_translated = "/refdata/db.fasta"

    assert expected_translated == modified_json["refdata.file1"]
    assert expected_not_translated == modified_json["refdata.db"]


def test_refdata_in_different_form(refdata_inputs_missing_slash, refdata_wdl):
    inputs_json = os.path.join(refdata_inputs_missing_slash, "inputs.json")
    text_file = os.path.join(refdata_inputs_missing_slash, "file1.txt")
    wdl_file = os.path.join(refdata_wdl, "refdata.wdl")
    import jaws_client.workflow

    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, wdl_file)
    modified_json = inputs.prepend_paths_to_json("/remote/uploads/NERSC/staging")
    expected = {
        "refdata.file1": "/remote/uploads/NERSC/staging" + text_file,
        "refdata.db": "/refdata",
    }
    dict_comparison(modified_json, expected)


def test_refdata_in_inputs_json(refdata_inputs, monkeypatch, refdata_wdl):
    inputs_json = os.path.join(refdata_inputs, "inputs.json")
    wdl_file = os.path.join(refdata_wdl, "refdata.wdl")
    import jaws_client.workflow

    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, wdl_file)
    assert "/refdata/db.fasta" in inputs.src_file_inputs


def test_rel_path_in_input_files(relative_inputs):
    """test_json = {
        "fileX": "./fileX.txt",
        "fileY": "../fileY.txt",
        "fileZ": "/home/profx/fileZ.txt"
    }"""
    basedir = relative_inputs
    inputs_json = os.path.join(basedir, "test", "relative.json")
    wdl_file = os.path.join(basedir, "test", "relative.wdl")
    import jaws_client.workflow

    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, wdl_file)
    for path in inputs.src_file_inputs:
        assert path.startswith(basedir)


def test_struct(struct_inputs):
    root_dir = struct_inputs
    inputs_json = os.path.join(root_dir, "test", "struct.json")
    wdl_file = os.path.join(root_dir, "test", "struct.wdl")
    import jaws_client.workflow

    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, wdl_file)
    assert len(inputs.inputs_json["test_struct.product_list"]) == 3
    for path in inputs.src_file_inputs:
        assert path.startswith(root_dir)
    assert inputs.inputs_json["test_struct.product_list"][0]["name"] == "Apple"

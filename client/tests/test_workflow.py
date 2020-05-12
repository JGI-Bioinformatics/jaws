import pytest
import os
import shutil
import uuid

import jaws_client.config
import jaws_client.workflow


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


def test_create_destination_json(configuration, dap_seq_example):

    root_dir = dap_seq_example
    staging_dir = os.path.join(root_dir, "jaws_central", "staging")

    expected = {
        "jgi_dap_leo.adapters": f"{staging_dir}/global/projectb/sandbox/gaag/bbtools/data/adapters.fa",
        "jgi_dap_leo.genome_fasta": f"{staging_dir}/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417.fasta",  # noqa
        "jgi_dap_leo.bt2index_dir": f"{staging_dir}/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417_bt2index",  # noqa
        "jgi_dap_leo.bt2index_name": "PsimiaeWCS417",
        "jgi_dap_leo.effgsize": 6169071,
        "jgi_dap_leo.genes_gff": f"{staging_dir}/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417.gff",  # noqa
        "jgi_dap_leo.bgmodel": f"{staging_dir}/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417.bgmodel",  # noqa
        "jgi_dap_leo.outdir": ".",
        "jgi_dap_leo.expt_raw_fastqs": [
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TATCAGC-TATCAGC.fastq.gz"
        ],  # noqa
        "jgi_dap_leo.ctl_raw_fastqs": [],
        "jgi_dap_leo.library_names_map": {
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TATCAGC-TATCAGC.fastq.gz": "CTTZN",  # noqa
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.GGTTGAT-GGTTGAT.fastq.gz": "CTUGZ",  # noqa
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TTGCTGG-TTGCTGG.fastq.gz": "CTUHA",  # noqa
        },
        "jgi_dap_leo.sample_names_map": {
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TATCAGC-TATCAGC.fastq.gz": "TF4",
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.GGTTGAT-GGTTGAT.fastq.gz": "negCtl2",  # noqa
            f"{staging_dir}/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TTGCTGG-TTGCTGG.fastq.gz": "negCtl1",  # noqa
        },
        "jgi_dap_leo.expt_bam": f"{staging_dir}/global/projectb/scratch/jaws/jfroula/leo_dap/CTTZN_TF4.bam",
        "jgi_dap_leo.expt_bai": f"{staging_dir}/global/projectb/scratch/jaws/jfroula/leo_dap/CTTZN_TF4.bam",
    }

    json_filepath = os.path.join(root_dir, "test.json")
    uuid = "12345"
    wf_inputs = jaws_client.workflow.WorkflowInputs(json_filepath, uuid)
    actual = wf_inputs.prepend_paths_to_json(staging_dir)
    assert actual
    dict_comparison(expected, actual.inputs_json)


def test_src_json_inputs(configuration, inputs_json):
    uuid = "1234"
    inputs = jaws_client.workflow.WorkflowInputs(inputs_json, uuid)
    expected = ["path/to/file1", "path/to/file2"]

    assert 2 == len(inputs.src_file_inputs)

    for expect in expected:
        assert expect in inputs.src_file_inputs


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_validation(configuration, simple_wdl_example):
    wdl = jaws_client.workflow.WdlFile(os.path.join(simple_wdl_example, "align.wdl"), "1234")
    wdl.validate()


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_wdl_subworkflows(configuration, subworkflows_example):
    basedir = subworkflows_example
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


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed"
)
def test_calculate_wdl_max_ram_gb(configuration, dap_seq_example):
    wdl = jaws_client.workflow.WdlFile(os.path.join(dap_seq_example, "test.wdl"), "1234")
    assert 5 == wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_calculate_wdl_max_ram_gb_with_subworkflows(
    configuration, subworkflows_example
):
    wdl = jaws_client.workflow.WdlFile(os.path.join(subworkflows_example, "main.wdl"), "1234")
    assert 0 == wdl.max_ram_gb


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_appropriate_staging_dir_for_all_wdls(configuration, subworkflows_example):
    basedir = subworkflows_example
    staging = os.path.join(basedir, "staging")
    os.mkdir(staging)

    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "main.wdl"), "1234")

    new_wdl_path = os.path.join(staging, wdl.name)
    wdl.write_to(new_wdl_path)

    assert staging not in wdl.file_location
    assert os.path.exists(os.path.join(staging, wdl.name))

    zip_path = os.path.join(staging, wdl.submission_id)
    os.mkdir(zip_path)

    for sub in wdl.subworkflows:
        new_path = os.path.join(zip_path, sub.name)
        sub.write_to(new_path)
        assert os.path.exists(new_path)


def test_remove_invalid_backend(wdl_with_invalid_backend):
    wdl = jaws_client.workflow.WdlFile(wdl_with_invalid_backend, "1234")
    wdl_with_backend_removed = wdl.sanitized_wdl()
    for line in wdl_with_backend_removed.contents:
        assert "backend" not in line


def test_move_input_files_to_destination(configuration, sample_workflow):
    inputs = os.path.join(sample_workflow, "workflow", "sample.json")
    staging = os.path.join(sample_workflow, "staging")
    inputs_json = jaws_client.workflow.WorkflowInputs(inputs, uuid.uuid4())
    jaws_client.workflow.move_input_files(inputs_json, os.path.join(staging, "NERSC"))


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_zipping_up_of_subworkflow_files(configuration, subworkflows_example):
    basedir = subworkflows_example
    staging = os.path.join(basedir, "staging")
    wdl = jaws_client.workflow.WdlFile(os.path.join(basedir, "main.wdl"), "1234")
    staged_wdl, zip_file = jaws_client.workflow.compress_wdls(
        wdl, compressed_path=staging
    )
    assert os.path.exists(staged_wdl)
    assert os.path.exists(zip_file)


def test_manifest_file(staged_files):
    stage_dir, dest_dir = staged_files
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
    wdl_file = os.path.join(subworkflows_example, "main.wdl")
    wdl = jaws_client.workflow.WdlFile(wdl_file, submission_id)
    zip_path = os.path.join(subworkflows_example, "zip_directory")
    zip_wdl, _ = jaws_client.workflow.compress_wdls(wdl, compressed_path=zip_path)
    assert os.path.basename(zip_wdl).strip(".wdl") == submission_id

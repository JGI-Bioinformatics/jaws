import pytest
import os
import jaws_client.config


# flake8: noqa

@pytest.fixture()
def configuration(tmp_path):
    config_path = tmp_path / "jaws_client.ini"
    user_config_path = tmp_path / "jaws_user.ini"
    os.environ["JAWS_CLIENT_CONFIG"] = config_path.as_posix()
    os.environ["JAWS_USER_CONFIG"] = user_config_path.as_posix()
    os.environ["JAWS_TZ"] = "US/PACIFIC"

    globus_basedir = tmp_path / "globus_basedir"
    staging_dir = globus_basedir / "staging"
    globus_basedir.mkdir()
    staging_dir.mkdir()

    contents = """
[JAWS]
site_id = CORI
url = http://localhost:5001/api/v2
womtool_jar =
uploads_subdir = {0}/globus/staging/uploads
staging_dir = {0}/globus/staging/users
data_repo_basedir = {0}/globus/data-repository-dev
shared_endpoint_group = genome
default_container = debian:latest
[GLOBUS]
client_id =
endpoint_id =
host_path = {0}/globus
""".format(tmp_path.as_posix())
    config_path.write_text(contents)

    user_contents = """
[USER]
token = "xasdasdasfasdasdasfas"
"""
    user_config_path.write_text(user_contents)

    return (config_path.as_posix(), user_config_path.as_posix())


@pytest.fixture
def wdl_path(tmp_path):
    wdls = tmp_path / "wdls"
    wdls.mkdir()
    return wdls


@pytest.fixture()
def input_file(wdl_path):
    wdl_dir = wdl_path
    inputs = wdl_dir / "test.json"
    path = wdl_dir.as_posix()
    contents = """
{
    "test.file1": "%s/test.fasta"
}"""
    inputs.write_text(contents % path)
    test_wdl = wdl_dir / "test.wdl"
    test_wdl.write_text("""
workflow test {
  File file1
  call hello_world
}
task hello_world {
  String hi = "Hello world"

  command {
    echo ${hi}
  }
}
    """)
    test_file = wdl_dir / "test.fasta"
    test_file.write_text("""
>test
GATTACA
    """)
    return path


@pytest.fixture
def no_runtime_file(wdl_path):
    wdl_dir = wdl_path
    no_runtime_wdl = wdl_dir / "test.wdl"
    contents = """
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
}"""
    no_runtime_wdl.write_text(contents)
    return no_runtime_wdl


@pytest.fixture
def wdl_with_invalid_backend(wdl_path):
    wdl_dir = wdl_path
    good_runtime = wdl_dir / "backend.wdl"
    contents = """
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
  runtime {
    cpu: "1"
    backend: "AWS"
  }
}"""
    good_runtime.write_text(contents)
    return good_runtime


@pytest.fixture
def bad_runtime_file(wdl_path):
    wdl_dir = wdl_path
    bad_runtime = wdl_dir / "test.wdl"
    contents = """
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
  runtime {
    cpus: "1"
  }
}"""
    bad_runtime.write_text(contents)
    return bad_runtime


@pytest.fixture()
def manifest_file(tmp_path):
    p = tmp_path / "staging"
    try:
        p.mkdir()
    except FileExistsError:
        pass
        # Exists already so we just do nothing

    manifest = p / "manifest.tsv"
    contents = """
src\tdest\tinode_type
    """
    manifest.write_text(contents)
    return manifest.as_posix()


@pytest.fixture()
def simple_wdl_example(tmp_path):
    bbtools = tmp_path / "align.wdl"
    inputs = tmp_path / "inputs.json"

    wdl_contents = """
workflow bbtools {
    File reads
    File ref

    call alignment {
       input: fastq=reads,
              fasta=ref
    }
    call samtools {
       input: sam=alignment.sam
   }
}

task alignment {
    File fastq
    File fasta

    command <<<
        shifterimg pull jfroula/bbtools:1.2.1 && \
        shifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=${fastq} ref=${fasta} out=test.sam
    >>>
    output {
       File sam = "test.sam"
    }
}

task samtools {
    File sam

    command {
       shifter --image=jfroula/bbtools:1.2.1 samtools view -b -F0x4 ${sam} | shifter --image=jfroula/bbtools:1.2.1 samtools sort - > test.sorted.bam
    }
    output {
       File bam = "test.sorted.bam"
    }
        #runtime {
        #  docker: "jfroula/bbtools:1.2.1"
        #}
}
"""
    inputs_content = """
{
    "bbtools.reads": "/global/dna/shared/data/jfroula/JAWS/data/5min_reads.fq",
    "bbtools.ref": "/global/dna/shared/data/jfroula/JAWS/data/5min_ref.fasta"
}
    """
    bbtools.write_text(wdl_contents)
    inputs.write_text(inputs_content)
    return tmp_path.as_posix()


@pytest.fixture()
def dap_seq_example(tmp_path):
    inputs = tmp_path / "test.json"
    wdl_file = tmp_path / "test.wdl"
    wdl_contents = """
workflow jgi_dap_leo {

    ### from inputs.json:

    Array[File] expt_raw_fastqs
    Array[File] ctl_raw_fastqs
    File adapters
    File genome_fasta
    File bt2index_dir
    String bt2index_name
    Int effgsize
    File genes_gff
    File bgmodel
    Map[File, String] library_names_map
    Map[File, String] sample_names_map
    File outdir
    File expt_bam
    File expt_bai
    File? merged_bam
    File? merged_bai

    ### defined here:

    Int threads = 4
    #1.875 GB per thread on cori haswell config


### Process negative controls in parallel, then merge bams

### Process experimental samples in parallel
### using merged negative control bam as background in findPeaks

    scatter (raw_fastq in expt_raw_fastqs) {

        String basename = library_names_map[raw_fastq] + "_" + sample_names_map[raw_fastq]

        call findPeaks {
            input:  expt_bam=expt_bam,
                    expt_bai=expt_bai,
                    ctl_bam=merged_bam, # optional
                    ctl_bai=merged_bai, # optional
                    basename=basename,
                    effgsize=effgsize
        }


        call motifInputs {
            input:  peaks_narrow=findPeaks.peaks_narrow,
                    basename=basename,
                    genome_fasta=genome_fasta
        }

        call findMotifs {
            input:  summit_seqs=motifInputs.summit_seqs,
                    peak_seqs=motifInputs.peak_seqs,
                    bgmodel=bgmodel,
                    genome_fasta=genome_fasta
        }
    }
}

### Task definitions

# macs2 to call peaks, using optional control.bam for background
task findPeaks {
    File expt_bam
    File expt_bai
    File? ctl_bam
    File? ctl_bai
    String basename
    Int effgsize
    Int threads = 1
    Int memory_gb = 5

    runtime {
        cluster: "cori"
        time: "01:00:00"
        cpu: 1
        memory: "0G"
    }
    command {
        shifter --image=leobaumgart/dap_py2:2.0 find_peaks.sh \
        -i ${expt_bam} \
        ${"-c" + ctl_bam} \
        -n ${basename} \
        -e ${effgsize}
    }
    output {
        File peaks_narrow = "${basename}_peaks.narrowPeak"
        File summits_bed = "${basename}_summits.bed"
        Int numpeaks = read_int("${basename}_numpeaks.txt")
        File macs2_stats = "macs2_stats.txt"
    }
}

# Write fasta files corresponding to the sequences in peak regions
# for inputs in the motif caller (meme).
# Call the same script twice, once to get entire peak region seqs,
# and once to get only seqs the summits +/- 30bp
task motifInputs {
    File peaks_narrow
    String basename
    File genome_fasta
    Int threads = 1
    Int memory_gb = 5

    runtime {
        cluster: "cori"
        time: "01:00:00"
        cpu: 1
        memory: "0G"
    }
    command {
        shifter --image=leobaumgart/dap_py3:2.0 /bin/bash -c \
        " \
        narrowPeak_to_fasta.py \
        -narrowPeak ${peaks_narrow} \
        -out ${basename}_summits.fasta \
        -ref ${genome_fasta} \
        -maxpeaks 100 \
        -extend 30 \
        -weights \
        && \
        narrowPeak_to_fasta.py \
        -narrowPeak ${peaks_narrow} \
        -out ${basename}_peaks.fasta \
        -ref ${genome_fasta} \
        -maxpeaks 100 \
        -extend all \
        -weights \
        "
    }
    output {
        File summit_seqs="${basename}_summits.fasta"
        File peak_seqs="${basename}_peaks.fasta"
    }
}

# meme (from meme-suite) to find motifs in peak regions
# fimo to map the found motifs back to the genome
task findMotifs {
    File summit_seqs
    File peak_seqs
    File bgmodel
    File genome_fasta
    Int threads = 1
    Int memory_gb = 5

    runtime {
        cluster: "denovo"
        time: "00:30:00"
        cpu: 1
        memory: "5G"
        poolname: "test"
        poolsize: 1
    }
    command {
        shifter --image=leobaumgart/dap_py2:2.0 find_motifs.sh \
        -s ${summit_seqs} \
        -p ${peak_seqs} \
        -b ${bgmodel} \
        -r ${genome_fasta}
    }
    output {
        File motifs_summits="meme_summits"
        File motifs_summits_pal="meme_summits_pal"
        File motifs_peaks="meme_peaks"
        File motifs_peaks_pal="meme_peaks_pal"
    }
}

    """
    inputs_contents = """{
    "jgi_dap_leo.adapters": "/global/projectb/sandbox/gaag/bbtools/data/adapters.fa",
    "jgi_dap_leo.genome_fasta": "/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417.fasta",
    "jgi_dap_leo.bt2index_dir": "/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417_bt2index",
    "jgi_dap_leo.bt2index_name": "PsimiaeWCS417",
    "jgi_dap_leo.effgsize": 6169071,
    "jgi_dap_leo.genes_gff": "/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417.gff",
    "jgi_dap_leo.bgmodel": "/global/projectb/sandbox/rnaseq/DAP/leo/cromwell_genomes/PsimiaeWCS417/PsimiaeWCS417.bgmodel",
    "jgi_dap_leo.outdir": ".",
    "jgi_dap_leo.expt_raw_fastqs": [
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TATCAGC-TATCAGC.fastq.gz"
    ],
    "jgi_dap_leo.ctl_raw_fastqs": [
    ],
    "jgi_dap_leo.library_names_map": {
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TATCAGC-TATCAGC.fastq.gz": "CTTZN",
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.GGTTGAT-GGTTGAT.fastq.gz": "CTUGZ",
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TTGCTGG-TTGCTGG.fastq.gz": "CTUHA"
    },
    "jgi_dap_leo.sample_names_map": {
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TATCAGC-TATCAGC.fastq.gz": "TF4",
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.GGTTGAT-GGTTGAT.fastq.gz": "negCtl2",
        "/global/dna/dm_archive/sdm/illumina/01/25/27/12527.1.262232.TTGCTGG-TTGCTGG.fastq.gz": "negCtl1"
    },
    "jgi_dap_leo.expt_bam": "/global/projectb/scratch/jaws/jfroula/leo_dap/CTTZN_TF4.bam",
    "jgi_dap_leo.expt_bai": "/global/projectb/scratch/jaws/jfroula/leo_dap/CTTZN_TF4.bam"
}
"""
    inputs.write_text(inputs_contents)
    wdl_file.write_text(wdl_contents)
    return tmp_path.as_posix()


@pytest.fixture()
def subworkflows_example(tmp_path):
    main = tmp_path / "main.wdl"
    sub1 = tmp_path / "sub1.wdl"
    sub2 = tmp_path / "sub2.wdl"
    inputs_json = tmp_path / "inputs.json"

    main_contents = """
import "sub1.wdl" as sub1
import "sub2.wdl" as sub2

workflow main_wdl {
    Boolean flag

    call run_preprocess { input: preprocess_input = "we are running preprocessing steps" }

    if( flag ){
        call sub1.sub1_workflow { input: sub1_input = "we are running sub-workflow 1" }
    }
    if( ! flag ){
        call sub2.sub2_workflow { input: sub2_input = "we are running sub-workflow 2"}
    }

    output {
        String? main_output_wdl1 = sub1_workflow.run_task1_output
        String? main_output_wdl2 = sub2_workflow.run_task2_output
    }
}

task run_preprocess {
    String preprocess_input

    command
    {
        echo ${preprocess_input}
        }
}
"""
    sub1_contents = """
    workflow sub1_workflow {
    String sub1_input

    call run_task1 { input: sub1_input = "we are running sub-workflow 1" }

    output { String run_task1_output = run_task1.status  }
}

task run_task1 {
    String sub1_input

    command
    {
        echo ${sub1_input}
        }
    output {
        String status = "task1 was run successfully"
    }

    runtime {
        docker: "jfroula/aligner-bbmap:1.1.9"
        poolname: "extrasmall"
        shared: 1
        node: 1
        nwpn: 1
        memory: "4G"
        time: "00:10:00"
    }
}
"""
    sub2_contents = """
    workflow sub2_workflow {
    String sub2_input

    call run_task2 { input: sub2_input = "we are running sub-workflow 2" }

    output { String run_task2_output = run_task2.status  }
}

task run_task2 {
    String sub2_input

    command
    {
        echo ${sub2_input}
        }
    output {
        String status = "task2 was run successfully"
    }

    runtime {
        docker: "jfroula/aligner-bbmap:1.1.9"
        poolname: "extrasmall"
        shared: 1
        node: 1
        nwpn: 1
        memory: "6G"
        time: "00:10:00"
    }
}
"""

    inputs_content = """
{
    "main_wdl.flag": false
}
"""

    main.write_text(main_contents)
    sub1.write_text(sub1_contents)
    sub2.write_text(sub2_contents)
    inputs_json.write_text(inputs_content)

    return tmp_path.as_posix()


@pytest.fixture()
def sample_workflow(tmp_path):
    """
    Creates a directory filled with files that are used for submitting
    to JAWS Central via JAWS Client. It contains a WDL, an inputs.json file
    and also a bunch of txt files.

    The directory is created in a temp path with the following hierarchy

    /workflow
       ----- sample.wdl
       ----- sample.json
       ----- simple.txt
       ----- file_0.txt
       ----- file_1.txt
       ----- file_2.txt
             ...
    /staging
    /globus_basedir

    Also note that in temp path are the staging sub directory and globus base
    directory where files are rsync'd over.

    :param tmp_path:
    :return:
    """
    workflow_dir = tmp_path / "workflow"
    workflow_dir.mkdir()
    wdl_file = workflow_dir / "sample.wdl"
    inputs_file = workflow_dir / "sample.json"

    simple_file = workflow_dir / "simple.txt"
    file_0_map = workflow_dir / "file_0.txt"
    file_1_map = workflow_dir / "file_1.txt"

    file_other_map = workflow_dir / "file_2.txt"
    other_file_map = workflow_dir / "file_3.txt"

    array_file_0 = workflow_dir / "file_arr_0.txt"
    array_file_1 = workflow_dir / "file_arr_1.txt"

    simple_file.write_text("This is a simple file")

    # We create another file since we have a file file mapping
    file_0_map.write_text("This is the key to file file map")
    file_1_map.write_text("This is the value of the file file map")

    file_other_map.write_text("This is file to other map")
    other_file_map.write_text("This is an other to file map")

    array_file_0.write_text("First entry in array of files.")
    array_file_1.write_text("Second entry in array of files.")

    wdl_contents = """
workflow simple_workflow {
    File simple
    Array[File] arr_files
    Map[File, File] file_to_file_map
    Map[File, String] file_to_other_map
    Map[String, File] other_to_file_map

    call print {
        input:  simple=simple,
                arr_files=arr_files,
                file_to_file_map=file_to_file_map,
                file_to_other_map=file_to_other_map,
                other_to_file_map=other_to_file_map
    }
}

task print {
    File simple
    Array[File] arr_files
    Map[File, File] file_to_file_map
    Map[File, String] file_to_other_map
    Map[String, File] other_to_file_map

    command {
        echo ${simple} ${file_to_file_map} ${file_to_other_map} \
        ${other_to_file_map} ${arr_files}
    }
    output {
        String status = "task completed"
    }
}
""" # noqa
    inputs_contents = """{{
  "simple_workflow.other_to_file_map": {{"{0}": "{1}"}},
  "simple_workflow.file_to_file_map":  {{"{2}": "{3}"}},
  "simple_workflow.file_to_other_map": {{"{4}": "{5}"}},
  "simple_workflow.simple": "{6}",
  "simple_workflow.arr_files": ["{7}", "{8}"]
}}
""".format("STRING",
           other_file_map.as_posix(),
           file_0_map.as_posix(),
           file_1_map.as_posix(),
           file_other_map.as_posix(),
           "STRING",
           simple_file.as_posix(),
           array_file_0.as_posix(),
           array_file_1.as_posix()
           )

    wdl_file.write_text(wdl_contents)
    inputs_file.write_text(inputs_contents)

    return tmp_path.as_posix()


@pytest.fixture()
def files_inputs(tmp_path):
    inputs = tmp_path / "files.json"
    wdl = tmp_path / "files.wdl"
    wdl_contents = """
version 1.0
workflow workflow_files {
    input {
        File one
        Array[File] two
        Array[Array[File]] three
        Pair[String, Array[File]] four
        Map[String, Int] five
        Map[String, Array[File]] six
        String seven
        File url
    }
    call task_files
}
task task_files {
    command <<<
        echo "Hello!"
    >>>
}
"""
    wdl.write_text(wdl_contents)
    contents = """{"workflow_files.two": ["/path/1", "2"], "workflow_files.six": {"foo": ["/path/3", "https://path4"], "bar": ["5", "6"]}, "workflow_files.seven": "abc", "workflow_files.one": "ftp://xyz", "workflow_files.three": [["/path/7", "/path/8"], ["9"]], "workflow_files.five": {"bob": 1234, "ali": 3457}, "workflow_files.four": {"Left": "elk", "Right": ["/path/10"]}, "workflow_files.url": "http://abf"}"""
    inputs.write_text(contents)
    return tmp_path.as_posix()


@pytest.fixture()
def staged_files(tmp_path):
    staging_dir = tmp_path / "staging"

    destination_dir = tmp_path / "jaws_site" / "staging" / "src_site"
    destination_dir.mkdir(parents=True)

    wdl_file = tmp_path / "example.wdl"
    inputs_json = tmp_path / "example.json"

    wdl_file.write_text("This is a WDL file")
    inputs_json.write_text("This is a JSON file")

    return staging_dir.as_posix(), destination_dir.as_posix()


@pytest.fixture()
def refdata_wdl(tmp_path):
    wdl = tmp_path / "refdata.wdl"
    wdl_contents = """
workflow refdata {
    File file1
    call runblastplus_sub
}

task runblastplus_sub {
    File ncbi_nt
    command <<<
        echo "hi!"
    >>>
}
"""
    wdl.write_text(wdl_contents)
    return tmp_path.as_posix()


@pytest.fixture()
def refdata_inputs(tmp_path):
    inputs = tmp_path / "inputs.json"
    text_file = tmp_path / "file1.txt"
    text_file.write_text("This is a file.")

    contents = """{{
    "refdata.file1": "{0}",
    "refdata.runblastplus_sub.ncbi_nt": "/refdata/nt"
}}
""".format(text_file)

    inputs.write_text(contents)
    return tmp_path.as_posix()


@pytest.fixture()
def relative_inputs(tmp_path):
    rel_dir = tmp_path / "test"
    rel_dir.mkdir()
    inputs = rel_dir / "relative.json"
    fileX = rel_dir / "fileX.txt"
    fileX.write_text("This is fileX.")
    fileY = tmp_path / "fileY.txt"
    fileY.write_text("This is fileY.")
    fileZ = tmp_path / "fileZ.txt"
    fileZ.write_text("This is fileZ.")
    contents = """{{
    "relative.fileX": "./fileX.txt",
    "relative.fileY": "../fileY.txt",
    "relative.fileZ": "{0}"
}}
""".format(fileZ)

    inputs.write_text(contents)
    wdl = rel_dir / "relative.wdl"
    wdl_contents = """
workflow relative {
    File fileX
    File fileY
    File fileZ
    call task_relative
}

task task_relative {
    command <<<
        echo "hi!"
    >>>
}
"""
    wdl.write_text(wdl_contents)
    return tmp_path.as_posix()


@pytest.fixture()
def struct_inputs(tmp_path):
    struct_dir = tmp_path / "test"
    struct_dir.mkdir()
    inputs = struct_dir / "struct.json"
    apple = struct_dir / "apple.txt"
    apple.write_text("Apple starts with an A.")
    brown = tmp_path / "brown.txt"
    brown.write_text("Brown starts with a B.")
    crown = tmp_path / "crown.txt"
    crown.write_text("Crown starts with a C.")
    contents = """{{
    "test_struct.product_list": [
    {
        "name" : "Apple",
        "batch" : {
            "Left": 2234,
            "Right" : 2020
        },
        "locations": ["{apple}"],
        "info": {
            "manufacture": "India",
            "distribution": "India"
        }
    },
    {
        "name" : "Brown",
        "batch" : {
            "Left": 9876,
            "Right" : 2022
        },
        "locations": ["{brown}"],
        "info": {
            "manufacture": "Germany",
            "distribution": "Spain"
        }
    },
    {
        "name" : "Crown",
        "batch" : {
            "Left": 4506,
            "Right" : 2019
        },
        "locations": ["{crown}"],
        "info": {
            "manufacture": "Spain",
            "distribution": "France"
        }
    }]
}}
""".format(apple=apple, brown=brown, crown=crown)

    inputs.write_text(contents)
    wdl = struct_dir / "struct.wdl"
    wdl_contents = """
version 1.0

struct Product {
    String name
    Pair[Int, Int] batch
    Map[String, String] info
}

workflow test_struct {
    input {
        Array[Product] product_list
    }

    scatter (product in product_list) {
        call write_product_params {
            input:
                product=product
        }
    }
}

task write_product_params {
    input {
        Product product
    }
    command <<<
        echo {product.name}
    >>>
}
"""
    wdl.write_text(wdl_contents)
    return tmp_path.as_posix()


@pytest.fixture()
def refdata_inputs_missing_slash(tmp_path):
    inputs = tmp_path / "inputs.json"
    text_file = tmp_path / "file1.txt"
    text_file.write_text("This is a file.")

    contents = """{{
    "refdata.file1": "{0}",
    "refdata.runblastplus_sub.ncbi_nt": "/refdata"
}}
""".format(text_file)

    inputs.write_text(contents)
    return tmp_path.as_posix()


@pytest.fixture()
def incorrect_wdl(tmp_path):
    wdl = tmp_path / "main.wdl"
    contents = """workflow bbtools {
    ;File reads
    File ref

    call alignment {
       input: fastq=reads,
              fasta=ref
    }
    call samtools {
       input: sam=alignment.sam
   }
}

task alignment {
    File fastq
    File fasta

    command {
        bbmap.sh in=${fastq} ref=${fasta} out=test.sam
    }
    output {
       File sam = "test.sam"
    }
}

task samtools {
    File sam

    command {
       samtools view -b -F0x4 ${sam} | samtools sort - > test.sorted.bam
    }
    output {
       File bam = "test.sorted.bam"
    }
}
"""
    wdl.write_text(contents)
    return wdl


@pytest.fixture()
def no_subworkflows_present(tmp_path):
    """
    Fixture where only the main.wdl is present
    but not the specified workflow
    """
    wdl = tmp_path / "main.wdl"
    contents = """import "alignment.wdl" as align

workflow main_wdl {
    File fastq
    File reference

    # this task calls the sub-workflow named bbmap_shard_wf which
    # is the alignment.wdl.
    # It's output is "merged_bam_file"
    call align.bbmap_shard_wf {
           input: reads = fastq,
                  reference = reference
    }
    call bam_stats {
           input: merged_bam = bbmap_shard_wf.merged_bam_file
    }
}

task bam_stats {
    String merged_bam

    command {
        reformat.sh in=${merged_bam} out=stdout.fq | \
        bbstats.sh in=stdin.fq out=stats
    }

    output {
        File alignment_stats = "stats"
    }

    runtime {
        docker: "jfroula/aligner-bbmap:1.1.9"
        poolname: "extrasmall"
        shared: 1
        node: 1
        nwpn: 1
        memory: "5G"
        time: "00:10:00"
    }
}

"""
    wdl.write_text(contents)
    return wdl


@pytest.fixture()
def output_example(tmp_path):
    run_dir = tmp_path / "run1"
    run_dir.mkdir()
    task_dir = run_dir / "task1"
    task_dir.mkdir()
    task_inputs = task_dir / "inputs"
    task_inputs.mkdir()
    task_infile = task_inputs / "infile"
    infile_contents = "EXAMPLE INPUT"
    task_infile.write_text(infile_contents)
    task_execution = task_dir / "execution"
    task_execution.mkdir()
    task_outfile = task_execution / "stdout"
    outfile_contents = "EXAMPLE OUTPUT"
    task_outfile.write_text(outfile_contents)
    return tmp_path.as_posix()


@pytest.fixture()
def deprecated_mem_example(tmp_path):
    wdl_file = tmp_path / "deprecated_mem.wdl"
    wdl_contents = """
workflow fq_count {
    File fastq_file
    call count_seqs { input: infile = fastq_file }
    output {
        File outfile = count_seqs.outfile
    }
}

task count_seqs {
    File infile
    command <<<
        wc -l ${infile} | perl -ne 'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, " sequences\\n"} else {print STDERR "Invalid Fastq file\\n"}' > num_seqs.txt
    >>>
    output {
        File outfile = "num_seqs.txt"
    }
    runtime {
        mem: "1G"
        time: "00:10:00"
    }
}
    """
    wdl_file.write_text(wdl_contents)
    return tmp_path.as_posix()

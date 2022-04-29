import jaws_condor.add


def test_calculate_node_needed():
    """
    $ condor_q -pr ../fmt_nobatch_id.cpf
    190     1024.0 KB    1024.0 B  1
    191     1024.0 KB    1024.0 B  1
    192      390.6 GB    1024.0 B  4
    193       97.7 GB    1024.0 B  10
    194       97.7 GB    1024.0 B  10
    195       97.7 GB    1024.0 B  10
    196       97.7 GB    1024.0 B  10

    ==> idle_list
    [
    [],
    [[190, 0.0009765625, 1024.0, 1], [191, 0.0009765625, 1024.0, 1], [193, 97.7, 1024.0, 10], [194, 97.7, 1024.0, 10]], [194, 97.7, 1024.0, 10]],  # noqa
    [],
    [[192, 390.6, 1024.0, 4]]
    ]

    ==> total number of nodes needed to sbatch
    max(sum(mem_requested) / each_max_mem, sum(cpu_requested) / each_max_cpu) +  1

    ==> final sbatch numbers
    [0, 2, 0, 1]
    """
    job_list = [[], [[190, 0.0009765625, 1024.0, 1], [191, 0.0009765625, 1024.0, 1], [193, 97.7, 1024.0, 10], [194, 97.7, 1024.0, 10], [195, 97.7, 1024.0, 10], [196, 97.7, 1024.0, 10]], [], [[192, 390.6, 1024.0, 4]]]  # noqa
    r_type = 1  # MEDIUM
    mem_range = ["0-0", "0-364", "0-0", "364-1480"]
    cpu_range = [36, 36, 36, 36]
    summed = [sum(x) for x in zip(*job_list[r_type])]
    assert [1159, 390.801953125, 6144.0, 42] == summed

    ret = jaws_condor.add.calculate_node_needed(job_list[r_type], mem_range[r_type], cpu_range[r_type])
    assert ret == 2

    r_type = 3  # XLARGE
    ret = jaws_condor.add.calculate_node_needed(job_list[r_type], mem_range[r_type], cpu_range[r_type])
    assert ret == 1

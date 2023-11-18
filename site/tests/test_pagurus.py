import argparse
import json
import os
import signal
import sys

import psutil
import pytest
from jaws_site import pagurus

pids = ["1", "2", "3"]


@pytest.fixture
def mock_psutil(monkeypatch):
    def mock_pids(*args):
        return pids

    class MockProcess:
        def __init__(self):
            self.data = {
                "ppid": 123,
                "name": "somename",
                "num_threads": 1,
                "memory_percent": 50,
                "num_fds": 456,
                "cwd": "mypath",
            }

        @staticmethod
        def username():
            return "myuser"

        @staticmethod
        def name():
            return "myname"

        def as_dict(self):
            return self.data

    def mock_process(*args):
        return MockProcess()

    monkeypatch.setattr(psutil, "pids", mock_pids)
    monkeypatch.setattr(psutil, "Process", mock_process)

    return MockProcess()


def test_FileWriter_csv(tmp_path):
    print("Output file path", tmp_path)
    header = ["name1", "name2", "name3"]
    outfile = tmp_path / "test.csv"

    print(f"Output file={outfile}")

    fw = pagurus.FileWriter(outfile=outfile, header=header, write_header=True)
    assert fw.header == header
    assert fw.write_header is True
    assert fw.outfile == outfile
    assert outfile.is_file()
    fw.write(0, 1, 2)
    fw.close()
    with open(outfile) as f:
        contents = f.read()
        print(contents)
        lines = contents.splitlines()
        assert lines[0] == ",".join(header)
        assert lines[1] == "0,1,2"
    fw.outfile.unlink()


def test_FileWriter_csv_envvar(tmp_path):
    header = ["name1", "name2", "name3"]
    outfile = tmp_path / "test.csv"

    print(f"Output file={outfile}")

    # Clear out environment variables we will be testing with
    env = {"testytesty": "testytoasty"}
    fw = pagurus.FileWriter(outfile=outfile, header=header, write_header=True, env=env)

    assert fw.header == header + list(env.keys())
    assert fw.write_header is True
    assert fw.outfile == outfile
    assert outfile.is_file()

    # ignores test and test2 because it's outside of the header
    fw.write(0, 1, 2)

    fw.close()
    with open(outfile) as f:
        contents = f.read()
        print(contents)
        lines = contents.splitlines()
        assert lines[0] == ",".join(header + list(env.keys()))
        assert lines[1] == "0,1,2,testytoasty"
    fw.outfile.unlink()


def test_FileWriter_json(tmp_path):
    header = ["name1", "name2", "name3"]
    outfile = tmp_path / "test.json"

    print("Output file path", tmp_path)

    fw = pagurus.FileWriter(
        outfile=outfile, header=header, write_header=True, jsonout=True
    )
    assert fw.header == header
    assert fw.write_header is False
    assert fw.outfile == outfile
    assert outfile.is_file()
    fw.write(0, 1, 2)
    fw.close()
    with open(outfile) as f:
        contents = f.read()
        print(contents)
        lines = contents.splitlines()
        jsonout = json.dumps(dict(zip(header, [0, 1, 2])))
        assert lines[0] == jsonout
    fw.outfile.unlink()


def test_FileWriter_json_envvar(tmp_path):
    header = ["name1", "name2", "name3"]
    outfile = tmp_path / "test.json"
    env = {"testytesty": "testytoasty"}

    fw = pagurus.FileWriter(
        outfile=outfile,
        header=header,
        write_header=True,
        jsonout=True,
        env=env,
    )
    header = header + list(env.keys())
    assert fw.header == header
    assert fw.write_header is False
    assert fw.outfile == outfile
    assert outfile.is_file()
    fw.write(0, 1, 2)
    fw.close()

    with open(str(fw.outfile), "r") as fh:
        obs_data = json.load(fh)

    exp_data = {"name1": 0, "name2": 1, "name3": 2, "testytesty": "testytoasty"}
    assert obs_data == exp_data

    fw.outfile.unlink()


def test_GracefulKiller_exit_gracefully(tmp_path, monkeypatch):
    class MockSignals:
        def __init__(self, name=None):
            self.name = name

        def __iter__(self):
            return iter([111, 222, 333])

    def mock_signal(*args):
        return

    def mock_exit(*args):
        return

    # monkeypatch signal
    signal.Signals = MockSignals()
    signal.SIGKILL = 444
    signal.SIGSTOP = 555

    monkeypatch.setattr(signal, "signal", mock_signal)
    monkeypatch.setattr(sys, "exit", mock_exit)

    running_dir = tmp_path / "running"
    done_dir = tmp_path / "done"

    running_dir.mkdir(exist_ok=True)
    done_dir.mkdir(exist_ok=True)
    outfile = running_dir / "test.csv"
    moved_outfile = done_dir / "test.csv"

    outfile.touch()
    print(outfile, type(outfile))
    fname = os.path.basename(outfile)
    gk = pagurus.GracefulKiller(running_dir, done_dir, fname)

    signal.Signals = MockSignals
    gk.exit_gracefully("somesig", None)
    assert moved_outfile.is_file()
    assert gk.kill_now is True

    if moved_outfile.is_file():
        moved_outfile.unlink()


def test_get_all_user_procs(mock_psutil, monkeypatch):
    def mock_getpid():
        for pid in pids:
            yield pid

    monkeypatch.setattr(os, "getpid", mock_getpid)

    procs = pagurus.get_all_user_procs("myuser")
    assert procs == [int(p) for p in pids]

    procs = pagurus.get_all_user_procs("otheruser")
    assert procs == []


def test_get_iocounters():
    class MockIOcounter:
        def __init__(self):
            self.read_count = 123
            self.write_count = 456
            self.read_chars = 111
            self.write_chars = 222

    results = pagurus.get_iocounters({"io_counters": MockIOcounter()})
    assert results == (123, 456, 111, 222)


def test_get_meminfo():
    class MockMeminfo:
        def __init__(self):
            self.rss = 123
            self.vms = 456

    results = pagurus.get_meminfo({"memory_info": MockMeminfo()})
    assert results == (123, 456)


def test_get_cputimes():
    class MockCputimes:
        def __init__(self):
            self.user = 123
            self.system = 456
            self.children_user = 111
            self.children_system = 222
            self.iowait = 333
            self.iowait = 444
            self.cpu_num = 555
            self.idle = 666

    inputs = {
        "cpu_times": MockCputimes(),
        "cpu_num": 555,
    }
    results = pagurus.get_cputimes(inputs)
    assert results == (555, 123, 456, 444, 222, 111, 666)

    results = pagurus.get_cputimes({})
    assert results == ("nan", "nan", "nan", "nan", "nan", "nan", "nan")


def test_cmd_data():
    inputs = {"cmdline": "abc,def,ghi"}
    results = pagurus.cmd_data(inputs)
    assert results == "a|b|c|||d|e|f|||g|h|i"

    results = pagurus.cmd_data({"cmdline": ""})
    assert results == "nan"


def test_runner(mock_psutil, tmp_path, monkeypatch):
    class ExitWhileLoopException(Exception):
        pass

    class MockGracefulKiller:
        def __init__(self, *args, **kwargs):
            self.kill_now = False

    def mock_get_all_user_procs(*args, **kwargs):
        return [123]

    def mock_get_cpu_times(*args, **kwargs):
        return [1, 2, 3, 4, 5, 6, 7]

    def mock_get_meminfo(*args, **kwargs):
        return [11, 22]

    def mock_get_get_iocounters(*args, **kwargs):
        return [111, 222, 333, 444]

    def mock_cmd_data(*args, **kwargs):
        return "mycmd"

    pagurus.GracefulKiller = MockGracefulKiller
    pagurus.get_all_user_procs = mock_get_all_user_procs
    pagurus.get_cputimes = mock_get_cpu_times
    pagurus.get_meminfo = mock_get_meminfo
    pagurus.get_iocounters = mock_get_get_iocounters
    pagurus.cmd_data = mock_cmd_data

    global ct
    ct = 0
    # Using sleep function to exit loop inside runner

    def mock_sleep(*args):
        global ct
        ct += 1
        if ct == 2:
            raise ExitWhileLoopException

    pagurus.sleep = mock_sleep

    try:
        pagurus.runner(path=str(tmp_path), move=True, rolling=0)
    except ExitWhileLoopException:
        assert True

    exp_vals = [
        "@timestamp,pid,ppid,name,num_threads,cpu_num,cpu_user,cpu_system,cpu_iowait,cpu_children_system,cpu_children_user,idle,mem_rss,mem_vms,memory_percent,num_fds,read_count,write_count,read_chars,write_chars,cmdline,current_dir\n",  # noqa
        "01-01-2023 10:10:10,123,123,somename,1,1,2,3,4,5,6,7,11,22,50,456,111,222,333,444,mycmd,mypath\n",  # noqa
    ]
    stats_file = os.path.join(tmp_path, "running", "stats.csv")
    with open(stats_file, "r") as fh:
        obs_vals = fh.readlines()

    assert isinstance(obs_vals, list)

    # Replace obs_val date entry
    n = obs_vals[1]
    n = ",".join(v for idx, v in enumerate(n.split(",")) if idx != 0)
    obs_vals[1] = f"01-01-2023 10:10:10,{n}"

    assert obs_vals == exp_vals


def test_main(monkeypatch):
    class MockArgumentParser:
        def __init__(self):
            self.args = {
                "rolling": 1,
                "path": "mypath",
                "outfile": "myfile",
                "rate": 2,
                "user": "someuser",
                "no_header": False,
                "move": False,
                "json": False,
                "debug": False,
                "envvar": {},
            }
            for k in self.args:
                setattr(self, k, self.args[k])

        def add_argument(*args, **kwargs):
            return

        def parse_args(self, *args, **kwargs):
            return self

    def mock_runner(*args, **kwargs):
        return 123

    def mock_argumentparser():
        return MockArgumentParser()

    monkeypatch.setattr(argparse, "ArgumentParser", mock_argumentparser)
    pagurus.runner = mock_runner
    pagurus.main()
    assert True

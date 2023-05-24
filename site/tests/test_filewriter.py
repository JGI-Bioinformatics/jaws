import pytest
import json
import os

from jaws_site.filewriter import FileWriter


def test_FileWriter_csv(tmp_path):
    print("Output file path", tmp_path)
    header = ['name1', "name2", "name3"]
    outfile = tmp_path/f"test.csv"
    fw = FileWriter(outfile=outfile, header=header, write_header=True)
    assert fw.header == header
    assert fw.write_header == True
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
    print("Output file path", tmp_path)
    header = ['name1', "name2", "name3"]
    outfile = tmp_path/f"test.csv"

    environ = {'testytesty': 'test',
               'testytoasty':'test2'}
    fw = FileWriter(outfile=outfile, header=header, write_header=True, env=environ)
    header2 = header + ["testytesty","testytoasty"]

    assert fw.header == header2
    assert fw.write_header == True
    assert fw.outfile == outfile
    assert outfile.is_file()
    fw.write(0, 1, 2, "test", "test2")
    fw.close()
    with open(outfile) as f:
        contents = f.read()
        print(contents)
        lines = contents.splitlines()
        assert lines[0] == ",".join(header2)
        assert lines[1] == '0,1,2,test,test2'
    fw.outfile.unlink()


def test_FileWriter_json(tmp_path):
    print("Output file path", tmp_path)
    header = ['name1',"name2","name3"]
    outfile = tmp_path/f"test.json"
    fw = FileWriter(outfile=outfile, header=header, write_header=True, jsonout=True)
    assert fw.header == header
    assert fw.write_header == False
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
    print("Output file path", tmp_path)
    header = ['name1',"name2","name3"]
    outfile = tmp_path/f"test.json"

    environ = {'testytesty': 'test',
               'testytoasty':'test2'}
    fw = FileWriter(outfile=outfile, header=header, write_header=True, jsonout=True, env=environ)
    header2 = header + ["testytesty", "testytoasty"]

    assert fw.header == header2
    assert fw.write_header == False
    assert fw.outfile == outfile
    assert outfile.is_file()
    fw.write(0, 1, 2,"test", "test2")
    fw.close()
    with open(outfile) as f:
        contents = f.read()
        print(contents)
        lines = contents.splitlines()
        jsonout = json.dumps(dict(zip(header2, [0, 1, 2, "test", "test2"])))
        assert lines[0] == jsonout
    fw.outfile.unlink()


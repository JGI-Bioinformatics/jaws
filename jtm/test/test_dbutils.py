# test_run.py
# This need to change to use MOCK
from jaws_jtm.lib.dbutils import DbSqlMy


def test_db_connection():
    db = DbSqlMy(db="test")
    assert db is not None
    db.close()


def test_db_creation():
    try:
        db = DbSqlMy(db="test")
        db.ddl("CREATE DATABASE IF NOT EXISTS test")
        db.close()
    except Exception:
        assert False
    assert True


def test_execute(capsys):
    try:
        db = DbSqlMy(db="test")
        assert "test" in [i[0] for i in db.selectAll("show databases;")]
        # sys.stderr.write(db.selectAll("show databases;"))
        # captured = capsys.readouterr()
        # assert captured.out == "hello\n"
        db.close()
    except Exception:
        assert False
    assert True

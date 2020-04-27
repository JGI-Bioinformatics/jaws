#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
#
import sys
import atexit
import time

from jaws_jtm.common import logger


class SqlField:
    """This class defines SQL table field"""
    def __init__(self, name=None, type=None, null=True):
        self.name = name
        self.type = type
        self.null = null

    def nullSql(self):
        if self.null:
            return ""
        else:
            return "NOT NULL"

    def tableDef(self):
        return " ".join([(attr if attr else "") for attr in (self.name, self.type, self.nullSql())])


class SqlTable:
    """This class defines SQL table"""
    @classmethod
    def fromDb(klass, db, name):
        """Generate SqlTable instance by querying an actual table in the database.
        @param db DbSql instance
        @param name Table name
        @return new SqlTable instance
        """
        # mnemonic names for indexes into cursor.description list
        I_NAME = 0
        I_TYPE = 1
        I_SIZE = 3

        descr = db.getTableDescr(name)
        fields = [SqlField(name=f[I_NAME], type="%s(%s)" % (f[I_TYPE], f[I_SIZE])) for f in descr]
        return klass(name=name, fields=fields)

    def __init__(self, name, fields):
        self.name = name
        self.fields = fields

    def createSql(self):
        """Return SQL DDL string that constructs this table"""
        return """
        create table %s
        (
        %s
        )
        """ % (self.name, ",\n".join([field.tableDef() for field in self.fields]))

    def insertSql(self):
        """Return SQL DML string that can be passed to Python DB-API cursor.executemany() method"""
        return """
        insert into %s
        (
        %s
        )
        values
        (
        %s
        )
        """ % (self.name, ",\n".join([field.name for field in self.fields]), ",\n".join(["?" for field in self.fields]))


class DbSql(object):
    """
    Wrapper around DB-API Connection class with some convenience methods.
    @todo Convert this to using SQL Alchemy. SQL Alchemy imposed too big abstraction penalty in the past,
    but this might not be the case anymore assuming carefully following its best use practices.
    We will still likely to need our bulk loading methods.
    SQAlchemy also supports a very limited set of backends, and implementing a new backend
    interface seems rather tedious where.
    """

    # Default SQL field definition - a fall-back field type to create
    # when nothing more specific is provided by the user
    def_field = SqlField(name="fld", type="char(40)")

    def __init__(self):
        import mysql.connector as dbmod
        self.dbmod = dbmod
        atexit.register(dbClose, dbObj=self)
        self.debug = 0

        # This will be lazy constructed because the needed DBAPI module is
        # only available after descendant class __init__ is called

        self.numpyTypeMap = None

    def dbapi(self):
        return self.dbmod

    def ddlIgnoreErr(self, *l, **kw):
        curs = self.cursor()
        try:
            curs.execute(*l, **kw)
        except Exception as msg:
            logger.debug((l, kw))
            logger.debug(msg)
        curs.close()

    def ddl(self, sql, ifDialect=None, dropList=tuple(), ignoreError=False, **kw):
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            for dbobj in dropList:
                words = [x.strip().lower() for x in dbobj.strip().split()]
                if words[0] == 'index':
                    assert words[2] == 'on'
                    self.dropIndex(words[1], words[3], rawName=True)
                else:
                    try:
                        dropSql = "drop %s" % (dbobj,)
                        if self.debug > 0:
                            logger.debug(dropSql)
                        curs.execute(dropSql)
                    except Exception:
                        pass
            watch = SqlWatch(sql, self.debug)
            try:
                curs.execute(sql, **kw)
            except Exception as msg:
                if ignoreError:
                    logger.debug((sql, kw))
                    logger.debug(msg)
                else:
                    sys.stderr.write("%s %s" % (sql, kw))
                    raise
            watch()
            curs.close()

    def execute(self, sql, ifDialect=None, debug=False, **kw):
        if debug:
            logger.debug("SQL: {}".format(sql))
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            watch = SqlWatch(sql, self.debug)
            try:
                curs.execute(sql, **kw)
            except Exception:
                sys.stderr.write("%s %s" % (sql, kw))
                raise
            watch()
            return curs
        else:
            return None

    def executemany(self, sql, data, ifDialect=None, **kw):
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            watch = SqlWatch(sql, self.debug)
            try:
                curs.executemany(sql, data, **kw)
            except Exception:
                sys.stderr.write("%s %s" % (sql, kw))
                raise
            watch()
            curs.close()

    def getTableDescr(self, name):
        """
        Return table description as seen by DB API module
        """
        curs = self.execute("select * from %s limit 1" % (name,))
        descr = curs.description
        curs.close()
        return descr

    def dumpCursor(self, cursor):
        pass

    def executeAndPrint(self, sql, **kw):
        curs = self.execute(sql, **kw)
        if curs is not None:
            logger.debug((curs.fetchall()))
            curs.close()

    def executeAndAssertEmpty(self, sql, message=None, **kw):
        self.executeAndAssert(sql, tuple(), message=message, **kw)

    def executeAndAssertZero(self, sql, message=None, **kw):
        self.executeAndAssert(sql, ((0,),), message=message, **kw)

    def selectAll(self, sql, debug=False, **kw):
        """
        Convenience method that does for select statement and execute+fetchall in one step.
        Use for statements with small result set.
        @param sql SQL SELECT statement
        @return result of cursor.fetchall (sequence of tuples)
        """
        if debug:
            logger.debug("SQL: {}".format(sql))
        curs = self.execute(sql, **kw)
        ret = curs.fetchall()
        curs.close()
        return ret

    def selectScalar(self, sql, debug=False, **kw):
        """
        Execute sql that must return a single row with a single column and return result as scalar value.
        """
        if debug:
            logger.debug("SQL: {}".format(sql))
        ret = self.selectAll(sql=sql, **kw)
        assert len(ret) == 1 and len(ret[0]) == 1, "Non-scalar value obtained in 'selectScalar()'"
        ret = ret[0][0]
        return ret

    def selectAsNx1Dict(self, sql, debug=False, **kw):
        """
        Execute sql that must return two columns with Nx1 relation and return result as dict(first->second).
        """
        if debug:
            logger.debug("SQL: {}".format(sql))
        ret = self.selectAll(sql=sql, **kw)
        assert len(ret) == 0 or len(ret[0]) == 2, "Result set must be two columns"
        dret = dict(ret)
        assert len(ret) == len(dret), "Result set must be Nx1 relation. Multi-valued keys were found."
        return dret

    def selectAs1Col(self, sql, debug=False, **kw):
        """
        Execute sql that must return one column and return result as 1D sequence.
        """
        if debug:
            logger.debug("SQL: {}".format(sql))
        ret = self.selectAll(sql=sql, **kw)
        assert len(ret) == 0 or len(ret[0]) == 1, "Result set must have one column"
        return [row[0] for row in ret]

    def dropTables(self, names):
        for name in names:
            self.dropTable(name)

    def dropTable(self, name):
        try:
            self.ddl("drop table " + name)
        except Exception:
            pass

    def dropIndex(self, name, table, rawName=False):
        """
        We always create and drop indices as tablename_indexname because in some DBMS (MonetDB)
        index names should be globally unique.
        """
        try:
            if rawName:
                self.ddl("drop index %s on %s" % (name, table))
            else:
                self.ddl("drop index ix_%s_%s on %s" % (table, name, table))
        except Exception:
            pass

    def createIndices(self, names, table):
        """
        We always create and drop indices as tablename_indexname because in some DBMS (MonetDB)
        index names should be globally unique.
        This can be specialized for MySQL and others that support ALTER TABLE ... ADD INDEX ... ADD INDEX
        """
        for name in names:
            self.dropIndex(name, table)
            self.ddl("create index ix_%s_%s on %s ( %s )" % (table, name, table, name))

    def connection(self):
        return self.con

    def cursor(self):
        return self.con.cursor()

    def commit(self):
        if hasattr(self, 'con'):
            self.con.commit()

    def reconnect(self):
        self.close()
        self.open()

    def open(self):
        pass

    def close(self):
        pass

    def dialectMatch(self, dialect):
        return dialect is None

    def analyze(self, table):
        self.ddl("ANALYZE TABLE " + table, ifDialect="mysql")

    def createTableAs(self, name, select, indices=None):
        """
        Save the results of SQL SELECT as a new table.
        This abstracts "create ... as ..." operation from minor differences in SQL dialects.
        For example, MonetDB Feb2008 required 'with data' at the end, MySQL 5 and SQLite do
        not recognize 'with data'
        Override in derived classes if necessary.

        :param name: name of table to (re-)create - existing table will be replaced
        :param select: select SQL select statement
        :param indices: indices if present, will be passed to createIndices() method
        :return:
        """

        self.ddl("""
        CREATE TABLE %s AS
        """ % (name,) + select + """
        """,
                 dropList=["table %s" % (name,)])
        if self.debug:
            curs = self.execute("select count(*) from %s" % (name,))
            if curs is not None:
                logger.debug(("%s rows created in table %s" % (curs.fetchall(), name)))
                curs.close()
        if indices is not None and len(indices) > 0:
            self.createIndices(table=name, **indices)
            self.ddl("analyze table %s" % name, ifDialect="mysql")


class DbSqlLite(DbSql):
    """
    Derivative of DbSql specific for SQLite DB engine
    """
    def_field = SqlField(name="fld", type="text")

    def __init__(self, db_path, str_type=str, dry_run=False, strategy="default", **kw):
        """
        Constructor.

        :param strategy: A string either of exclusive_unsafe|default; exclusive_unsafe
        attempts to maximize speed for use patterns where the database is created from
        scratch by a single process and will be re-created if the application fails,
        so we can turn off transactions (by switching off journal creation) and lock
        the database file exclusively.
        """
        # use the standard Python module (python >=2.5):
        import sqlite3 as dbmod
        kw = kw.copy()
        self.dbmod = dbmod
        DbSql.__init__(self)
        self.str_type = str_type
        self.db_path = db_path
        self.dry_run = dry_run
        if strategy not in ("default", "exclusive_unsafe"):
            raise ValueError("Unknown value of 'strategy': %s" % (strategy,))
        if strategy == "exclusive_unsafe":
            kw.setdefault("isolation_level", "EXCLUSIVE")
        self.con = self.dbmod.connect(self.db_path, **kw)
        self.con.text_factory = self.str_type
        self.execute("PRAGMA cache_size = 1000000")  # 1GB
        if strategy == "exclusive_unsafe":
            self.execute("PRAGMA journal_mode = OFF")
            self.execute("PRAGMA synchronous = OFF")
            self.execute("PRAGMA temp_store = MEMORY")

    def close(self):
        self.commit()
        if hasattr(self, 'con'):
            self.con.close()
            delattr(self, 'con')

    def dialectMatch(self, dialect):
        return dialect is None or dialect == "sqlite"

    def dropIndex(self, name, table, rawName=False):
        """
        Specialization for SQLite because SQLite does not understand 'on <table_name>' clause.
        """
        if rawName:
            self.ddl("drop index if exists %s" % (name,))
        else:
            self.ddl("drop index if exists ix_%s_%s" % (table, name))

    def createIndices(self, table, names=None, primary=None, compounds=None, attrib={}):
        """
        This is a specialization for SQLite, which does not support ALTER TABLE ... ADD INDEX ... ADD INDEX.

        :param attrib: optional index attributes. Currently supported is 'unique', e.g. attrib={'id':{'unique':True}}
        :note primary cannot be created in sqlite outside of create table statement, so we fake it with a unique index,
        which is still not equal to 'primary' constraint because we cannot add 'not null' constraint.
        """
        drop_names = []
        sql_tpl = "CREATE %%s INDEX %%s ON %s (%%s)" % (table,)
        if primary is not None:
            if names is None:
                names = []
            else:
                names = list(names)
            names.append(primary)
            attrib[primary] = {"unique": True}
        add_index = []
        if names is not None:
            for name in names:
                drop_names.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                add_index.append(sql_tpl % (unique, "ix_%s_%s" % (table, name), name))
        if compounds is not None:
            for name in list(compounds.keys()):
                drop_names.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                add_index.append(sql_tpl % (unique, "ix_%s_%s" % (table, name), compounds[name]))
        for name in drop_names:
            self.dropIndex(name, table)
        for sql in add_index:
            self.ddl(sql)


class DbSqlMysql(DbSql):
    """
    Derivative of DbSql specific for MySQL DB engine
    """

    def __init__(self, config=None, dry_run=False, **kw):
        import mysql.connector as dbmod
        self.dbmod = dbmod

        if config:
            self.config = config
            self.MYSQL_HOST = self.config.configparser.get("MYSQL", "host")
            self.MYSQL_USER = self.config.configparser.get("MYSQL", "user")
            self.MYSQL_PW = self.config.configparser.get("MYSQL", "password")
            self.MYSQL_PORT = self.config.configparser.getint("MYSQL", "port")
            self.MYSQL_DB = self.config.configparser.get("MYSQL", "db")

        DbSql.__init__(self)
        self.open(**kw)

        # TODO: handle situation due to long periods of computations w/o SQL calls:
        #  Exception _mysql_exceptions.OperationalError: (2006, 'MySQL server has gone away')
        #  or
        #  _mysql_exceptions.OperationalError: (2013, 'Lost connection to MySQL server during query')

    def open(self, **kw):
        self.close()
        if "conn" in kw and kw['conn'] is not None:
            self.con = kw['conn']
        else:
            config = {
                "user": self.MYSQL_USER,
                "password": self.MYSQL_PW,
                "host": self.MYSQL_HOST,
                "port": "%d" % self.MYSQL_PORT,
                "database": self.MYSQL_DB}
            self.con = self.dbmod.connect(**config)

    def close(self):
        self.commit()
        if hasattr(self, 'con'):
            self.con.close()
            del self.con

    def dropTable(self, name):
        self.ddl("drop table if exists " + name)

    def dialectMatch(self, dialect):
        return dialect is None or dialect == "mysql"

    def createIndices(self, table, names=None, primary=None, compounds=None, attrib={}):
        """
        This is a specialization for MySQL, which supports ALTER TABLE ... ADD INDEX ... ADD INDEX.

        :param attrib: optional index attributes. Currently supported is 'unique', e.g. attrib={'id':{'unique':True}}
        """
        drop_names = []
        sql = "ALTER TABLE %s " % (table,)
        comma = ""
        if primary is not None:
            sql = sql + "%s\nADD PRIMARY KEY (%s)" % (comma, primary)
            comma = ","
        if names is not None:
            add_index = []
            for name in names:
                drop_names.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                add_index.append("ADD %s INDEX ix_%s_%s (%s)" % (unique, table, name, name))
            if len(add_index) > 0:
                sql = sql + "%s\n" % (comma,) + ",\n".join(add_index)
                comma = ","
        if compounds is not None:
            add_index = []
            for name in list(compounds.keys()):
                drop_names.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                add_index.append("ADD %s INDEX ix_%s_%s (%s)" % (unique, table, name, compounds[name]))
            if len(add_index) > 0:
                sql = sql + "%s\n" % (comma,) + ",\n".join(add_index)
                comma = ","
        for name in drop_names:
            self.dropIndex(name, table)
        self.ddl(sql)


class SqlWatch:
    def __init__(self, sql, debug):
        self.debug = debug
        if self.debug > 0:
            # time.clock() seems to be broken on SuSe 10 x86_64 Python 2.4
            # - it always returns the same value
            self.start = time.time()
            logger.debug(sql)

    def __call__(self):
        if self.debug > 0:
            finish = time.time()
            logger.debug(("SQL finished in %.3f sec" % (finish-self.start)))
            self.start = finish


def dbClose(dbObj):
    dbObj.close()

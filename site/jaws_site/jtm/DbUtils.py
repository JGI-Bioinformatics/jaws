#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)

import sys
import atexit
# import mysql.connector.pooling
import time
import numpy

from Config import *

## Classes that describe definitions of SQL fields and table
## These are used primarily to construct DDL statements.
# class SqlResultAssertionError(StandardError):
#
#     def __init__(self, receive, expect, message=None):
#         self.receive = receive
#         self.expect = expect
#         self.message = message
#
#     def __str__(self):
#         return "SqlResultAssertionError\nMessage:%s\nReceived:\n%s\n\nExpected:\n%s\n" % (
#         self.message, self.receive, self.expect)


class SqlField:
    """This class defines SQL table field"""

    def __init__(self,name=None,type=None,null=True):
        self.name = name
        self.type = type
        self.null = null

    def nullSql(self):
        if self.null:
            return ""
        else:
            return "NOT NULL"

    def tableDef(self):
        return " ".join( [ (attr if attr else "") for attr in (self.name,self.type, self.nullSql()) ] )


class SqlTable:
    """This class defines SQL table"""

    @classmethod
    def fromDb(klass, db, name):
        """Generate SqlTable instance by querying an actual table in the database.
        @param db DbSql instance
        @param name Table name
        @return new SqlTable instance
        """
        ## mnemonic names for indexes into cursor.description list

        I_NAME = 0
        I_TYPE = 1
        I_SIZE = 3
        I_PREC = 4
        I_SCALE = 5

        descr = db.getTableDescr(name)
        fields = [SqlField(name=f[I_NAME], type="%s(%s)" % (f[I_TYPE], f[I_SIZE])) \
                  for f in descr]
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

    ## Default SQL field definition - a fall-back field type to create
    ## when nothing more specific is provided by the user
    defField = SqlField(name="fld", type="char(40)")

    def __init__(self):
        import MySQLdb as dbmod
        self.dbmod = dbmod
        atexit.register(dbClose, dbObj=self)
        self.debug = 0

        # This will be lazy constructed because the needed DBAPI module is
        # only available after descendant class __init__ is called

        self.numpyTypeMap = None

    def dbapi(self):
        return self.dbmod

    def getNumpyTypeMap(self):
        if self.numpyTypeMap is None:
            self.numpyTypeMap = SqlNumpyTypeMap(self.dbapi())
        return self.numpyTypeMap

    def ddlIgnoreErr(self, *l, **kw):
        curs = self.cursor()
        try:
            curs.execute(*l, **kw)
        except StandardError as msg:
            print(l, kw)
            print(msg)
            pass
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
                            print(dropSql)
                        curs.execute(dropSql)
                    except StandardError:
                        pass
            watch = SqlWatch(sql, self.debug)
            try:
                curs.execute(sql, **kw)
            except StandardError as msg:
                if ignoreError:
                    print(sql, kw)
                    print(msg)
                else:
                    sys.stderr.write("%s %s" % (sql, kw))
                    raise
            watch()
            curs.close()

    def execute(self, sql, ifDialect=None, **kw):
        # if not sql.strip().startswith("drop"):
        #    return
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            watch = SqlWatch(sql, self.debug)
            try:
                curs.execute(sql, **kw)
            except StandardError:
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
            except StandardError:
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
            print(curs.fetchall())
            curs.close()

    # def executeAndAssert(self, sql, expect, message=None, **kw):
    #     if self.debug > 0:
    #         print("Asserting that query will return this: %s" % (expect,))
    #     curs = self.execute(sql, **kw)
    #     if curs is None:
    #         raise SqlResultAssertionError(receive=None, expect=expect, message=message)
    #     rows = curs.fetchall()
    #     if len(rows) != len(expect):
    #         raise SqlResultAssertionError(receive=rows, expect=expect, message=message)
    #     for rowReceive, rowExpect in zip(rows, expect):
    #         if list(rowReceive) != list(rowExpect):
    #             raise SqlResultAssertionError(receive=rows, expect=expect, message=message)
    #     curs.close()

    def executeAndAssertEmpty(self, sql, message=None, **kw):
        self.executeAndAssert(sql, tuple(), message=message, **kw)

    def executeAndAssertZero(self, sql, message=None, **kw):
        self.executeAndAssert(sql, ((0,),), message=message, **kw)

    def selectAll(self, sql, **kw):
        """
        Convenience method that does for select statement and execute+fetchall in one step.
        Use for statements with small result set.
        @param sql SQL SELECT statement
        @return result of cursor.fetchall (sequence of tuples)
        """
        curs = self.execute(sql, **kw)
        ret = curs.fetchall()
        curs.close()
        return ret

    def selectScalar(self, sql, **kw):
        """
        Execute sql that must return a single row with a single column and return result as scalar value.
        """
        ret = self.selectAll(sql=sql, **kw)
        assert len(ret) == 1 and len(ret[0]) == 1, "Non-scalar value obtained in 'selectScalar()'"
        ret = ret[0][0]
        return ret

    def selectAsNx1Dict(self, sql, **kw):
        """
        Execute sql that must return two columns with Nx1 relation and return result as dict(first->second).
        """
        ret = self.selectAll(sql=sql, **kw)
        assert len(ret) == 0 or len(ret[0]) == 2, "Result set must be two columns"
        dret = dict(ret)
        assert len(ret) == len(dret), "Result set must be Nx1 relation. Multi-valued keys were found."
        return dret

    def selectAs1Col(self, sql, **kw):
        """
        Execute sql that must return one column and return result as 1D sequence.
        """
        ret = self.selectAll(sql=sql, **kw)
        assert len(ret) == 0 or len(ret[0]) == 1, "Result set must have one column"
        return [row[0] for row in ret]

    def dropTables(self, names):
        for name in names:
            self.dropTable(name)

    def dropTable(self, name):
        try:
            self.ddl("drop table " + name)
        except StandardError:
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
        except StandardError:
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
                print("%s rows created in table %s" % (curs.fetchall(), name))
                curs.close()
        if indices is not None and len(indices) > 0:
            self.createIndices(table=name, **indices)
            self.ddl("analyze table %s" % name, ifDialect="mysql")

    # def makeBulkInserterFile(self, **kw):
    #     kw = copy(kw)
    #     if hasattr(kw["table"], "name"):
    #         kw["table"] = kw["table"].name
    #     kw.setdefault("tmpDir", self.tmpDir)
    #     return BulkInserterFile(self, **kw)

    # def makeBulkInserter(self, *l, **kw):
    #     return BulkInserter(self, *l, **kw)

    # def makeBulkReader(self, *l, **kw):
    #     return BulkReader(db=self, *l, **kw)

    # def saveRecords(self, records, table, indices=None):
    #     """
    #     Create a new table and save records into it.
    #
    #     :param records: iterable with each element been a sequence of field values itself
    #     :param table: instance of SqlTable description class
    #     :param indices:
    #     :return:
    #     """
    #
    #     self.ddl(table.createSql(), dropList=["table " + table.name])
    #     inserter = self.makeBulkInserterFile(table=table.name)
    #     for rec in records:
    #         inserter(rec)
    #     inserter.flush()
    #     if indices is not None and len(indices) > 0:
    #         self.createIndices(table=table.name, **indices)
    #         self.ddl("analyze table %s" % (table.name,), ifDialect="mysql")

    # def createTableFromKeyVal(self, name, records, keyName, keyType, valName, valType, valNull=False, indices=None):
    #     """
    #     Create a table from in-memory key-value records.
    #
    #     :param name: The name of the new table
    #     :param records: Either mapping type or iterable of pair sequences
    #     :param keyName: name of key field to create
    #     :param keyType: SQL type of key field
    #     :param valName: name of value field to create
    #     :param valType: SQL type of value field
    #     :param valNull: Logical flag to allow or not NULL in value field. Only important if
    #     updates to the table are expected afterwards
    #     :param indices: if not None, should be a dictionary with arguments to createIndices,
    #     except the 'table' argument, which will be taken from 'name'. "primary" argument
    #     will be set to keyName.
    #     :return:
    #     """
    #
    #     table = SqlTable(
    #         name,
    #         fields=[
    #             SqlField(
    #                 keyName, type=keyType, null=False
    #             ),
    #             SqlField(
    #                 valName, type=valType, null=valNull
    #             )
    #         ]
    #     )
    #     if isinstance(records, collections.Mapping):
    #         records = records.items()
    #     if indices is None:
    #         indices = dict()
    #     else:
    #         indices = indices.copy()
    #     indices.setdefault("primary", keyName)
    #     self.saveRecords(records, table, indices=indices)

    # def createTableFromArray(self, name, arr, withData=True, returnInserter=False, indices=None):
    #     """
    #     Create a table that reflects the fields of Numpy record array.
    #
    #     :return BulkInserter: object or None
    #     :param name: The name of the new table
    #     :param arr: Numpy array to use as template
    #     :param withData: if True, also load data from array
    #     :param returnInserter: if True, return a BulkInserter object;
    #     The caller then is resposible for closing the inserter object.
    #     :param indices: if not None, should be a dictionary with arguments to createIndices,
    #     except the 'table' argument, which will be taken from 'name'.
    #     All fields are constrained as NOT NULL, as Numpy does not have NULL values,
    #     and NOT NULL constraint speeds up queries.
    #     """
    #     numMap = self.getNumpyTypeMap()
    #     flds = numMap.sqlFromDtype(dtype=arr.dtype)
    #     table = SqlTable(name=name, fields=[SqlField(name=f[0], type=f[1], null=False) for f in flds])
    #     self.ddl(table.createSql(), dropList=["table " + table.name])
    #     inserter = None
    #     if withData or returnInserter:
    #         inserter = self.makeBulkInserterFile(table=table.name)
    #     if withData:
    #         for rec in arr:
    #             inserter(rec)
    #         if not returnInserter:
    #             inserter.close()
    #     if indices is not None and len(indices) > 0:
    #         self.createIndices(table=name, **indices)
    #         self.ddl("analyze table %s" % name, ifDialect="mysql")
    #     return inserter

    # def createTableFromCsv(self, name, csvFile, fieldsMap={}, defField=None, hasHeader=False,
    #                        dialect="excel-tab", dialect_options={}, indices=None, preProc=None, hdrPreProc=None):
    #     """
    #     Create and fill a table from CSV file.
    #     The intention is to provide a single command to load a CSV file into SQL table
    #     where reasonable default values can be generated for all options.
    #
    #     :param name: The name of the new table
    #     :param csv: Either a file name, in which case csv.reader(openCompresed(),dialect=dialect)
    #     will be used to open the file, or it should be an existing file stream, or csv.reader object
    #     (in the latter case dialect and dialect_options parameters are ignored).
    #     :param fieldsMap: A dictionary that either maps field position to SqlField instances
    #     (if hasHeader is False) or maps field names to SqlField instances (if hasHeader is True).
    #     In the latter case, each SqlField instance can replace the name from the header, or leave the
    #     name as None, in which case the name supplied by the header will be used.
    #     :param defField: SqlField instance to generate default SQL definitions for fields not
    #     mapped by fieldsMap. If hasHeader is False, defField.name will be used as a
    #     prefix, such that the field name becomes prefix_xxx where xxx is the absolute field
    #     position in CSV row. Otherwise, defField.name is ignored, and only the other attributes are
    #     used to define default field type etc.
    #     :param hasHeader: tells whether to treat the first row of CSV file as header
    #     :param dialect" Dialect string defined in csv module
    #     :param dialect_options: Passed to csv.reader(**dialect_options)
    #     :param indices: if not None, should be a dictionary with arguments to createIndices,
    #     except the 'table' argument, which will be taken from 'name'
    #     :param preProc: If specified, this should be a method that will be applied to each
    #     row returned by csv.reader. The method must return a sequence (possibly empty) of new
    #     rows, which will be inserted into SQL table instead of the original row. Their size and
    #     types must match the original row. The preProc must have this signature:
    #     preProc(row,fields,nameToInd) where row is returned by csv.reader.next(); fields is a list
    #     of SqlField objects matching the row fields; nameToInd is a dictionary mapping field names
    #     to indexes of fields in the row. The last two parameters allow preProc's code to access row
    #     elements by field names, e.g. row[nameToInd["seqid"]]. The preProc parameter addresses a common
    #     use case where the input file is very large but we need to load into the SQL DB only a small
    #     subset of it for which a simple filter condition exists such as set membership. It also covers
    #     simple manipulation of input data such as various string substitutions.
    #     Example:
    #     preProc = lambda row,fields,nameToInd,idSet=set(1,2,3): \
    #             ( row, ) if row[namesToInd["seqId"]] in idSet else (,)
    #     createTableFromCsv(...,preProc=preProc)
    #     """
    #     if preProc is None:
    #         preProc = lambda row, fields, nameToInd: (row,)
    #     if hdrPreProc is None:
    #         hdrPreProc = lambda row: row
    #     if defField is None:
    #         defField = copy(self.defField)
    #     if isinstance(csvFile, str):
    #         closeCsv = True
    #         csvFileInp = openCompressed(csvFile, "r")
    #         csvFile = csv.reader(csvFileInp, dialect=dialect, **dialect_options)
    #     else:
    #         closeCsv = False
    #         # Here we need to figure out if csvFile is just a file stream, or
    #         # a CSV reader already. Unfortunately, csv module does not specify
    #         # any common base class for CSV readers, so we have to rely on
    #         # tests for attribute presence
    #         if not (hasattr(csvFile, "dialect") and hasattr(csvFile, "line_num")):
    #             # this is NOT a result of calling csv.reader() or compatible interface,
    #             # so we assume it to be a file stream object, and call csv.reader
    #             # on it to create a CSV reader
    #             csvFile = csv.reader(csvFile, dialect=dialect, **dialect_options)
    #
    #     if not hasHeader:
    #         firstRow = csvFile.next()
    #         nFields = len(firstRow)
    #     else:
    #         firstRow = None
    #         nFields = 0
    #         hdr = hdrPreProc(csvFile.next())
    #         # now fieldsMap is assumed to map names to SqlField, and we convert it
    #         # to positional map, checking first that all mapped fields are present
    #         # in the header
    #         flds_hdr = [f.strip() for f in hdr]
    #         flds_hdr_set = set(flds_hdr)
    #         for fldname in fieldsMap:
    #             if not fldname in flds_hdr_set:
    #                 raise ValueError, "fieldsMap argument contains field name that is " + \
    #                                   "not found in the CSV header: %s" % (fldname,)
    #         fieldsMapPos = {}
    #         for (ifld, fldname) in enumerate(flds_hdr):
    #             if fldname in fieldsMap:
    #                 fieldsMapPos[ifld] = fieldsMap[fldname]
    #                 # if not already set, the field name is taken from the header
    #                 if fieldsMapPos[ifld].name is None:
    #                     fieldsMapPos[ifld].name = fldname
    #             else:
    #                 flddef = copy(defField)
    #                 flddef.name = fldname
    #                 fieldsMapPos[ifld] = flddef
    #         fieldsMap = fieldsMapPos
    #         # with header, all fields get defined in fieldsMap
    #         defField = None
    #
    #     fields = buildSqlFields(fieldsMap=fieldsMap, nFields=nFields, defField=defField)
    #     table = SqlTable(name=name, fields=fields)
    #     self.ddl(table.createSql(), dropList=["table " + table.name])
    #     inserter = self.makeBulkInserterFile(table=table)
    #     nameToFieldInd = dict(((field.name, iField) for (iField, field) \
    #                            in enumerate(fields)))
    #     if firstRow is not None:
    #         newRows = preProc(firstRow, fields, nameToFieldInd)
    #         for newRow in newRows:
    #             inserter(newRow)
    #     for row in csvFile:
    #         newRows = preProc(row, fields, nameToFieldInd)
    #         for newRow in newRows:
    #             inserter(newRow)
    #     inserter.close()
    #     if indices is not None and len(indices) > 0:
    #         self.createIndices(table=name, **indices)
    #         self.ddl("analyze table %s" % name, ifDialect="mysql")
    #     if closeCsv:
    #         csvFileInp.close()

    # def selectAsArray(self, sql):
    #     """
    #     Execute SQL select and return the entire result set as Numpy record array.
    #     """
    #     reader = self.makeBulkReader(sql=sql)
    #     ret = reader.allAsArray()
    #     reader.close()
    #     return ret

    # def exportAsCsv(self,
    #                 sql,
    #                 out,
    #                 mode="w",
    #                 withHeader=True,
    #                 bufLen=100000,
    #                 dialect="excel-tab",
    #                 dialect_options={"lineterminator": "\n"},
    #                 comment=None,
    #                 sqlAsComment=False,
    #                 commentEscape='#',
    #                 epilog=None):
    #     """
    #     Excecute SQL and export the result as CSV file.
    #
    #     :param sql: SQL select statement to export results of
    #     :param out: Either file name, or file stream object, or CSV writer object
    #     :param mode: mode parameter to open file if 'out' is a string, otherwise ignored
    #     :param withHeader: If True, write the field names as the header
    #     :param bufLen: Size (in number of records) of the internal memory buffer
    #     used when moving SQL result set into the output file
    #     :param dialect: Dialect string defined in csv module
    #     :param dialect_options: Passed to csv.writer(**dialect_options)
    #     :param comment: if not None, this string will be printed at the top
    #     :param sqlAsComment: if True, will print sql statement as an extra comment line
    #     :param commentEscape: this string will be inserted at the start of every
    #     comment line
    #     :note To output any comments, out should not be a csv.writer instance
    #     :note We set the default lineterminator to Linux style '\n', as opposed to
    #     Python's default of Windows style '\r\n'
    #     """
    #     if isinstance(out, str):
    #         out = openCompressed(out, mode)
    #         doClose = True
    #         w = csv.writer(out, dialect=dialect, **dialect_options)
    #     else:
    #         doClose = False
    #         if not (hasattr(out, "dialect") and hasattr(out, "writerows")):
    #             # this is NOT a result of calling csv.writer() or compatible interface,
    #             # so we assume it to be a file stream object, and call csv.writer()
    #             # on it to create a CSV writer
    #             w = csv.writer(out, dialect=dialect, **dialect_options)
    #         else:
    #             if comment is not None or sqlAsComment:
    #                 raise ValueError("Illegal to write comment lines into csv.writer")
    #             w = out
    #
    #     reader = self.makeBulkReader(sql=sql, format="list", bufLen=bufLen)
    #     chunks = reader.chunks()
    #     names = reader.fieldNames()
    #     if sqlAsComment:
    #         if comment is None:
    #             comment = sql
    #         else:
    #             comment = comment + '\n' + sql
    #     if comment is not None:
    #         for line in ((commentEscape + l + '\n') for l in comment.split('\n')):
    #             out.write(line)
    #     if withHeader:
    #         w.writerow(names)
    #     for rows in chunks:
    #         w.writerows(rows)
    #     reader.close()
    #     if epilog is not None:
    #         out.write(epilog)
    #     if doClose:
    #         out.close()

    # def exportAsPivotCsv(self, sql,
    #                      rowField, colField, valField,
    #                      out,
    #                      mode="w",
    #                      withHeader=True,
    #                      rowFieldOut=None,
    #                      colFieldOrderBy=None,
    #                      restval=0,
    #                      bufLen=100000,
    #                      dialect="excel-tab",
    #                      dialect_options={"lineterminator": "\n"},
    #                      comment=None,
    #                      sqlAsComment=False,
    #                      commentEscape='#',
    #                      epilog=None,
    #                      valFormatStr="%s"):
    #     """
    #     Excecute SQL and export the result as CSV file.
    #
    #     :param sql: SQL select statement to export results of
    #     :param out: Either file name, or file stream object, or CSV writer object
    #     :param mode: mode parameter to open file if 'out' is a string, otherwise ignored
    #     :param restval: What to write for missing values in each cell (default is 0)
    #     :param withHeader: If True, write the field names as the header
    #     :param bufLen: Size (in number of records) of the internal memory buffer
    #     used when moving SQL result set into the output file
    #     :param dialect: Dialect string defined in csv module
    #     :param dialect_options: Passed to csv.writer(**dialect_options)
    #     :param comment: if not None, this string will be printed at the top
    #     :param sqlAsComment: if True, will print sql statement as an extra comment line
    #     :param commentEscape: this string will be inserted at the start of every
    #     comment line
    #     :note To output any comments, out should not be a csv.writer instance
    #     :note We set the default lineterminator to Linux style '\n', as opposed to
    #     Python's default of Windows style '\r\n'
    #     """
    #     if colFieldOrderBy is None:
    #         colFieldOrderBy = colField
    #     cols = self.selectAs1Col("""
    #     select
    #         distinct %s
    #     from
    #         ( %s ) a
    #     order by
    #         %s
    #     """ % (colField, sql, colFieldOrderBy))
    #     if rowFieldOut is None:
    #         rowFieldOut = rowField
    #     assert rowFieldOut.strip().lower() not in [c.strip().lower() for c in cols], \
    #         "%s name for row ID header name conflicts with one of the pivot column names" % (rowFieldOut,)
    #     names = [rowFieldOut] + cols
    #     if isinstance(out, str):
    #         out = openCompressed(out, mode)
    #         doClose = True
    #         w = csv.DictWriter(out, fieldnames=names, restval=restval, dialect=dialect, **dialect_options)
    #     else:
    #         doClose = False
    #         if not (hasattr(out, "dialect") and hasattr(out, "writerows")):
    #             # this is NOT a result of calling csv.writer() or compatible interface,
    #             # so we assume it to be a file stream object, and call csv.writer()
    #             # on it to create a CSV writer
    #             w = csv.DictWriter(out, fieldnames=names, restval=restval, dialect=dialect, **dialect_options)
    #         else:
    #             if comment is not None or sqlAsComment:
    #                 raise ValueError("Illegal to write comment lines into csv.writer")
    #             w = out
    #
    #     if sqlAsComment:
    #         if comment is None:
    #             comment = sql
    #         else:
    #             comment = comment + '\n' + sql
    #     if comment is not None:
    #         for line in ((commentEscape + l + '\n') for l in comment.split('\n')):
    #             out.write(line)
    #     if withHeader:
    #         w.writerow(dict([(name, name) for name in names]))
    #     reader = self.makeBulkReader(sql=sql, format="list", bufLen=bufLen)
    #     rdrCols = reader.fieldNames(format="dict")
    #     # convert KeyError exception into more descriptive messages because naming mismatch between SQL
    #     # output fields and requested key-value fields can be a fairly typical user error
    #     try:
    #         rowFieldInd = rdrCols[rowField]
    #     except KeyError:
    #         raise ValueError("Row key field name is not found in SQL cursor recordset: %s" % (rowField,))
    #     try:
    #         colFieldInd = rdrCols[colField]
    #     except KeyError:
    #         raise ValueError("Column key field name is not found in SQL cursor recordset: %s" % (colField,))
    #     try:
    #         valFieldInd = rdrCols[valField]
    #     except KeyError:
    #         raise ValueError("Value field name is not found in SQL cursor recordset: %s" % (valField,))
    #     for rowKey, rows in it.groupby(reader.rows(), lambda r: r[rowFieldInd]):
    #         rowOut = dict(((row[colFieldInd], valFormatStr % (row[valFieldInd],)) for row in rows))
    #         rowOut[rowFieldOut] = rowKey
    #         w.writerow(rowOut)
    #     reader.close()
    #     if epilog is not None:
    #         out.write(epilog)
    #     if doClose:
    #         out.close()


class DbSqlLite(DbSql):
    """
    Derivative of DbSql specific for SQLite DB engine
    """
    
    defField = SqlField(name="fld", type="text")
    
    def __init__(self,dbpath,strType=str,dryRun=False,strategy="default",**kw):
        """
        Constructor.

        :param strategy: A string either of exclusive_unsafe|default; exclusive_unsafe
        attempts to maximize speed for use patterns where the database is created from
        scratch by a single process and will be re-created if the application fails, 
        so we can turn off transactions (by switching off journal creation) and lock
        the database file exclusively.
        """
        #use the standard Python module (python >=2.5):
        import sqlite3 as dbmod
        kw = kw.copy()
        self.dbmod = dbmod
        DbSql.__init__(self)
        self.strType = strType
        self.dbpath = dbpath
        self.dryRun = dryRun
        if strategy not in ("default", "exclusive_unsafe"):
            raise ValueError("Unknown value of 'strategy': %s" % (strategy,))
        if strategy == "exclusive_unsafe":
            kw.setdefault("isolation_level", "EXCLUSIVE")
        self.con = self.dbmod.connect(self.dbpath, **kw)
        self.con.text_factory = self.strType
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
    
    # def makeBulkInserterFile(self,**kw):
    #     #@todo There is a lot of room for optimizing bulk insertion
    #     #in SQLite. E.g. this
    #     #$echo "create table mytable ( col1 int, col2 int);" | sqlite3 foo.sqlite
    #     #$echo ".import demotab.txt mytable"  | sqlite3 foo.sqlite
    #     #or other methods that use true sql insert but wrap many of them in a
    #     #single transction as well as play with memory settings, described e.g.
    #     #here: http://stackoverflow.com/questions/364017/faster-bulk-inserts-in-sqlite3
    #     k = {}
    #     if kw.has_key("bufLen"):
    #         k["bufLen"] = kw["bufLen"]
    #     k["table"] = kw["table"]
    #     return self.makeBulkInserter(**k)
   
    def dropIndex(self, name, table, rawName=False):
        """
        Specialization for SQLite because SQLite does not understand 'on <table_name>' clause.
        """
        if rawName:
            self.ddl("drop index if exists %s" % (name,))
        else:
            self.ddl("drop index if exists ix_%s_%s" % (table,name))

    def createIndices(self, table, names=None, primary=None, compounds=None, attrib={}):
        """
        This is a specialization for SQLite, which does not support ALTER TABLE ... ADD INDEX ... ADD INDEX.

        :param attrib: optional index attributes. Currently supported is 'unique', e.g. attrib={'id':{'unique':True}}
        :note primary cannot be created in sqlite outside of create table statement, so we fake it with a unique index,
        which is still not equal to 'primary' constraint because we cannot add 'not null' constraint.
        """
        dropNames = []
        sqlTpl = "CREATE %%s INDEX %%s ON %s (%%s)" % (table,)
        if primary is not None:
            if names is None:
                names = []
            else:
                names = list(names)
            names.append(primary)
            attrib[primary] = {"unique":True}
        addindex = []
        if names is not None:
            for name in names:
                dropNames.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                addindex.append(sqlTpl % (unique,"ix_%s_%s" % (table,name),name))
        if compounds is not None:
            for name in compounds.keys():
                dropNames.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                addindex.append(sqlTpl % (unique,"ix_%s_%s" % (table,name),compounds[name]))
        for name in dropNames:
            self.dropIndex(name,table)
        for sql in addindex:
            self.ddl(sql)


class DbSqlMy(DbSql):
    """
    Derivative of DbSql specific for MySQL DB engine
    """
    
    def __init__(self, dryRun=False, **kw):
        import MySQLdb as dbmod
        self.dbmod = dbmod
        DbSql.__init__(self)
        self.open(**kw)
        # self.open_pool(**kw)
        #self.dbmod.server_init(("phyla","--defaults-file=my.cnf"),('server',))


        ## TODO: handle situation due to long periods of computations w/o SQL calls:
        ## Exception _mysql_exceptions.OperationalError: (2006, 'MySQL server has gone away')
        ## or
        ## _mysql_exceptions.OperationalError: (2013, 'Lost connection to MySQL server during query')

    def open(self, **kw):
        import MySQLdb.cursors
        self.close()
        if "conn" in kw and kw['conn'] is not None:
            self.con = kw['conn']
        else:
            self.con = self.dbmod.connect(port=kw.get("port", MYSQL_PORT),
                                          #unix_socket=kw.get("unix_socket", "/tmp/sulsj.mysql.sock"),
                                          host=kw.get("host", MYSQL_HOST),
                                          db=kw.get("db", MYSQL_DB),
                                          user=MYSQL_USER,
                                          passwd=MYSQL_PW,
                                          read_default_file="my.cnf",
                                          read_default_group="client")
                                          # ssl_disabled='True')
                                          ## SSCursor fetches results row-by-row from the server,
                                          ## although I did not see any difference in overall speed
                                          #cursorclass=MySQLdb.cursors.SSCursor)

    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()
            del self.con
        #self.dbmod.server_end()

    def open_pool(self, **kw):
        res = {}
        res["host"] = MYSQL_HOST
        res["port"] = MYSQL_PORT
        res["user"] = MYSQL_USER
        res["password"] = MYSQL_PW
        res["database"] = MYSQL_DB
        self.dbconfig = res
        self.pool = self.create_pool(pool_name=MYSQL_POOLNAME, pool_size=MYSQL_POOLSIZE)
        self.con = self.pool.get_connection()

    def create_pool(self, pool_name=MYSQL_POOLNAME, pool_size=MYSQL_POOLSIZE):
        """
        Create a connection pool, after created, the request of connecting
        MySQL could get a connection from this pool instead of request to
        create a connection.

        :param pool_name: the name of pool, default is "jtm"
        :param pool_size: the size of pool, default is 3
        :return: connection pool
        """
        pool = mysql.connector.pooling.MySQLConnectionPool(
            pool_name=pool_name,
            pool_size=pool_size,
            pool_reset_session=True,
            **self.dbconfig)
        return pool

    def get_pool_conn(self):
        return self.pool.get_connection()

    # def close_pool_conn(self, conn, cursor):
    #     cursor.close()
    #     conn.close()

    def dropTable(self,name):
        self.ddl("drop table if exists " + name)

    def dialectMatch(self,dialect):
        return dialect is None or dialect == "mysql"

    # def makeBulkInserterFile(self,**kw):
    #     kw = copy(kw)
    #     if hasattr(kw["table"],"name"):
    #         kw["table"] = kw["table"].name
    #     return BulkInserterFileMy(self,**kw)

    def createIndices(self, table, names=None, primary=None, compounds=None, attrib={}):
        """
        This is a specialization for MySQL, which supports ALTER TABLE ... ADD INDEX ... ADD INDEX.

        :param attrib: optional index attributes. Currently supported is 'unique', e.g. attrib={'id':{'unique':True}}
        """
        dropNames = []
        sql = "ALTER TABLE %s " % (table,)
        comma = ""
        if primary is not None:
            sql = sql + "%s\nADD PRIMARY KEY (%s)" % (comma,primary)
            comma = ","
        if names is not None:
            addindex = []
            for name in names:
                dropNames.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                addindex.append("ADD %s INDEX ix_%s_%s (%s)" % (unique,table,name,name))
            if len(addindex) > 0:
                sql = sql + "%s\n" % (comma,) + ",\n".join(addindex)
                comma = ","
        if compounds is not None:
            addindex = []
            for name in compounds.keys():
                dropNames.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                    pass
                addindex.append("ADD %s INDEX ix_%s_%s (%s)" % (unique,table,name,compounds[name]))
            if len(addindex) > 0:
                sql = sql + "%s\n" % (comma,) + ",\n".join(addindex)
                comma = ","
        for name in dropNames:
            self.dropIndex(name,table)
        self.ddl(sql)

    # def exportToFile(self,sql1,sql2,fileName,fieldsTerm='|',linesTerm=r'\n'):
    #     rmf(fileName)
    #     self.ddl("""
    #     %s
    #     INTO    OUTFILE '%s'
    #     FIELDS TERMINATED BY '%s'
    #     LINES TERMINATED BY '%s'
    #     %s
    #     """ % (sql1,fileName,fieldsTerm,linesTerm,sql2))

    # def exportToStream(self,sql1,sql2,fieldsTerm='|',linesTerm=r'\n'):
    #     (fobj,bulkFile) = makeTmpFile(dir=self.tmpDir,prefix="sqlBulkExp_",suffix=".csv",createParents=True)
    #     fobj.close()
    #     self.exportToFile(sql1=sql1,sql2=sql2,fileName=bulkFile,fieldsTerm=fieldsTerm,linesTerm=linesTerm)
    #     fobj = open(bulkFile,'r',2**20)
    #     os.remove(bulkFile) # will cause exception on Windows.
    #     return fobj

class SqlWatch:
    def __init__(self, sql, debug):
        self.debug = debug
        if self.debug > 0:
            ##time.clock() seems to be broken on SuSe 10 x86_64 Python 2.4
            ##- it always returns the same value
            self.start = time.time()
            print(sql)

    def __call__(self):
        if self.debug > 0:
            finish = time.time()
            print("SQL finished in %.3f sec" % (finish-self.start))
            self.start = finish

def dbClose(dbObj):
    #print "Running 'atexit()' handler"
    dbObj.close()


class DbSqlMyPool(object):
    def __init__(self, **kw):
        self.open_pool(**kw)

    def open_pool(self, **kw):
        res = {}
        res["host"] = MYSQL_HOST
        res["port"] = MYSQL_PORT
        res["user"] = MYSQL_USER
        res["password"] = MYSQL_PW
        res["database"] = MYSQL_DB
        self.dbconfig = res
        self.pool = self.create_pool(pool_name=MYSQL_POOLNAME, pool_size=MYSQL_POOLSIZE)
        # self.con = self.pool.get_connection()

    def create_pool(self, pool_name=MYSQL_POOLNAME, pool_size=MYSQL_POOLSIZE):
        """
        Create a connection pool, after created, the request of connecting
        MySQL could get a connection from this pool instead of request to
        create a connection.

        :param pool_name: the name of pool, default is "jtm"
        :param pool_size: the size of pool, default = 5
        :return: connection pool
        """
        pool = mysql.connector.pooling.MySQLConnectionPool(
            pool_name=pool_name,
            pool_size=pool_size,
            pool_reset_session=True,
            **self.dbconfig)
        return pool

    def get_pool_conn(self):
        return self.pool.get_connection()

    def close_pool_conn(self):
        # cursor.close()
        # conn.close()
        # self.pool.close()  # Todo: how to close pool explicitly?
        pass


class SqlNumpyTypeMap:
    """Map between SQL data types and NumPy data types.
    Precision in digits vs number of bytes for integer NUMBER types are taken
    from MySQL reference (but should be universal)
    http://dev.mysql.com/doc/refman/5.0/en/numeric-types.html
    and from SQL code:
    db.ddl("create table tmp_types (f_bool bool, f_tinyint tinyint, f_smallint smallint, f_med mediumint, f_int int,"+
    " f_big bigint, f_char char(5), f_float float)",dropList=["table tmp_types"])
    curs = db.execute("select * from tmp_types limit 1")
    print curs.description
    """

    intDigitsBytes = numpy.array([('bool', 1, 1, 1),
                                  ('tinyint', 4, 1, 1),
                                  ('smallint', 6, 2, 2),
                                  ('mediumint', 9, 3, 4),
                                  ('int', 11, 4, 4),
                                  ('bigint', 20, 8, 8)],
                                  dtype=[('typeSql', 'S15'), ('digits', 'i4'), ('bytesSql', 'i4'), ('bytesNpy', 'i4')])

    def __init__(self, dbapi):
        """@param dbapi - Python DBAPI module instance."""
        self.dbapi = dbapi

    def close(self):
        self.dbapi = None

    def digitsToBytes(self, dig):
        return self.intDigitsBytes['bytesNpy'][numpy.digitize([dig - 0.0001], self.intDigitsBytes['digits'])[0]]

    def dtype(self, descr):
        """Take cursor.description object and return NumPy dtype suitable for record array construction.
        dtype is returned in a 'list of tuples' representation e.g. [('id','int8'),('taxid','int4')]
        SQL is case insensitive, so we convert field names to lower case when constructing dtype.
        @bug 'smallint' in MySQL DBAPI has typecode 2, which is not recognized as NUMBER by that DBAPI.
        You should cast 'smallint' fields to some other integer type in your SQL SELECT statement."""

        ## mnemonic names for indexes into cursor.description list

        I_NAME = 0
        I_TYPE = 1
        I_SIZE = 3
        I_PREC = 4
        I_SCALE = 5

        dbapi = self.dbapi
        dt = []
        for fld in descr:
            npy_t = None
            fld_t = fld[I_TYPE]
            if fld_t is None:
                npy_t = 'O'
            elif fld_t == dbapi.NUMBER:
                if fld[I_SCALE] > 0:
                    npy_t = 'f8'
                else:
                    npy_t = 'i%s' % (self.digitsToBytes(fld[I_PREC]),)
            elif fld_t == dbapi.STRING:
                npy_t = 'S%s' % (fld[I_SIZE],)
            elif fld_t == dbapi.DATE or fld_t == dbapi.DATETIME:
                ## DATE is string for numpy
                npy_t = 'S%s' % (fld[I_SIZE],)
            else:
                # just try 'double' a a last resort. MySQL dbapi lacks NUMERIC (DECIMAL) constant
                npy_t = 'f8'
                # raise TypeError("Unsuported SQL DBAPI typecode: %s. Field is %s" % (fld_t,fld))
            dt.append((fld[I_NAME].lower(), npy_t))
        return dt

    def sqlFromDtype(self, dtype):
        sql = []
        for fName in dtype.names:
            fDtype = dtype.fields[fName][0]
            kind = fDtype.kind
            itemsize = fDtype.itemsize
            if kind == 'b':
                spec = 'bool'
            elif kind in 'iu':
                if itemsize <= 4:
                    spec = 'integer'
                else:
                    spec = 'bigint'
            elif kind == 'f':
                spec = 'float(%s)' % itemsize
            elif kind in 'cS':
                spec = 'char(%s)' % itemsize
            else:
                raise ValueError("Unsupported Numpy datatype for SQL conversion: %s" % fDtype)
            sql.append((fName, spec))
            # pdb.set_trace()
        return sql
import click
import sqlite3
import sys

# prints all table names
# cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
# print(cursor.fetchall())
# tables = ['workflow','task','try','node','block','status','resource']

# gets status table schema
# query =  "pragma table_info({})".format('status')
# table_info = cursor.execute(query).fetchall()
# prints:
# [(0, 'task_id', 'INTEGER', 1, None, 1),
#  (1, 'task_status_name', 'TEXT', 1, None, 3),
#  (2, 'timestamp', 'DATETIME', 1, None, 4),
#  (3, 'run_id', 'TEXT', 1, None, 2),
#  (4, 'try_id', 'INTEGER', 1, None, 0)]

# workflow table schema
# query =  "pragma table_info({})".format('workflow')
# table_info = cursor.execute(query).fetchall()
# print(table_info)
# [(0, 'run_id', 'TEXT', 1, None, 1),
#  (1, 'workflow_name', 'TEXT', 0, None, 0),
#  (2, 'workflow_version', 'TEXT', 0, None, 0),
#  (3, 'time_began', 'DATETIME', 1, None, 0),
#  (4, 'time_completed', 'DATETIME', 0, None, 0),
#  (5, 'host', 'TEXT', 1, None, 0),
#  (6, 'user', 'TEXT', 1, None, 0),
#  (7, 'rundir', 'TEXT', 1, None, 0),
#  (8, 'tasks_failed_count', 'INTEGER', 1, None, 0),
#  (9, 'tasks_completed_count', 'INTEGER', 1, None, 0)]


@click.command()
@click.option("-t", "--taskid", help="Parsl task ID")
@click.option("-d", "--db", default="monitoring.db", help="Path to Parsl monitoring DB")
def read_db(taskid, db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    # get most recent (current) parsl run id
    query = 'SELECT run_id FROM WORKFLOW ORDER BY time_began DESC LIMIT 1'
    cursor.execute(query)
    r = cursor.fetchone()
    runid = "\"" + r[0] + "\""  # formatting for SQL query

    taskid = str(taskid)
    query = 'SELECT task_status_name FROM status '
    query += 'WHERE run_id = {} AND task_id = {} '.format(runid, taskid)
    query += 'ORDER BY timestamp DESC LIMIT 1'
    cursor.execute(query)
    r = cursor.fetchone()
    if r is None:
        cursor.close()
        conn.close()
        raise ValueError("Task {} not found in Parsl monitoring DB.".format(taskid))
    print("Task status: " + r[0])

    cursor.close()
    conn.close()

    # setting exit code for cromwell check-alive
    # if still running, exit code 0
    # else (finished or failed), exit code 1
    if r[0] == 'exec_done' or r[0] == 'failed':
        sys.exit(1)


if __name__ == "__main__":
    read_db()

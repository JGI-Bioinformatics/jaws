"""
This class provides:
 * a database of Tasks' workers
 * a database of Tasks' resources used and accounting data
 * a priority queue for Tasks
 * estimation of queue-wait times for a Task
"""

class Tasks

def __init__():
    """
    Init Tasks object.
    """

def add_task():
    """
    Create new task.
    """

def retrieve_by_id():
    """
    Retrieve a Task by id.
    """

def retrieve_next():
    """
    Get the highest priority Task matching given criteria.
    """

def estimated_wait():
    """
    Given a task id, determine it's place in the queue and the estimated wait time.
    The estimated wait is the total max time requested of all higher-priority tasks,
    however tasks with higher priority may be submitted anytime before the task begins.
    """

def update():
    """
    Update a task's parameters, which may affect it's priority.
    """

def delete():
    """
    Delete a task.
    """

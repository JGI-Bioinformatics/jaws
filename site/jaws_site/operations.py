from jaws_site import tasks

operations = {
    "task_log": {
        "required_params": ["cromwell_run_id", "task_dir", "status", "timestamp"],
        "function": tasks.save_task_log,
    },
}

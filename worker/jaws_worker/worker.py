import subprocess
import psutil

def run(cmd, max_min=):
    try:
        proc = subprocess.run(cmd.split(), timeout=max_sec)
    except subprocess.TimeoutExpired as error:
        print(f"Maximum time exceeded.")
        kill(pid)

def kill(parent_pid):
    parent = psutil.Process(parent_pid)
    for child in parent.children(recursive=True):  # or parent.children() for recursive=False
        child.kill()
    parent.kill()

    # WRITE "rc" FILE SO CROMWELL DETECTS FAILED TASK


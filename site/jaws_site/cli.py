#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import os
import click

from jaws_site import database, rpc_server, jawsd, log, models
from jaws_site.dispatch import dispatch
from jaws_site.config import jaws_config


log_file = os.environ["JAWS_SITE_LOG"] if "JAWS_SITE_LOG" in os.environ else "./jaws-site.log"
logger = log.setup_logger(__package__, log_file)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def serve():
    """JAWS-Site"""
    pass


@serve.command()
def server():
    """Start JAWS-Site RPC-server."""
    app = rpc_server.RpcServer()
    app.start_server()


@serve.command()
def daemon():
    """Start JAWS-Site daemon."""
    daemon = jawsd.JAWSd(db)
    daemon.start_daemon()


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def dispatcher():
    """JAWS-Site RPC functions"""
    pass


@dispatcher.command()
def server_status():
    """Check Cromwell server status."""
    params = {}
    result = dispatch("server_status", params)
    print(result)


@dispatcher.command()
@click.argument("cromwell_id")
def run_metadata(cromwell_id: str):
    """Get Cromwell metadata for a run."""
    params = {"cromwell_id": cromwell_id}
    result = dispatch("run_metadata", params)
    print(result)


@dispatcher.command()
@click.argument("cromwell_id")
def task_status(cromwell_id: str):
    """Get status of a run."""
    params = {"cromwell_id": cromwell_id}
    result = dispatch("task_status", params)
    print(result)


@dispatcher.command()
@click.argument("cromwell_id")
def task_ids(cromwell_id: str):
    """Retrieve the task IDs for a run."""
    params = {"cromwell_id": cromwell_id}
    result = dispatch("task_ids", params)
    print(result)


@dispatcher.command()
@click.argument("cromwell_id")
def cancel_run(cromwell_id: str):
    """Cancel a run."""
    params = {"cromwell_id": cromwell_id}
    result = dispatch("cancel_run", params)
    print(result)


@dispatcher.command()
@click.argument("cromwell_id")
def run_logs(cromwell_id: str):
    """Get the logs for a run."""
    params = {"cromwell_id": cromwell_id}
    result = dispatch("run_logs", params)
    print(result)


@dispatcher.command()
@click.argument("cromwell_id")
def failure_logs(cromwell_id: str):
    """Get the logs for failed tasks."""
    params = {"cromwell_id": cromwell_id}
    result = dispatch("failure_logs", params)
    print(result)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def jaws_daemon():
    """JAWS-Site RPC functions"""
    pass


@jaws_daemon.command()
def check_runs():
    """Check if any runs need action."""
    conf = jaws_config.conf
    daemon = jawsd.JAWSd(conf)
    daemon.check_runs()


@jaws_daemon.command()
@click.argument("user_id")
@click.argument("submission_id")
@click.argument("dest_endpoint")
@click.argument("dest_dir")
def submit_run(user_id: str, submission_id: str, dest_endpoint: str, dest_dir: str):
    """Submit a new run."""
    conf = jaws_config.conf
    daemon = jawsd.JAWSd(conf)
    user = daemon.session.query(models.User).get(user_id).one_or_none()
    if user is None:
        raise jawsd.DataError(f"User does not exist: {user_id}")
    run = models.Run(
        user_id=user_id,
        status="staged",
        submission_id=submission_id,
        upload_task_id="NA",
        globus_endpoint=dest_endpoint,
        outdir=dest_dir
    )
    daemon.session.add(run)
    daemon.session.commit()
    daemon.submit_run(run)


cli = click.CommandCollection(sources=[serve, dispatcher, jaws_daemon])


if __name__ == "__main__":
    db = database.JawsDb(jaws_config)
    cli()

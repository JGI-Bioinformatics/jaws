#!/usr/bin/env python3
import argparse
import bz2
import gzip
from math import ceil

from time import sleep
from datetime import datetime
import os
import logging
import sys
import signal
import json

from typing import Dict, List
from pathlib import Path

# This script is a copy of the pagurus script in github: https://github.com/tylern4/pagurus


# Exit if psutil is not installed
try:
    import psutil
except ImportError:
    print("Error: Install psutil to get statistics!", file=sys.stderr)
    exit(0)


class FileWriter:
    compressed = False

    def __init__(self, outfile,
                 header: List[str] = [""],
                 write_header: bool = True,
                 rolling: bool = False,
                 jsonout: bool = False,
                 env: Dict = {}) -> None:
        self.extensions = {
            'gz': 'csv.gz',
            'bz2': 'csv.bz2',
            'csv': 'csv'
        }

        self.header: List[str] = header
        self.number: int = 0
        self.write_header: bool = write_header
        self.rolling: bool = rolling
        self.env: Dict = env

        # Create an appropriate formatting function
        if jsonout:
            # Formatter function that outputs a dictionary in JSON
            if header == [""]:
                raise Exception("header cannot be blank for JSON output")
            if write_header:
                logging.debug("forcing write_header to false due to --json flag")
                self.write_header = False

            # Adds envs to the end of the dict for fmt_writer
            def fmt(*args):
                temp = dict(zip(self.header, args))
                temp.update(env)
                return "{}\n".format(json.dumps(temp))
            self.fmt_func = lambda *args: fmt(*args)

        else:
            # Make formatter function based on number of metrics in header
            fmt = ",".join(["{}" for _ in range(len(self.header))])
            fmt_writer = fmt + "\n"
            self.fmt_func = lambda *args: fmt_writer.format(*args)
        self.outfile: Path = outfile
        self.next_file()

    def write(self, *args):
        self.output_file.write(self.fmt_func(*args))

    def flush(self):
        self.output_file.flush()

    def close(self):
        logging.info(f"Closing {self.outfile}")
        self.output_file.close()

    def _open_file(self):
        extenstion = self.outfile.name.split('.')[-1]
        if extenstion == 'gz':
            self.compressed = True
            logging.info(f"Writing gzip pagurus gzip file {self.outfile}")
            self.output_file = gzip.open(self.outfile, 'wt')
        elif extenstion == 'bz2':
            self.compressed = True
            logging.info(f"Writing bz2 pagurus bz2 file {self.outfile}")
            self.output_file = bz2.open(self.outfile, 'wt')
        else:
            logging.info(f"Writing pagurus file {self.outfile}")
            self.output_file = open(self.outfile, "w")

    def _write_header(self):
        self.output_file.write(self.fmt_func(*self.header))
        self.output_file.flush()

    def _renamer(self):
        pass

    def next_file(self):
        if self.number != 0:
            self.close()

        if self.rolling:
            # Split name into it's parts
            name_split = self.outfile.as_posix().split('.')
            # Get the right extention
            ext = self.extensions[name_split[-1]]
            # Add in the file number into the name
            new_outfile_name = f"{name_split[0]}.{self.number}.{ext}"
            # Replace the path
            self.outfile = Path(new_outfile_name)

        self._open_file()
        if self.write_header:
            self._write_header()
        self.number += 1


class GracefulKiller:
    """
    Kills the process graefully when it gets a signal from the OS.
    https://stackoverflow.com/questions/18499497/how-to-process-sigterm-signal-gracefully
    """
    kill_now = False

    def __init__(self, old_prefix=None, new_prefix=None, filename=None):
        # Catch for all signals and exit gracefully
        good_sigs = set(signal.Signals) - {signal.SIGKILL, signal.SIGSTOP}
        for sig in good_sigs:
            signal.signal(sig, self.exit_gracefully)

        self.filename = filename
        # Use pathlib features
        self.old_prefix = Path(f"{old_prefix}")
        self.running_file = self.old_prefix / f"{self.filename}"

        self.new_prefix = new_prefix

    def exit_gracefully(self, sig, frame):
        if self.filename != None:
            self.running_file.rename(f"{self.new_prefix}/{self.filename}")
            logging.info(f"Killed with {signal.Signals(sig).name}")
            logging.info(
                f"Moving {self.running_file} to {self.new_prefix}/{self.filename}")

        logging.info("Killing Process")
        self.kill_now = True
        sys.exit(0)


def get_all_user_procs(username=None) -> List[int]:
    # Get all processes running on the system
    total_procs = psutil.pids()
    # Make a list to return the processes to return
    user_procs = []

    # If username is not given get it from the os
    if username is None:
        try:
            username = os.getlogin()
        except OSError:
            user = "runner"

    total_procs = [int(p) for p in total_procs]
    for p in total_procs:
        # In case the process stoped between
        try:
            user = psutil.Process(p).username()
            name = psutil.Process(p).name()
        except psutil.NoSuchProcess:
            continue

        # Add user processes to list to return
        if username == user and p != os.getpid() and name != 'slurm_script':
            user_procs.append(p)

    return user_procs


def get_iocounters(pData: Dict):
    if 'io_counters' in pData and pData['io_counters'] is not None:
        read_count = pData['io_counters'].read_count
        write_count = pData['io_counters'].write_count
        read_chars = pData['io_counters'].read_chars
        write_chars = pData['io_counters'].write_chars
    else:
        read_count = "nan"
        write_count = "nan"
        read_chars = "nan"
        write_chars = "nan"
    return read_count, write_count, read_chars, write_chars


def get_meminfo(pData):
    rss = pData['memory_info'].rss
    vms = pData['memory_info'].vms

    return rss, vms


def get_cputimes(pData):
    # pcputimes(user=0.05, system=0.02, children_user=0.0, children_system=0.0, iowait=0.0)
    if 'cpu_times' in pData and pData['cpu_times'] is not None:
        user = pData['cpu_times'].user
        system = pData['cpu_times'].system
        try:
            children_user = pData['cpu_times'].children_user
            children_system = pData['cpu_times'].children_system
        except:
            children_system = "nan"
            children_user = "nan"

        try:
            iowait = pData['cpu_times'].iowait
            cpu_num = pData['cpu_num']
        except:
            iowait = "nan"
            cpu_num = "nan"

        try:
            idle = pData['cpu_times'].idle
        except:
            idle = "nan"

    else:
        user = "nan"
        system = "nan"
        iowait = "nan"
        cpu_num = "nan"
        idle = "nan"
        children_system = "nan"
        children_user = "nan"

    return cpu_num, user, system, iowait, children_system, children_user, idle


# @lru_cache
def cmd_data(pData):
    """
    Gets data from the command line arguments and serilaizes it for csv
    """
    cmd = pData['cmdline']
    if len(cmd) == 0:
        return "nan"
    cmd = [c.replace(",", "|") for c in cmd]

    return "|".join(cmd)


def runner(
        path: str = ".", filename: str = "stats.csv",
        pole_rate: float = 0.1, username: str = "",
        write_header: bool = True,
        move: bool = False, rolling: int = 0,
        json: bool = False, env: Dict = {}):
    """
    Runs while your executable is still running and logs info
    about running process to the output file, defaulting to CSV format
    unless the --json flag is set
    Args:
        outfile (str, optional): output filename. Defaults to "stats.csv".
        poleRate (float, optional): Time to sleep before getting new data. Defaults to 0.1.
        username (str, optional): User name to look for when getting statistics.
    """

    # Sets moving files if rolling is true
    if move and rolling > 0:
        move = False

    sleep(1)

    out_dir = Path(f"{path}")
    out_dir.mkdir(exist_ok=True)

    if move:
        killer = GracefulKiller(old_prefix=f"{path}/running",
                                new_prefix=f"{path}/done",
                                filename=filename)
        running = out_dir/"running"
        running.mkdir(exist_ok=True)
        (out_dir/"done").mkdir(exist_ok=True)
        outfile = running/f"{filename}"
    else:
        killer = GracefulKiller()
        outfile = out_dir/f"{filename}"

    header = ["@timestamp", "pid", "ppid", "name", "num_threads", "cpu_num",
              "cpu_user", "cpu_system", "cpu_iowait",
              "cpu_children_system", "cpu_children_user", "idle",
              "mem_rss", "mem_vms", 'memory_percent',
              "num_fds", "read_count", "write_count", "read_chars",
              "write_chars", "cmdline", "current_dir"]

    stats_file = FileWriter(outfile=outfile, header=header,
                            write_header=write_header,
                            rolling=True if rolling > 0 else False,
                            jsonout=json, env=env)
    itteration = 0
    # Keep pulling data from the process while it's running
    while not killer.kill_now:
        user_procs = get_all_user_procs(username=username)
        for proc_num in user_procs:
            try:
                proc = psutil.Process(proc_num)
                pData = proc.as_dict()

                # Add new line to the file with relevant data
                stats = [datetime.now().strftime("%m-%d-%Y %H:%M:%S.%f"),
                         proc_num,
                         pData['ppid'],
                         pData['name'],
                         pData['num_threads'],
                         *get_cputimes(pData),
                         *get_meminfo(pData),
                         pData['memory_percent'],
                         pData['num_fds'],
                         *get_iocounters(pData),
                         cmd_data(pData),
                         pData['cwd']]

                stats_file.write(*stats)

            except psutil.NoSuchProcess as e:
                # Comes when a process is killed between getting the number and getting the data
                pass
            except AttributeError as e:
                # logging.debug(f'Error ({type(e).__name__}): {e}')
                pass
            except TypeError as e:
                # logging.debug(f'Error ({type(e).__name__}): {e}')
                pass
            except Exception as e:
                logging.error(f'Error ({type(e).__name__}): {e}')
                pass
                # Breaks out of just the loop and not the function

        # Write and Sleep for a number of seconds before going to the next loop
        stats_file.flush()
        sleep(pole_rate)

        if rolling > 0:
            itteration += 1
            if itteration % ceil(rolling/pole_rate) == 0:
                stats_file.next_file()

    # Finally we close the file.
    stats_file.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", type=str,
                        help="File name for csv.",
                        default="stats.csv")
    parser.add_argument("-p", "--path", type=str,
                        help="Path to put csv file.",
                        default=".")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="Run with debugging info.",
                        default=False)
    parser.add_argument("-r", "--rate", type=float,
                        help="Polling rate for process.", default=0.1)
    parser.add_argument("-u", "--user", type=str,
                        help="Username to get stats for.", default=None)
    parser.add_argument("-noh", "--no-header",
                        help="Turn off writting the header.", default=True, action='store_false')
    parser.add_argument("-mv", "--move",
                        help="Moves file from 'running' to 'done' directories", default=False, action='store_true')
    parser.add_argument("-l", "--rolling", type=int, help="Time to roll file over to number to file name in ~minutes.",
                        default=0)
    parser.add_argument("--json", default=False, action="store_true", help="Output JSON strings instead of CSV lines")
    parser.add_argument("--envvar", action="append", default=[],
                        help="add environment var to output (can be specified multiple times)")

    args = parser.parse_args()

    # Turn on logging if in debug mode
    if args.debug:
        logging.basicConfig(
            format='%(asctime)s %(levelname)s ==> %(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.FATAL)

    logging.debug("Running in debug mode")

    # 10 gives ~a minute for swapping files
    # so rolling*10 should be okay for minutes (on average)
    rolling = args.rolling * 10

    # Get's the environment variables once and places them into a dict
    env = {ev: os.getenv(ev) for ev in args.envvar}

    # Start the recorder
    runner(path=args.path, filename=args.outfile, pole_rate=args.rate,
           username=args.user, write_header=args.no_header, move=args.move,
           rolling=rolling, json=args.json, env=env)


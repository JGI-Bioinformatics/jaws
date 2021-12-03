#!/usr/bin/env python3
import argparse

from time import sleep
from datetime import datetime
import os
import logging
import sys
import signal
from typing import List

# Exit if psutil is not installed
try:
    import psutil
except ImportError:
    print("Error: Install psutil to get statistics!", file=sys.stderr)
    exit(0)


class GracefulKiller:
    """
    Kills the process graefully when it gets a signal from the OS.
    https://stackoverflow.com/questions/18499497/how-to-process-sigterm-signal-gracefully
    """
    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, *args):
        logging.info("Killing Process\n\n")
        self.kill_now = True


def get_all_user_procs(username=None) -> List[int]:
    # Get all processes running on the system
    total_procs = psutil.pids()
    # Make a list to return the processes to return
    user_procs = []

    # If username is not given get it from the os
    if username is None:
        username = os.getlogin()

    total_procs = [int(p) for p in total_procs]
    for p in total_procs:
        # In case the process stoped between
        try:
            user = psutil.Process(p).username()
        except psutil.NoSuchProcess:
            continue

        # Add user processes to list to return
        if username == user:
            user_procs.append(p)

    return user_procs


def get_iocounters(pData):
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
        children_user = pData['cpu_times'].children_user
        children_system = pData['cpu_times'].children_system
        iowait = pData['cpu_times'].iowait
    else:
        user = "nan"
        system = "nan"
        children_user = "nan"
        children_system = "nan"
        iowait = "nan"

    return user, system, children_user, children_system, iowait


def runner(outfile: str = "stats.csv", poleRate: float = 0.1, username: str = ""):
    """
    Runs while your executable is still running and logs info
    about running process to a csv file.

    Args:
        outfile (str, optional): output filename. Defaults to "stats.csv".
        poleRate (float, optional): Time to sleep before getting new data. Defaults to 0.1.
        username (str, optional): User name to look for when getting statistics.
    """

    sleep(2)

    stats_file = open(outfile, "w")

    header = ["datetime", "pid", "username", "name", "num_threads",
              "cpu_percent", "cpu_t_user", "cpu_t_system", "cpu_num",
              "user", "system", "children_user", "children_system",
              "iowait", "mem_rss", "mem_vms", 'memory_percent',
              'num_fds', "read_count", "write_count", "read_chars",
              "write_chars"]
    # Make formater based on number of metrics in heased
    ftm_writer = ",".join(["{}" for _ in range(len(header))])
    ftm_writer = ftm_writer + "\n"

    stats_file.write(ftm_writer.format(*header))
    stats_file.flush()
    # Keep pulling data from the process while it's running
    killer = GracefulKiller()
    while not killer.kill_now:
        user_procs = get_all_user_procs(username=username)
        for proc_num in user_procs:
            try:
                proc = psutil.Process(proc_num)
                pData = proc.as_dict()

                # Add new line to the file with relevant data
                stats_file.write(ftm_writer.format(
                    datetime.now().strftime("%m-%d-%Y %H:%M:%S.%f"),
                    proc_num,
                    pData['username'],
                    pData['name'],
                    pData['num_threads'],
                    pData['cpu_percent'],
                    pData['cpu_times'].user,
                    pData['cpu_times'].system,
                    pData['cpu_num'] if 'cpu_times' in pData else "nan",
                    *get_cputimes(pData),
                    *get_meminfo(pData),
                    pData['memory_percent'],
                    pData['num_fds'],
                    *get_iocounters(pData)
                ))

            except psutil.NoSuchProcess as e:
                print('Error:', e)
                continue
            except Exception as e:
                print('Error:', e)
                # Breaks out of just the loop and not the function

        # Write and Sleep for a number of seconds before going to the next loop
        stats_file.flush()
        sleep(poleRate)

    # Finally we close the file.
    stats_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tag", type=str,
                        help="Tags the process and gives to the statistcs csv file.",
                        default="stats")
    parser.add_argument("-o", "--outfile", type=str,
                        help="Full path to output file name for csv.",
                        default=None)
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

    args = parser.parse_args()

    # Turn on logging if in debug mode
    if args.debug:
        logging.basicConfig(
            format='%(asctime)s %(levelname)s ==> %(message)s', level=logging.INFO)
    else:
        logging.basicConfig(level=logging.FATAL)

    # Get current time
    nowtime = datetime.now().strftime("%m-%d-%Y-%H.%M.%S")

    # Make out output stats csv name
    if args.outfile is not None:
        outfile = args.outfile
    elif args.tag is not None:
        outfile = f'{args.path}/{args.tag}_{nowtime}.csv'
    else:
        outfile = f'{args.path}/stats_{nowtime}.csv'

    logging.info(f'Saving csv to {outfile}')

    # Start the recorder
    runner(outfile=outfile, poleRate=args.rate, username=args.user)

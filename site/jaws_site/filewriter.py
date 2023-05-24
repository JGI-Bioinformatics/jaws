"""
Filewriter class used by Nick Tyler's Pagarus utility
"""

import logging
import bz2
import gzip
import json
from pathlib import Path
from typing import Dict, List


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

        self.header: List[str] = header + env.keys()
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
            self.fmt_func = lambda *args: fmt_writer.format(*args + env.values())
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

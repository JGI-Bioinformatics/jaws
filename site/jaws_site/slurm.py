"""
This module provides a collection of functions called by rpc_operations.py.
"""


#import requests
import logging
#import os
#import glob
#import json
#import io
#from datetime import datetime
#from dateutil import parser
#from collections import deque
#import boto3
#import botocore


def status(logger):
        """
        Check the status of all the sites by using the "sbatch --test-only" command which returns an estimatated 
        start time for a job. The jobs are asking for the same resources as condor asks for when creating slurm pools.  
        Each site has a combination of  small, medium, large, xlarge pools.
        """

        logger.debug(f"Getting status of ")
        
		return {"sites": "cmds"}


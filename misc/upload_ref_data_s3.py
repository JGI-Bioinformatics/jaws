#!/usr/bin/env python3

import sys
import os
from typing import List, Tuple
import boto3
from botocore.exceptions import ClientError



def get_paths(src_dir: str, dest_dir: str) -> Tuple[List[str], List[str]]:
    """
    :param src_dir: Root dir of transfer
    :param dest_dir: Root destination dir
    :return: list of src and destination file paths
    :rtype: list
    """
    src_dir = os.path.abspath(src_dir)
    src_paths = []
    dest_paths = []
    for root, dirs, files in os.walk(src_dir):
        for file in files:
            src_file = os.path.join(root, file)
            dest_file = os.path.normpath(os.path.join(dest_dir, os.path.relpath(root, src_dir), file))
            src_paths.append(src_file)
            dest_paths.append(dest_file)
    return src_paths, dest_paths

def transfer(aws_access_key_id, aws_secret_access_key, aws_region_name, s3_bucket, src_dir, dest_dir):
    s3_session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name
    )
    s3_obj = s3_session.resource('s3')
    bucket_obj = s3_obj.Bucket(s3_bucket)

    src_paths, dest_paths = get_paths(src_dir, dest_dir)
    for src_path, dest_path in zip(src_paths, dest_paths):
        print(f"XFER: {src_path} -> {dest_path}")
        with open(src_path, "rb") as fh:
            bucket_obj.upload_fileobj(fh, dest_path)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Required parameters: src_dir dest_dir")
    src_dir = sys.argv[1]
    dest_dir = sys.argv[2]
    print(f"Upload {src_dir} -> {dest_dir}")
    aws_access_key_id = os.environ["AWS_ACCESS_KEY_ID"]
    aws_secret_access_key = os.environ["AWS_SECRET_ACCESS_KEY"]
    aws_s3_bucket = os.environ["AWS_S3_BUCKET"]
    aws_region_name = os.environ["AWS_REGION_NAME"]
    if not aws_access_key_id:
        sys.exit("AWS_ACCESS_KEY_ID env var required")
    if not aws_secret_access_key:
        sys.exit("AWS_SECRET_ACCESS_KEY env var required")
    if not aws_s3_bucket:
        sys.exit("AWS_S3_BUCKET env var required")
    if not aws_region_name:
        sys.exit("AWS_REGION_NAME env var required")
    transfer(aws_access_key_id, aws_secret_access_key, aws_region_name, aws_s3_bucket, src_dir, dest_dir)

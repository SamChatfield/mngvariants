"""General S3 functionality.
"""
import os

import boto3
import botocore.exceptions
from botocore.handlers import disable_signing


class S3Connection:
    def __init__(self, url, bucket_name, aws_config_file_path=None, debug=False):
        if aws_config_file_path:
            os.environ['AWS_CONFIG_FILE'] = str(aws_config_file_path)

        if debug:
            boto3.set_stream_logger('')

        self.s3_obj = boto3.resource('s3', endpoint_url=url)
        self.s3_obj.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
        self.bucket = self.s3_obj.Bucket(bucket_name)

    def is_valid_path(self, path):
        """Check if the given path exists in S3."""
        try:
            self.bucket.Object('{}'.format(path)).load()
            return True
        except botocore.exceptions.ClientError as err:
            if err.response['Error']['Code'] in ['403', '404']:
                return False
            raise

    def download_file(self, s3_path, out_path):
        """Download the file at s3_path to the out_path location locally"""
        # TODO: Check if file exists in S3 and throw exception if not
        file_s3_obj = self.bucket.Object('{}'.format(s3_path))

        if out_path.exists():
            raise IOError('File already exists at {}'.format(out_path))
        else:
            out_path_str = str(out_path.resolve())
            file_s3_obj.download_file(out_path_str)

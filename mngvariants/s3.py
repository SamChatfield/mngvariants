"""General S3 functionality.
"""
import boto3
import botocore.exceptions
from botocore.client import Config


class S3Connection:
    def __init__(self, url, bucket_name, s3_config, debug=False):
        if debug:
            boto3.set_stream_logger('')

        config_obj = Config(s3=s3_config)

        self.s3_obj = boto3.resource('s3', endpoint_url=url, config=config_obj)
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
            out_path_str = str(out_path.expanduser().resolve())
            file_s3_obj.download_file(out_path_str)

    def upload_file(self, local_path, s3_path, acl='public-read'):
        """Upload the file at local_path to S3 at s3_path"""
        file_s3_obj = self.bucket.Object('{}'.format(s3_path))

        local_path_str = str(local_path.expanduser().resolve())
        file_s3_obj.upload_file(local_path_str, ExtraArgs={'ACL': acl})

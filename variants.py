import argparse
from pathlib import Path
import requests
import boto3
from arg_parser import arg_parser
from botocore.handlers import disable_signing
from datetime import datetime

# boto3.set_stream_logger('')
s3 = boto3.resource('s3', endpoint_url='***REMOVED***')
s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
s3_bucket = s3.Bucket('microbesng')
# bucket.download_file('b0e1596842_20180911_Forde1/reads.zip', 'reads.zip')

# Base URL for LIMS
lims_base_url = '***REMOVED***'
# Base URL for S3
s3_base_url = 'http://microbesng.s3.climb.ac.uk'

# S3 / Boto3 Setup
# boto3.set_stream_logger('')
# Tell Boto3 to use the local aws config file which gives region and path style
aws_config_file_path = Path(sys.path[0], '.aws', 'config')
os.environ['AWS_CONFIG_FILE'] = str(aws_config_file_path)
s3 = boto3.resource('s3', endpoint_url='***REMOVED***')
s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
s3_bucket = s3.Bucket('microbesng')

def is_valid_results_path(results_path):
    """Check if the provided results path is present in S3 by sending a HEAD request."""
    results_res = requests.head('{}/{}/data.html'.format(s3_base_url, results_path))
    return results_res.status_code == 200

def get_results_path(uuid):
    """Get the results path of a project by querying LIMS using the UUID."""
    # Get the JSON representation of the project from LIMS
    project_res = requests.get(
        '{}/layout/project_api/uuid={}.json'.format(lims_base_url, uuid),
        params={ 'RFMkey': '***REMOVED***' }
    )
    # Raise an exception if response was HTTP error code
    project_res.raise_for_status()
    # Extract the results path from the response
    results_path = project_res.json()['data'][0]['results_path']

    # Check that the results path actually exists, if it doesn't raise an Exception
    if not is_valid_results_path(results_path):
        raise Exception('The returned results path was not valid')

    return results_path

def download_reads(project_dir, results_path):
    """Download the reads zip for this project from S3"""
    reads_zip_path = project_dir / 'reads.zip'
    reads_s3_obj = s3_bucket.Object('{}/reads.zip'.format(results_path))
    if reads_zip_path.is_file():
        print('Reads already downloaded for this project')
        # TODO: Check if the local reads are up to date with the S3 reads (non-trivial)
    else:
        print('Downloading reads...')
        path_str = str(reads_zip_path.resolve())
        reads_s3_obj.download_file(path_str)
        print('Download complete')
    return reads_zip_path

def unzip_samples():
    """Unzip only the required samples from the reads zip"""
    pass

def main(args):
    # Get the S3 results path from the LIMS
    results_path = get_results_path(args.uuid)
    print('Results Path: {}'.format(results_path))

    # Create a directory inside workspace/projects/ for this project by UUID
    project_dir = args.workspace / 'projects' / args.uuid
    try:
        project_dir.mkdir(parents=True)
        print('Project directory {} created'.format(project_dir))
    except FileExistsError:
        print('Project Directory {} already exists'.format(project_dir))
    
    # Download the reads zip for this project from S3
    reads_zip_path = download_reads(project_dir, results_path)

if __name__ == '__main__':
    # Parse the command line arguments against the valid arguments defined in arg_parser.py
    args = arg_parser().parse_args()
    print('UUID: {}, type: {}'.format(args.uuid, type(args.uuid)))
    print('Samples: {}, type: {}'.format(args.samples, type(args.samples)))
    print('Reference: {}, type: {}'.format(args.reference, type(args.reference)))
    print('Workspace: {}, type: {}'.format(args.workspace, type(args.workspace)))
    main(args)

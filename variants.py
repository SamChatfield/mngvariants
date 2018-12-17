import argparse
import os
import re
import shutil
import sys
import urllib.request
import zipfile
from itertools import chain
from pathlib import Path

import boto3
import botocore.exceptions
import pandas as pd
import requests
from botocore.handlers import disable_signing

from arg_parser import arg_parser

# Base URL for LIMS
lims_base_url = '***REMOVED***'
# Base URL for S3
s3_base_url = 'http://microbesng.s3.climb.ac.uk'
# RefSeq Bacteria Assembly Summary URL
refseq_bacteria_assembly_summary_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'

# S3 / Boto3 Setup
# boto3.set_stream_logger('')
# Tell Boto3 to use the local aws config file which gives region and path style
aws_config_file_path = Path(sys.path[0], '.aws', 'config')
os.environ['AWS_CONFIG_FILE'] = str(aws_config_file_path)
s3 = boto3.resource('s3', endpoint_url='***REMOVED***')
s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
s3_bucket = s3.Bucket('microbesng')

def is_valid_results_path(results_path):
    """Check if the provided results path exists in S3."""
    try:
        s3_bucket.Object('{}/data.html'.format(results_path)).load()
        return True
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] in ['403', '404']:
            return False
        else:
            raise

def get_results_path(uuid):
    """Get the results path of a project by querying LIMS using the UUID."""
    # Get the JSON representation of the project from LIMS
    project_res = requests.get(
        '{}/layout/project_api/uuid={}.json'.format(lims_base_url, uuid),
        params={ 'RFMkey': '***REMOVED***' },
        timeout=10
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
    """Download the reads zip for this project from S3."""
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

def unzip_samples(project_dir, reads_zip_path, samples):
    """Unzip only the required samples' forward reads and reverse reads from the zip."""
    # Create the reads directory if it doesn't already exist
    reads_dir = project_dir / 'reads'
    try:
        reads_dir.mkdir()
    except FileExistsError:
        print('Reads directory {} already exists'.format(reads_dir))
    
    # Compute the list of sample files to extract from the zip, two for each sample
    sample_filenames = list(chain.from_iterable([
        ('{}_1_trimmed.fastq.gz'.format(s), '{}_2_trimmed.fastq.gz'.format(s)) for s in samples
    ]))
    
    # Open the reads zip file
    with zipfile.ZipFile(reads_zip_path) as reads_zip:
        # Compute the paths within the zip of the samples to extract
        zip_filepaths = []
        for sample_filename in sample_filenames:
            sample_filepaths = [Path(filepath) for filepath in reads_zip.namelist() if sample_filename in filepath]
            zip_filepaths += sample_filepaths
        # Extract all of the files identified, except those that have already been extracted
        print('Extracting:')
        for zip_filepath in zip_filepaths:
            reads_local_file = reads_dir / zip_filepath.name
            if reads_local_file.is_file():
                print('Sample {} already extracted, skipping'.format(reads_local_file.name))
            else:
                with open(reads_local_file, 'wb') as local_file:
                    print('Zip {}\n-> Local {}'.format(zip_filepath, reads_local_file))
                    local_file.write(reads_zip.read(str(zip_filepath)))

def ftp_to_https(ftp_url):
    """Convert a URL from FTP to HTTPS protocol."""
    https_url = re.sub(r'^ftp:\/\/(.+)$', 'https://\\1', ftp_url)
    print('FTP {}\n-> HTTPS {}'.format(ftp_url, https_url))
    return https_url

def get_refseq_url(reference):
    """Get the RefSeq directory HTTPS URL for the given reference by reading the bacteria assembly summary."""
    print('Getting RefSeq URL for reference "{}"...'.format(reference))
    assembly_summary_data = pd.read_table(refseq_bacteria_assembly_summary_url, header=1, index_col=0, dtype={
        'relation_to_type_material': str
    })
    assembly_summary_data.rename_axis('assembly_accession', axis='index', inplace=True)
    try:
        reference_data = assembly_summary_data.loc[reference]
    except KeyError:
        raise Exception('Error: reference not found in RefSeq bacteria assembly_summary')
    refseq_dir_ftp = reference_data['ftp_path']
    refseq_dir_https = ftp_to_https(refseq_dir_ftp)
    print('RefSeq dir for this reference: {}'.format(refseq_dir_https))
    return refseq_dir_https

def download_file(url, local_path):
    """Download file from the URL to the local path."""
    print('Download: {} -> {}'.format(url, local_path))
    with urllib.request.urlopen(url) as response, open(local_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)

def add_reference_to_config(config_file, reference, refseq_url):
    """Add the snpEff.config entries for the reference."""
    assembly_report_data = pd.read_table(
        '{}/{}_assembly_report.txt'.format(refseq_url, refseq_url.split('/')[-1]),
        comment='#',
        names=[
            'Sequence-Name',
            'Sequence-Role',
            'Assigned-Molecule',
            'Assigned-Molecule-Location/Type',
            'GenBank-Accn',
            'Relationship',
            'RefSeq-Accn',
            'Assembly-Unit',
            'Sequence-Length',
            'UCSC-style-name'
        ]
    )
    # Extract the list of accessions without duplicates and in order
    accessions = list(dict.fromkeys(assembly_report_data['RefSeq-Accn']))
    # Write the lines for this reference with these accessions to snpEff.config
    with open(config_file, 'a') as cfg:
        new_lines = [
            '',
            '{0}.genome : {0}'.format(reference),
            '\t{}.chromosome : {}'.format(reference, ', '.join(accessions)),
        ]
        new_lines += ['\t{}.{}.codonTable : Bacterial_and_Plant_Plastid'.format(reference, accession) for accession in accessions]
        new_lines = ['{}\n'.format(l) for l in new_lines]
        cfg.writelines(new_lines)

def get_reference(workspace, reference):
    """Get the reference genome sequence and annotations from refseq if we don't already have them."""
    print('Get reference: {}'.format(reference))

    config_file = workspace / 'snpEff.config'
    # Check that the snpEff.config file exists
    if not config_file.is_file():
        raise Exception('Workspace directory does not contain snpEff.config')

    reference_dir = workspace / 'references' / reference
    sequences_file = reference_dir / 'sequences.fa.gz'
    genes_file = reference_dir / 'genes.gff.gz'

    # Create the directory for this reference in workspace/snpEff/ if it doesn't already exist
    try:
        reference_dir.mkdir()
    except FileExistsError:
        print('Reference directory {} already exists'.format(reference_dir))

    # Check if sequences.fa.gz and genes.gff.gz are already downloaded
    (download_sequences, download_genes) = (not sequences_file.is_file(), not genes_file.is_file())

    # Check if snpEff.config already contains the configuration for this reference
    regex = r'^{}.genome'.format(re.escape(reference))
    with open(config_file) as cfg:
        modify_config = not re.search(regex, cfg.read(), flags=re.MULTILINE)

    # Only get the RefSeq URL if we absolutely have to because it's slow
    if download_sequences or download_genes or modify_config:
        # Get the RefSeq directory URL for the reference
        refseq_url = get_refseq_url(reference)

        # Download sequences.fa.gz
        if download_sequences:
            sequences_url = '{}/{}_genomic.fna.gz'.format(refseq_url, refseq_url.split('/')[-1])
            download_file(sequences_url, sequences_file)
        # Download genes.gff.gz
        if download_genes:
            genes_url = '{}/{}_genomic.gff.gz'.format(refseq_url, refseq_url.split('/')[-1])
            download_file(genes_url, genes_file)
        
        # Add reference to snpEff.config
        if modify_config:
            print('Adding snpEff.config entry for reference {}'.format(reference))
            add_reference_to_config(config_file, reference, refseq_url)
    
    return reference_dir

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

    # Unzip required samples to reads directory
    unzip_samples(project_dir, reads_zip_path, args.samples)

    # Get reference genome and annotations from RefSeq
    ref_dir = get_reference(args.workspace, args.reference)

    # TODO: index reference (can run this on gz and it gives equivalent indices to ungz - tested by md5 hash)
    # TODO: check if alignment then works and is equivalent whether run against gz indicies or ungz indicies
    # TODO: do alignment as in gobwa script

if __name__ == '__main__':
    # Parse the command line arguments against the valid arguments defined in arg_parser.py
    args = arg_parser().parse_args()
    print('UUID: {}, type: {}'.format(args.uuid, type(args.uuid)))
    print('Samples: {}, type: {}'.format(args.samples, type(args.samples)))
    print('Reference: {}, type: {}'.format(args.reference, type(args.reference)))
    print('Workspace: {}, type: {}'.format(args.workspace, type(args.workspace)))
    main(args)

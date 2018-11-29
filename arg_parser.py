import argparse
from pathlib import Path

def convert_dir(dir_string):
    path = Path(dir_string).expanduser().resolve()
    if not path.is_dir():
        raise argparse.ArgumentTypeError('Workspace directory {} not found'.format(dir_string))
    return path

def arg_parser():
    # Create the parser with the help text to display in -h
    parser = argparse.ArgumentParser(
        description='Variant calling from MicrobesNG LIMS and S3 against references from refseq'
    )
    # Parse UUID
    parser.add_argument(
        '--uuid',
        type=str,
        required=True,
        help='project UUID'
    )
    # Parse list of samples
    parser.add_argument(
        '--samples',
        nargs='+',
        type=str,
        required=True,
        help='list of samples on which to perform variant calling'
    )
    # Parse reference genome
    parser.add_argument(
        '--reference',
        type=str,
        required=True,
        help='reference genome against which to perform variant calling'
    )
    # Parse the workspace directory
    parser.add_argument(
        '--workspace',
        default='~/variant_ws',
        type=convert_dir,
        help='directory to use as a workspace for variant calling, default ~/variant_ws'
    )
    return parser

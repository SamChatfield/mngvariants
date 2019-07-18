import argparse
import logging
from pathlib import Path

DEFAULT_WORKSPACE = '~/variant_ws/'


def convert_dir(dir_string):
    path = Path(dir_string).expanduser().resolve()
    if not path.is_dir() and dir_string == DEFAULT_WORKSPACE:
        path.mkdir()
    elif not path.is_dir():
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
    # Parse workspace directory
    parser.add_argument(
        '--workspace',
        default=DEFAULT_WORKSPACE,
        type=convert_dir,
        help='directory to use as a workspace for variant calling, default ~/variant_ws'
    )
    # Verbose flag - set log level to INFO
    parser.add_argument(
        '-v', '--verbose',
        action="store_const", dest="loglevel", const=logging.INFO,
        default=logging.WARNING,
        help="Increase verbosity of output"
    )
    # Debug flag - set log level to DEBUG
    parser.add_argument(
        '-d', '--debug',
        action="store_const", dest="loglevel", const=logging.DEBUG,
        help="Print debug information"
    )
    return parser


def parse():
    return arg_parser().parse_args()

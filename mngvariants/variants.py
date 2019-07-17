import csv
import gzip
import json
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from itertools import chain
from pathlib import Path

import pandas as pd
import pkg_resources
import requests
from dotenv import load_dotenv

from . import s3, util

# Load environment variables from .env file
load_dotenv(dotenv_path=Path(Path(pkg_resources.resource_filename(__name__, '.env'))))

# Setup console logging
LOG_FORMAT = '%(asctime)s [%(levelname)s] : %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter(LOG_FORMAT))
logger.addHandler(console_handler)

# Base URL for LIMS
LIMS_BASE_URL = os.environ['LIMS_BASE_URL']
# URL and bucket for S3
S3_URL = os.environ['S3_URL']
S3_BUCKET_NAME = os.environ['S3_BUCKET']
# RefSeq Bacteria Assembly Summary URL
REFSEQ_ASSEMBLY_SUMMARY_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'

# S3 configuration
S3_CONFIG = {
    'addressing_style': 'virtual'
}
# Create S3Connection object to handle all S3 interactions
S3_CONN = s3.S3Connection(S3_URL, S3_BUCKET_NAME, S3_CONFIG)

# Name of the variant calling results zip file
RESULTS_ZIP_NAME = 'variants_new'


def get_results_path(uuid):
    """Get the results path of a project by querying LIMS using the UUID."""
    # Get the JSON representation of the project from LIMS
    try:
        project_res = requests.get(
            '{}/layout/project_api/uuid={}.json'.format(LIMS_BASE_URL, uuid),
            params={'RFMkey': os.environ['LIMS_RESTFM_KEY']},
            timeout=10
        )
    except requests.exceptions.ConnectTimeout:
        print('Error: Could not connect to LIMS', file=sys.stderr)
        sys.exit(1)

    # Raise an exception if response was HTTP error code
    project_res.raise_for_status()

    # Extract the results path from the response
    results_path = project_res.json()['data'][0]['results_path']

    # Check that the results path actually exists, if it doesn't raise an Exception
    if not S3_CONN.is_valid_path('{}/data.html'.format(results_path)):
        raise Exception('The returned results path was not valid')

    return results_path


def download_reads(project_dir, results_path):
    """Download the reads zip for this project from S3."""
    reads_s3_path = '{}/reads.zip'.format(results_path)
    reads_zip_path = project_dir / 'reads.zip'

    if reads_zip_path.is_file():
        print('Reads already downloaded for this project, skipping')
        # TODO: Check if local reads are up to date with S3 reads (non-trivial)
    else:
        print('Downloading reads...')
        S3_CONN.download_file(reads_s3_path, reads_zip_path)
        print('Download reads complete')

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
    return reads_dir


def get_refseq_url(reference):
    """Get the RefSeq directory HTTPS URL for the given reference by reading the bacteria assembly
    summary.
    """
    print('Getting RefSeq URL for reference "{}"...'.format(reference))
    assembly_summary_data = pd.read_table(REFSEQ_ASSEMBLY_SUMMARY_URL, header=1, index_col=0, dtype={
        'relation_to_type_material': str
    })
    assembly_summary_data.rename_axis('assembly_accession', axis='index', inplace=True)
    try:
        reference_data = assembly_summary_data.loc[reference]
    except KeyError:
        raise Exception('Error: reference not found in RefSeq bacteria assembly_summary')
    refseq_dir_ftp = reference_data['ftp_path']
    refseq_dir_https = util.ftp_to_https(refseq_dir_ftp)
    print('RefSeq dir for this reference: {}'.format(refseq_dir_https))
    return refseq_dir_https


def get_accessions(refseq_url):
    """Get a list of the accessions for the reference from the given refseq_url's
    assembly_report.txt file. Return the accessions with the version suffixes removed.
    """
    print('Getting accessions from {}...'.format(refseq_url))
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

    # Extract the list of accessions from the assembly report without duplicates and in order
    accessions = list(dict.fromkeys(assembly_report_data['RefSeq-Accn']))
    # Strip off the version suffixes
    accessions_stripped = [accn.split('.')[0] for accn in accessions]

    return accessions_stripped


def download_file(url, local_path):
    """Download file from the URL to the local path."""
    print('Download: {} -> {}'.format(url, local_path))
    with urllib.request.urlopen(url) as response, open(local_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)


def extract_file(in_gzip_path, out_file_path):
    """Extract gzip file at the path in_gzip_path to the path at out_file_path."""
    print('Extract: {} -> {}'.format(in_gzip_path, out_file_path))
    with gzip.open(in_gzip_path, 'rb') as in_file:
        with open(out_file_path, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)


def rewrite_accessions(file_path, accessions):
    """Rewrite instances of the accessions in the file given by file_path without the version
    suffix.
    """
    regex = r'({})\.\d+'.format('|'.join(accessions))
    with open(file_path, 'r+') as f:
        updated_contents = re.sub(regex, '\\1', f.read())
        f.seek(0)
        f.write(updated_contents)
        f.truncate()


def add_reference_to_config(config_file, reference, refseq_url, accessions):
    """Add the snpEff.config entries for the reference."""
    print('Adding snpEff.config entry for reference {}...'.format(reference))
    # Write the lines for this reference with the accessions to snpEff.config
    with open(config_file, 'a') as cfg:
        new_lines = [
            '',
            '{0}.genome : {0}'.format(reference),
            '\t{}.chromosome : {}'.format(reference, ', '.join(accessions)),
        ]
        new_lines += [
            '\t{}.{}.codonTable : Bacterial_and_Plant_Plastid'.format(reference, accn)
            for accn
            in accessions
        ]
        new_lines = ['{}\n'.format(l) for l in new_lines]
        cfg.writelines(new_lines)


def build_snpeff_database(config_file, references_dir, reference):
    print('Building SnpEff database...')
    subprocess.call([
        'snpEff',
        'build',
        '-c {}'.format(config_file.resolve()),
        '-dataDir {}'.format(references_dir.resolve()),
        '-genbank',
        '-v',
        '{}'.format(reference)
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def get_reference(workspace, reference):
    """Get the reference genome sequence and annotations from refseq if we don't already have them.
    """
    print('Get reference {}...'.format(reference))

    config_file = workspace / 'snpEff.config'
    # Check that the snpEff.config file exists
    if not config_file.is_file():
        raise FileNotFoundError('Workspace directory does not contain snpEff.config')

    references_dir = workspace / 'references'
    reference_dir = references_dir / reference
    sequences_gzip = reference_dir / 'sequences.fa.gz'
    sequences_file = reference_dir / 'sequences.fa'
    genes_gzip = reference_dir / 'genes.gbk.gz'
    genes_file = reference_dir / 'genes.gbk'

    # Create the directory for this reference in workspace/references/ if it doesn't already exist
    try:
        reference_dir.mkdir(parents=True)
    except FileExistsError:
        print('Reference directory {} already exists'.format(reference_dir))

    # Check if sequences.fa.gz and genes.gbk.gz are already downloaded
    (download_sequences, download_genes) = (not sequences_gzip.is_file(), not genes_gzip.is_file())
    # Check if the sequences and/or genes need to be extracted again
    (extract_sequences, extract_genes) = (
        download_sequences or not sequences_file.is_file(),
        download_genes or not genes_file.is_file()
    )

    # Check if snpEff.config already contains the configuration for this reference
    regex = r'^{}.genome'.format(re.escape(reference))
    with open(config_file) as cfg:
        modify_config = not re.search(regex, cfg.read(), flags=re.MULTILINE)

    # Only get the RefSeq URL if we absolutely have to because it's slow
    if download_sequences or download_genes or modify_config:
        # Get the RefSeq directory URL for the reference
        refseq_url = get_refseq_url(reference)
        accessions = get_accessions(refseq_url)

        # Download sequences.fa.gz
        if download_sequences:
            sequences_url = '{}/{}_genomic.fna.gz'.format(refseq_url, refseq_url.split('/')[-1])
            download_file(sequences_url, sequences_gzip)
        # Extract sequences.fa.gz to sequences.fa and rewrite accessions without version suffixes
        if extract_sequences:
            extract_file(sequences_gzip, sequences_file)
            rewrite_accessions(sequences_file, accessions)

        # Download genes.gbff.gz as genes.gbk.gz
        if download_genes:
            genes_url = '{}/{}_genomic.gbff.gz'.format(refseq_url, refseq_url.split('/')[-1])
            download_file(genes_url, genes_gzip)
        # Extract genes.gbk.gz to sequences.gbk and rewrite accessions without version suffixes
        if extract_genes:
            extract_file(genes_gzip, genes_file)
            rewrite_accessions(genes_file, accessions)

        # Add reference to snpEff.config
        if modify_config:
            add_reference_to_config(config_file, reference, refseq_url, accessions)

        build_snpeff_database(config_file, references_dir, reference)

    return (sequences_file, genes_file)


def index_sequences(sequences_file):
    print('Indexing sequences file {}...'.format(sequences_file))
    subprocess.call(['bwa', 'index', '{}'.format(sequences_file)], stderr=subprocess.DEVNULL)


def align(project_dir, reads_dir, sequences_file, samples):
    print('Aligning reads to reference...')
    cpu_count = multiprocessing.cpu_count()
    alignment_files = []

    for sample in samples:
        print('Aligning sample {}'.format(sample))
        fwd_read = reads_dir / '{}_1_trimmed.fastq.gz'.format(sample)
        rev_read = reads_dir / '{}_2_trimmed.fastq.gz'.format(sample)
        sorted_sample_file = project_dir / '{}.sorted.bam'.format(sample)
        alignment_files.append(sorted_sample_file)

        bwa_mem = subprocess.Popen([
            'bwa',
            'mem',
            '-t{}'.format(cpu_count),
            '{}'.format(sequences_file),
            '{}'.format(fwd_read),
            '{}'.format(rev_read),
        ], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        samtools_view = subprocess.Popen([
            'samtools',
            'view',
            '-Shu',
            '-'
        ], stdin=bwa_mem.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        bwa_mem.stdout.close()

        samtools_sort = subprocess.Popen([
            'samtools',
            'sort',
            '-',
            '-o',
            '{}'.format(sorted_sample_file.resolve())
        ], stdin=samtools_view.stdout, stderr=subprocess.DEVNULL)
        samtools_view.stdout.close()

        samtools_sort.communicate()
    return alignment_files


def index_align(alignment_files):
    print('Indexing alignment files...')
    procs = []
    for alignment_file in alignment_files:
        proc = subprocess.Popen([
            'samtools',
            'index',
            '{}'.format(alignment_file)
        ])
        procs.append(proc)
    for proc in procs:
        proc.communicate()


def no_coverage(project_dir, alignment_files):
    print('Finding regions of no coverage...')
    nocov_files = []
    procs = []

    for alignment_file in alignment_files:
        sample_name = alignment_file.name.rsplit('.', 2)[0]
        bed_file = project_dir / '{}_nocov.bed'.format(sample_name)
        nocov_files.append(bed_file)

        with open(bed_file, 'w') as out_file:
            bedtools_genomecov = subprocess.Popen([
                'bedtools',
                'genomecov',
                '-bga',
                '-ibam',
                '{}'.format(alignment_file)
            ], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

            awk = subprocess.Popen([
                'awk',
                '{ if ($4==0) print }'
            ], stdin=bedtools_genomecov.stdout, stdout=out_file, stderr=subprocess.DEVNULL)
            bedtools_genomecov.stdout.close()

            procs.append(awk)
    for proc in procs:
        proc.communicate()

    return nocov_files


def generate_mpileup(project_dir, sequences_file, alignment_files):
    print('Generating mpileup file...')
    mpileup_file = project_dir / 'data.mpileup'

    alignment_file_paths = [filepath.resolve() for filepath in alignment_files]

    with open(mpileup_file, 'w') as out_file:
        subprocess.call([
            'samtools',
            'mpileup',
            '-f',
            '{}'.format(sequences_file)
        ] + alignment_file_paths, stdout=out_file, stderr=subprocess.DEVNULL)

    return mpileup_file


def write_sample_list(project_dir, samples):
    sample_list_file = project_dir / 'samples_list.txt'
    with open(sample_list_file, 'w') as f:
        for sample in samples:
            f.write('{}\n'.format(sample))
    return sample_list_file


def varscan_cmd(min_coverage, min_var_freq, p_value, mpileup_file, sample_list_file):
    return [
        'varscan',
        'mpileup2cns',
        '{}'.format(mpileup_file),
        '--min-coverage {}'.format(min_coverage),
        '--min-var-freq {}'.format(min_var_freq),
        '--p-value {}'.format(p_value),
        '--output-vcf',
        '--variants 1',
        '--vcf-sample-list {}'.format(sample_list_file.resolve())
    ]


def variant_calling(project_dir, sample_list_file, mpileup_file):
    print('Performing variant calling...')
    spec_file = project_dir / 'spec_variants.vcf'
    sens_file = project_dir / 'sens_variants.vcf'

    with open(spec_file, 'w') as spec_out, open(sens_file, 'w') as sens_out:
        # TODO: Define these in a config file?
        spec_cmd = varscan_cmd(3, 0.1, 0.05, mpileup_file, sample_list_file)
        sens_cmd = varscan_cmd(10, 0.9, 0.05, mpileup_file, sample_list_file)
        spec_varscan = subprocess.Popen(spec_cmd, stdout=spec_out, stderr=subprocess.DEVNULL)
        sens_varscan = subprocess.Popen(sens_cmd, stdout=sens_out, stderr=subprocess.DEVNULL)
        spec_varscan.communicate()
        sens_varscan.communicate()

    return (spec_file, sens_file)


def snpeff_cmd(workspace_dir, reference, in_file):
    config_file = workspace_dir / 'snpEff.config'
    references_dir = workspace_dir / 'references'
    return [
        'snpEff',
        'eff',
        '-v',
        '-c {}'.format(config_file.resolve()),
        '-dataDir {}'.format(references_dir.resolve()),
        '-o vcf',
        '-no-downstream',
        '-no-upstream',
        '-classic',
        '-noStats',
        '{}'.format(reference),
        '{}'.format(in_file.resolve())
    ]


def snpeff(workspace_dir, project_dir, reference, spec_file, sens_file):
    print('Running SnpEff...')
    annotated_spec_file = project_dir / 'spec_variants_annotated.vcf'
    annotated_sens_file = project_dir / 'sens_variants_annotated.vcf'

    with open(annotated_spec_file, 'w') as spec_out, open(annotated_sens_file, 'w') as sens_out:
        spec_cmd = snpeff_cmd(workspace_dir, reference, spec_file)
        sens_cmd = snpeff_cmd(workspace_dir, reference, sens_file)
        spec_snpeff = subprocess.Popen(spec_cmd, stdout=spec_out, stderr=subprocess.DEVNULL)
        sens_snpeff = subprocess.Popen(sens_cmd, stdout=sens_out, stderr=subprocess.DEVNULL)
        spec_snpeff.communicate()
        sens_snpeff.communicate()

    return (annotated_spec_file, annotated_sens_file)


def create_tsv(project_dir, annotated_spec_file, annotated_sens_file):
    print('Converting VCF to TSV...')
    spec_txt_file = project_dir / 'spec_variants_annotated.txt'
    sens_txt_file = project_dir / 'sens_variants_annotated.txt'

    # TODO: Find a better way of doing this (rewrite vcf2tab?)
    vcf2tab_path = Path(__file__).parent.resolve() / 'vcf2tab.py'
    vcf2tab_cmd = ['python', str(vcf2tab_path)]

    with open(annotated_spec_file) as spec_in, open(spec_txt_file, 'w') as spec_out:
        spec_vcf2tab = subprocess.Popen(vcf2tab_cmd, stdin=spec_in, stdout=spec_out)
        spec_vcf2tab.communicate()

    with open(annotated_sens_file) as sens_in, open(sens_txt_file, 'w') as sens_out:
        sens_vcf2tab = subprocess.Popen(vcf2tab_cmd, stdin=sens_in, stdout=sens_out)
        sens_vcf2tab.communicate()

    return (spec_txt_file, sens_txt_file)


def create_json(project_dir, spec_txt_file, sens_txt_file):
    print('Converting TSV to JSON...')
    spec_json_file = project_dir / 'spec_variants_annotated.json'
    sens_json_file = project_dir / 'sens_variants_annotated.json'

    with open(spec_txt_file) as spec_in, open(spec_json_file, 'w') as spec_out:
        spec_reader = csv.DictReader(spec_in, delimiter='\t')
        json.dump([row for row in spec_reader], spec_out, indent=4)

    with open(sens_txt_file) as sens_in, open(sens_json_file, 'w') as sens_out:
        sens_reader = csv.DictReader(sens_in, delimiter='\t')
        json.dump([row for row in sens_reader], sens_out, indent=4)

    return (spec_json_file, sens_json_file)


def package_results(project_dir, sequences_file, genes_file, samples):
    print('Packaging results as zip...')
    results_zip = project_dir / '{}.zip'.format(RESULTS_ZIP_NAME)

    # Reference fasta and gff
    filepaths = [sequences_file, genes_file]
    # Sample bam, bam.bai and bed files
    filepaths += list(chain.from_iterable([(
        project_dir / '{}.sorted.bam'.format(s),
        project_dir / '{}.sorted.bam.bai'.format(s),
        project_dir / '{}_nocov.bed'.format(s)
    ) for s in samples]))
    # Variant VCF, TXT and JSON files
    filepaths += list(project_dir.glob('*variants_annotated*'))
    # README
    filepaths.append(Path(pkg_resources.resource_filename(__name__, 'data/README.txt')))

    with zipfile.ZipFile(results_zip, 'w', zipfile.ZIP_DEFLATED) as pz:
        for fp in filepaths:
            pz.write(fp, '{}/{}'.format(RESULTS_ZIP_NAME, fp.name))

    return results_zip


def upload_results(results_path, results_zip, spec_json_file, sens_json_file):
    print('Uploading results...')
    filepaths = [spec_json_file, sens_json_file, results_zip]
    for fp in filepaths:
        s3_path = '{}/{}'.format(results_path, fp.name)
        print('Uploading {} to S3 at {}'.format(fp.name, s3_path))
        S3_CONN.upload_file(fp, s3_path)


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
        print('Project directory {} already exists'.format(project_dir))

    # Setup file logging to project directory
    log_dir = project_dir / 'logs'
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / 'mngvariants.log'
    if log_file.is_file():
        log_file.unlink()
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter(LOG_FORMAT))
    logger.addHandler(file_handler)

    # Download the reads zip for this project from S3
    reads_zip_path = download_reads(project_dir, results_path)

    # Unzip required samples to reads directory
    reads_dir = unzip_samples(project_dir, reads_zip_path, args.samples)

    # Get reference genome and annotations from RefSeq
    (sequences_file, genes_file) = get_reference(args.workspace, args.reference)

    # Index the reference sequences file using bwa index
    index_sequences(sequences_file)

    # Align the sample reads to the reference using bwa mem and sort the bam file
    alignment_files = align(project_dir, reads_dir, sequences_file, args.samples)

    # Index alignment files
    index_align(alignment_files)

    # Compute regions of no coverage for each sample outputting bed files
    no_coverage(project_dir, alignment_files)

    # Generate mpileup
    mpileup_file = generate_mpileup(project_dir, sequences_file, alignment_files)

    # Write sample list file as required by VarScan's vcf-sample-list option
    sample_list_file = write_sample_list(project_dir, args.samples)

    # Perform variant calling
    (spec_file, sens_file) = variant_calling(project_dir, sample_list_file, mpileup_file)

    # Call SnpEff to annotate variants
    (annotated_spec_file, annotated_sens_file) = snpeff(args.workspace, project_dir, args.reference, spec_file, sens_file)

    # Create the TSV representation of the annotated VCFs
    (spec_txt_file, sens_txt_file) = create_tsv(project_dir, annotated_spec_file, annotated_sens_file)

    # Create JSON DataTable file
    (spec_json_file, sens_json_file) = create_json(project_dir, spec_txt_file, sens_txt_file)

    # Package up the results into a zip to be delivered to the customer
    results_zip = package_results(project_dir, sequences_file, genes_file, args.samples)

    # Upload the results zip and JSON DataTable files to S3
    upload_results(results_path, results_zip, spec_json_file, sens_json_file)

    print('\nDone')

import gzip
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
import requests
from dotenv import load_dotenv

from . import s3, util

# Load environment variables from .env file
load_dotenv(dotenv_path=Path(__file__).parent.resolve() / '.env')

# Base URL for LIMS
LIMS_BASE_URL = os.environ['LIMS_BASE_URL']
# URL and bucket for S3
S3_URL = os.environ['S3_URL']
S3_BUCKET_NAME = os.environ['S3_BUCKET']
# RefSeq Bacteria Assembly Summary URL
REFSEQ_ASSEMBLY_SUMMARY_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'

# Create S3Connection object to handle all S3 interactions
AWS_CONFIG_FILE_PATH = Path(__file__).parent.resolve() / '.aws' / 'config'
S3_CONN = s3.S3Connection(S3_URL, S3_BUCKET_NAME, AWS_CONFIG_FILE_PATH)

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

def download_file(url, local_path):
    """Download file from the URL to the local path."""
    print('Download: {} -> {}'.format(url, local_path))
    with urllib.request.urlopen(url) as response, open(local_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)

def add_reference_to_config(config_file, reference, refseq_url):
    """Add the snpEff.config entries for the reference."""
    print('Adding snpEff.config entry for reference {}...'.format(reference))
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
        new_lines += [
            '\t{}.{}.codonTable : Bacterial_and_Plant_Plastid'.format(reference, accession)
            for accession
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
        '-gff3',
        '-v',
        '{}'.format(reference)
    ], stderr=subprocess.DEVNULL)

def get_reference(workspace, reference):
    """Get the reference genome sequence and annotations from refseq if we don't already have them.
    """
    print('Get reference {}...'.format(reference))

    config_file = workspace / 'snpEff.config'
    # Check that the snpEff.config file exists
    if not config_file.is_file():
        raise Exception('Workspace directory does not contain snpEff.config')

    references_dir = workspace / 'references'
    reference_dir = references_dir / reference
    sequences_file = reference_dir / 'sequences.fa.gz'
    genes_file = reference_dir / 'genes.gff.gz'

    # Create the directory for this reference in workspace/snpEff/ if it doesn't already exist
    try:
        reference_dir.mkdir(parents=True)
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
            add_reference_to_config(config_file, reference, refseq_url)

        build_snpeff_database(config_file, references_dir, reference)

    return reference_dir

def extract_reference(reference_directory):
    """Extract the reference sequence and annotations and return the decompressed files."""
    print('Extracting reference files...')
    sequences_out = reference_directory / 'sequences.fa'
    genes_out = reference_directory / 'genes.gff'

    with gzip.open(reference_directory / 'sequences.fa.gz', 'rb') as seq_in:
        with open(sequences_out, 'wb') as seq_out:
            shutil.copyfileobj(seq_in, seq_out)

    with gzip.open(reference_directory / 'genes.gff.gz', 'rb') as gen_in:
        with open(genes_out, 'wb') as gen_out:
            shutil.copyfileobj(gen_in, gen_out)

    return sequences_out

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

        if sorted_sample_file.is_file():
            print('Sample {} already aligned, skipping'.format(sample))
        else:
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

def generate_mpileup(project_dir, sequences_file, alignment_files):
    print('Generating mpileup file...')
    mpileup_file = project_dir / 'data.mpileup'

    if mpileup_file.is_file():
        print('mpileup already exists, skipping')
    else:
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

    if spec_file.is_file() and sens_file.is_file():
        print('Variant VCFs already exist, skipping')
    else:
        with open(spec_file, 'w') as spec_out, open(sens_file, 'w') as sens_out:
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

    if annotated_spec_file.is_file() and annotated_sens_file.is_file():
        print('Annotated variant VCFs already exist, skipping')
    else:
        with open(annotated_spec_file, 'w') as spec_out, open(annotated_sens_file, 'w') as sens_out:
            spec_cmd = snpeff_cmd(workspace_dir, reference, spec_file)
            sens_cmd = snpeff_cmd(workspace_dir, reference, sens_file)
            spec_snpeff = subprocess.Popen(spec_cmd, stdout=spec_out, stderr=subprocess.DEVNULL)
            sens_snpeff = subprocess.Popen(sens_cmd, stdout=sens_out, stderr=subprocess.DEVNULL)
            spec_snpeff.communicate()
            sens_snpeff.communicate()

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

    # Download the reads zip for this project from S3
    reads_zip_path = download_reads(project_dir, results_path)

    # Unzip required samples to reads directory
    reads_dir = unzip_samples(project_dir, reads_zip_path, args.samples)

    # Get reference genome and annotations from RefSeq
    ref_dir = get_reference(args.workspace, args.reference)

    # Extract the reference files and reads
    sequences_file = extract_reference(ref_dir)

    # Index the reference sequences file using bwa index
    index_sequences(sequences_file)

    # Align the sample reads to the reference using bwa mem and sort the bam file
    alignment_files = align(project_dir, reads_dir, sequences_file, args.samples)

    # Index alignment files
    index_align(alignment_files)

    # Generate mpileup
    mpileup_file = generate_mpileup(project_dir, sequences_file, alignment_files)

    # Write sample list file as required by VarScan's vcf-sample-list option
    sample_list_file = write_sample_list(project_dir, args.samples)

    # Perform variant calling
    (spec_file, sens_file) = variant_calling(project_dir, sample_list_file, mpileup_file)

    # Call SnpEff to annotate variants
    snpeff(args.workspace, project_dir, args.reference, spec_file, sens_file)

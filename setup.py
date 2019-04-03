from setuptools import find_packages, setup

# Package metadata
NAME = 'mngvariants'
VERSION = '0.1.0'
DESCRIPTION = 'MicrobesNG variant calling.'
URL = 'https://github.com/SamChatfield/mngvariants'
EMAIL = 'samchatfield97@gmail.com'
AUTHOR = 'Sam Chatfield'
REQUIRES_PYTHON = '>=3.6.0'

# Required dependencies
REQUIRED = [
    'requests==2.20.1',
    'boto3==1.9.53',
    'pandas==0.23.4',
    'python-dotenv==0.10.1'
]

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(),
    package_data={'': ['.env', '.aws/*']},
    entry_points={
        'console_scripts': [
            'mngvariants=mngvariants.cli:main'
        ]
    },
    install_requires=REQUIRED
)

import subprocess
import setuptools
from setuptools import setup

subprocess.check_call('make')

setup(
    name='healcov',
    version='0.1',
    author='Pier Fiedorowicz',
    author_email='pierfied@email.arizona.edu',
    description='healcov',
    long_description='',
    packages=setuptools.find_packages(),
    package_data={'healcov': ['*']},
    zip_safe=False,
)
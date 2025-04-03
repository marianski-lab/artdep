from setuptools import setup, find_packages
import subprocess
import os

print(os.path.abspath("."))

version_check = subprocess.run(['./version.sh'], stdout=subprocess.PIPE)
version_num = version_check.stdout.decode('utf-8')

setup(
    name='matmacore',
    version=version_num,
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'colormaps ',
        'networkx'
    ],
)


from setuptools import setup, find_packages
import subprocess

version_check = subprocess.run(['version.sh'], stdout=subprocess.PIPE).stdout.decode('utf-8')
version_num = version_check.stdout

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


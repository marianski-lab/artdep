from setuptools import setup, find_packages
import subprocess

result = subprocess.run(['curl -s https://api.github.com/repos/marianski-lab/artdep/tags | jq -r first(.[].name | select(test("^v?[0-9]")))'], stdout=subprocess.PIPE)
version_num = result.stdout.decode('utf-8').strip()

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

from setuptools import setup, find_packages
import requests

response = requests.get("https://api.github.com/repos/marianski-lab/artdep/tags")
version_num = response.json()[0]['name']

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

from setuptools import setup, find_packages

setup(
    name='artdep',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'colormaps ',
        'networkx'
    ],
)


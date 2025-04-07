from setuptools import setup, find_packages
import dynamic_versioning

version_num = '0.1.10'

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
    cmdclass={
        "egg_info": dynamic_versioning.DynamicVersioningEggInfo
    }
)


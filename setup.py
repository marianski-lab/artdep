from setuptools import setup, find_packages
import dynamic_versioning

setup(
    name='matmacore',
    packages=find_packages(),
    version="0.1.11",
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


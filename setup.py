from setuptools import setup, find_packages

setup(
    name='matmacore',
    packages=find_packages(),
    version='0.1.14',
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


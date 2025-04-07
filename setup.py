from setuptools import setup, find_packages
import dynamic_versioning

version_num = '0.1.10'

setup(
    cmdclass={
        "egg_info": dynamic_versioning.DynamicVersioningEggInfo
    }
)


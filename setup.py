from setuptools import setup, find_packages
from codecs import open
from os import path


# Get the directory where this current file is saved
here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

req_path = path.join(here, 'requirements.txt')
with open(req_path, "r") as f:
    install_reqs = f.read().strip()
    install_reqs = install_reqs.split("\n")

setup(
    name='csa-atrophy',
    version='0.1.0',
    python_requires='>=3.6',
    description='Evaluate the sensitivity of atrophy detection with SCT.',
    author='NeuroPoly Lab, Polytechnique Montreal',
    author_email='neuropoly@googlegroups.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='',
    install_requires=install_reqs,
    entry_points={
        'console_scripts': [
            'affine_transfo=affine_transfo:main',
            'affine_rescale=affine_rescale:main',
            'yaml_parser=yaml_parser:main',
            'csa_rescale_stat=csa_rescale_stat:main',
            'manual_correction=manual_correction:main',
            'get_pmj=get_pmj:main'
        ],
    },
)

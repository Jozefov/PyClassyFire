# setup.py

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='PyClassyFire',
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='A Python client for the ClassyFire API for large-scale chemical compound classification.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Jozefov/PyClassyFire',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'requests>=2.32.3',
        'click>=8.1.7',
        'rdkit>=2024.3.5',
        'tqdm>=4.66.5',
    ],
    entry_points={
        'console_scripts': [
            'classyfire=pyclassyfire.cli:main',
        ],
    },
)

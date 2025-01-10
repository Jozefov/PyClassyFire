from setuptools import setup, find_packages
import pathlib

HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text(encoding="utf-8")

setup(
    name='PyClassyFire',
    version='0.1.0',
    author='Filip Jozefov',
    author_email='your.email@example.com',
    description='A Python client for the ClassyFire API for large-scale chemical compound classification.',
    long_description=README,
    long_description_content_type="text/markdown",
    url='https://github.com/Jozefov/PyClassyFire',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    python_requires='>=3.6',
    install_requires=[
        'requests>=2.32.3',
        'click>=8.1.7',
        'tqdm>=4.66.5',
        'rdkit>=2024.3.5',
        'pandas>=2.2.3',
    ],
    entry_points={
        'console_scripts': [
            'pyclassyfire=pyclassyfire.cli:cli',
        ],
    },
    project_urls={
        'Source': 'https://github.com/Jozefov/PyClassyFire',
    },
)
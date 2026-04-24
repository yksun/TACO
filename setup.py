"""Setup script for TACO - Telomere-Aware Contig Optimization."""
from setuptools import setup, find_packages

setup(
    name="taco-genome",
    version="1.3.0",
    description="Telomere-Aware Contig Optimization pipeline",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Yukun Sun",
    author_email="ysun@fieldmuseum.org",
    url="https://github.com/yksun/TACO",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "taco=taco.__main__:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

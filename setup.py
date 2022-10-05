import setuptools
import os

with open("README.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

# requirements = []
# with open("requirements.txt", "r") as requirements_file:
#     for req in (line.strip() for line in requirements_file):
#         requirements.append(req)

setuptools.setup(
    name="ukbb_common",
    version="0.1.3",
    author="",
    author_email="",
    description="Common functions for UK Biobank Data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Nealelab/ukbb_common",
    project_urls={
        "Source Code": "https://github.com/Nealelab/ukbb_common",
        "Issues": "https://github.com/Nealelab/ukbb_common/issues",
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.6",
    package_dir={"ukbb_common": "."},
    packages=['ukbb_common', 'ukbb_common.resources', 'ukbb_common.utils'],
    install_requires=['hail', 'batch'],
)

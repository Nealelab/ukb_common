import setuptools


with open("README.md", "r", encoding="utf-8") as readme_file:
    long_description = readme_file.read()

install_requires = []
with open("requirements.txt", "r") as requirements_file:
    for req in (line.strip() for line in requirements_file):
        if req != "hail":
            install_requires.append(req)

setuptools.setup(
    name="ukb_common",
    version="0.1.0",
    author="",
    author_email="",
    description="Common functions for UK Biobank Data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Nealelab/ukb_common",
    packages=setuptools.find_namespace_packages(include=["ukb_common.*"]),
    project_urls={
        "Documentation": "https://broadinstitute.github.io/ukb_common/",
        "Source Code": "https://github.com/broadinstitute/ukb_common",
        "Issues": "https://github.com/broadinstitute/ukb_common/issues",
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.6",
    install_requires=install_requires,
)

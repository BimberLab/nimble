import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="immunogenotyper",
    version="0.1.0",
    author="Sebastian Benjamin",
    author_email="benjamse@ohsu.edu",
    description="A tool for generating calls on bulk-seq and 10x RNA genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/devsebb/ImmunoGenotyperLibraryifier",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

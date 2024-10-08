import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="nimble",
    version="0.9.0",
    author="Sebastian Benjamin",
    author_email="benjamse@ohsu.edu",
    description="A tool for generating calls on bulk-seq and 10x RNA genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BimberLab/nimble",
    packages=setuptools.find_packages(),
    package_data={'nimble': ['templates/*.html']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["requests", "biopython", "pandas", "distro", "pysam", "seaborn", "jinja2", "matplotlib"],
    python_requires=">=3.9",
)

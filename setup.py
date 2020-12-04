import setuptools

#Get long description from readme file
with open("README.md", 'r') as readme_file:
    long_description = readme_file.read()

setuptools.setup(
    name="inequipy",
    version="1.0.8",
    author="Tom Logan, Mitchell Anderson",
    author_email="tom.logan@canterbury.ac.nz",
    description="Inequality Metrics contains functions for the Kolm-Pollak, Atkinson and Gini Index aproaches to calculating inequality of a distribution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/urutau-nz/inequipy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'scipy',
        'numpy'
    ],
    python_requires='>=3',
)

#There are other things which can be modified

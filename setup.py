import setuptools

#Get long description from readme file
with open("README.md", 'r') as readme_file:
    long_description = readme_file.read()

setuptools.setup(
    name="InequalityMetrics",
    version="1.0.0",
    author="Tom Logan",
    author_email="micae@example.com",
    description="Inequality Metrics contains KMG functions for calculating inquality distributions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://kcna.kp/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

#There are other things which can be modified
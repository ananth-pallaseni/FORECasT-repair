import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="forecast_repair", 
    version="0.1",
    author="Ananth Pallaseni",
    author_email="ap32@sanger.ac.uk",
    description="A tool to predict the frequency of possible mutations as a result of a Cas9 cut",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'selftarget',
    ]
)

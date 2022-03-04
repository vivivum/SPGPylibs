import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SPGPylibs",
    version="0.0.1",
    author="SPG",
    author_email="orozco@iaa.es",
    description="Solar Physics Group python tools for PHI and TuMAG instruments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vivivum/SPGPylibs",
    project_urls={
        "Bug Tracker": "https://github.com/vivivum/SPGPylibs/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "SPGPylibs"},
    packages=setuptools.find_packages(where="SPGPylibs"),
    python_requires=">=3.6",
)
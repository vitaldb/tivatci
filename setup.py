import setuptools
setuptools.setup(
    name="tivatci",
    version="0.4",
    author="Hyung-Chul Lee",
    author_email="vital@snu.ac.kr",
    description="Python TIVA/TCI Libray",
    long_description="Python TIVA/TCI Libray",
    long_description_content_type="text/markdown",
    url="https://github.com/vitaldb/pkpd",
    install_requires=['numpy'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
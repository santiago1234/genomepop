from setuptools import setup, find_packages

setup(
    name='genomepop',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'scikit-allel'
    ],
    author='Santiago G. Medina Muñoz',
    author_email='santiago.medina@cinvestav.mx',
    description='A package for analyzing whole genome data'
)

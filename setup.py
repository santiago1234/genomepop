from setuptools import setup, find_packages

setup(
    name='genomepop',
    version='0.0',
    packages=find_packages(),
    install_requires=[
        'click',
        'pandas',
        'numpy',
        'scikit-allel'
    ],
    entry_points='''
    [console_scripts]
    hello=genomepop.cltools:hello
    gnomix2tracts=genomepop.cltools:gnomix2tracts
    ''',
    author='Santiago G. Medina Muñoz',
    author_email='santiago.medina@cinvestav.mx',
    url='https://github.com/santiago1234/genomepop',
    description='A package for analyzing whole genome data'
)

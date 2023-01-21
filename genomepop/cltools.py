"""
cltools
========
This submodule contains command-line tools that can be used to interact with the module.

Command line tools available:
    - tool1: This command line tool does XYZ.
        usage: tool1 arg1 arg2
    - tool2: This command line tool does ABC.
        usage: tool2 arg1 arg2 --flag1 --flag2
    - tool3: This command line tool does DEF.
        usage: tool3 arg1 --flag1

"""
import click

from . import localancestry


@click.command()
@click.option('--msp',
              '-m',
              help='Path to Genomix msp file (predicted local ancestry)',
              required=True)
def gnomix2tracts(msp):
    """
    Convert Gnomix output (msp file)
    to a bed file that can be used to run Tracts
    analysis.
    """
    click.echo('converting to tracts ...')


# Another cli tool
@click.command()
def hello():
    print('Hi!')

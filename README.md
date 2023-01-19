# genomepop

genomepop is a package for performing population genetics analysis on whole genome sequences.

## Installation

To install `genomepop`, you can use `pip`:

```bash
pip install
```

`genomepop` requires Python 3.6 or later and the following dependencies:
- `numpy`
- `scipy`
- `pandas`
- `biopython`

## Usage
Here is an example of how to use `genomepop` to perform a basic population genetics analysis:


```from genomepop import PopGen
genome_file = "example.fasta"
popgen = PopGen(genome_file)
popgen.calculate_fst()
popgen.plot_pca()
```

## Documentation
For more information on how to use `genomepop`, please refer to the [documentation](https://genomepop.readthedocs.io/).

## Citation
If you use `genomepop` in your research, please cite:

```@article{popgenanalyzer_citation,
title={genomepop: A package for population genetics analysis},
author={Author, A.},
journal={BMC Bioinformatics},
volume={20},
pages={1--10},
year={2019}
}
```

## Contribution
If you would like to contribute to `genomepop`, please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for information on how to do so.

# Contributing to *thermofeel*

Please report [bug reports](https://github.com/ecmwf-projects/thermofeel/issues)
or [pull requests](https://github.com/ecmwf-projects/thermofeel/issues) on
[GitHub](https://github.com/ecmwf-projects/thermofeel).

We want your feedback, please e-mail: <user-services@ecmwf.int>

The package is installed from PyPI with:

```bash
pip install thermofeel
```

How you could use thermofeel:

- Test out the code in this library with data from different models and in
  different structures
- Compare the methods with other outputs
- Plot outputs, examples are provided in the examples directory
- Go further than our examples and let us know how it goes

## Development setup

This project is pure Python (runtime dependency: `numpy`) and uses a top-level
`Makefile` as the single entry point. Run `make help` to see every target.

```bash
make venv     # create a local .venv (uv) with the package + test deps
make all      # the full gate: check + test + lint
make fmt      # apply isort + black formatting
make docs     # build the MkDocs documentation site
```

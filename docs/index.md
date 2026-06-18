# Welcome to thermofeel's documentation

**thermofeel** is a Python library developed by the European Centre for
Medium-Range Weather Forecasts (ECMWF) to integrate thermal (heat and cold)
indices into operational systems and to maintain them over time.

It allows users to:

- integrate these methods into their Numerical Weather Prediction (NWP)
  services, and
- calculate thermal indices from their own data.

This supports transparency and reproducibility, and lets users expand the scope
of these methods.

## Contents

- [Overview](guide/overview.md) — what thermofeel computes
- **Indices** — one guide page per thermal index
- **Supporting quantities** — relative humidity, vapour pressure, mean radiant
  temperature
- [API reference](api.md) — generated from the source docstrings
- [Examples](examples.md) — runnable usage

## Install

```bash
pip install thermofeel
```

The core depends only on `numpy`.

## License

*thermofeel* is available under the open source
[Apache License](http://www.apache.org/licenses/LICENSE-2.0.html). In applying
this licence, ECMWF does not waive the privileges and immunities granted to it
by virtue of its status as an intergovernmental organisation nor does it submit
to any jurisdiction.

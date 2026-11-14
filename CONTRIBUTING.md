# Contributing to `ggtree-py`

Thank you for considering contributing to `ggtree-py`!

## Development Setup

```bash
git clone <repo-url>
cd ggtree.py
python3 -m pip install -e ".[dev]"
```

Only `numpy>=1.22` is required at runtime (Python `>= 3.10`). Installing `rsvg-convert` (from `librsvg`) on your `PATH` enables high-resolution vector-to-PNG rasterization in `ggsave("out.png")`; without it, `ggtree` falls back to its built-in pure-NumPy PNG rasterizer.

## Running Tests and Examples

Run the full unit test suite (all 49 tests covering ported R manual examples and `testthat` suites):

```bash
python3 -s -m unittest discover -s tests -v
```

Generate all 14 composite example figures (SVG + PNG) in `examples/output/`:

```bash
python3 -s examples/run_all_examples.py
```

## Building Distributions

```bash
uv build
```

## Submitting Pull Requests

- Keep core runtime dependencies restricted to the Python standard library and `numpy`.
- Add or update tests in `tests/test_ggtree.py` or `tests/test_examples.py`.
- Update the relevant documentation under `docs/` when adding or modifying public APIs.
- Follow the [Contributor Code of Conduct](CONDUCT.md).

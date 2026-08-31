# Docs starter for pyCOMET

## Install

```bash
pip install "py-comet[docs]"
```

## Preview

```bash
mkdocs serve
```

## Publish (GitHub Pages)

```bash
mkdocs gh-deploy
```

> `mkdocstrings` imports `comet` to generate the API reference, so COMET must be
> importable in the same environment. Installing the `docs` extra as above
> covers that; from a source checkout use `pip install -e ".[docs]"`.

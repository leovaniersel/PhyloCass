# Releasing to PyPI

The package is ready to publish: metadata is complete, both artifacts build and
pass `twine check`, and each has been installed into a clean environment and
exercised. What is left needs a PyPI account, which is why it is not automated
away here.

## Two things that used to block this, now settled

**Source visibility.** A PyPI release publishes the sdist, which contains the
whole source tree; there is no private mode. `leovaniersel/PhyloCass` is public,
so this is moot.

**Repository links.** `pyproject.toml` points Homepage, Repository and Issues at
`github.com/leovaniersel/PhyloCass`. Those resolve, now that the repository is
public.

The one prerequisite left is registering a trusted publisher on PyPI — see
below. It has to exist *before* the release fires, because the workflow
authenticates against it.

## Test it on TestPyPI first

TestPyPI is a separate index with separate accounts, and it is the cheap way to
see the rendered page before it is permanent. Note that a version number on
PyPI can never be reused, even after deleting a release.

```bash
python -m build
python -m twine check dist/*
python -m twine upload --repository testpypi dist/*

# then, in a clean environment (phylozoo comes from real PyPI):
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ phylocass
```

## Publish

Either upload directly:

```bash
python -m build
python -m twine check dist/*
python -m twine upload dist/*
```

with a PyPI API token as the password and `__token__` as the username — or use
the workflow below, which needs no token on your machine at all.

## Trusted publishing via GitHub Actions (recommended)

[`.github/workflows/publish.yml`](.github/workflows/publish.yml) builds, tests
and uploads whenever you publish a GitHub release. It uses PyPI's trusted
publishing, so there is no API token to create, store or rotate.

One-time setup, on <https://pypi.org/manage/account/publishing/>. Because
`phylocass` does not exist on PyPI yet, this goes under **"Add a new pending
publisher"** — a project-scoped publisher can only be added to a project that
already exists:

| field | value |
| --- | --- |
| PyPI project name | `phylocass` |
| Owner | `leovaniersel` |
| Repository name | `PhyloCass` |
| Workflow name | `publish.yml` |
| Environment name | `pypi` |

The `pypi` GitHub environment already exists. Then release:

```bash
git tag v0.3.0 && git push origin v0.3.0     # already done for v0.3.0
gh release create v0.3.0 --generate-notes    # or the web UI, or the REST API
```

`gh` is not installed on the author's machine; the web UI at
`github.com/leovaniersel/PhyloCass/releases/new?tag=v0.3.0` and the REST API
(`POST /repos/{owner}/{repo}/releases`) both work just as well.

## Version numbers

The version lives in two places and the tests check they agree
(`test_package_metadata_matches_dunder_version`):

- `pyproject.toml` → `version`
- `src/phylocass/__init__.py` → `__version__`

Bump both, run `pytest`, then tag.

## Checklist

- [ ] `pytest` passes
- [ ] version bumped in both places
- [ ] `python -m build && python -m twine check dist/*`
- [ ] installed the wheel *and* the sdist into clean environments and run
      `phylocass --version` plus one real input
- [ ] decided about repository visibility and the `[project.urls]` links

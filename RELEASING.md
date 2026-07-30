# Releasing to PyPI

The package is ready to publish: metadata is complete, both artifacts build and
pass `twine check`, and each has been installed into a clean environment and
exercised. What is left needs a PyPI account, which is why it is not automated
away here.

## Before the first release, decide two things

**1. Publishing makes the source public.** `leovaniersel/PhyloCass` is a private
repository, but a PyPI release publishes the sdist — which contains the whole
source tree — to anyone who asks for it. There is no private mode. If the code
should stay private, do not publish; install from the repository instead:

```bash
pip install git+ssh://git@github.com/leovaniersel/PhyloCass.git
```

**2. The metadata links to the repository.** `pyproject.toml` points Homepage,
Repository and Issues at `github.com/leovaniersel/PhyloCass`. While the repo is
private those links 404 for everyone reading the PyPI page. Either make the
repository public alongside the release, or drop the `[project.urls]` table.

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

One-time setup, on <https://pypi.org/manage/account/publishing/>:

| field | value |
| --- | --- |
| PyPI project name | `phylocass` |
| Owner | `leovaniersel` |
| Repository name | `PhyloCass` |
| Workflow name | `publish.yml` |
| Environment name | `pypi` |

Then create the `pypi` environment in the repository settings, and release:

```bash
git tag v0.3.0
git push origin v0.3.0
gh release create v0.3.0 --generate-notes    # or use the web UI
```

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

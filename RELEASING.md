# Releasing to PyPI

<https://pypi.org/project/phylocass/> — published since v0.3.0.

Releases are automated. Publishing a GitHub release runs the tests on Python
3.10–3.13, builds the sdist and wheel, checks them, and uploads to PyPI via
trusted publishing. No API token exists on anyone's machine.

## Cutting a release

```bash
# 1. bump the version in both places (see below), then
pytest

# 2. tag
git tag v0.4.0
git push origin v0.4.0

# 3. publish a GitHub release for that tag, which fires the workflow
gh release create v0.4.0 --generate-notes
```

`gh` is not installed on the author's machine. The web UI
(`github.com/leovaniersel/PhyloCass/releases/new?tag=v0.4.0`) and the REST API
(`POST /repos/leovaniersel/PhyloCass/releases` with `{"tag_name": "v0.4.0"}`)
do the same job.

Then watch <https://github.com/leovaniersel/PhyloCass/actions>. A failing
`test` job skips `build`, and `publish` never runs — so a broken release
cannot reach PyPI.

## Trusted publishing (already configured)

[`.github/workflows/publish.yml`](.github/workflows/publish.yml) authenticates
to PyPI with a short-lived OIDC token instead of a stored secret. The publisher
registered on <https://pypi.org/manage/account/publishing/> is:

| field | value |
| --- | --- |
| PyPI project | `phylocass` |
| Owner | `leovaniersel` |
| Repository | `PhyloCass` |
| Workflow | `publish.yml` |
| Environment | `pypi` |

The matching `pypi` environment exists in the repository settings. A brand-new
project needs a *pending* publisher; once it has been published to, PyPI
converts that into a project-scoped one automatically.

## Version numbers

The version lives in two places, and a test checks they agree
(`test_package_metadata_matches_dunder_version`):

- `pyproject.toml` → `version`
- `src/phylocass/__init__.py` → `__version__`

A version number on PyPI can never be reused, even after deleting a release, so
bump rather than retry.

## Rehearsing on TestPyPI

Rarely needed now the pipeline is proven, but useful before a release that
changes packaging. TestPyPI is a separate index with separate accounts:

```bash
python -m build
python -m twine check dist/*
python -m twine upload --repository testpypi dist/*

# phylozoo still comes from the real index
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ phylocass
```

## Checklist

- [ ] `pytest` passes
- [ ] version bumped in both places
- [ ] release notes say what changed
- [ ] after the run: `pip install phylocass` in a clean environment, then
      `phylocass --version` and one real input

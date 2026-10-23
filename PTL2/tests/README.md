# PTL2 unit tests

These tests check that the PTL2 scripts make the right decisions on small, hand-built examples. Right now they cover [`Scripts/contamination.py`](../Scripts/contamination.py).

- They do **not** run MAFFT, IQ-TREE, FastTree, trimAl or Guidance, so you don't need any of those installed.
- They use tiny trees (6–8 tips) and write only to temporary folders. Your data and output folders are never touched.
- The whole suite runs in a few seconds on a laptop.

## Running the tests

All commands are run from the top-level `EukPhylo/` folder.

**1. Create the environment (once).** This installs Python 3.11, Biopython, ete3, tqdm and pytest into a separate conda environment called `eukphylo-dev`:

```bash
conda env create -f dependencies.yml
```

**2. Activate it (every new terminal):**

```bash
conda activate eukphylo-dev
```

**3. Run the tests:**

```bash
pytest PTL2/tests
```

## Reading the output

A normal run ends like this:

```
36 passed, 1 xfailed in 3.10s
```

| Result | Meaning |
|---|---|
| `passed` | The script did what the test expected. |
| `xfailed` | An **expected** failure. The test describes a known bug that hasn't been fixed yet, so it is supposed to fail. This is not a problem. |
| `failed` | The script's behavior doesn't match the test. pytest prints which test failed, the line, and the expected vs. actual values. |
| `XPASS(strict)` on an `xfail` test | A known bug has been fixed. Remove the `@pytest.mark.xfail` line above that test. |
| `error` | The test couldn't run at all, usually because of the environment (see Troubleshooting). |

To see the name of every test as it runs, add `-v`:

```bash
pytest PTL2/tests -v
```

## Running only some tests

```bash
pytest PTL2/tests -k get_sisters                                                # tests with "get_sisters" in the name
pytest PTL2/tests/test_contamination.py::test_get_sisters_subsister_rule        # one test
pytest PTL2/tests -x                                                            # stop at the first failure
```

To see which lines of `contamination.py` the tests don't reach:

```bash
pytest PTL2/tests --cov=contamination --cov-report=term-missing
```

## What is tested

| Function | Scenarios |
|---|---|
| `get_newick` | Newick and Nexus files, quotes/backslashes in names, a file with no tree |
| `reroot` | rooting on prokaryotes, falling back to Opisthokonts, no qualifying clade |
| `get_subtrees` (clade-based removal) | a target sequence outside a clean clade is removed, non-targets are kept, `num_contams` as a count and a proportion, `min_target_presence`, `required_taxa`, target taxa from a file, the exceptions file, the rules file, missing input files |
| `get_sisters` (sister-based removal) | removal when sister to a contaminant, long branches kept, mixed sisters kept, co-contaminants, the sub-sister rule, a missing co-contaminants file, single-taxon trees (known bug, `xfail`) |
| `write_new_preguidance` | only kept sequences are written, excluded taxa at 2/5/10-character codes, zero or multiple matching sequence files, a missing exclude file |
| `cl_fasttree` | builds the right FastTree command, skips non-FASTA files (FastTree itself is mocked) |

Not yet covered: `run()`, `cl_mafft`, `cl_iqtree` and `cl_iqtree_fast`.

## Writing a new test

1. Add a function whose name starts with `test_` to `test_contamination.py`, or to a new file named `test_<script>.py` for another script. pytest finds them automatically.
2. Use the helpers at the top of `test_contamination.py`:
   - `make_args(...)` builds the parameter object the scripts expect, with sensible defaults. Override only what your test is about, e.g. `make_args(num_contams=1)`.
   - `write_file(tmp_path, name, content)` writes a tree or rules file into a temporary folder and returns its path.
3. Name tips like real EukPhylo sequences: the first 10 characters are the taxon code (e.g. `Sr_ci_AAAA_1`), and the first 2 are the major clade (`Ba`, `Za`, `Op`, `Pl`, `Am`, `Ex`, `Sr`).
4. Keep trees small enough that you can work out the right answer by hand, and say what the tree looks like in a comment.
5. If you find a bug and aren't fixing it in the same change, write the test for the correct behavior and mark it as an expected failure:

   ```python
   @pytest.mark.xfail(strict=True, reason='Known bug: <what goes wrong>')
   ```

## Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: No module named 'Bio'` (or `ete3`, `pytest`) | The environment isn't active. Run `conda activate eukphylo-dev`. |
| `ERROR: file or directory not found: PTL2/tests` | You're not in the top-level `EukPhylo/` folder. `cd` there first. |
| ete3 fails to install or import | ete3 doesn't support the newest Python versions. Use the Python version pinned in `dependencies.yml` (3.11). |
| VS Code underlines imports with "could not be resolved" | Select the environment as the interpreter: Cmd/Ctrl+Shift+P → "Python: Select Interpreter" → `eukphylo-dev`. |

import os, sys

# The pipeline scripts import each other as top-level modules (e.g. `import guidance`),
# so Scripts/ has to be on the path for the tests to import them.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Scripts'))

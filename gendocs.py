import os
os.system("sphinx-apidoc -f -o docs tivatci/__init__.py")
os.system("sphinx-build docs docs/_build")

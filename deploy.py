import os
import shutil

PROJECT_NAME = 'tivatci'

# generate docs
os.system(f'sphinx-apidoc -f -o docs {PROJECT_NAME}/__init__.py')
os.system("sphinx-build docs docs/_build")

# generate wheel
f = open("setup.py", "rt")
ver = ''
for line in f.readlines():
    if line.find('version=') > 0:
        tabs = line.strip().split('=')
        if len(tabs) >= 2:
            ver = tabs[1].strip('",')
if not ver:
    print('version not found in setup.py')
    quit()
os.system('python setup.py bdist_wheel')

# upload to pypi
os.system('twine upload dist\\' + PROJECT_NAME + '-' + ver + '-py3-none-any.whl')

# remove cache dirs
print('removing cache dirs', end='...', flush=True)
shutil.rmtree('build')
shutil.rmtree(PROJECT_NAME + '.egg-info')
shutil.rmtree('dist')
print('done')
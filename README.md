# Seismic data processing

I would recommend installing this in a `conda` environment.

This can be installed by:

1. clone repo
2. `cd` to top-level directory (i.e. with setup.py in it)
3. run  `pip install .` (or `pip install -e .`)
4. run `pip install -r requirements.txt`

Hopefully then the notebooks in `notebooks/` should work. There may still be an issue with
certain matplotlib plotting libraries (one of which, if I recall, requires installing in a virtualenv or conda environment)

# Creating a conda environment

## Links
* Install conda: [link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
* Create environment: [link](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

## Example

* Once you've installed `conda` you should be able to run:

```bash
# create env
conda create --name myenv
# install pip because
# not everything is conda installable
conda install pip
# activate env
conda activate myenv
# install requirements
pip install -r requirements.txt
# install seispy
pip install -e .
```


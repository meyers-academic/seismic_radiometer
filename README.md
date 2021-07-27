# A brief note

Coming from the gravitational-wave community, much of this code is written using code developed for the purpose of detecting and characterizing gravitational-waves. What that means is that, while many of the signal processing methods are similar, the packages used are non-standard in the seismic community. That means that, for now, if you'd like to use this code, please contact me (pmeyers279@gmail.com) and I will be happy to help.

I am hoping to include an interface that allows one to pull data using obspy and run this code. Presently, though, this code is written mainly to read "Gravitational-Wave Frames" (gwf), which is one data format in which we stored our data from the Homestake seismometer array.

A link to the related paper is here:

Also, I recognize that many of the names of classes, etc. are probably confusing. I'm working towards clearing up ambiguities.


# Seismic radiometer code

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
# activate env
conda activate myenv
# install pip because
# not everything is conda installable
conda install pip
# install requirements
pip install -r requirements.txt
# install seispy
pip install -e .
```

## Running jupyter notebooks

There are quite a few example Jupyter notebooks. You can run these through your browser by running (from this directory):

```bash
conda activate myenv
cd notebooks/
jupyter notebook
```

It should then bring up a browser window with a file tree with all of the available notebooks.

## General notes:

* Currently, we have included the locations for the homestake seismometer array in `dict` format in `seispy/station/station.py`. The specifications are units of meters, where the first value is EW, second value is NS, and third value is either depth (if coords=='depth', which is standard), or elevation (if coords=='utc'). Units are in meters. Accessing these are fairly straightforward, and examples are shown in the notebooks.

## FAQs:

# Seismic data processing

## Installation

### OSX [assuming you have root access]

**We're installing into a virtualenv here...**

Non-python dependencies are best done with macports

`sudo port install lal lalframe py27-lal py27-lalframe`

It's recommended that you create a virtualenv to install this you can install this with. I keep my virtualenvs in an "opt" directory in my home directory:

```
sudo port select --set virtualenv virtualenv-2.7
virtualenv ~/opt/seispy/
```

You'll need to add the lal swig python bindings to your python path for your virtualenv. Find where lal and lalframe are installed with:

`port contents py27-lal` and `port contents py27-lalframe`

For me it's here:

`/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lalframe`

and

`/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lal`

Now add this to the end of the activation script for your virtualenv:

`echo 'PYTHONPATH=$PYTHONPATH:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lalframe:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lal' >> ~/opt/seispy/bin/activate`

Now you can activate your virtualenv and it should be able to find lal properly.


### Linux

I'm working on it...

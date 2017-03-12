# Seismic data processing

## Installation

### OSX [assuming you have root access]

**We're installing into a virtualenv here...**

Non-python dependencies are best done with macports

`sudo port install lal lalframe py27-lal py27-lalframe`

It's recommended that you create a virtualenv to install in. I keep my virtualenvs in an "opt" directory in my home directory:

```
sudo port select --set virtualenv virtualenv-2.7
virtualenv ~/opt/seispy/
```


You'll need to add the python lal and lalframe bindings explicitly to the PYTHONPATH when you activate your virtualenv. This can be done by looking for where they're installed with these commands:

`port contents py27-lal` and `port contents py27-lalframe`

For me it's here:

`/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lalframe`

and

`/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lal`

You can then add a command to append them to your `PYTHONPATH` at the end of the activation script for your virtualenv:

`echo 'PYTHONPATH=$PYTHONPATH:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lalframe:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/lal' >> ~/opt/seispy/bin/activate`

Now you can source your virtualenv and hopefully things will work. You can test whether your installation is working by moving to the top level directory and running

```
python setup.py test
```

### Linux

I'm working on it...

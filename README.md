# ComputationalPhysics

This repository contains examples and assessments for the Computational Physics
(CP) module of the physics honours course 2021 (PHY4000W) at the University of
Cape Town.

## Usage

You can run the examples on your own laptop, on cloud service such as
https://mybinder.org/ or https://cocalc.com, or on the HEP cluster of the
physics department.

### On your own laptop

If you do not want to risk messing up your main Python istallation, you can set
up a [virtual environment]:
```
python3 -m venv venv
```

You have to activate the environment every time you want to use it:
```
. venv/bin/activate
```
If the command was successful, your shell prompt will (probably) change and
display the active virtual environment.

You can install the dependencies required to run the examples in this
repository with the provided requirements.txt file:
```
pip install -r requirements.txt
```
This step can be repeated when some dependencies are missing. If this does not
fix the problem, please report this as a bug.

With this setup, you can start the Jupyter Notebook:
```
jupyter notebook
```

Once you have done the setup, it will be sufficient to activate the environment
```
. venv/bin/activate
jupyter notebook
```

[virtual environment]: https://docs.python.org/3/library/venv.html

### On the HEP Cluster

The HEP groups in the physics department run a computer cluster where you can
run Jupyter Notebooks by logging into a JupyterHub with your UCT credentials:

https://hep02.phy.uct.ac.za/jupyter

At the moment, you will have to manually upload/download the examples to the
hub. This is a new service, so please give it a try and shout if something does
not work.

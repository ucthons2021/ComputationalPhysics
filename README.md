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

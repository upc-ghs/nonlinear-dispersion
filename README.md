## nonlinear-dispersion
```python``` routines to perform simple simulations of nonlinear dispersion with the [MODFLOW-API](https://www.usgs.gov/publications/modflow-application-programming-interface-simulationcontrol-and-software). 

### Setup
Install project requirements with the instruction 
 
```
pip install -r requirements.txt
```

The file ```src/config.py``` contains some global configurations pointing to the folder where the libraries of the MODFLOW-API is located and the output folder for simulations.

### Run 
The folder ```src``` stores self-contained routines for configuring and running models. All routines write the simulation data to the folder ```sims``` (defined in ```config.py```), which is not tracked by the versioning system. Routines include a set of simple arguments for controlling the execution. For example, the following example will write and run the routine ```flopymf6fujita.py```:

```
python flopymf6fujita.py --write --run --force
```

## Resources
- [MODFLOW-API](https://www.usgs.gov/publications/modflow-application-programming-interface-simulationcontrol-and-software) 

- [xmipy](https://github.com/Deltares/xmipy) 

## nonlinear-dispersion
```python``` routines to perform simple simulations of nonlinear dispersion with the [MODFLOW-API](https://www.usgs.gov/publications/modflow-application-programming-interface-simulationcontrol-and-software). 

### Setup
Install project requirements with the instruction 
 
```
pip install -r requirements.txt
```

The file ```src/config.py``` contains some global configurations indicating the folder where the libraries of the MODFLOW-API are located and the output folder for simulations.

**Note**: the Linux and Windows libraries are also included locally in the ```lib``` folder. Remember to select the appropiate for your operating system before running models. Libraries can also be downloaded [here](https://github.com/MODFLOW-USGS/executables)

### Run 
The folder ```src``` stores self-contained routines for configuring and running models. All routines write the simulation data to the folder ```sims``` (defined in ```config.py```), which is not tracked by the versioning system. 

Routines include a set of simple arguments for controlling the execution. For example, the following example will write and run the routine ```flopymf6fujita.py```:

```
python flopymf6fujita.py --write --run --force
```

## License
MIT License

## References
[Pérez-Illanes R., Saaltink M. W., Fernàndez-Garcia, D., 2024, Nonlinear Formulation of Multicomponent Reactive Transport With Species-Specific Dispersion Properties, Water Resources Research, doi:10.1029/2023WR036358](https://doi.org/10.1029/2023WR036358)


## Resources
- [MODFLOW-API](https://www.usgs.gov/publications/modflow-application-programming-interface-simulationcontrol-and-software) 

- [xmipy](https://github.com/Deltares/xmipy) 

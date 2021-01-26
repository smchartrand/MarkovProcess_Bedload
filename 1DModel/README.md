# Reduced Complexity Model of Bedload Transport for Uniform Particle Diameters

## Requirements

This model requires:

- Python3
- Pip3

Installing Python requirements:

```bash
 pip3 install -r requirements.txt
```

## Important Terminology in Model

- **Entrainment**:
- **Entrainment Event**:
- **Hop**: The distance a particle moves during entrainment
- **Desired Hop**: The distance a 
- **Bed Particle**: A particle embedded in the stream bed. Does not move.
- **Model Particle**: A particle subject to entrainment events. Moves.
- **Event Particle**: A model particle which has been selected for entrainment.
- **Stream**: Consists of Bed and Model particles.

## Running the Model

### Running from command line

1. Set current directory to **`RCM_BedloadTransport/model`**:

```bash
 cd /path/to/RCM_BedloadTransport/model
```

2. Run the model:
    - Run the model _with_ stream plots:

    ```bash
    python3 run.py 
    ```

   - Run the model _without_ stream plots:

    ```bash
    python3 run.py --no-plots
    ```

### Running in Spyder

1. Open **`run.py`** and **`parameters.py`** in Spyder 

2. Navigate to **`run.py`** and execute using the `'Run file (F5)'` button, or through the Console pane:

```bash
runfile('run.py')
```

_Note_: to use the above command, your console must be located in the model directory

If you wish to view the plots in the plot pane, ensure the call to `plot_stream` at the end of **`run.py`** uses the parameter ```to_file=False```. The function call should look like this:

```python
    ml.plot_stream(iteration, bed_particles, model_particles, pm.x_max, 10, 
                   available_vertices, to_file=False)
```

If `to_file` is not set to `False` then change it so it is.

## Setting Parameters

The model comes with a parameters file, `parameters.py`, which users can edit. The file is located in the **`model`** directory:

```
project
│   README.md
└───model
│   │   parameters.py
│   │   run.py
│   │   ...
```

Comments indicate what values have been tested for each paramter and within what range. Users are welcome to enter parameters beyond and below these bounds, _but results are untested_ and they could crash the model.

## Results

You can view the plots in the **`plots`** directory:
```
RCM_BedloadTransport
│   README.md
└───model
│   │   ...
└───plots <--- here!
│   │   ...
```

If you chose to run the model with plotting, then this directory will contain plots of the stream during each iteration as well as the final flux histogram. If you ran the model without plotting then the flux histogram should be the only plot produced.

**_Note_**: If doing another multiple runs, it is strongly suggested that you move the plots, rename them, or delete them. The model will overwrite any files whose name already exists.
# Wild plotter SVG driver
Simple SVG driver for a Wild TA-10 drawing/cutting plotter. For the purposes of Technologia Incognita Wild TA-10 plotter (https://wiki.techinc.nl/Wild_TA-10).

## Build instructions
Use `make` to build. The binary is called `wild_driver`. Only things required to build is Make and a C++17-compatible version of gcc.

## Usage instructions

Basic example:
```
$ ./wild_driver input.svg output.wild
parsing input file input.svg ...
input image of size 264.10mm by 264.10mm
transformed extents are:
    x: (0.00mm, 264.10mm)
    y: (0.00mm, 264.10mm)
rasterizing plot...
reordering paths...
saving visualization as vis.svg
using the pen
the tool is lowered
number of paths rendered: 72
done!
stats:
    points count: 77921
    total path length: 52022.41mm
    total free travel length: 274.42mm
```

More advanced example:
```
$ ./wild_driver input.svg --scale 0.2 --no_dry_run --translate 50 0 --box --min_step_size 1.5 -o output.wild
parsing input file input.svg ...
input image of size 264.10mm by 264.10mm
transformed extents are:
    x: (50.00mm, 102.82mm)
    y: (0.00mm, 52.82mm)
drawing the box outline of the plot only
reordering paths...
saving visualization as vis.svg
using the pen
the tool is lowered
number of paths rendered: 1
done!
stats:
    points count: 5
    total path length: 211.28mm
    total free travel length: 50.00mm
```

All options are documented in the `--help`:
```
$ ./wild_driver --help
Wild TA-10 plotter SVG driver
Usage: ./wild_driver [OPTIONS] input [output]

Positionals:
  input TEXT:FILE REQUIRED              Input svg file.
  output TEXT=out.wild                  Output plotter commands file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input svg file.
  -o,--output TEXT=out.wild             Output plotter commands file.
  --scale FLOAT:POSITIVE=1              Scale of the plot.
  --pre_translate [FLOAT,FLOAT]=0 0     Translate before scaling (ie. to center the svg file).
  --translate [FLOAT,FLOAT]=0 0         Translate after scaling (ie. to move the plot on the work surface).
  --box                                 Draw the box around the extents of the plot instead. Forces --dry_run unless otherwise specified.
  --dry_run,--no_dry_run{false}         Dry run (have the tool lifted when drawing). Off by default.
  --cut                                 Select cutting mode (initialize the cutting head).
  --lift_angle INT=20                   For cutting only. Threshold (in degrees) after which cutting head will be lifted before changing directions. 0 disables this feature.
  --speed INT:INT in [1 - 37]           Plotting/cutting speed in x*8mm/s.
  --no_order_paths{false}               Optimize ordering of rendered paths to save on free (tool up) travel length. Warning: may slow down for over 5000 paths rendered (O(n^2) complexity).
  --min_step_size FLOAT:POSITIVE=0.5    Minimal step size for the rasterization process in mm. Smaller values will increase resolution, but may slow down the plotting speed due to rs232 transfer rate limits.
  --hatch                               Hatch the inside of (filled in) svg shapes
  --hatch_density FLOAT:POSITIVE=2      Distance (in mm) between two consecutive hatch lines.
  --no_visualize{false}                 Visualize the path as an svg file.
```
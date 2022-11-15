# Ray Tracer (v0.2)

Physics simulation ray-tracing algorithm



## To-do list
- Initialise `Ray2D` WITHOUT freq and then `propagate()` and `populate(frequency_array)`  => NO MORE RAYPACK GROUPING BY FREQ!!!!
- Clean up physics/
- Use `Ray2D.calc_z` to generate heatmap (avoids same-ray energy overlap and energy cell overcrowding)
- ? Update `Ray2D` to v3.3 with support for an array of frequencies
- ? Rename `Simulation2D` to `Source2D` => `Simulation2D` can be plotted and takes an env and multiple `Source2D` objects
- Actually a Raypack is a source

## Method
### Base assumptions
* Ray path is frequency-independent (but rebounds aren't)

### Rebound calculation
Stuff

### Ray focusing
* Initial sweep
* Ray refining

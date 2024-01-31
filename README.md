# Ray Tracer (v0.2)

Physics simulation ray-tracing algorithm



## To-do list
- Clean up physics/
- Use `Ray2D.calc_z` to generate heatmap (avoids same-ray energy overlap and energy cell overcrowding)
- ? Update `Ray2D` to v3.3 with support for an array of frequencies
- ? Rename `Simulation2D` to `Source2D` => `Simulation2D` can be plotted and takes an env and multiple `Source2D` objects
- Actually a Raypack is a source

## Method

### Sound theory
http://www.sengpielaudio.com/calculator-FactorRatioLevelDecibel.htm
Sound rays…

### Base assumptions
* Ray path is frequency-independent (but rebounds aren't)

### Rebound calculation
Stuff
TODO: Write equation…
TODO: refactor code (./src)

### Ray focusing
* Initial sweep
* Ray refining

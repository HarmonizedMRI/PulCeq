# Convert Pulseq file to Pulserver format

## Requirements

* A block must contain at least one of the following Pulseq events: rf, gradient, delay, adc

For GE:
* RF raster time must be multiple of 2us 
* Gradient raster time must be multiple of 4us

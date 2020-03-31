IJ plugin for 

- neuron image patch extraction from the existing seto of mosaics using their annotations (.log files)
- possibility to rescale the outputs: annots are 500x500, outputs can be rescaled into custom square size
- visualization of the annotated neuron patches (Plugins>ViewLog) for inspection

Place 
m01.tif, m01.log,
m02.tif, m02.log,
...
m08.tif, m08.log
in the same directory (first argument of the plugin)
the rest of the arguments are the size of the output patch and the overlap (separate anntation overlap % to sample positive and negative patch), I guess 50-50 is what was used in the paper. Can use 60% positive / 40% negative to be more strict and make clearer separation.
package advantra.critpoint;

import advantra.feature.CircHaarFeat;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class TestCircHaarFeat implements PlugIn {

	public void run(String arg0) {
		int M = 16;
		CircHaarFeat c = new CircHaarFeat(M);
		c.createFeatures();
		c.showFeatures(512).show();
		
	}


}

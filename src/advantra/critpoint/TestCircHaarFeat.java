package advantra.critpoint;

import advantra.feature.CircHaarFeat;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class TestCircHaarFeat implements PlugIn {

	public void run(String arg0) {
		
		int M = 16;
		CircHaarFeat featsPi = new CircHaarFeat(M);
		featsPi.createFeatures();
		ImagePlus im = featsPi.showFeatures();
		im.show();
		
		// zoom several times
	    for (int i = 0; i < 7; i++) {
	    	im.getCanvas().zoomIn(0, 0);
	    }
	    
	    ImagePlus im1 = featsPi.showFeature(0);
	    im1.show();
		
	 // zoom several times
	    for (int i = 0; i < 7; i++) {
	    	im1.getCanvas().zoomIn(0, 0);
	    }
	    
	    
	}


}

package advantra.critpoint;

import advantra.feature.CircHaarFeat;
import ij.plugin.PlugIn;

public class TestCircHaarFeat implements PlugIn {

	public void run(String arg0) {
		int M = 16;
		CircHaarFeat c = new CircHaarFeat(M);
		c.createFeatures();
		c.showFeatures(512).show();
		
	}


}

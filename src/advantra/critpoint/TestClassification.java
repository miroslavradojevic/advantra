package advantra.critpoint;

import ij.ImagePlus;
import ij.plugin.PlugIn;

public class TestClassification  implements PlugIn {

	String folder_name = System.getProperty("user.dir");
	// to extract the features again
	ImagePlus 	train_img;
	double[][] 	train_loc;
	int[]		train_cls;
	// extraction parameters (standard deviations of the Gaussan used for smoothing/scaling)
	double 	sigma_1 = 2.0;
	double 	sigma_2 = 2.0;
	int		nr 		= 1;
	
	
	public void run(String arg0) {
		
	}

}

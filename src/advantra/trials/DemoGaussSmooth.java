package advantra.trials;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import imagescience.feature.Smoother;
import imagescience.image.FloatImage;

//this plug in is intended as a demo for imagescience
// gaussian derivative function 

public class DemoGaussSmooth implements PlugIn {
	
	ImagePlus im;
	ImageProcessor imp;
	
	public void run(String arg0) {
	
		IJ.run("Close All");
		im = IJ.openImage();
		im.show();
		
		// WRONG
		// FloatImage im_float = new FloatImage(im);// this wouldn't work directly because of the type of the image that was read
		
		// INSTEAD - CONVERT IT TO FLOAT IMAGEPRCESSOR AND USE THATH ONE TO INITIALIZE FLOATIMAGE
		imp = im.getProcessor().convertToFloat();
		FloatImage im_float = new FloatImage(new ImagePlus("Processed image", imp));
		
		Smoother sm = new Smoother(); // class to do the smoothing
		
		///// generic dialog /////
		String menu_title = "Demo Gaussian filtering";
		
		float scale = 2f;
		
		GenericDialog gd = new GenericDialog(menu_title, IJ.getInstance());
		
		gd.addMessage("Scale");
		gd.addNumericField( "gaussian smoothing scale :", scale, 2, 5, ">0.0" );
		
		gd.showDialog();
		
		if (gd.wasCanceled()) {
			return;
		}
		
		scale 		= (float)gd.getNextNumber();
		
		// scaling - result is stored in im_float
		sm.gauss(im_float, scale);
		im_float.imageplus().show();
		
	}

}

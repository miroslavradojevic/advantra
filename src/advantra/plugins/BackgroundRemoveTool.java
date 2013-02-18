package advantra.plugins;

import advantra.tools.OtsuBinarisation;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import imagescience.feature.Smoother;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import net.sf.ij.im3d.Util;
import net.sf.ij.im3d.morphology.Morpho;

public class BackgroundRemoveTool implements PlugInFilter  {
	
	ImagePlus 	img;
	ImagePlus 	out;
	String 		method;
	double 		scale;

	public void run(ImageProcessor arg0) {
		
		GenericDialog gd = new GenericDialog("Remove background", IJ.getInstance());
		
		String[] methods = new String[3];
		methods[0] = "Method_0";
		methods[1] = "Method_1";
		methods[2] = "Method_2";
		gd.addChoice("method", methods, methods[0]);
		//gd.addNumericField("Choose method:", method, 0, 3, "");
		
		gd.addMessage("0: orig-gaussan(scale)");
		gd.addMessage("1: binarize(otsu)+dilate+erode");
		gd.addMessage("2: binarize(otsu)+erode+dilate");
		
		gd.addNumericField("Scale (method 0 only)", scale, 1, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		method 		= (String)gd.getNextChoice(); //getNextNumber();
		scale		= (double)gd.getNextNumber();
		
		switch (method) {
		case "Method_0":
			Image I 						= new FloatImage(Image.wrap(img));
			Dimensions dims 				= I.dimensions();

			Image Ig 						= new FloatImage(dims);
			double[] aI 					= new double[dims.x];
			double[] aIg 					= new double[dims.x];
			
			Smoother sm = new Smoother();
			Ig = sm.gauss(I.duplicate(), scale); 
			
			Coordinates coords 	= new Coordinates();
			for (coords.z=0; coords.z<dims.z; ++coords.z) {
				
				I.axes(Axes.X);	Ig.axes(Axes.X);
				
				for (coords.y=0; coords.y<dims.y; ++coords.y) {
					I.get(coords, aI);	Ig.get(coords, aIg);
					for (int k = 0; k < aIg.length; k++) {
						aIg[k] = aI[k] - aIg[k];
					}
					Ig.set(coords, aIg);
				}

			}
			
			out = Ig.imageplus();
			out.setTitle("Original-GaussianAtScale"+String.valueOf(scale));
			out.show();
			
			break;

		case "Method_1":
			OtsuBinarisation otsu = new OtsuBinarisation(img);
			ImagePlus img_binarized = otsu.run();
			
			Morpho morpho = new Morpho();
			ImageStack after_dilate = Util.duplicateEmpty(img_binarized.getStack());
			ImageStack after_erode = Util.duplicateEmpty(img_binarized.getStack());
			
			morpho.dilate(img_binarized.getStack(), after_dilate);
			morpho.erode(after_dilate, after_erode);
			
			out = new ImagePlus("otsu+dilate+erode", after_erode);
			out.show();
			
			break;

			

		default:
			IJ.showMessage("Wrong method selected.");
			break;
		}
		
	}

	public int setup(String arg0, ImagePlus img) {
		
		this.img = img;
		
		// reset calibration before going further
		Calibration cal = new Calibration(img);
		cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1.0; cal.setUnit("pixel");
		img.setCalibration(cal);
		
		return DOES_ALL+NO_CHANGES;
	}

}

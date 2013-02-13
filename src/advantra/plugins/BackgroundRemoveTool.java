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
	int 		method;
	double 		scale;

	public void run(ImageProcessor arg0) {
		
		GenericDialog gd = new GenericDialog("    b     ackground", IJ.getInstance());
		gd.addMessage("0: orig-gaussan(scale)");
		gd.addMessage("1: binarize(otsu)+dilate+erode");
		gd.addNumericField("Choose method:", method, 0, 3, "");
		gd.addNumericField("Scale (method 0 only)", scale, 1, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		method 		= (int)gd.getNextNumber();
		scale		= (double)gd.getNextNumber();
		
		System.out.println("loaded method: "+method+"  scale: "+scale);
		
		switch (method) {
		case 0:
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

		case 1:
			OtsuBinarisation otsu = new OtsuBinarisation(img);
			ImagePlus img_binarized = otsu.run();
			System.out.println("binarized!"+(img_binarized.getType() != ImagePlus.GRAY8));
			img_binarized.show();
			
			Morpho morpho = new Morpho();
			System.out.println("made class");
			// = new ImageStack(img.getWidth(), img.getHeight(), img.getStackSize());
			ImageStack dest = Util.duplicateEmpty(img_binarized.getStack());
			morpho.dilate(img_binarized.getStack(), dest);
			System.out.println("done.");
			out = new ImagePlus("otsu+dilate+erode", dest);
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

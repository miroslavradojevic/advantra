package advantra.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.feature.Differentiator;
import imagescience.image.ByteImage;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class DifferentialStructure implements PlugInFilter {

	public static void main(String[] args){
		
	}
	
	
	/*
	 * imagej part
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	
	ImagePlus img;
	
	double 	sigma_1 = 1.0;
	double 	sigma_2 = 1.0;
	int		nr 		= 1;
	
	boolean Lw_enable, Lww_enable, Lvv_enable;
	
	Image[] Lx, Ly, Lxx, Lyy, Lxy;
	Image[] Lz, Lzz, Lxz, Lyz; // 3d image
	
	public void run(ImageProcessor arg0) {
		
		// dialog to enter input values
		GenericDialog gd = new GenericDialog("DifferentialStructure", IJ.getInstance());

		gd.addMessage("Choose scales:");
		
		gd.addNumericField( "sigma start:", sigma_1, 0, 5, "pix" );
		gd.addNumericField( "sigma end  :", sigma_2, 0, 5, "pix" );	// default value, number of zeros
		gd.addNumericField( "number of scales : ", nr,  0, 5, "");
		
		gd.addMessage("Choose descriptors:");
		gd.addMessage("1st order gauge");
		
		gd.addCheckbox("Lw", true);
		gd.addCheckbox("Lww", true);
		//gd.addCheckbox("Lv", true);
		gd.addCheckbox("Lvv", true);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		sigma_1 = (double)gd.getNextNumber();
		sigma_2 = (double)gd.getNextNumber();
		nr 		= (int)gd.getNextNumber();
		Lw_enable = (boolean)gd.getNextBoolean();
		Lww_enable = (boolean)gd.getNextBoolean();
		Lvv_enable = (boolean)gd.getNextBoolean();
		
		double[] sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = sigma_1+i*((sigma_2-sigma_1)/(nr-1));
			System.out.println("scale : "+sc[i]);
		}
		
		System.out.println("calculate derivatives... ");
		IJ.run(img, "32-bit", "");
		Image image = Image.wrap(img);
		
		image.imageplus().show();
		
		Differentiator df = new Differentiator();
		
		/*
		 * 1st derivatives
		 */
		Lx 		= new Image[nr];
		Ly 		= new Image[nr];
		Lxy 	= new Image[nr];
		
		for (int i = 0; i < nr; i++) {
			Lx[i] 	= df.run(image, sc[i], 1, 0, 0); Lx[i].imageplus().show();
			Ly[i] 	= df.run(image, sc[i], 0, 1, 0); Ly[i].imageplus().show();
			Lxy[i] 	= df.run(image, sc[i], 1, 1, 0); Lxy[i].imageplus().show();
		}
		
		System.out.println("done.");
		
	}

	public int setup(String arg0, ImagePlus im) {
		this.img = im;
		return DOES_8G+NO_CHANGES;
	}

	
	
}

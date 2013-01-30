package advantra.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.feature.Differentiator;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class DifferentialStructure implements PlugInFilter {

	public static void main(String[] args){
		
	}
	
	
	/*
	 * imagej part
	 * 
	 */
	
	ImagePlus img;
	
	double 	sigma_1 = 1.0;
	double 	sigma_2 = 1.0;
	int		nr 		= 1;
	
	boolean Lw_enable, Lww_enable, Lvv_enable;
	
	Image[] Lx, Ly, Lxx, Lyy, Lxy;
	Image[] Lz, Lzz, Lxz, Lyz; // 3D image
	
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
		
		sigma_1 	= (double)gd.getNextNumber();
		sigma_2 	= (double)gd.getNextNumber();
		nr 			= (int)gd.getNextNumber();
		Lw_enable 	= (boolean)gd.getNextBoolean();
		Lww_enable 	= (boolean)gd.getNextBoolean();
		Lvv_enable 	= (boolean)gd.getNextBoolean();
		
		double[] sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
			System.out.println("scale : "+sc[i]);
		}
		
		Image input_image = new FloatImage(Image.wrap(img));
		//input_image.imageplus().show();
		
		final Dimensions dims = input_image.dimensions();
		Differentiator df = new Differentiator();
		if (dims.z == 1) { // 2D case
			// calculate derivatives
			Lx 		= new Image[nr];
			Ly 		= new Image[nr];
			Lxx 	= new Image[nr];
			Lyy		= new Image[nr];
			Lxy 	= new Image[nr];
			
			for (int i = 0; i < nr; i++) {
				Lx[i] 	= df.run(input_image.duplicate(), sc[i], 1, 0, 0); 
				Ly[i] 	= df.run(input_image.duplicate(), sc[i], 0, 1, 0);
				Lxx[i] 	= df.run(input_image.duplicate(), sc[i], 2, 0, 0);
				Lyy[i] 	= df.run(input_image.duplicate(), sc[i], 0, 2, 0);
				Lxy[i] 	= df.run(input_image.duplicate(), sc[i], 1, 1, 0); 
				if (Lw_enable) {
					// gradient
					//double[] Lw = new double[];
				}
				
			}
			
		}
		else{ //3D case
			// calculate derivatives
			Lx 		= new Image[nr];
			Ly 		= new Image[nr];
			Lz 		= new Image[nr];
			Lxx 	= new Image[nr];
			Lyy		= new Image[nr];
			Lzz		= new Image[nr];
			Lxy 	= new Image[nr];
			Lxz 	= new Image[nr];
			Lyz 	= new Image[nr];
			
			for (int i = 0; i < nr; i++) {
				Lx[i] 	= df.run(input_image.duplicate(), sc[i], 1, 0, 0); 
				Ly[i] 	= df.run(input_image.duplicate(), sc[i], 0, 1, 0);
				Lz[i] 	= df.run(input_image.duplicate(), sc[i], 0, 0, 1);
				Lxx[i] 	= df.run(input_image.duplicate(), sc[i], 2, 0, 0);
				Lyy[i] 	= df.run(input_image.duplicate(), sc[i], 0, 2, 0);
				Lzz[i] 	= df.run(input_image.duplicate(), sc[i], 0, 0, 2);
				Lxy[i] 	= df.run(input_image.duplicate(), sc[i], 1, 1, 0);
				Lxz[i] 	= df.run(input_image.duplicate(), sc[i], 1, 0, 1);
				Lyz[i] 	= df.run(input_image.duplicate(), sc[i], 0, 1, 1);
			}
			
		}
		
		System.out.println("done calculating derivatives.");
		
		
		
	}

	public int setup(String arg0, ImagePlus im) {
		this.img = im;
		return DOES_8G+NO_CHANGES;
	}

	
	
}

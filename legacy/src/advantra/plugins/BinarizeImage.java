package advantra.plugins;

import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.io.File;

import advantra.tools.OtsuBinarisation;

public class BinarizeImage implements PlugInFilter  {
	
	ImagePlus img;
	
	public static void main(String[] args){
		
		if(args.length!=1) {
			System.out.println("usage: set image path or name as argument");
			return;
		}
		
		if((new File(args[0])).exists() || (new File(args[0])).isDirectory()){
			return;
		}
		
		ImagePlus img = new ImagePlus((new File(args[0]).getAbsolutePath()));
		
		OtsuBinarisation otsu = new OtsuBinarisation(img);
		ImagePlus out = otsu.run();
		if(out.getStack().getSize()>1) {
			(new FileSaver(out)).saveAsTiffStack(out.getTitle()+"_binarised.tif");
		}
		else {
			(new FileSaver(out)).saveAsTiff(out.getTitle()+"_binarised.tif");
		}
		
	}

	public void run(ImageProcessor imp) {
		OtsuBinarisation otsu = new OtsuBinarisation(img);
		ImagePlus out = otsu.run();
		out.show();
	}

	public int setup(String arg0, ImagePlus img) {
		this.img = img;
		return DOES_8G+NO_CHANGES;
	}

}

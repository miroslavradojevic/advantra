package advantra.trials;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import imagescience.image.Axes;
import imagescience.image.ColorImage;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.Image;

public class DemoImageToArray implements PlugIn {
	
	public void run(String arg){
		
		IJ.run("Close All");
		ImagePlus imp = IJ.openImage("");
		System.out.println("Image type...");
		if(imp.getType()==ImagePlus.GRAY16){
			System.out.println("GRAY16");
		}
		else if(imp.getType()==ImagePlus.GRAY8){
			System.out.println("GRAY8");
		}
		else if(imp.getType()==ImagePlus.GRAY32){
			System.out.println("GRAY32");
		}
		else if(imp.getType()==ImagePlus.COLOR_RGB){
			System.out.println("RGB");
		}
		else{
			System.out.println("SOME OTHER TYPE...");
		}
		
		imp.show();
		
		// represent it in arrays using imagescience... 
		Image img = Image.wrap(imp);
		Coordinates cin = new Coordinates(); // use default coords
		img.axes(Axes.X + Axes.Y);
		Dimensions dims = img.dimensions();
		double[][] imageZ = new double[dims.y][dims.x];
		img.get(cin, imageZ);
		// and the image is in imageZ array
		
		// save it as red image - fill the value of grey level as red channel
		Image outimg = Image.create(new Dimensions(dims.x, dims.y, 3), "imagescience.image.ColorImage");
		cin.z = 0;
		outimg.axes(Axes.X + Axes.Y);
		((ColorImage)outimg).component(ColorImage.RED);
		outimg.set(cin, imageZ);
		outimg.imageplus().show();
		
		
	}

	public void main(String[] args){
		System.out.println("main()");
	}

	
}

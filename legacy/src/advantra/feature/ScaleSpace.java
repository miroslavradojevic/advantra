package advantra.feature;

import imagescience.feature.Smoother;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

/*
 * this class will extract image in scale-space
 */

public class ScaleSpace {
	
	double Gamma = 1.0;

	public ScaleSpace(){}
	
//	public Vector<Image> extract(Image input2D, double[] sigmas){
//		
//		// normalized derivatives (derivatives are blob in this case)
//		Vector<Image> scsp = new Vector<Image>(sigmas.length);
//		
//		Smoother sm = new Smoother();
//		for (int i = 0; i < sigmas.length; i++) {
//			scsp.add(sm.gauss(input2D.duplicate(), sigmas[i]));
//		}
//		
//		//Image Scl = new FloatImage(input2D.dimensions()); // Lap.axes(Axes.X);
//		
//		return scsp;
//		
//	}
	
	public Image extractAsStack(Image input2D, double[] sigmas){
		
		Dimensions in_dims = input2D.dimensions();
		Image scale_space = new FloatImage(new Dimensions(in_dims.x, in_dims.y, sigmas.length));
		scale_space.axes(Axes.X);
		
		double[] aIn 		= new double[in_dims.x];
		double[] aScale_space 	= new double[in_dims.x];
		
		Smoother sm = new Smoother();
		Image in = null;
		
		Coordinates coords 	= new Coordinates();
		
		for (int i = 0; i < sigmas.length; ++i) {
			
			in = sm.gauss(input2D.duplicate(), sigmas[i]);
			in.axes(Axes.X);
			
			for (coords.y = 0; coords.y < in_dims.y; ++coords.y) {
				
				in.get(coords, aIn);
				
				for (int x = 0; x < in_dims.x; ++x) {
					aScale_space[x] = aIn[x];
				}
				coords.z = i;
				scale_space.set(coords, aScale_space);
				coords.z = 0;
				
			}
			
		}
		
		return scale_space;
		
	}
	
}

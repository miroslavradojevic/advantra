package advantra.feature;

import imagescience.feature.Hessian;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import java.util.Vector;

public class Ballness2D {
	
	public static Vector<Image> calculateArray(Image input2D, double[] sigmas) {
		
		Vector<Image> out = new Vector<Image>(sigmas.length);
		
		Dimensions dim = input2D.dimensions();
		double[] aL1 = new double[dim.x];
		double[] aBness = new double[dim.x];
		
		for (int i = 0; i < sigmas.length; i++) {
			
			Image Bness = new FloatImage(input2D.dimensions()); Bness.axes(Axes.X);
			
			Hessian hs = new Hessian();
			Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
			L.get(1).axes(Axes.X);
			
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				L.get(1).get(coords,aL1);
				for (int x=0; x<dim.x; ++x){
					aBness[x] = Math.abs(aL1[x]);// * aL2[x];
				}
				Bness.set(coords, aBness);
			}
			
			out.add(Bness);
			
		}
		return out;
	}
	
	public static Image calculateImg(Image input2D, double[] sigmas){
	
		Dimensions dim = input2D.dimensions();
		
		Image out = new FloatImage(new Dimensions(dim.x, dim.y, sigmas.length));
		out.axes(Axes.X);
		
		double[] aIn 			= new double[dim.x];
		double[] aScale_space 	= new double[dim.x];
		
		Hessian hs 			= new Hessian();
		
		Coordinates coords 	= new Coordinates();
		
		for (int i = 0; i < sigmas.length; ++i) {
			
			Vector<Image> L = hs.run(input2D, sigmas[i], true);
			L.get(1).axes(Axes.X); // smaller abs value
			
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				L.get(1).get(coords,aIn);
				for (int x = 0; x < dim.x; ++x) {
					aScale_space[x] = aIn[x];
				}
				
				coords.z = i;
				out.set(coords, aScale_space);
				coords.z = 0;
				
			}
			
		}
		
		return out;
	
	}
	
	

}

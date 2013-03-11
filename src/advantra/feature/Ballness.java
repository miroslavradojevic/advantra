package advantra.feature;

import imagescience.feature.Hessian;
import imagescience.feature.Smoother;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import java.util.Vector;

/*
 * this class will extract |l1| at different scales
 */

public class Ballness {
	
	public Ballness(){}
	
//	public Vector<Image> run(Image input2D, double[] sigmas) {
//		
//		Vector<Image> out = new Vector<Image>(sigmas.length);
//		
//		Dimensions dim = input2D.dimensions();
//		double[] aL1 = new double[dim.x];
//		//double[] aL2 = new double[dim.x];
//		double[] aBness = new double[dim.x];
//		
//		for (int i = 0; i < sigmas.length; i++) {
//			
//			Image Bness = new FloatImage(input2D.dimensions()); Bness.axes(Axes.X);
//			
//			Hessian hs = new Hessian();
//			Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
//			//L.get(0).axes(Axes.X);
//			L.get(1).axes(Axes.X);
//			
//			Coordinates coords 	= new Coordinates();
//			for (coords.y=0; coords.y<dim.y; ++coords.y) {
//				//L.get(0).get(coords,aL1);
//				L.get(1).get(coords,aL1);
//				for (int x=0; x<dim.x; ++x){
//					aBness[x] = Math.abs(aL1[x]);// * aL2[x];
//				}
//				Bness.set(coords, aBness);
//			}
//			
//			out.add(Bness);
//			
//		}
//		return out;
//	}
	
	public Image extractAsStack(Image input2D, double[] sigmas){
	
		Dimensions in_dims = input2D.dimensions();
		Image bness = new FloatImage(new Dimensions(in_dims.x, in_dims.y, sigmas.length));
		bness.axes(Axes.X);
		
		double[] aIn 		= new double[in_dims.x];
		double[] aScale_space 	= new double[in_dims.x];
		
		Hessian hs 			= new Hessian();
		Vector<Image> in 	= new Vector<Image>();
		
		Coordinates coords 	= new Coordinates();
		
		for (int i = 0; i < sigmas.length; ++i) {
			
			in = hs.run(input2D, sigmas[i], true); //sm.gauss(input2D.duplicate(), sigmas[i]);
			
			for (int x = 0; x < in_dims.x; ++x) {
				aScale_space[x] = aIn[x];
			}
			
		}
		
		return bness;
	
	}

}

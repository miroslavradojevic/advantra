package advantra.feature;

import java.util.Vector;

import imagescience.feature.Differentiator;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class Laplacian2D {
	
	public static Vector<Image> calculateArray(Image input2D, double[] sigmas) {
		
		Vector<Image> out = new Vector<Image>(sigmas.length);
		
		Image Lxx = null;
		Image Lyy = null;
		
		Dimensions dim = input2D.dimensions();
		double[] aLxx = new double[dim.x];
		double[] aLyy = new double[dim.x];
		double[] aLap = new double[dim.x];
		
		for (int i = 0; i < sigmas.length; i++) {
			
			Image Lap = new FloatImage(input2D.dimensions()); Lap.axes(Axes.X);
			
			Differentiator df = new Differentiator();
			Lxx = df.run(input2D.duplicate(), sigmas[i], 2, 0, 0); Lxx.axes(Axes.X);
			Lyy = df.run(input2D.duplicate(), sigmas[i], 0, 2, 0); Lyy.axes(Axes.X);
			
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				Lxx.get(coords,aLxx);
				Lyy.get(coords,aLyy);
				for (int x=0; x<dim.x; ++x){
					aLap[x] = aLxx[x] + aLyy[x];
				}
				Lap.set(coords, aLap);
			}
			
			out.add(Lap);
			
		}
		
		return out;
		
	}
	
	public static Image calculateImg(Image input2D, double[] sigmas){
		
		Image Lxx = null;
		Image Lyy = null;
		
		Dimensions dim = input2D.dimensions();
		
		Image out = new FloatImage(new Dimensions(dim.x, dim.y, sigmas.length));
		out.axes(Axes.X);
		
		double[] aLxx = new double[dim.x];
		double[] aLyy = new double[dim.x];
		double[] aLap = new double[dim.x];
		
		Coordinates coords 	= new Coordinates();
		
		for (int i = 0; i < sigmas.length; i++) {
			
			Differentiator df = new Differentiator();
			Lxx = df.run(input2D.duplicate(), sigmas[i], 2, 0, 0); Lxx.axes(Axes.X);
			Lyy = df.run(input2D.duplicate(), sigmas[i], 0, 2, 0); Lyy.axes(Axes.X);
			
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				Lxx.get(coords,aLxx);
				Lyy.get(coords,aLyy);
				for (int x=0; x<dim.x; ++x){
					aLap[x] = aLxx[x] + aLyy[x];
				}
				
				coords.z = i;
				out.set(coords, aLap);
				coords.z = 0;
			
			}
			
		}
		
		return out;
		
	
		
	}
	
}

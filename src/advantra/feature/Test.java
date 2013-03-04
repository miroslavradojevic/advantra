package advantra.feature;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class Test implements PlugIn {

	public void run(String arg0) {
		
//		double[][] mat = new double[][]
//				{{ 5, 	-4, 	 2},
//				 {-4, 	 5, 	 2},
//				 { 2,  	 2, 	-1}
//				};
			
			MyHessianFlanagan 	myhf = new MyHessianFlanagan();
//			MyHessianJama		myhj = new MyHessianJama();
//			
//			System.out.println("FLANAGAN SAYS:");
//			
//			double[] a = myhf.eigL(mat);
//			System.out.println(Utils.arrayToString(a));
//			
//			double[][] aa = myhf.eigV(mat);
//			System.out.println(Utils.arrayToString(aa));
//			
//			System.out.println("JAMA SAYS:");
//			
//			a = myhf.eigL(mat);
//			System.out.println(Utils.arrayToString(a));
//			
//			aa = myhf.eigV(mat);
//			System.out.println(Utils.arrayToString(aa));
			
			ImagePlus ip = IJ.openImage();
			ImageStack is = new ImageStack(ip.getWidth(), ip.getHeight());
			for (int i = 0; i < ip.getStackSize(); i++) {
				is.addSlice(ip.getStack().getProcessor(i+1).convertToFloat());
			}
			
			ImagePlus ip_float = new ImagePlus("converted", is);
			ip_float.show();
			Image img = new FloatImage(ip_float);
			
			long t1, t2;
			
			int dimensionality = (img.dimensions().z==1)? 2 : 3 ; 
			
			t1 = System.currentTimeMillis();
			Vector<Image> outputsFL = myhf.eigs(img.duplicate(), 2.0, true);
			t2 = System.currentTimeMillis();
			System.out.println("calculated in: "+((t2-t1)/1000f)+" sec.");
			
			for (int i = 0; i < 2*dimensionality; i++) {
				outputsFL.get(i).imageplus().show();
				
			}
			
			// save them as csv files if 2d
			
			if(dimensionality==2){
				
				saveAsCsv(outputsFL.get(0), "l2.csv");
				saveAsCsv(outputsFL.get(1), "l1.csv");
				saveAsCsv(outputsFL.get(2), "v11.csv");
				saveAsCsv(outputsFL.get(3), "v12.csv");
				
			}
			
			
	}
	
	public static void saveAsCsv(Image img, String path){
		
		FileWriter fw;
		Dimensions img_dims = img.dimensions();
		double[] aImg = new double[img_dims.x];
		img.axes(Axes.X);
		
		try{
			fw = new FileWriter(path);
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<img_dims.y; ++coords.y) {
				System.out.println("get y="+coords.y);
				img.get(coords, aImg);
				for (int x=0; x<img_dims.x; ++x){
					fw.write(aImg[x]+"");
					if(x==img_dims.x-1){
						fw.write("\n");
					}
					else{
						fw.write(", ");
					}
				}
			}
		} 
		catch(IOException exIO){
			System.out.printf("Couldn't write/open the file "+path+" ");
		}
		
		System.out.println("file saved!");
		
	}
	
}

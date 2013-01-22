package advantra.trials;

import advantra.shapes.Cylinder;
import advantra.shapes.Sphere;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.plugin.PlugIn;

public class DrawShapes  implements PlugIn {

	// this code draws the  defined geometric shapes
	public void run(String arg0) {
		
		IJ.run("Close All");
		int w = 128;
		int h = 128;
		int l = 60;
		
		ImagePlus img = NewImage.createByteImage("Shapes", w, h, l, NewImage.FILL_BLACK);
		img.show();
		
		// convert it to the format for processing - integer values - use byte[] to mediate
		byte[] read_layer;// = new byte[w*h];
		int[][]  img_stack = new int[l][w*h];
		
		for(int i = 1; i <= l; i++){
			read_layer 	= (byte[])img.getStack().getPixels(i);
			for (int j = 0; j < w*h; j++) {
				img_stack[i-1][j] = read_layer[j] & 0xff;
			}
		}
		
		// sphere
		Sphere sph = new Sphere(0, 0, 0, 20);
		// to store extracted voxels
		int[][] roi_coord 	= new int[3][Sphere.numberOfVoxInSphere(20)]; 
		int[]   roi_vals  	= new int[Sphere.numberOfVoxInSphere(20)];
		
		int count = sph.extractVox(img_stack, h, w, l, roi_coord, roi_vals);
		
		//set values from the sphere to 255
		for (int i = 0; i < count; i++) {
			int row = roi_coord[0][i]; 
			int col = roi_coord[1][i]; 
			int lay = roi_coord[2][i];
			img_stack[lay][col+row*w] = 255;
		}
		// another sphere
		sph = new Sphere(64, 64, 30, 10);
		int[][] roi_1_coord = new int[3][Sphere.numberOfVoxInSphere(10)];
		int[]	roi_1_vals	= new int[Sphere.numberOfVoxInSphere(10)];
		
		
		count = sph.extractVox(img_stack, h, w, l, roi_1_coord, roi_1_vals);
		System.out.println("coords: "+roi_1_coord.length+"x"+roi_1_coord[0].length);
		System.out.println("values: "+roi_1_vals.length);
		for (int i = 0; i < count; i++) {
			int row = roi_1_coord[0][i]; 
			int col = roi_1_coord[1][i]; 
			int lay = roi_1_coord[2][i];
			img_stack[lay][col+row*w] = 255;
		}
		// cylinder (center, radius, height, orientation)
		int R=10, H=35;
		Cylinder cyl = new Cylinder(40, 40, 40, R, H, 1, 1, 1);
		int allocate = Sphere.numberOfVoxInSphere((int)Math.ceil(Math.sqrt(R*R+H*H)));
		
		int[][] roi_2_coords	= new int[3][allocate];
		// the rest is actually not necessary here - just allocated to be able to 
		// execute the function
		int[]		roi_2_vals		= new int[allocate];
		double[][] 	roi_2_xy		= new double[2][allocate];
		double[]   	roi_2_dists		= new double[allocate];
		
		count = cyl.extractVox(
				img_stack, h, w, l, 
				roi_2_coords, 
				roi_2_vals,
				roi_2_xy,
				roi_2_dists); 
		
		System.out.println("count (in cylinder):" + count);
		
		for (int i = 0; i < count; i++) {
			int row = roi_2_coords[0][i]; 
			int col = roi_2_coords[1][i]; 
			int lay = roi_2_coords[2][i];
			img_stack[lay][col+row*w] = 255;
		}
		
		// convert it back to byte
		for (int i = 1; i <= l; i++) {
			read_layer = (byte[])img.getStack().getPixels(i);
			for (int j = 0; j < w*h; j++) {
				read_layer[j] = (byte)(img_stack[i-1][j]  & 0xff);
			}
		}
		
		// show it
		img.show();
		
		// save it
		IJ.saveAs(img, "Tiff", "output.tif");
		System.out.println("saved.");
		
		
	}
	
}

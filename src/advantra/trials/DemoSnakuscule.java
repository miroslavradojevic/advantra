package advantra.trials;

//import snake2D.Snake2D;
//import snake2D.Snake2DKeeper;
//import snake2D.Snake2DNode;
//import snake2D.Snake2DScale;
//
//import java.awt.geom.Point2D;

//import advantra.snakuscule.Snake2D;
import advantra.snakuscule.Snake2DKeeper;
//import advantra.snakuscule.Snake2DNode;
//import advantra.snakuscule.Snake2DScale;
import advantra.snakuscule.MySnakuscule;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.plugin.PlugIn;

public class DemoSnakuscule implements PlugIn {
	
	// this code tests the Snakuscules - the aim is to fit the snake (curve)
	// on the bright blob (neurite cross-section)
	
	//create the image with gaussian shaped blob at certain position
	public void run(String arg0) {
		
		IJ.run("Close All");
		
		int w = 128;
		int h = 128;
		
		ImagePlus img = NewImage.createByteImage("Gaussian blob", w, h, 1, NewImage.FILL_BLACK);
		
		int[]  img_array = new int[w*h];
		byte[] read_layer 	= (byte[])img.getStack().getPixels(1);
		
		
		for (int i = 0; i < w*h; i++) {
			img_array[i] = read_layer[i] & 0xff;
		}
		
		// make the gaussian blob
		int center_col = (int) (0.3 * w);
		int center_row = (int) (0.7 * h);
		
		System.out.println("center_row:  "+center_row);
		System.out.println("center_col:  "+center_col);
		
		double a = 255;
		double c = 5;
		
		for (int i = 0; i < w*h; i++) {
			
			int col = i % w;
			int row = i / w;
			double dist = Math.sqrt(Math.pow(center_row-row, 2) + Math.pow(center_col-col, 2));
			img_array[i] = (int) (a*Math.exp(-(dist*dist)/(2*c*c)));
			
		}
		// convert it back to byte and show
		read_layer = (byte[])img.getStack().getPixels(1);
		for (int i = 0; i < w*h; i++) {
			read_layer[i] = (byte)(img_array[i]  & 0xff);
		}
		
		img.show();
		
		// try to cover it with the snake
		// HOW IT'S USED
		double[] snake_nodes = new double[4]; // [node_1_x, node_1_y, node_2_x, node_2_y] 
		MySnakuscule snk 		= new MySnakuscule(img.getProcessor()); // it is actually derived from Snake2D
		Snake2DKeeper keeper 	= new Snake2DKeeper();
		keeper.justOptimize(snk, snake_nodes); // result will be stored in node_ variables
		// print the results
		System.out.println("Detected nodes:");
		System.out.println("nodes    : "+snake_nodes[0]+", "+snake_nodes[1]+", "+snake_nodes[2]+", "+snake_nodes[3]);
		double radius = 0.5 * Math.sqrt(Math.pow(snake_nodes[0]-snake_nodes[2], 2)+Math.pow(snake_nodes[1]-snake_nodes[3], 2));
		System.out.println("radius   : "+radius);
		System.out.println("center   : "+ ((snake_nodes[0]+snake_nodes[2])/2) +" real was ("+center_col+"), "+ ((snake_nodes[1]+snake_nodes[3])/2)+" real was ("+center_row+")");
		
		// make it a function
		
		
	}

}
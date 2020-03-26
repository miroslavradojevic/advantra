package advantra.snakuscule;

import ij.IJ;
import ij.process.ImageProcessor;

public class SnakusculeTools {

		public static double[] extractCircularBlob(ImageProcessor ip){
			
			double[] output = new double[3];
			
			double[] snake_nodes = new double[4]; // [node_1_x, node_1_y, node_2_x, node_2_y] 
			
			MySnakuscule snk 		= new MySnakuscule(ip);
			
			Snake2DKeeper keeper 	= new Snake2DKeeper();
			
			keeper.justOptimize(snk, snake_nodes); // result will be stored in node_ variables
			
			IJ.log("snake nodes after optimization: "+snake_nodes[0]+" <> "+snake_nodes[1]+" <> "+snake_nodes[2]+" <> "+snake_nodes[3]);
			
			// output[0] 	-----> radius
			output[0] = 0.5 * Math.sqrt(Math.pow(snake_nodes[0]-snake_nodes[2], 2)+Math.pow(snake_nodes[1]-snake_nodes[3], 2));
			
			// output[1..2] -----> center in image coordinates (row, col)
			output[1] = (snake_nodes[0]+snake_nodes[2])/2; // this one corresponds to column 
			output[2] = (snake_nodes[1]+snake_nodes[3])/2; // corresponds to row
			
			return output;
		}
		
}

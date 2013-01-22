package advantra.general;

//import ij.ImagePlus;
//import ij.gui.NewImage;

public class ArrayHandling {
	// TODO: change the name of this class because it works with indexing as well 
	
	/*
	 * contains some array or indexing operations 
	 */
	
	// convert index (of a 2d matrix element) to [row, col]
	public static int[] index2sub_2d(int index, int width){
		
		// width  is necessary to do the mapping
		int[] sub = new int[2];
		sub[0] = index / width;	// row
		sub[1] = index % width;	// col
		
		//return sub
		return sub;
		
	}
	
	public static int 	index2row_2d(int index, int width){
		
		// width  is necessary to do the mapping

		//return sub
		return index / width;	// row
		
	}
	
	public static int 	index2col_2d(int index, int width){
		
		// width  is necessary to do the mapping

		//return sub
		return index % width;	// col
		
	}
	
	/* 
	 * angle indexing
	 */
	
	public static enum IdxMode {
		LAST_INCLUDED, 
		LAST_EXCLUDED;
	}
	
	public static int value2index(
			double 		value_to_index,
			IdxMode 	indexing_mode, 
			double 		start_interval, 
			double 		end_interval, 
			int 		N){// zero indexing
		
		if(value_to_index<start_interval || value_to_index>end_interval){
			
			System.err.println("ArrayHandling:interval2index():value that's indexed is not within the set interval");
			System.exit(1);
			
		}
		
		int idx = 0;
		double a = value_to_index-start_interval;
		double A = end_interval-start_interval;
		
		// idx ~ a/A * N
		
		switch(indexing_mode) {
		
		case LAST_EXCLUDED:		
								idx = (int) Math.floor( a / ( A / N) );
								break;
		case LAST_INCLUDED:		
								idx = (int) Math.ceil( a / ( A /(2*(N-1))) / 2);					
								break;
		default:				
								System.err.println("ArrayHandling:value2index():This type of indexing is not supported.");
								System.exit(1);
								break;
		
		}
		
		return idx;
		
	}
	
	public static double index2value(
			int 		input_index,
			IdxMode 	indexing_mode, 
			double 		start_interval, 
			double 		end_interval, 
			int 		N){
		
		double out_val = 0;
		
//		double a = value_to_index-start_interval;
		double A = end_interval-start_interval;
		
		// idx ~ a/A * N
		
		switch(indexing_mode) {
		
		case LAST_EXCLUDED:		
								out_val = start_interval + input_index * (A/N);//Math.floor( a / ( A / N) );
								break;
		case LAST_INCLUDED:		
								out_val = start_interval + input_index * (A/(N-1));//(int) Math.ceil( a / ( A /(2*(N-1))) / 2);					
								break;
		default:				
								System.err.println("ArrayHandling:index2value():This type of indexing is not supported.");
								System.exit(1);
								break;
		
		}
		
		return out_val;
				
		
	}
	
	// convert 2d matrix [row, col] to element index
	public static int sub2index_2d(int row, int col, int width){
		
		// width  is necessary to do the mapping
		int index = row*width + col;
		
		return index; // return index
		
	}
	
	
	// reshape 1d float array to 2d array
	public static float[][] reshape(float[] input_array, int row, int col){
		
		float[][] output_array = new float[row][col];
		
		if(row*col != input_array.length){
			System.out.println("Array sizes do not match!");
			System.exit(-1);
		}
		
		for (int i = 0; i < input_array.length; i++) {
			output_array[i/col][i%row] = input_array[i];
		}
		
		return output_array;
		
	}
	
	// print 2d array on standard output
	public static void print2DArray(double[][] array_to_print){
		for (int i = 0; i < array_to_print.length; i++) {
			for (int j = 0; j < array_to_print[0].length; j++) {
				if(j==array_to_print[0].length-1){
					System.out.format("%9.5f ",array_to_print[i][j]);//at the end
				}
				else{
					System.out.printf("%9.5f, ",array_to_print[i][j]);
				}
			}
			System.out.println();
		}
	}
	
	
	// print 1d array on standard output
	public static void print1DArray(double[] array_to_print){
		for (int i = 0; i < array_to_print.length; i++) {
			
				if(i==array_to_print.length-1){
					System.out.printf("%9.5f \n",array_to_print[i]);//at the end
				}
				else{
					System.out.printf("%9.5f, ",array_to_print[i]);
				}
			
		}
	}
	// print 1d array on standard output
	public static void print1DArray(int[] array_to_print){
		for (int i = 0; i < array_to_print.length; i++) {
			
				if(i==array_to_print.length-1){
					System.out.printf("%5d \n",array_to_print[i]);//at the end
				}
				else{
					System.out.printf("%5d, ",array_to_print[i]);
				}
			
		}
	}
	
	// print 2d array on standard output
	public static void print2DArray(int[][] array_to_print){
		for (int i = 0; i < array_to_print.length; i++) {
			for (int j = 0; j < array_to_print[0].length; j++) {
				if(j==array_to_print[0].length-1){
					System.out.printf("%5d \n",array_to_print[i][j]);//at the end
				}
				else{
					System.out.printf("%5d, ",array_to_print[i][j]);
				}
			}
		}
	}
	
	public static int[] round(double[] in_array){
		int[] out_array = new int[in_array.length];
		for (int i = 0; i < in_array.length; i++) {
			out_array[i] = (int)Math.round(in_array[i]);
		}
		return out_array;
	}

	public static int[][] round(double[][] in_array){
		int[][] out_array = new int[in_array.length][in_array[0].length];
		for (int i = 0; i < in_array.length; i++) {
			for (int j = 0; j < in_array[0].length; j++) {
				out_array[i][j] = (int)Math.round(in_array[i][j]);
			}
		}
		return out_array;
	}
	
	public static double[] linspace(double radius_start, double radius_step, double radius_end){
		int radius_cnt = 0;
		
		for (double r = radius_start; r <= radius_end; r+=radius_step) {
			radius_cnt++;
		}
		
		double[] radiuses = new double[radius_cnt];
		
		for (int i = 0; i < radius_cnt; i++) {
			radiuses[i] = radius_start+i*radius_step;
		}
		
		return radiuses;
	}
	
}

package advantra.filter;

import advantra.general.ArrayHandling;

public class Template {
	
	public static float[] radius(int rows, int cols){

		// modify limits
		int dim2 = (cols%2==1)? (cols-1) : cols;
		int dim1 = (rows%2==1)? (rows-1) : rows;
		
		float[] radius 	= new float[dim1*dim2];
		
		for (int x = 0; x < cols; x++) {
			for (int y = 0; y < rows; y++) {
				
				double x_norm = (x - dim2/2 + x/dim2)/(double)dim2;
				double y_norm = (y - dim1/2 + y/dim1)/(double)dim1;
				
				float radius_value 	= 
						(float) Math.sqrt(Math.pow(x_norm, 2) + Math.pow(y_norm, 2)); 	// radius

				// position the values in array: y represents row, x represents column
				if ((y<dim1/2)&&(x<dim2/2)) {
					radius[ArrayHandling.sub2index_2d( y+(dim1+1)/2 , x+(dim2+1)/2, dim2)] 	
							= radius_value;// 1st -> 3rd quadrant (necessary later for fft)
				}
				else if ((y<dim1/2)&&(x>=dim2/2)) { 
					radius[ArrayHandling.sub2index_2d( y+(dim1+1)/2 , x-dim2/2, dim2)] 		
							= radius_value;// 2nd -> 4th quadrant
				}
				else if ((y>=dim1/2)&&(x<dim2/2)) {
					radius[ArrayHandling.sub2index_2d( y-dim1/2 , x+(dim2+1)/2, dim2)] 		
							= radius_value;// 4th -> 2nd
				}
				else if ((y>=dim1/2)&&(x>=dim2/2)){
					radius[ArrayHandling.sub2index_2d( y-dim1/2 , x-dim2/2, dim2)] 			
							= radius_value; // 3rd -> 1st
				}
			}
		}

		//radius[0] = 1;
		
		return radius;
	}
	
	public static float[] theta(int rows, int cols){
		
		// modify limits
		int dim2 = (cols%2==1)? (cols-1) : cols;
		int dim1 = (rows%2==1)? (rows-1) : rows;
		
		float[] theta 	= new float[dim1*dim2];
		
		for (int x = 0; x < cols; x++) {
			for (int y = 0; y < rows; y++) {
				
				double x_norm = (x - dim2/2 + x/dim2)/(double)dim2;
				double y_norm = (y - dim1/2 + y/dim1)/(double)dim1;
				
				float theta_value 	=
						(float) Math.atan2(-y_norm, x_norm); 					// theta
				
				// position the values in array: y represents row, x represents column
				if ((y<dim1/2)&&(x<dim2/2)) {
					theta[ArrayHandling.sub2index_2d( y+(dim1+1)/2 , x+(dim2+1)/2, dim2)]	
							= theta_value;
				}
				else if ((y<dim1/2)&&(x>=dim2/2)) { 
					theta[ArrayHandling.sub2index_2d( y+(dim1+1)/2 , x-dim2/2, dim2)]		
							= theta_value;
				}
				else if ((y>=dim1/2)&&(x<dim2/2)) {
					theta[ArrayHandling.sub2index_2d( y-dim1/2 , x+(dim2+1)/2, dim2)]		
							= theta_value;
				}
				else if ((y>=dim1/2)&&(x>=dim2/2)){
					theta[ArrayHandling.sub2index_2d( y-dim1/2 , x-dim2/2, dim2)] 			
							= theta_value;
				}
			}
		}
		return theta;
	}
	
}

package advantra.processing;

import flanagan.math.Matrix;

public class EigenCalculations {

	public static void getDominantEigVec3D(Matrix covariance_matrix, String ordering, double[] orient, double[] eigValsSort){
		
		double[]   eigVals = new double[3];
		double[][] eigVecs = new double[3][3];
		
		//double ratio = 0;
		
		if(orient.length != 3){
			System.err.printf("EigenCalculations::getDominantEigVec3D(): \noutput array orient should have length 3!");
			return ;
		}
		
		if(eigValsSort.length != 3){
			System.err.printf("EigenCalculations::getDominantEigVec3D(): \neigValsSort array should have length 3!");
			return ;
		}
		
		eigVals = covariance_matrix.getEigenValues();
		eigVecs = covariance_matrix.getEigenVectorsAsColumns();
		
		// take abs for all of them
		eigVals[0] = Math.abs(eigVals[0]);
		eigVals[1] = Math.abs(eigVals[1]);
		eigVals[2] = Math.abs(eigVals[2]);
		
		if(ordering.equals("highest")){
			// corresponding to highest absolute eigen value 
			
			if(eigVals[1]>=eigVals[0] && eigVals[1]>=eigVals[2]){
				// eigVals[1] has highest abs
				orient[0] = eigVecs[0][1]; 
				orient[1] = eigVecs[1][1];  
				orient[2] = eigVecs[2][1];
				
				eigValsSort[0] = eigVals[1];
				
				if(eigVals[0]>eigVals[2]){
					eigValsSort[1] = eigVals[0]; 
					eigValsSort[2] = eigVals[2];
				}
				else{
					eigValsSort[1] = eigVals[2]; 
					eigValsSort[2] = eigVals[0];
				}
				
			}
			else if(eigVals[2]>=eigVals[0] && eigVals[2]>=eigVals[1]){
				// eigVals[2] has highest abs
				orient[0] = eigVecs[0][2]; 
				orient[1] = eigVecs[1][2];  
				orient[2] = eigVecs[2][2];
				
				eigValsSort[0] = eigVals[2];
				
				if(eigVals[0]>eigVals[1]){
					eigValsSort[1] = eigVals[0]; 
					eigValsSort[2] = eigVals[1];
				}
				else{
					eigValsSort[1] = eigVals[1]; 
					eigValsSort[2] = eigVals[0];
				}
				
			}
			else{
				// eigVals[0] has highest abs
				orient[0] = eigVecs[0][0]; 
				orient[1] = eigVecs[1][0];  
				orient[2] = eigVecs[2][0];
				
				eigValsSort[0] = eigVals[0];
				
				if(eigVals[1]>eigVals[2]){
					eigValsSort[1] = eigVals[1]; 
					eigValsSort[2] = eigVals[2];
				}
				else{
					eigValsSort[1] = eigVals[2]; 
					eigValsSort[2] = eigVals[1];
				}
				
			}
			
		}
		else if(ordering.equals("lowest")){
			// corresponding to lowest abs eigen value
			
			if(eigVals[1]<=eigVals[0] && eigVals[1]<=eigVals[2]){
				// eigVals[1] has lowest abs
				orient[0] = eigVecs[0][1]; 
				orient[1] = eigVecs[1][1];  
				orient[2] = eigVecs[2][1];
				
				eigValsSort[0] = eigVals[1];
				
				if(eigVals[0]<eigVals[2]){
					eigValsSort[1] = eigVals[0]; 
					eigValsSort[2] = eigVals[2];
				}
				else{
					eigValsSort[1] = eigVals[2]; 
					eigValsSort[2] = eigVals[0];
				}				
				
			}
			else if(eigVals[2]<=eigVals[0] && eigVals[2]<=eigVals[1]){
				// eigVals[2] has lowest abs
				orient[0] = eigVecs[0][2]; 
				orient[1] = eigVecs[1][2];  
				orient[2] = eigVecs[2][2];
				
				eigValsSort[0] = eigVals[2];
				
				if(eigVals[0]<eigVals[1]){
					eigValsSort[1] = eigVals[0]; 
					eigValsSort[2] = eigVals[1];
				}
				else{
					eigValsSort[1] = eigVals[1]; 
					eigValsSort[2] = eigVals[0];
				}
				
			}
			else{
				// eigVals[0] has lowest abs
				orient[0] = eigVecs[0][0]; 
				orient[1] = eigVecs[1][0];  
				orient[2] = eigVecs[2][0];
				
				eigValsSort[0] = eigVals[0];
				
				if(eigVals[1]<eigVals[2]){
					eigValsSort[1] = eigVals[1]; 
					eigValsSort[2] = eigVals[2];
				}
				else{
					eigValsSort[1] = eigVals[2]; 
					eigValsSort[2] = eigVals[1];
				}
				
			}
			
		}
		else{
			System.err.println("EigenCalculations::getDominantEigVec3D(): \nordering has to be either 'highest' or 'lowest'.");
			return ;
		}
		
	}
	
}

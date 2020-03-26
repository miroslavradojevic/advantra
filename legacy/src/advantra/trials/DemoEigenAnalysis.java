package advantra.trials;

//import java.io.BufferedReader;
//import java.io.IOException;
//import java.io.InputStreamReader;

import advantra.general.FlanaganTools;
import flanagan.math.Matrix;
import ij.plugin.PlugIn;

public class DemoEigenAnalysis implements PlugIn  {
	
	public void run(String arg0) {
		
		/*
		System.out.print("Enter 2x2 matrix values : ");
		
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

	    String values_string = null;

	      try {
	    	  values_string = br.readLine();
	      } catch (IOException ioe) {
	         System.out.println("IO error trying to read your name!");
	         System.exit(1);
	      }
	      
	    String[] after_splitting = values_string.split(",\\s*"); 
	    int number_of_values = after_splitting.length;
	    double[] values = new double[number_of_values];
	    
	    for (int i = 0; i < number_of_values; i++) {
				
	    	values[i] = Double.valueOf(after_splitting[i]);
			
		}
	      
		Matrix aa = new Matrix(2,2);
		double[][] array22 = {{3.2, 1.8},{1.1, 4.5}};
		if(number_of_values>=4){
			array22[0][0] = values[0];
			array22[0][1] = values[1];
			array22[1][0] = values[2];
			array22[1][1] = values[3];
		}
		
		aa.setSubMatrix(0, 0, array22);
		System.out.println("Input Matrix: ");
		FlanaganTools.printMatrix(aa);
		
		////// eig val 
		double[] eigVal = new double[2];
		eigVal = aa.getEigenValues();
		System.out.println("Eigen values: ");
		for (int i = 0; i < eigVal.length; i++) {
			System.out.print(eigVal[i]+",  ");
		}
		System.out.println("");
		
		
		////// eig vec
		double[][] eigVec = new double[2][2];
		eigVec = aa.getEigenVectorsAsColumns();
		System.out.println("Eigen vectors (cols): ");
		for (int i = 0; i < eigVec.length; i++) {
			for (int j = 0; j < eigVec[0].length; j++) {
				System.out.print(eigVec[i][j]+",  ");
			}
			System.out.println("");
		}
		System.out.println("");
		
		*/
		///////////////////////////////////
		System.out.println("just a stupid test...");
		
		double [] 		eigVals = {0.0, 0.0};
//		double [][] 	reset22	= {{0.0, 0.0}, {0.0, 0.0}};
		System.out.println("Matrix: ");
//		Matrix bb 		= new Matrix(2,2);
		Matrix bb_add 	= new Matrix(2,2);
		
		bb_add.setElement(0, 0, -0.1901);
		bb_add.setElement(0, 1,  0.0992);
		bb_add.setElement(1, 0,  0.0992);
		bb_add.setElement(1, 1, -0.0880);
		FlanaganTools.printMatrix(bb_add);
		eigVals = bb_add.getEigenValues();
		System.out.println("eig. values: "+eigVals[0]+", "+eigVals[1]);
		
//		bb.plusEquals(bb_add);
//		
//		FlanaganTools.printMatrix(bb);
//		eigVals = bb.getSortedEigenValues();
//		System.out.println("eig. values: "+eigVals[0]+", "+eigVals[1]);
//		
//		System.out.println("Redefine Matrix: ");
//		
//		bb.setSubMatrix(0, 0, reset22);
		
		bb_add.setElement(0, 0, -0.4089);
		bb_add.setElement(0, 1,  0.1215);
		bb_add.setElement(1, 0,  0.1215);
		bb_add.setElement(1, 1, -0.3601);
		FlanaganTools.printMatrix(bb_add);
		eigVals = bb_add.getEigenValues();
		System.out.println("eig. values: "+eigVals[0]+", "+eigVals[1]);
		
		
//		System.out.println("bb:");
//		FlanaganTools.printMatrix(bb);
//		System.out.println("bb_add:");
//		FlanaganTools.printMatrix(bb_add);
//		bb.plusEquals(bb_add);
//		System.out.println("bb:");
//		FlanaganTools.printMatrix(bb);
//		eigVals = bb.getSortedEigenValues();
//		System.out.println("eig. values: "+eigVals[0]+", "+eigVals[1]);
		
		System.out.println("over...");
		
	}

	
	
}

package advantra.general;

import flanagan.math.Matrix;

public class FlanaganTools {

	public static void printMatrix(Matrix a) {
		
		int rows = a.getNrow();
		int cols = a.getNcol();
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				System.out.printf("%5.4f  ", a.getElement(i, j));
			}
			System.out.println();
		}
		System.out.println();
	}
	
}

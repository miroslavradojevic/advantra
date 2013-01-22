package advantra.processing;

public class Moments {
	
	// calculates the moments for voxels/images taking the 3D coordinates and the intensities
	// TODO: do it for sphere and give moments directly as output
	public static double extract_moments_3D(
			int[][] roi_coord, 
			int[] roi_vals, 
			int elements_in_region, 
			double[] moments){
		
		// moments[0] = centroid(x) 
		// moments[1] = centroid(y) 
		// moments[2] = centroid(z)
		// moments[3] = central_moment(x,y) 
		// moments[4] = central_moment(x,z) 
		// moments[5] = central_moment(y,z)
		// moments[6] = central_moment(x,x) 
		// moments[7] = central_moment(y,y)
		// moments[8] = central_moment(z,z)
		double M000  = 0;
		double M100  = 0;  double M010 = 0;  double M001 = 0;
		double M110  = 0;  double M101 = 0;  double M011 = 0;
		double M200  = 0;  double M020 = 0;  double M002 = 0;
		
		// 1st order... - M100, M010, M001...  2nd order mu110, mu101, mu011, mu200 etc...
		// image_region has to be 4xN where rows 1..3 are coords and row 4 is intensity
		
		// check if the size is ok
		if(roi_coord.length!=3){
			System.err.println("Moments:extract_moments_3D():\n " +
					"Matrix used to store 3d coordinates has to have 3 rows.");
			System.exit(-1);
		}
		
		if(roi_coord[0].length<=0){
			System.err.println("Moments:extract_moments_3D():\n " +
					"Matrix used to store 3d coordinates has to have at least one column.");
			System.exit(-1);
		}
		
		if(roi_coord[0].length != roi_vals.length){
			System.err.println("Moments:extract_moments_3D():\n " +
					"Output matrices for coordinates and grey-level intensistes have to have the same length.");
			System.exit(-1);
		}
		
		for (int i = 0; i < elements_in_region; i++) {
			
			M000 += roi_vals[i];
			
			M100 += roi_coord[0][i] * roi_vals[i];
			M010 += roi_coord[1][i] * roi_vals[i];
			M001 += roi_coord[2][i] * roi_vals[i];
			
			M110 += roi_coord[0][i] * roi_coord[1][i] * roi_vals[i];
			M101 += roi_coord[0][i] * roi_coord[2][i] * roi_vals[i];
			M011 += roi_coord[1][i] * roi_coord[2][i] * roi_vals[i];
			M200 += roi_coord[0][i] * roi_coord[0][i] * roi_vals[i];
			M020 += roi_coord[1][i] * roi_coord[1][i] * roi_vals[i];
			M002 += roi_coord[2][i] * roi_coord[2][i] * roi_vals[i];
			
		}
		
		if(M000>0){
			// first centroids
			moments[0] = M100 / M000; // centroid(x)
			moments[1] = M010 / M000; // centroid(y)
			moments[2] = M001 / M000; // centroid(z)
			// second central moments
			moments[3] = M110/M000 - moments[0]*moments[1]; // central_moment(x,y) 
			moments[4] = M101/M000 - moments[0]*moments[2]; // central_moment(x,z)
			moments[5] = M011/M000 - moments[1]*moments[2]; // central_moment(y,z)
			
			moments[6] = M200/M000 - moments[0]*moments[0]; // central_moment(x,x) 
			moments[7] = M020/M000 - moments[1]*moments[1]; // central_moment(y,y)
			moments[8] = M002/M000 - moments[2]*moments[2]; // central_moment(z,z)
		}
		else{
			System.err.println("all the intensities inside roi are zero! Moments will not be changed");
		}
		
		return M000;
		
	}

	public static void extract_moments_2D(
			double[][] image_coordinates, 
			int[] image_values, 
			int elements_in_region, 
			double[] moments){
		
		// moments[0] = centroid(x) 
		// moments[1] = centroid(y) 
		// moments[2] = central_moment(x,x) 
		// moments[3] = central_moment(y,y)
		// moments[4] = central_moment(x,y) 
		
		double M00  = 0;
		double M10  = 0;  
		double M01  = 0; 
		double M11  = 0;  
		double M20  = 0;  
		double M02 	= 0;
		
		// 1st order... - M10, M01,...  2nd order mu11, mu20, mu02 etc...
		// image_coordinates has to be 2xN where rows 1..2 are coords
		// image_values      has to be 1xN and contains image intensities
		
		if(image_coordinates.length!=2){
			System.err.println("Moments:extract_moments_2D():\n Matrix used as coordinates input for 2D centroid extraction has to have 2 rows.");
			System.exit(1);
		}
		
		if(moments.length!=5){
			System.err.println("Moments:extract_moments_2D():\n Vector for storing moments has to be 5 dimensional (5 moments extracted in 2D).");
			System.exit(1);
		}
		
		for (int i = 0; i < elements_in_region; i++) {
			
			M00 += image_values[i];
			M10 += image_values[i] * image_coordinates[0][i];
			M01 += image_values[i] * image_coordinates[1][i];
			M11 += image_values[i] * image_coordinates[0][i] * image_coordinates[1][i];
			M20 += image_values[i] * image_coordinates[0][i] * image_coordinates[0][i];
			M02 += image_values[i] * image_coordinates[1][i] * image_coordinates[1][i];
			
		}
		// first centroids
		moments[0] = M10 / M00; // centroid(x)
		moments[1] = M01 / M00; // centroid(y)
		// second central moments
		moments[2] = M11/M00 - moments[0]*moments[1]; // central_moment(x,y) 
		moments[3] = M20/M00 - moments[0]*moments[0]; // central_moment(x,x) 
		moments[4] = M02/M00 - moments[1]*moments[1]; // central_moment(y,y)
		
	}
	
	
}

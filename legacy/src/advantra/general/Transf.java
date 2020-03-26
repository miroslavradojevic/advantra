package advantra.general;

import flanagan.math.Matrix;

public class Transf {
	
	/*
	 * general calculations/transformations within coordinate systems
	 */

	public static void cartesian(
			double ax, 
			double ay, 
			double az, 
			double[] b, 
			double[] c
			){
		// take the 3d orientation [ax, ay, az] and calculate the values of the other two orthogonal vectors of the cartesian coordinate system
		// defined by [ax, ay, az]
		// outputs are b: [bx, by, bz] and c: [cx, cy, cz]
		double norm = Math.sqrt(Math.pow(ax, 2)+Math.pow(ay, 2)+Math.pow(az, 2));
		
		double yaw 		= 	Math.atan2(ay/norm, ax/norm);
		double pitch 	= 	-Math.atan2(az/norm, Math.sqrt(Math.pow(ax/norm, 2)+Math.pow(ay/norm, 2))) - Math.PI/2;// - Math.PI/2;
		
		double[][] rotZ = new double [4][4];// transformation matrices
		double[][] rotY = new double [4][4];
		// double[][] trans = new double [4][4]; // if necessary
		
		rotZ[0][0] = Math.cos(yaw); 	rotZ[0][1] = -Math.sin(yaw); 	rotZ[0][2] = 0; 				rotZ[0][3] = 0;
		rotZ[1][0] = Math.sin(yaw); 	rotZ[1][1] =  Math.cos(yaw); 	rotZ[1][2] = 0; 				rotZ[1][3] = 0;
		rotZ[2][0] = 0; 				rotZ[2][1] =  0;			 	rotZ[2][2] = 1; 				rotZ[2][3] = 0;
		rotZ[3][0] = 0; 				rotZ[3][1] =  0;			 	rotZ[3][2] = 0; 				rotZ[3][3] = 1;
		
		rotY[0][0] = Math.cos(pitch);	rotY[0][1] =  0; 				rotY[0][2] = Math.sin(pitch); 	rotY[0][3] = 0;
		rotY[1][0] = 0;				 	rotY[1][1] =  1; 				rotY[1][2] = 0; 				rotY[1][3] = 0;
		rotY[2][0] = -Math.sin(pitch); 	rotY[2][1] =  0;			 	rotY[2][2] = Math.cos(pitch);	rotY[2][3] = 0;
		rotY[3][0] = 0; 				rotY[3][1] =  0;				rotY[3][2] = 0; 				rotY[3][3] = 1;
		

		
		Matrix RtZ = new Matrix(rotZ);
		Matrix RtY = new Matrix(rotY);
		Matrix x_unit   = new Matrix(4,1); x_unit.setElement(0, 0, 1); x_unit.setElement(3, 0, 1);
		Matrix y_unit   = new Matrix(4,1); y_unit.setElement(1, 0, 1); y_unit.setElement(3, 0, 1);
		
		Matrix transformationX = RtZ.times(RtY.times(x_unit));
		Matrix transformationY = RtZ.times(RtY.times(y_unit));
		
		b[0] = transformationX.getElementCopy(0, 0); 
		b[1] = transformationX.getElementCopy(1, 0); 
		b[2] = transformationX.getElementCopy(2, 0); 
		
		c[0] = transformationY.getElementCopy(0, 0);
		c[1] = transformationY.getElementCopy(1, 0);
		c[2] = transformationY.getElementCopy(2, 0);
		
	}
	
	public static void rotate3xN(
			double ax,
			double ay,
			double az,
			final double[][] vecs_to_rotate // 3xN
			){
		/*
		 * will rotate all 3d points in the same manner as z was rotated from [0 0 1] to [ax ay az]
		 * input_vecs, rotated_vecs are 3xN where N>=1
		 */
		
		if(vecs_to_rotate.length!=3){
			System.err.println("Sphere:rotate3xN(): Vectors that are rotated need to have length 3 (coordinates in 3d space).");
			return;
		}
		
		// flanagan matrices - for multiplication
		Matrix input_vecs_matrix 	= new Matrix(vecs_to_rotate.length+1, vecs_to_rotate[0].length, 1); // allocate it and fill with ones
		input_vecs_matrix.setSubMatrix(0, 0, vecs_to_rotate);
		
		// define the transformation wrt to [ax, ay, az]
		double yaw 		= 	Math.atan2(ay, ax); 							// azimuth, theta
		double pitch 	= 	Math.acos(az / Math.sqrt(ax*ax+ay*ay+az*az)); 	// inclination, phi
		
		Matrix RtZ = new Matrix(4,4,0);
		Matrix RtY = new Matrix(4,4,0);
		
		// set elements for the transformation
		//-----------------------------------------------------------------------------------------------------------------------------------
		RtZ.setElement(0, 0, Math.cos(yaw));		RtZ.setElement(0, 1, -Math.sin(yaw));	
		RtZ.setElement(1, 0, Math.sin(yaw));		RtZ.setElement(1, 1,  Math.cos(yaw));
																							RtZ.setElement(2, 2,  1);
																																	RtZ.setElement(3, 3,  1);
		//-----------------------------------------------------------------------------------------------------------------------------------																											
		RtY.setElement(0, 0, Math.cos(pitch));												RtY.setElement(0, 2, Math.sin(pitch));	
													RtY.setElement(1, 1, 1);	
		RtY.setElement(2, 0, -Math.sin(pitch));												RtY.setElement(2, 2,  Math.cos(pitch));
																																	RtY.setElement(3, 3,  1);	
		//-----------------------------------------------------------------------------------------------------------------------------------																													
		
		Matrix output_vecs_matrix = RtZ.times(RtY.times(input_vecs_matrix));																															
		// set the values to ones obtained after rotation... take the submatrix of the obtained one
		for (int i = 0; i < vecs_to_rotate.length; i++) {
			for (int j = 0; j < vecs_to_rotate[0].length; j++) {
				vecs_to_rotate[i][j] = output_vecs_matrix.getElement(i, j); //, 2, number_of_input_vecs-1).getArrayCopy();
			}
		}
		
	}

	public static void rotateNx3(
			double ax,
			double ay,
			double az,
			final double[][] vecs_to_rotate // Nx3
			){
		/*
		 * will rotate all 3d points in the same manner as z was rotated from [0 0 1] to [ax ay az]
		 * input_vecs, rotated_vecs are 3xN where N>=1
		 */
		
		// flanagan matrices - for multiplication
		Matrix input_vecs_matrix 	= new Matrix(vecs_to_rotate.length, 4, 1); // allocate it and fill with ones
		input_vecs_matrix.setSubMatrix(0, 0, vecs_to_rotate);
		
		// define the transformation wrt to [ax, ay, az]
		double yaw 		= 	Math.atan2(ay, ax); 							// azimuth, theta
		double pitch 	= 	Math.acos(az / Math.sqrt(ax*ax+ay*ay+az*az)); 	// inclination, phi
		
		Matrix RtZ = new Matrix(4,4,0);
		Matrix RtY = new Matrix(4,4,0);
		
		// set elements for the transformation
		//-----------------------------------------------------------------------------------------------------------------------------------
		RtZ.setElement(0, 0, Math.cos(yaw));		RtZ.setElement(0, 1, -Math.sin(yaw));	
		RtZ.setElement(1, 0, Math.sin(yaw));		RtZ.setElement(1, 1,  Math.cos(yaw));
																							RtZ.setElement(2, 2,  1);
																																	RtZ.setElement(3, 3,  1);
		//-----------------------------------------------------------------------------------------------------------------------------------																											
		RtY.setElement(0, 0, Math.cos(pitch));												RtY.setElement(0, 2, Math.sin(pitch));	
													RtY.setElement(1, 1, 1);	
		RtY.setElement(2, 0, -Math.sin(pitch));												RtY.setElement(2, 2,  Math.cos(pitch));
																																	RtY.setElement(3, 3,  1);	
		//-----------------------------------------------------------------------------------------------------------------------------------																													
		
		Matrix tf = RtZ.times(RtY).transpose();			
		
		Matrix output_vecs_matrix = input_vecs_matrix.times(tf);
		
		for (int i = 0; i < vecs_to_rotate.length; i++) {
			for (int j = 0; j < 3; j++) {
				vecs_to_rotate[i][j] = output_vecs_matrix.getElement(i, j); 
			}
		}
		
	}

	
	public static void rotate(
			double inclination,
			double azimuth,
			final double[][] vecs_to_rotate // 3xN
			){
		/*
		 * will rotate all 3d points in the same manner as z was rotated from [0 0 1] being azimuth=all inclination=0
		 * input_vecs, rotated_vecs are 3xN where N>=1
		 */
		
		Matrix input_vecs_matrix 	= new Matrix(3+1, vecs_to_rotate[0].length, 1); // allocate it and fill with ones
		input_vecs_matrix.setSubMatrix(0, 0, vecs_to_rotate);
		
		// define the transformation wrt to [ax, ay, az]
		//double yaw 		= 	azimuth; 
		//double pitch 	= 	inclination; 	
		
		Matrix RtZ = new Matrix(4,4,0);
		Matrix RtY = new Matrix(4,4,0);
		
		// set elements for the transformation
		//-----------------------------------------------------------------------------------------------------------------------------------
		RtZ.setElement(0, 0, Math.cos(azimuth));		RtZ.setElement(0, 1, -Math.sin(azimuth));	
		RtZ.setElement(1, 0, Math.sin(azimuth));		RtZ.setElement(1, 1,  Math.cos(azimuth));
																							RtZ.setElement(2, 2,  1);
																																	RtZ.setElement(3, 3,  1);
		//-----------------------------------------------------------------------------------------------------------------------------------																											
		RtY.setElement(0, 0, Math.cos(inclination));												RtY.setElement(0, 2, Math.sin(inclination));	
													RtY.setElement(1, 1, 1);	
		RtY.setElement(2, 0, -Math.sin(inclination));												RtY.setElement(2, 2,  Math.cos(inclination));
																																	RtY.setElement(3, 3,  1);	
		//-----------------------------------------------------------------------------------------------------------------------------------																													
		
		Matrix output_vecs_matrix = RtZ.times(RtY.times(input_vecs_matrix));																															
		// set the values to ones obtained after rotation... take the (3xN) submatrix of the obtained (4xN) matrix
		for (int i = 0; i < vecs_to_rotate.length; i++) {
			for (int j = 0; j < vecs_to_rotate[0].length; j++) {
				vecs_to_rotate[i][j] = output_vecs_matrix.getElement(i, j); //, 2, number_of_input_vecs-1).getArrayCopy();
			}
		}
		
	}
	
	public static void rotateNx3(
			double inclination,
			double azimuth,
			double[][] vec_to_rotate // Nx3
			){
		
		Matrix input_vecs_matrix 	= new Matrix(vec_to_rotate.length, 3+1, 1); // allocate it and fill with ones
		input_vecs_matrix.setSubMatrix(0, 0, vec_to_rotate);
		
		// define the transformation wrt to [ax, ay, az]
		
		Matrix RtZ = new Matrix(4,4,0);
		Matrix RtY = new Matrix(4,4,0);
		
		// set elements for the transformation
		//-----------------------------------------------------------------------------------------------------------------------------------
		RtZ.setElement(0, 0, Math.cos(azimuth));		RtZ.setElement(0, 1, -Math.sin(azimuth));	
		RtZ.setElement(1, 0, Math.sin(azimuth));		RtZ.setElement(1, 1,  Math.cos(azimuth));
																							RtZ.setElement(2, 2,  1);
																																	RtZ.setElement(3, 3,  1);
		//-----------------------------------------------------------------------------------------------------------------------------------																											
		RtY.setElement(0, 0, Math.cos(inclination));												RtY.setElement(0, 2, Math.sin(inclination));	
													RtY.setElement(1, 1, 1);	
		RtY.setElement(2, 0, -Math.sin(inclination));												RtY.setElement(2, 2,  Math.cos(inclination));
																																	RtY.setElement(3, 3,  1);	
		//-----------------------------------------------------------------------------------------------------------------------------------																													
		
		Matrix tf = RtZ.times(RtY).transpose();																															
																																	
		Matrix output_vecs_matrix = input_vecs_matrix.times(tf);																															
		
		for (int i = 0; i < vec_to_rotate.length; i++) {
			for (int j = 0; j < 3; j++) {
				vec_to_rotate[i][j] = output_vecs_matrix.getElement(i, j); 
			}
		}
		
	}

	public static void rotate(double ax, double ay, double az, double[] vec_to_rotate){
		
		// define the transformation wrt to [ax, ay, az]
		double theta 		= 	Math.atan2(ay, ax); 								// azimuth
		double phi 			= 	Math.acos(az / Math.sqrt(ax*ax+ay*ay+az*az)); 		// inclination
		
		double[] out = new double[3];
		
		out[0] = 
				Math.cos(theta)*Math.cos(phi) 	* vec_to_rotate[0]
				-Math.sin(theta)				* vec_to_rotate[1]
				+Math.cos(theta)*Math.sin(phi) 	* vec_to_rotate[2];
		out[1] = 
				Math.sin(theta)*Math.cos(phi) 	* vec_to_rotate[0]
				+Math.cos(theta) 				* vec_to_rotate[1]
				+Math.sin(theta)*Math.sin(phi) 	* vec_to_rotate[2];
		out[2] = 
				-Math.sin(phi) 					* vec_to_rotate[0]
				+Math.cos(phi) 					* vec_to_rotate[2];
		
		vec_to_rotate[0] = out[0];
		vec_to_rotate[1] = out[1];
		vec_to_rotate[2] = out[2];
		
		
	}

	public static void shift(double ax, double ay, double az, double[] vec_to_shift){
		vec_to_shift[0] += ax;
		vec_to_shift[1] += ay;
		vec_to_shift[2] += az;
	}
	
	/*
	 * some vector calculations
	 */
	
	public static void 		thisVectorScale	(double constant, 	double[] vector){ 	// result is stored in vector!!!
		for (int i = 0; i < vector.length; i++) {
			vector[i] *= constant;
		}
	}
	
	public static double[] 	vectorScale	(double constant, 	double[] vector){		// result is in output (newly allocated)
		
		if(vector.length<=0){
			System.err.println("Transf:vectorScale():input vector has to have length > 0");
			System.exit(1);
		}
		
		double[] out = new double[vector.length];
		
		for (int i = 0; i < vector.length; i++) {
			out[i] = constant * vector[i];
		}
		
		return out;
	}
	
	public static void 		subtractVects(double[] a, double[] b, double[] a_minus_b){
		
		if	(
				a.length == b.length &&
				a.length == a_minus_b.length){
			// then all three are equal length
			for (int i = 0; i < a.length; i++) {
				a_minus_b[i] = a[i] - b[i];
			}
		}
		else{
			System.err.println("Vectors must have the same length to be subtracted!");
			System.exit(1);
		}
	}
	
	public static void 		subtractVectsAndStoreInFirst(double[] a, double[] b){
		
		if	(
				a.length == b.length){
			// then all three are equal length
			for (int i = 0; i < a.length; i++) {
				a[i] = a[i] - b[i];
			}
		}
		else{
			System.err.println("Vectors must have the same length to be subtracted!");
			System.exit(1);
		}
	}
	
	public static double 	vectorNorm(double[] vector){
		
		double norm = 0;
		
		for (int i = 0; i < vector.length; i++) {
			norm += vector[i] * vector[i];
		}
		
		return Math.sqrt(norm);
	}
	
	public static void 		vectorNormalize(double[] vector){
		
		double norm = 0;
		
		for (int i = 0; i < vector.length; i++) {
			norm += vector[i] * vector[i];
		}
		
		norm =  Math.sqrt(norm);
		
		for (int i = 0; i < vector.length; i++) {
			vector[i] = (norm>0)? (vector[i]/norm) : 0 ;
		}
	}	
	
	public static double 	dotProd(double[] a, double[] b){
		
		double dot_product = 0;
		
		if(a.length != b.length){
			System.err.println("Transf:dotProd():vectors have to have the same length!");
			System.exit(1);
		}
		else{
			
			for (int i = 0; i < a.length; i++) {
				dot_product += a[i] * b[i];
			}
		}
		
		return dot_product;
		
	}
	
	public static boolean 	isUnitLength(double[] vector){
		
		boolean isUnit = true;
		
		if(!(vectorNorm(vector)<1.001 && vectorNorm(vector)>0.999)){
			isUnit = false;
		}
		
		return isUnit;
	}
	
	public static double 	distance_point_to_line(double[] line_base_point, double[] line_unit_direction, double[] point){ 
		
		// line is in vector form (point+unit_direction)
		// general form
		double distance_from_line = 0;
		
		if		(
				line_base_point.length == point.length 					&&
				line_base_point.length == line_unit_direction.length 	
				){
			
			if(vectorNorm(line_unit_direction)< 0.99 || vectorNorm(line_unit_direction)> 1.01){
				System.err.println(
						"Transf:distance_point_to_line():Line's unit orientation vector must have unit length! \nCurrent length is "+
								vectorNorm(line_unit_direction));
				System.exit(1);
			}
			
			double[] subtract_ponts = new double[3]; 						// a - p define it...
			
			subtractVects(line_base_point, point, subtract_ponts); 			// a - p calculate it...
			
			subtractVectsAndStoreInFirst(subtract_ponts, vectorScale(dotProd(subtract_ponts, line_unit_direction), line_unit_direction));
			
			distance_from_line = vectorNorm(subtract_ponts);
		
		}
		else{
			System.err.println("Transf:distance_point_to_line():Vectors must have the same length in order for the distance to be calculated!");
			System.exit(1);
		}
		
		return distance_from_line;
	}
	
	public static double 	distance_point_to_line_3d(double[] line_base_point, double[] line_unit_direction, double[] point){ 
		
		// line is in vector form (point+unit_direction)

		double distance_from_line = 0;
		
		if		(
				line_base_point.length != 3 	||
				line_unit_direction.length != 3 ||
				point.length != 3
				){
			
			System.err.println(
					"Transf:distance_point_to_line_3d():Vectors must have length 3!");
			System.exit(1);
			
			
		}
			
		if(vectorNorm(line_unit_direction)< 0.99 || vectorNorm(line_unit_direction)> 1.01){
				System.err.println(
						"Transf:distance_point_to_line():Line's unit orientation vector must have unit length! \nCurrent length is "+
								vectorNorm(line_unit_direction));
				System.exit(1);
		}
			
		double[] distance_3d = new double[3]; 	
		
		// line_base_point - point
		distance_3d[0] = line_base_point[0] - point[0];
		distance_3d[1] = line_base_point[1] - point[1];
		distance_3d[2] = line_base_point[2] - point[2]; 	// it contains difference a-p
		
		double proj = 
				distance_3d[0] * line_unit_direction[0] +
				distance_3d[1] * line_unit_direction[1] +
				distance_3d[2] * line_unit_direction[2]; 	// dotproduct((a-p), n)
		
		distance_3d[0] = distance_3d[0] - proj * line_unit_direction[0];
		distance_3d[1] = distance_3d[1] - proj * line_unit_direction[1];
		distance_3d[2] = distance_3d[2] - proj * line_unit_direction[2]; // || (a-p) - ((a-p)*n) * n ||
 		
		distance_from_line = vectorNorm(distance_3d);
		
		return distance_from_line;
	}	

	public static double 	distance_point_to_line_3d(double[] line_base_point, double[] line_unit_direction, int[] point){ 
		
		// line is in vector form (point+unit_direction)

		double distance_from_line = 0;
		
		if		(
				line_base_point.length != 3 	||
				line_unit_direction.length != 3 ||
				point.length != 3
				){
			
			System.err.println(
					"Transf:distance_point_to_line_3d():Vectors must have length 3!");
			System.exit(1);
			
			
		}
			
		if(vectorNorm(line_unit_direction)< 0.99 || vectorNorm(line_unit_direction)> 1.01){
				System.err.println(
						"Transf:distance_point_to_line():Line's unit orientation vector must have unit length! \nCurrent length is "+
								vectorNorm(line_unit_direction));
				System.exit(1);
		}
			
		// http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
		
		double[] distance_3d = new double[3]; 	
		
		// line_base_point - point
		distance_3d[0] = line_base_point[0] - point[0];
		distance_3d[1] = line_base_point[1] - point[1];
		distance_3d[2] = line_base_point[2] - point[2]; 	// it contains difference a-p
		
		double proj = 
				distance_3d[0] * line_unit_direction[0] +
				distance_3d[1] * line_unit_direction[1] +
				distance_3d[2] * line_unit_direction[2]; 	// dotproduct((a-p), n)
		
		distance_3d[0] = distance_3d[0] - proj * line_unit_direction[0];
		distance_3d[1] = distance_3d[1] - proj * line_unit_direction[1];
		distance_3d[2] = distance_3d[2] - proj * line_unit_direction[2]; // || (a-p) - ((a-p)*n) * n ||
 		
		distance_from_line = vectorNorm(distance_3d);
		
		return distance_from_line;
	}

	public static double 	distance_point_to_line_3d(int[] line_base_point, double[] line_unit_direction, int[] point){ 
		
		// line is in vector form (point+unit_direction)

		double distance_from_line = 0;
		
		if		(
				line_base_point.length != 3 	||
				line_unit_direction.length != 3 ||
				point.length != 3
				){
			
			System.err.println(
					"Transf:distance_point_to_line_3d():Vectors must have length 3!");
			System.exit(1);
			
			
		}
			
		if(vectorNorm(line_unit_direction)< 0.99 || vectorNorm(line_unit_direction)> 1.01){
				System.err.println(
						"Transf:distance_point_to_line():Line's unit orientation vector must have unit length! \nCurrent length is "+
								vectorNorm(line_unit_direction));
				System.exit(1);
		}
			
		double[] distance_3d = new double[3]; 	
		
		// line_base_point - point
		distance_3d[0] = line_base_point[0] - point[0];
		distance_3d[1] = line_base_point[1] - point[1];
		distance_3d[2] = line_base_point[2] - point[2]; 	// it contains difference a-p
		
		double proj = 
				distance_3d[0] * line_unit_direction[0] +
				distance_3d[1] * line_unit_direction[1] +
				distance_3d[2] * line_unit_direction[2]; 	// dotproduct((a-p), n)
		
		distance_3d[0] = distance_3d[0] - proj * line_unit_direction[0];
		distance_3d[1] = distance_3d[1] - proj * line_unit_direction[1];
		distance_3d[2] = distance_3d[2] - proj * line_unit_direction[2]; // || (a-p) - ((a-p)*n) * n ||
 		
		distance_from_line = vectorNorm(distance_3d);
		
		return distance_from_line;
	}

	public static double 	distance_point_to_line_3d(int[] line_base_point, double[] line_unit_direction, double[] point){ 
		
		// line is in vector form (point+unit_direction)

		double distance_from_line = 0;
		
		if		(
				line_base_point.length != 3 	||
				line_unit_direction.length != 3 ||
				point.length != 3
				){
			
			System.err.println(
					"Transf:distance_point_to_line_3d():Vectors must have length 3!");
			System.exit(1);
			
			
		}
			
		if(vectorNorm(line_unit_direction)< 0.99 || vectorNorm(line_unit_direction)> 1.01){
				System.err.println(
						"Transf:distance_point_to_line():Line's unit orientation vector must have unit length! \nCurrent length is "+
								vectorNorm(line_unit_direction));
				System.exit(1);
		}
			
		double[] distance_3d = new double[3]; 	
		
		// line_base_point - point
		distance_3d[0] = line_base_point[0] - point[0];
		distance_3d[1] = line_base_point[1] - point[1];
		distance_3d[2] = line_base_point[2] - point[2]; 	// it contains difference a-p
		
		double proj = 
				distance_3d[0] * line_unit_direction[0] +
				distance_3d[1] * line_unit_direction[1] +
				distance_3d[2] * line_unit_direction[2]; 	// dotproduct((a-p), n)
		
		distance_3d[0] = distance_3d[0] - proj * line_unit_direction[0];
		distance_3d[1] = distance_3d[1] - proj * line_unit_direction[1];
		distance_3d[2] = distance_3d[2] - proj * line_unit_direction[2]; // || (a-p) - ((a-p)*n) * n ||
 		
		distance_from_line = vectorNorm(distance_3d);
		
		return distance_from_line;
	}

	/*
	 * angles
	 */
	
	public static double angle_between_vectors_3d_radians(double[] vector1, double[] vector2){
		if((vector1.length!=3) || (vector2.length!=3)){
			System.err.println("Transf:angle_radian_3d():\n both input vectors have to have length 3.");
			System.exit(1);
		}
		// TODO: there's no check whether the vectors are unit length
		// but we'll be using mostly those from Cylinder class, which are constrained to unit length
		double dot_prod = vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2];
		dot_prod = (dot_prod<=1)? dot_prod : 1.0 ; // so that Math.acos doesn't give NaN
		return Math.acos(dot_prod); // 0-pi value
		
	}
	
	/*
	 * wrapping angles
	 */
	
	// wraps it to (-pi/2, pi/2]
	public static double angle_wrap_pi_2(double input_angle){ 
		
		
		while(input_angle>Math.PI/2){
			input_angle -= Math.PI;
		}
		
		while(input_angle<=-Math.PI/2){
			input_angle += Math.PI;
		}
		
		return input_angle;
	}

	// wraps it to (-90, 90]
	public static double angle_wrap_90(double input_angle){ 
		
		while(input_angle>90){
			input_angle -= 180;
		}
		
		while(input_angle<=-90){
			input_angle += 180;
		}
		
		return input_angle;
	}
	
	//TODO: wrappings 0-pi, 0-180
	
	// wraps it to (-pi, pi]
	public static double angle_wrap_pi(double input_angle){ 
		
		while(input_angle>Math.PI){
			input_angle -= 2*Math.PI;
		}
		
		while(input_angle<=-Math.PI){
			input_angle += 2*Math.PI;
		}
		
		return input_angle;
	}

	// wraps it to (-180, 180]
	public static double angle_wrap_180(double input_angle){ 
		
		while(input_angle>180){
			input_angle -= 360;
		}
		
		while(input_angle<=-180){
			input_angle += 360;
		}
		
		return input_angle;
	}

	// wraps it to [0, 360)
	public static double angle_wrap_360(double input_angle){ 
		
		while(input_angle>=360){
			input_angle -= 360;
		}
		
		while(input_angle<0){
			input_angle += 360;
		}
		
		return input_angle;
	}

	// wraps it to [0, 2*pi)
	public static double angle_wrap_2_pi(double input_angle){ 
		
		while(input_angle>=2*Math.PI){
			input_angle -= 2*Math.PI;
		}
		
		while(input_angle<0){
			input_angle += 2*Math.PI;
		}
		
		return input_angle;
	}
	
	/* 
	 * spherical-cartesian system transformations theta = azimuth, phi = inclination, r = radius
	 */

	public static void 		cart2sph(double x_in, double y_in, double z_in, double[] r_phi_theta_out){
		
		if(r_phi_theta_out.length!=3){
			System.out.println("Transf:cart2sph():output length has to be 3.");
			System.exit(1);
		}
		// r
		r_phi_theta_out[0] = Math.sqrt(x_in*x_in+y_in*y_in+z_in*z_in); 
		// phi - polar (inclination) [0,pi]
		r_phi_theta_out[1] = z_in / r_phi_theta_out[0]; 
		r_phi_theta_out[1] = (Math.abs(r_phi_theta_out[1])>1.0)? 1.0 : r_phi_theta_out[1];
		r_phi_theta_out[1] = Math.acos(r_phi_theta_out[1]); 		// already wrapped
		// theta - azimuth
		r_phi_theta_out[2] = Math.atan2(y_in, x_in);  				// azimuth [-pi, pi]
		r_phi_theta_out[2] = angle_wrap_2_pi(r_phi_theta_out[2]); 	// azimuth [0, 2pi)
			
	}

	public static double 	cart2sph_r(double x_in, double y_in, double z_in) {
		double r = Math.sqrt(x_in * x_in + y_in * y_in + z_in * z_in);
		return r;
	}

	public static double 	cart2sph_phi(double x_in, double y_in, double z_in) {
		// phi - polar (inclination) [0,pi]
		double phi = z_in / Math.sqrt(x_in * x_in + y_in * y_in + z_in * z_in);
		phi = (Math.abs(phi) > 1.0) ? 1.0 : phi;
		phi = Math.acos(phi); // already wrapped
		return phi;
	}

	public static double 	cart2sph_theta(double x_in, double y_in, double z_in) {
		// theta - azimuth
		double theta = Math.atan2(y_in, x_in); // azimuth [-pi, pi]
		theta = angle_wrap_2_pi(theta); // azimuth [0, 2pi)
		return theta;
	}
	
	public static void 		sph2cart(double r_in, double phi_in, double theta_in,
			double[] x_y_z_out) {

		x_y_z_out[0] = r_in * Math.sin(phi_in) * Math.cos(theta_in); // x
		x_y_z_out[1] = r_in * Math.sin(phi_in) * Math.sin(theta_in); // y
		x_y_z_out[2] = r_in * Math.cos(phi_in);// z

	}

	public static double 	sph2cart_x(double r_in, double phi_in, double theta_in) {
		double x = r_in * Math.sin(phi_in) * Math.cos(theta_in);
		return x;
	}

	public static double 	sph2cart_y(double r_in, double phi_in, double theta_in) {
		double y = r_in * Math.sin(phi_in) * Math.sin(theta_in);
		return y;
	}

	public static double 	sph2cart_z(double r_in, double phi_in, double theta_in) {
		double z = r_in * Math.cos(phi_in);
		return z;
	}

	/*
	 * cylindrical<-->cartesian system transformations ro = radial distance, phi = azimuth, z = height
	 */

	public static void 			cart2cyl(double x_in, double y_in, double z_in, double[] ro_phi_z){
		ro_phi_z[0] = Math.sqrt(x_in*x_in + y_in*y_in);
		if(x_in==0 && y_in==0){
			ro_phi_z[1] = 0;
		}
		else if(x_in>=0){
			ro_phi_z[1] = Math.asin(y_in/ro_phi_z[0]);
		}
		else if(x_in<0){
			ro_phi_z[1] = -Math.asin(y_in/ro_phi_z[0])+Math.PI;
		}
		ro_phi_z[2] = z_in;
	}
	
	public static double 		cart2cyl_ro(double x_in, double y_in, double z_in){
		double ro = Math.sqrt(x_in*x_in + y_in*y_in);
		return ro;
	}
	
	public static double 		cart2cyl_phi(double x_in, double y_in, double z_in){
		double phi = 0;
		if(x_in==0 && y_in==0){
			phi = 0;
		}
		else if(x_in>=0){
			phi = Math.asin(y_in/Math.sqrt(x_in*x_in + y_in*y_in));
		}
		else if(x_in<0){
			phi = -Math.asin(y_in/Math.sqrt(x_in*x_in + y_in*y_in))+Math.PI;
		}
		return phi;
	}
	
	public static double 		cart2cyl_z(double x_in, double y_in, double z_in){	
		return z_in;
	}
	
	public static void 			cyl2cart(double ro, double phi, double z, double[] x_y_z){
		x_y_z[0] = ro*Math.cos(phi); // x
		x_y_z[1] = ro*Math.sin(phi); // y
		x_y_z[2] = z; // z
	}

	public static double 			cyl2cart_x(double ro, double phi, double z){
		return ro*Math.cos(phi); // x
	}
	
	public static double 			cyl2cart_y(double ro, double phi, double z){
		return ro*Math.sin(phi); // y
	}
	
	public static double 			cyl2cart_z(double ro, double phi, double z){
		return z; // z
	}
	
	public static Matrix sphChangeOrientation(double phi, double theta){//, double[] p_before, double[] p_after){
		
		Matrix rotTheta = new Matrix(
				new double[][]{
						{Math.cos(theta),	-Math.sin(theta),	0},
						{Math.sin(theta), 	 Math.cos(theta), 	0},
						{0,					0,					1}
				});
		
		Matrix rotPhi = new Matrix(
				new double[][]{
						{Math.cos(phi),		0,				Math.sin(phi)},
						{0, 	 			1, 				0},
						{-Math.sin(phi),	0,				Math.cos(phi)}
				});
		
		Matrix tf = rotTheta.times(rotPhi);
		
		return tf;
		
	}
	
	public static void sphChangeOrientation(double phi, double theta, double[] p_in, double[] p_out){
		
		Matrix tf = new Matrix(sphChangeOrientation(phi, theta));
		
		p_out = Matrix.rowMatrix(p_in).times(tf).getRowCopy(0); 
		
	}
	
}

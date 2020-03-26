package advantra.trace;

import advantra.shapes.Sphere;

public class DirectedPointsOnSemiSphere{
	
	int         d_point_number;
	double[][]  d_point_coordinate;
	double[][]  d_point_orientation;
	double[]  	d_point_cos_angle_wrt_axis; 
	Sphere		d_point_associated_sphere;
	
	
	public DirectedPointsOnSemiSphere(int Npoints, double[] center, double[] orientation, double radius){
		
		// doesn't make sense if the radius is <=1
		if(radius<=1){
			System.out.println("Radius of the semi-sphere cannot be <=1 in voxels.");
			System.exit(1);
		}
		
		// orientation cannot have all zeros
		if(isZeroVector(orientation, 0.0001)){
			System.out.printf("All elements in orientation vector are lower than %6.2f stopping...", 0.0001);
			System.exit(1);	
		}
		
		d_point_number = Npoints;
		
		// allocate
		d_point_coordinate			= new double[3][Npoints];
		d_point_orientation			= new double[3][Npoints];
		d_point_cos_angle_wrt_axis	= new double[Npoints];
		d_point_associated_sphere 	= new Sphere(center[0], center[1], center[2], radius);
		
		// fill the values in
		d_point_associated_sphere.generate3DSemiSpherePts(
				Npoints, 
				orientation[0], orientation[1], orientation[2], 
				d_point_coordinate, 
				d_point_orientation);
		
	}
	
	public void setPoints(double[] center, double[] orientation, double radius){
		
		// doesn't make sense if the radius is <=1
		if(radius<=1){
			System.out.println("Radius of the semi-sphere cannot be <=1 in voxels.");
			System.exit(1);
		}
		
		// orientation cannot have all zeros
		if(isZeroVector(orientation, 0.0001)){
			System.out.printf("All elements in orientation vector are lower than %6.2f stopping...", 0.0001);
			System.exit(1);	
		}
		
		// fill the values in
		d_point_associated_sphere.setSphere(center[0], center[1], center[2], radius);
		
		d_point_associated_sphere.generate3DSemiSpherePts(
				d_point_number, 
				orientation[0], orientation[1], orientation[2], 
				d_point_coordinate, 
				d_point_orientation);
		
	}
	
	public int getNpoints(){
		
		return d_point_number;
	}
	
	public boolean isZeroVector(double[] inputVec, double smallValue){
		boolean isZero = true;
		for (int i = 0; i < inputVec.length; i++) {
			if(inputVec[i]>smallValue){
				isZero = false;
				break;
			}
			
		}
		return isZero;
	}
	
	public double[] getPointCoordinate(int index){
		
		double[] pnt = new double[3];
		
		pnt[0] = d_point_coordinate[0][index];
		pnt[1] = d_point_coordinate[1][index];
		pnt[2] = d_point_coordinate[2][index];
		
		return pnt;
		
	}
	
	public double[] getPointOrientation(int index){
		
		double[] ort = new double[3];
		
		ort[0] = d_point_orientation[0][index];
		ort[1] = d_point_orientation[1][index];
		ort[2] = d_point_orientation[2][index];
		
		return ort;
		
	}

	public double getPointAngleDivergence(int index){
		
		return d_point_cos_angle_wrt_axis[index];
		
	}

	
	
}

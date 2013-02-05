package advantra.shapes;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import advantra.general.ArrayHandling;
import advantra.general.ColourTransf;
import advantra.general.ImageConversions;
import advantra.general.Transf;
import advantra.processing.IntensityCalc;

public class Cylinder extends RegionOfInterest {

	// designed for 3 dimensional cylinder data extracted from image	
	private double r;
	private double vx, vy, vz; //unit orientation vec
	private double h;
	
	public Cylinder(){
		// dummy constructor
		super(0,0,0);
		this.roi_type = RoiType.CYLINDER;
		this.r = 1;
		this.h = 1;
		this.vx = 1;
		this.vy = 0;
		this.vz = 0;
	}
	
	/*
	 * INITIALIZE CYLINDER
	 */
	
	public Cylinder(double x, double y, double z, double r, double h, double vx, double vy, double vz){
		super(x, y, z); // x,y,z will not be the center but the base point on one end cap, the other point will be determined with h*(vx, vy, vz)
		this.roi_type = RoiType.CYLINDER;
		this.r = r;
		this.vx = vx/Math.sqrt(vx*vx+vy*vy+vz*vz); // normalize
		this.vy = vy/Math.sqrt(vx*vx+vy*vy+vz*vz);
		this.vz = vz/Math.sqrt(vx*vx+vy*vy+vz*vz);
		this.h = h;
	}
	
	public Cylinder(double[] position_xyz, double r, double h, double[] orient_xyz){
		super(position_xyz[0], position_xyz[1], position_xyz[2]);
		this.roi_type = RoiType.CYLINDER;
		if((position_xyz.length!=3) || (orient_xyz.length!=3)){
			System.err.println("Cylinder:Cylinder(): \n" +
					"position_xyz and orient_xyz should have length 3!");
			System.exit(1);
		}
		this.r = r;
		this.h = h;
		double norm = Math.sqrt(Math.pow(orient_xyz[0],2)+Math.pow(orient_xyz[1],2)+Math.pow(orient_xyz[2],2));
		if(norm>0){
			this.vx = orient_xyz[0]/norm;
			this.vy = orient_xyz[1]/norm;
			this.vz = orient_xyz[2]/norm;
		}
		else{
			System.err.println("direction norm was not >0. Stoppping...");
			System.exit(1);
		}
		
	}
	
	/*
	 * SET CYLINDER
	 */
	
	public void setCylinder(double[] position_xyz, double r, double h, double[] orient_xyz){
		
		if((position_xyz.length!=3) || (orient_xyz.length!=3)){
			System.err.println("Cylinder:setCylinder(): \n" +
					"position_xyz and orient_xyz should have length 3!");
			System.exit(1);
		}
		
		this.x = position_xyz[0];
		this.y = position_xyz[1];
		this.z = position_xyz[2];
		
		this.r = r;
		this.h = h;
		
		this.vx = orient_xyz[0]/Math.sqrt(
				Math.pow(orient_xyz[0],2)+
				Math.pow(orient_xyz[1],2)+
				Math.pow(orient_xyz[2],2));
		this.vy = orient_xyz[1]/Math.sqrt(
				Math.pow(orient_xyz[0],2)+
				Math.pow(orient_xyz[1],2)+
				Math.pow(orient_xyz[2],2));
		this.vz = orient_xyz[2]/Math.sqrt(
				Math.pow(orient_xyz[0],2)+
				Math.pow(orient_xyz[1],2)+
				Math.pow(orient_xyz[2],2));
	}
		
	public void setCylinder(double x, double y, double z, double r, double h, double vx, double vy, double vz){
		this.x = x;
		this.y = y;
		this.z = z;
		this.r = r;
		this.h = h;
		this.vx = vx/Math.sqrt(vx*vx+vy*vy+vz*vz); // normalize
		this.vy = vy/Math.sqrt(vx*vx+vy*vy+vz*vz);
		this.vz = vz/Math.sqrt(vx*vx+vy*vy+vz*vz);
	}
	
	public void setCylinder(Cylinder input_cyl){
		this.x = input_cyl.x;
		this.y = input_cyl.y;
		this.z = input_cyl.z;
		
		this.r = input_cyl.r;
		this.h = input_cyl.h;		
		
		this.vx = input_cyl.vx;
		this.vy = input_cyl.vy;
		this.vz = input_cyl.vz;
		
	}
	
	/*
	 * RADIUS & HEIGHT
	 */

	public void setR(double r){
		this.r = r;
	}
	
	public double getR(){
		return this.r;
	}
	
	public void setH(double h){
		this.h = h;
	}
	
	public double getH(){
		return this.h;
	}
	
	/*
	 * ORIENTATION
	 */
	
	public void setV(double vx, double vy, double vz){
		this.vx = vx/Math.sqrt(vx*vx+vy*vy+vz*vz);
		this.vy = vy/Math.sqrt(vx*vx+vy*vy+vz*vz);
		this.vz = vz/Math.sqrt(vx*vx+vy*vy+vz*vz);
	}
	
	public void setV(double[] orientation){
		if(orientation.length!=3){
			System.err.println("Cylinder:setV(): \n" +
					"orientation vector has to have length 3.");
			System.exit(1);
		}
		this.vx = orientation[0]/Math.sqrt(
				Math.pow(orientation[0], 2)+
				Math.pow(orientation[1], 2)+
				Math.pow(orientation[2], 2));
		this.vy = orientation[1]/Math.sqrt(
				Math.pow(orientation[0], 2)+
				Math.pow(orientation[1], 2)+
				Math.pow(orientation[2], 2));
		this.vz = orientation[2]/Math.sqrt(
				Math.pow(orientation[0], 2)+
				Math.pow(orientation[1], 2)+
				Math.pow(orientation[2], 2));
	}
	
	public double[] getV(){
		double[] v_orient = new double[3];
		v_orient[0] = this.vx;
		v_orient[1] = this.vy;
		v_orient[2] = this.vz;
		return  v_orient;
	}

	public double getOrientX(){
		return this.vx;
	}

	public double getOrientY(){
		return this.vy;
	}

	public double getOrientZ(){
		return this.vz;
	}
	
	public double[] getPos(){
		double[] pos_xyz = new double[3];
		pos_xyz[0] = this.x;
		pos_xyz[1] = this.y;
		pos_xyz[2] = this.z;
		return  pos_xyz;
	}
	
	/*
	 * POSITION
	 */
	
	public void setPos(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public void setPos(double[] pos_xyz){
		if(pos_xyz.length!=3){
			System.err.println("Cylinder:setPos(): \n" +
					"position vector has to have length 3.");
			System.exit(1);
		}
		this.x = pos_xyz[0];
		this.y = pos_xyz[1];
		this.z = pos_xyz[2];
	}
	
	public void translatePos(double shift){
		this.x += shift;
		this.y += shift;
		this.z += shift;
	}
	
	/*	##################################
	 *  METHODS
	 * 	##################################
	 */	
	
	public int extractVox(
			final int[][] image_stack, 
			final int height, 
			final int width, 
			final int len, 
			final int[][] 		roi_coord, 
			final int[] 		roi_vals, 
			final double[][] 	roi_x_y_crosssec,
			final double[]		radial_dist
			){ 
		
		// TODO: try voxel sampling - maybe it is more efficient instead of taking all the pixels

		final int MAX_R = 20;
		final int MIN_H = 2;
		
		int dist = (int) Math.ceil(Math.sqrt(Math.pow(h/2, 2)+Math.pow(r, 2)));
		// TODO: add more conditions maybe
		if(
				roi_coord[0].length<Sphere.numberOfVoxInSphere(dist)
				){
			System.err.println("Cylinder:extractVox(): \n" +
					"'roi_coord' or 'roi_vals' matrices needs to have "+Sphere.numberOfVoxInSphere(dist)+"+ columns\n" +
							" to cover the cylinder with R="+r+" and H="+h+" "+
					"\nallocate some more space!");
			System.exit(1);
			//return;
		}

		if(
				roi_coord[0].length != roi_vals.length 				|| 
				roi_coord[0].length != roi_x_y_crosssec[0].length 	||
				roi_coord[0].length != radial_dist.length
				){
			System.err.println("Cylinder:extractVox(): \n" +
					"output matrices have irregular dimensions. They should have the same number of columns! \n" +
					"check allocated rows/columns!");
			System.exit(1);	
		}
		
		if(roi_coord.length!=3){
			System.err.println("Cylinder:extractVox(): \n'roi_coord' matrix needs to have 3 rows!");
			System.exit(1);			
		}
		
		if(roi_x_y_crosssec.length!=2){
			System.err.println("Cylinder:extractVox(): \n'roi_x_y_crosssec' matrix needs to have 2 rows!");
			System.exit(1);			
		}
		
		if(r>MAX_R){
			System.err.println("Cylinder:extractVox(): \nradius is at most "+MAX_R+" pixels...");
			System.exit(1);
		}
		
		if(h<MIN_H){
			System.err.println("Cylinder:extractVox(): \nH needs to be at least "+MIN_H+" pixels...");
			System.exit(1);
		}
		
		// center in pixel coordinates
		int[] basePoint = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range so that it doesn't collide with borders
		int startX = 0; if (basePoint[0]-dist > 0) startX = basePoint[0]-dist;
		int startY = 0; if (basePoint[1]-dist > 0) startY = basePoint[1]-dist;
		int startZ = 0; if (basePoint[2]-dist > 0) startZ = basePoint[2]-dist;
			
		int endX   = height;if (basePoint[0]+dist < height) 	endX = basePoint[0]+dist;
		int endY   = width; if (basePoint[1]+dist < width) 		endY = basePoint[1]+dist;
		int endZ   = len;  	if (basePoint[2]+dist < len) 		endZ = basePoint[2]+dist;

		int count=0;
		
		// extract the unit vects orthogonal to the direction "orient": "orient_x" and "orient_y"
		double[] orient_a = new double[3];
		double[] orient_b = new double[3];
		
		Transf.cartesian(vx, vy, vz, orient_a, orient_b); // calculate them
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					
					// x, y, z will be the coordinate, this.x, this.y, this.z is the center
					
					double pd_x = x - this.x;
					double pd_y = y - this.y;
					double pd_z = z - this.z;
					
//					double d_x = h * vx;
//					double d_y = h * vy;
//					double d_z = h * vz;
					
					double dot_prod = pd_x * vx + pd_y * vy + pd_z * vz;
					
					if(dot_prod < -0.5*h || dot_prod > 0.5*h){
						// don't add it at all
					}
					else
					{
						double dist_sq = (pd_x * pd_x + pd_y * pd_y + pd_z * pd_z) - (dot_prod * dot_prod);// squared distance from the cylinder axis
						if(dist_sq > r * r){
							
							// don't add it at all
							
						}
						else
						{
							// add it - it is in the cylinder
							// extract the positions
							roi_coord[0][count] = x;
							roi_coord[1][count] = y;
							roi_coord[2][count] = z;
							
							// extract the values
							roi_vals[count] 	= image_stack[z][y + width * x];
							
							// 2d cross section coords + distances
							// each will be dot product (cartesian coordinate) of the vector having cylinder center as origin with orientation orthogonals unit vectors
							// vector having cylinder center as origin: [ x-this.x, y-this.y, z-this.z ]
							// orthogonal unit vectors: orient_a and orient_b where angle(orient_a, v)=90, angle(orient_b, v)=90 and angle(orient_a, orient_b)=90 
							roi_x_y_crosssec[0][count] = (x-this.x)*orient_a[0] + (y-this.y)*orient_a[1] + (z-this.z)*orient_a[2] ; 
							roi_x_y_crosssec[1][count] = (x-this.x)*orient_b[0] + (y-this.y)*orient_b[1] + (z-this.z)*orient_b[2] ;
							
							// radial distance from the cross section center (in pixels)
							radial_dist[count] = 
									Math.sqrt(
											roi_x_y_crosssec[0][count]*roi_x_y_crosssec[0][count] + 
											roi_x_y_crosssec[1][count]*roi_x_y_crosssec[1][count]
													);
							
							count ++;
						}
					}
				}
			}
		}	

		return count;
		
	}
	
	public int extractVox(
			final ImagePlus 	input_image, 
			final int[][] 		roi_coord, 
			final int[] 		roi_vals, 
			final double[][] 	roi_x_y_crosssec 
			){ 
		
		// TODO: try voxel sampling - maybe it is more efficient instead of taking all the pixels

		final int MAX_R = 20;
		final int MIN_H = 2;
		
		int dist = (int) Math.ceil(Math.sqrt(Math.pow(h, 2)+Math.pow(r, 2)));
		// TODO: add more conditions maybe
		if(
				roi_coord[0].length<Sphere.numberOfVoxInSphere(dist)  || 
				roi_vals.length<Sphere.numberOfVoxInSphere(dist) ||
				roi_x_y_crosssec[0].length<Sphere.numberOfVoxInSphere(dist)
				){
			System.err.println("Cylinder:extractVox(): \n" +
					"matrices need to have "+Sphere.numberOfVoxInSphere(dist)+" & more columns\n" +
					" to cover the cylinder with R="+r+" and H="+h+" \n"+
					"allocate some more space!");
			System.exit(1);
			
		}

		if(
				roi_coord[0].length != roi_vals.length || 
				roi_vals.length 	!= roi_x_y_crosssec[0].length ||
				roi_coord[0].length != roi_x_y_crosssec[0].length){
			System.err.println("Cylinder:extractVox(): \n" +
					"matrices have irregular dimensions. They should have the same number of columns! \n" +
					"check allocated rows/columns!");
			System.exit(1);	
		}
		
		if(roi_coord.length!=3){
			System.err.println("Cylinder:extractVox(): \n" +
					"coordinates matrix needs to have 3 rows!");
			System.exit(1);			
		}
		
		if(roi_x_y_crosssec.length!=2){
			System.err.println("Cylinder:extractVox(): \n'roi_x_y_crosssec' matrix needs to have 2 rows!");
			System.exit(1);			
		}
		
		if(r>MAX_R){
			System.err.println("Cylinder:extractVox(): \nradius is at most "+MAX_R+" pixels...");
			System.exit(1);
		}
		
		if(h<MIN_H){
			System.err.println("Cylinder:extractVox(): \nh needs to be at least "+MIN_H+" pixels...");
			System.exit(1);
		}
		
		// center in pixel coordinates
		int[] basePoint = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range so that it doesn't collide with borders
		int startX = 0; if (basePoint[0]-dist > 0) startX = basePoint[0]-dist;
		int startY = 0; if (basePoint[1]-dist > 0) startY = basePoint[1]-dist;
		int startZ = 0; if (basePoint[2]-dist > 0) startZ = basePoint[2]-dist;
			
		int height = input_image.getStack().getHeight();
		int width  = input_image.getStack().getWidth();
		int len    = input_image.getStack().getSize();
		
		int endX   = height;if (basePoint[0]+dist < height) 	endX = basePoint[0]+dist;
		int endY   = width; if (basePoint[1]+dist < width) 		endY = basePoint[1]+dist;
		int endZ   = len;  	if (basePoint[2]+dist < len) 		endZ = basePoint[2]+dist;

		int count=0;
		
		// extract the unit vects orthogonal to the direction "orient": "orient_x" and "orient_y"
		double[] orient_a = new double[3];
		double[] orient_b = new double[3];
		
		Transf.cartesian(vx, vy, vz, orient_a, orient_b); // calculate them
		
		// array will be necessary to operate with the values 
		int[][] image_stack = ImageConversions.GraytoIntArray(input_image);
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					
					// x, y, z will be the coordinate, this.x, this.y, this.z is the center
					
					double pd_x = x - this.x;
					double pd_y = y - this.y;
					double pd_z = z - this.z;
					
					double d_x = h * vx;
					double d_y = h * vy;
					double d_z = h * vz;
					
					double dot_prod = pd_x * d_x + pd_y * d_y + pd_z * d_z;
					
					if(dot_prod < -0.5*h*h || dot_prod > 0.5*h*h){
						// don't add it at all
					}
					else
					{
						double dist_sq = (pd_x * pd_x + pd_y * pd_y + pd_z * pd_z) - (dot_prod * dot_prod) / (h * h);// squared distance from the cylinder axis
						if(dist_sq > r * r){
							// don't add it at all
							
						}
						else
						{
							// add it - it is in the cylinder
							// extract the positions
							roi_coord[0][count] = x;
							roi_coord[1][count] = y;
							roi_coord[2][count] = z;
							
							// extract the values
							roi_vals[count] 	= image_stack[z][y + width * x];
							
							// 2d cross section coords + distances
							// each will be dot product (cartesian coordinate) of the vector having cylinder center as origin with orientation orthogonals unit vectors
							// vector having cylinder center as origin: [ x-this.x, y-this.y, z-this.z ]
							// orthogonal unit vectors: orient_a and orient_b where angle(orient_a, v)=90, angle(orient_b, v)=90 and angle(orient_a, orient_b)=90 
							roi_x_y_crosssec[0][count] = (x-this.x)*orient_a[0] + (y-this.y)*orient_a[1] + (z-this.z)*orient_a[2] ; 
							roi_x_y_crosssec[1][count] = (x-this.x)*orient_b[0] + (y-this.y)*orient_b[1] + (z-this.z)*orient_b[2] ;
							
							// distance from the cross section
							//roi_dists[count] = Math.sqrt(roi_x_y_crosssec[0][count]*roi_x_y_crosssec[0][count] + roi_x_y_crosssec[1][count]*roi_x_y_crosssec[1][count]);
							
							count ++;
						}
					}
				}
			}
		}	

		return count;
		
	}

	public int extractVox(
			final ImageStack 	input_image_stack, 
			final int[][] 		roi_coord, 
			final int[] 		roi_vals, 
			final double[][] 	roi_x_y_crosssec 
			){ 
		
		final int MAX_R = 20;
		final int MIN_H = 2;
		
		int dist = (int) Math.ceil(Math.sqrt(Math.pow(h/2, 2)+Math.pow(r, 2)));
		if(
				roi_coord[0].length<Sphere.numberOfVoxInSphere(dist)  || 
				roi_vals.length<Sphere.numberOfVoxInSphere(dist) ||
				roi_x_y_crosssec[0].length<Sphere.numberOfVoxInSphere(dist)
				){
			System.err.println("Cylinder:extractVox(): \n" +
					"matrices need to have "+Sphere.numberOfVoxInSphere(dist)+" & more columns\n" +
					" to cover the cylinder with R="+r+" and H="+h+" \n"+
					"allocate some more space!");
			System.exit(1);
			
		}

		if(
				roi_coord[0].length != roi_vals.length || 
				roi_vals.length 	!= roi_x_y_crosssec[0].length ||
				roi_coord[0].length != roi_x_y_crosssec[0].length){
			System.err.println("Cylinder:extractVox(): \n" +
					"matrices have irregular dimensions. They should have the same number of columns! \n" +
					"check allocated rows/columns!");
			System.exit(1);	
		}
		
		if(roi_coord.length!=3){
			System.err.println("Cylinder:extractVox(): \n" +
					"coordinates matrix needs to have 3 rows!");
			System.exit(1);			
		}
		
		if(roi_x_y_crosssec.length!=2){
			System.err.println("Cylinder:extractVox(): \n'roi_x_y_crosssec' matrix needs to have 2 rows!");
			System.exit(1);			
		}
		
		if(r>MAX_R){
			System.err.println("Cylinder:extractVox(): \nradius is at most "+MAX_R+" pixels...");
			System.exit(1);
		}
		
		if(h<MIN_H){
			System.err.println("Cylinder:extractVox(): \nh needs to be at least "+MIN_H+" pixels...");
			System.exit(1);
		}
		
		// center in pixel coordinates
		int[] basePoint = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range so that it doesn't collide with borders
		int startX = 0; if (basePoint[0]-dist > 0) startX = basePoint[0]-dist;
		int startY = 0; if (basePoint[1]-dist > 0) startY = basePoint[1]-dist;
		int startZ = 0; if (basePoint[2]-dist > 0) startZ = basePoint[2]-dist;
			
		int height = input_image_stack.getHeight();
		int width  = input_image_stack.getWidth();
		int len    = input_image_stack.getSize();
		
		int endX   = height;if (basePoint[0]+dist < height) 	endX = basePoint[0]+dist;
		int endY   = width; if (basePoint[1]+dist < width) 		endY = basePoint[1]+dist;
		int endZ   = len;  	if (basePoint[2]+dist < len) 		endZ = basePoint[2]+dist;

		int count=0;
		
		// extract the unit vects orthogonal to the direction "orient": "orient_x" and "orient_y"
		double[] orient_a = new double[3];
		double[] orient_b = new double[3];
		
		Transf.cartesian(vx, vy, vz, orient_a, orient_b); // calculate them
		
		// array will be necessary to operate with the values 
		int[][] image_stack = ImageConversions.GraytoIntArray(input_image_stack);
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					
					// x, y, z will be the coordinate, this.x, this.y, this.z is the center
					
					double pd_x = x - this.x;
					double pd_y = y - this.y;
					double pd_z = z - this.z;
					
					double d_x = h * vx;
					double d_y = h * vy;
					double d_z = h * vz;
					
					double dot_prod = pd_x * d_x + pd_y * d_y + pd_z * d_z;
					
					if(dot_prod < -0.5*h*h || dot_prod > 0.5*h*h){
						// don't add it at all
					}
					else
					{
						double dist_sq = (pd_x * pd_x + pd_y * pd_y + pd_z * pd_z) - (dot_prod * dot_prod) / (h * h);// squared distance from the cylinder axis
						if(dist_sq > r * r){
							// don't add it at all
							
						}
						else
						{
							// add it - it is in the cylinder
							// extract the positions
							roi_coord[0][count] = x;
							roi_coord[1][count] = y;
							roi_coord[2][count] = z;
							
							// extract the values
							roi_vals[count] 	= image_stack[z][y + width * x];
							
							// 2d cross section coords + distances
							// each will be dot product (cartesian coordinate) of the vector having cylinder center as origin with orientation orthogonals unit vectors
							// vector having cylinder center as origin: [ x-this.x, y-this.y, z-this.z ]
							// orthogonal unit vectors: orient_a and orient_b where angle(orient_a, v)=90, angle(orient_b, v)=90 and angle(orient_a, orient_b)=90 
							roi_x_y_crosssec[0][count] = (x-this.x)*orient_a[0] + (y-this.y)*orient_a[1] + (z-this.z)*orient_a[2] ; 
							roi_x_y_crosssec[1][count] = (x-this.x)*orient_b[0] + (y-this.y)*orient_b[1] + (z-this.z)*orient_b[2] ;
							
							// distance from the cross section
							//roi_dists[count] = Math.sqrt(roi_x_y_crosssec[0][count]*roi_x_y_crosssec[0][count] + roi_x_y_crosssec[1][count]*roi_x_y_crosssec[1][count]);
							
							count ++;
						}
					}
				}
			}
		}	

		return count;
		
	}
	
	public int extractCoords(
			final ImagePlus input_image,
			final int[][] cyl_coord
			){
		
		// coordinates of voxels that make the cylinder
		
		final int MAX_R = 20;
		final int MIN_H = 2;
		
		int dist = (int) Math.ceil(Math.sqrt(Math.pow(h/2, 2)+Math.pow(r, 2)));
		// TODO: add more conditions maybe
		if(
				cyl_coord[0].length<Sphere.numberOfVoxInSphere(dist) ){
			System.err.println("Cylinder:extractCoords(): \n" +
					"output matrix needs to have "+Sphere.numberOfVoxInSphere(dist)+"+ & more columns\n" +
					" to cover the cylinder with R="+r+" and H="+h+" \n"+
					"allocate some more space!");
			System.exit(1);
		}

		if(cyl_coord.length!=3){
			System.err.println("Cylinder:extractCoords(): \n" +
					"coordinates matrix needs to have 3 rows!");
			System.exit(1);			
		}
		
		if(r>MAX_R){
			System.err.println("Cylinder:extractCoords(): \nradius is at most "+MAX_R+" pixels...");
			System.exit(1);
		}
		
		if(h<MIN_H){
			System.err.println("Cylinder:extractCoords(): \nh needs to be at least "+MIN_H+" pixels...");
			System.exit(1);
		}
		
		// center in pixel coordinates
		int[] basePoint = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range so that it doesn't collide with borders
		int startX = 0; if (basePoint[0]-dist > 0) startX = basePoint[0]-dist;
		int startY = 0; if (basePoint[1]-dist > 0) startY = basePoint[1]-dist;
		int startZ = 0; if (basePoint[2]-dist > 0) startZ = basePoint[2]-dist;
			
		int height = input_image.getStack().getHeight();
		int width  = input_image.getStack().getWidth();
		int len    = input_image.getStack().getSize();
		
		int endX   = height;if (basePoint[0]+dist < height) 	endX = basePoint[0]+dist;
		int endY   = width; if (basePoint[1]+dist < width) 		endY = basePoint[1]+dist;
		int endZ   = len;  	if (basePoint[2]+dist < len) 		endZ = basePoint[2]+dist;

		int count=0;
		
		// extract the unit vects orthogonal to the direction "orient": "orient_x" and "orient_y"
		double[] orient_a = new double[3];
		double[] orient_b = new double[3];
		
		Transf.cartesian(vx, vy, vz, orient_a, orient_b); // calculate them
		
		//int[][] image_stack = ArrayHandling.image_plus_stack_to_int_array(input_image); // array will be necessary to take the values out
		
		for(int x=startX; x<endX; x++){
			for(int y=startY; y<endY; y++){
				for(int z=startZ; z<endZ; z++){
					
					// x, y, z will be the coordinate, this.x, this.y, this.z is the center
					
					double pd_x = x - this.x;
					double pd_y = y - this.y;
					double pd_z = z - this.z;
					
//					double d_x = h * vx;
//					double d_y = h * vy;
//					double d_z = h * vz;
					
					double dot_prod = pd_x * vx + pd_y * vy + pd_z * vz;
					
					if(dot_prod < -0.5*h || dot_prod > 0.5*h){
						// don't add it at all
					}
					else
					{
						double dist_sq = (pd_x * pd_x + pd_y * pd_y + pd_z * pd_z) - (dot_prod * dot_prod);// squared distance from the cylinder axis
						if(dist_sq > r * r){
							// don't add it at all
							
						}
						else
						{
							// add it - it is in the cylinder
							// extract the positions
							cyl_coord[0][count] = x;
							cyl_coord[1][count] = y;
							cyl_coord[2][count] = z;
							
							// extract the values
							//roi_vals[count] 	= image_stack[z][y + width * x];
							
							// 2d cross section coords + distances
							// each will be dot product (cartesian coordinate) of the vector having cylinder center as origin with orientation orthogonals unit vectors
							// vector having cylinder center as origin: [ x-this.x, y-this.y, z-this.z ]
							// orthogonal unit vectors: orient_a and orient_b where angle(orient_a, v)=90, angle(orient_b, v)=90 and angle(orient_a, orient_b)=90 
							//roi_x_y_crosssec[0][count] = (x-this.x)*orient_a[0] + (y-this.y)*orient_a[1] + (z-this.z)*orient_a[2] ; 
							//roi_x_y_crosssec[1][count] = (x-this.x)*orient_b[0] + (y-this.y)*orient_b[1] + (z-this.z)*orient_b[2] ;
							
							// distance from the cross section
							//roi_dists[count] = Math.sqrt(roi_x_y_crosssec[0][count]*roi_x_y_crosssec[0][count] + roi_x_y_crosssec[1][count]*roi_x_y_crosssec[1][count]);
							
							count ++;
						}
					}
				}
			}
		}	

		return count;
	}

	/*	##################################
	 *  TRANSFORMATIONS
	 * 	##################################
	 */	
	
	// localCart 	means 2 values for position within the cylinder's cross-section and one for the height
	// localCyl  	means cylinder's ro, phi, z 
	// globalCart 	means coordinates within the image
	
	public double[] 	cylindrical2cartesian(double ro, double phi, double z){
		double[] x_y_z = new double[3];
		Transf.cyl2cart(ro, phi, z, x_y_z);
		// local cartesian 2 global cartesian
		Transf.rotate(vx, vy, vz, x_y_z);
		x_y_z[0] += this.x;
		x_y_z[1] += this.y;
		x_y_z[2] += this.z;
		return x_y_z;
	}
	
	public double[][] 	cylindrical2cartesian(double[][] ro_phi_z){
		double[][] x_y_z = new double[ro_phi_z.length][3];
		for (int i = 0; i < ro_phi_z.length; i++) {
			Transf.cyl2cart(ro_phi_z[i][0], ro_phi_z[i][1], ro_phi_z[i][2], x_y_z[i]);
		}
		// local cartesian 2 global cartesian
		Transf.rotateNx3(vx, vy, vz, x_y_z);
		for (int i = 0; i < ro_phi_z.length; i++) {
			x_y_z[i][0] += this.x;
			x_y_z[i][1] += this.y;
			x_y_z[i][2] += this.z;
		}
		return x_y_z;
	}
	
	public double[] 	cartesian2cylindrical(double x, double y, double z){
		double[] ro_phi_z = new double[3];
		Transf.cart2cyl(x-this.x, y-this.y, z-this.z, ro_phi_z);
		return ro_phi_z;
	}

	public double[][] 	cartesian2cylindrical(double[][] x_y_z){
		
		double[][] ro_phi_z = new double[x_y_z.length][3];
		for (int i = 0; i < x_y_z.length; i++) {
			Transf.cart2cyl(x_y_z[i][0]-this.x, x_y_z[i][1]-this.y, x_y_z[i][2]-this.z, ro_phi_z[i]);
		}
		return ro_phi_z;
	}
	
	public double[][] 	extractCylinderCoordinates(double dr, double dh){ // Nx3
		
		double[][] x_y_z = new double[numberOfCylinderPoints(dr, dh)][3];
		
		int count_in_cyl = 0;
		
		for (double ro = 0.0; ro <= this.r; ro+=dr) {
			
			int M = 1;
			if(ro>0){
				M = (int)Math.ceil((2*Math.PI)*(ro/dr));
			}
			
			for (int m = 1; m <= M; m++) {
				double phi = m*((2*Math.PI)/M);
				for (double z = -this.h/2; z <= this.h/2; z+=dh) {
					x_y_z[count_in_cyl] = cylindrical2cartesian(ro, phi, z);
					count_in_cyl++;
					
				}
			}
			
		}
		
		return x_y_z;
		
	}
	
	public double[][] 	extractCylinderCoordinatesAndValues(ImagePlus input_img, double dr, double dh){ // Nx3
		
		double[][] x_y_z_I = new double[numberOfCylinderPoints(dr, dh)][4]; // x(coord), y(coord), z(coord), I (intensity)
		
		IntensityCalc calc = new IntensityCalc(input_img.getStack());
		
		int count_in_cyl = 0;
		
		for (double ro = 0.0; ro <= this.r; ro+=dr) {
			
			int M = 1;
			if(ro>0){
				M = (int)Math.ceil((2*Math.PI)*(ro/dr));
			}
			
			for (int m = 1; m <= M; m++) {
				double phi = m*((2*Math.PI)/M);
				for (double z = -this.h/2; z <= this.h/2; z+=dh) {
					
					double[] loc 	= cylindrical2cartesian(ro, phi, z);
					float I 		= calc.interpolateAt_new((float)loc[0], (float)loc[1], (float)loc[2]);
					x_y_z_I[count_in_cyl][0] = loc[0];
					x_y_z_I[count_in_cyl][1] = loc[1];
					x_y_z_I[count_in_cyl][2] = loc[2];
					x_y_z_I[count_in_cyl][3] = I;
					count_in_cyl++;
					
				}
			}
			
		}
		
		return x_y_z_I;
		
	}
	
	public double[]		extractAvgInOut_new(IntensityCalc image_calc, double in_out_ratio, double dr, double dh){
		
		int 	cnt_in	= 0;
		int 	cnt_out	= 0;
		double 	avg_in 	= 0;
		double  avg_out	= 0;
		
		int extracted_samples = 0;

		for (double ro = 0.0; ro <= this.r; ro+=dr) {
			
			int M = 1;
			if(ro>0){
				M = (int)Math.ceil((2*Math.PI)*(ro/dr));
			}
			
			for (int m = 1; m <= M; m++) {
				
				double phi = m*((2*Math.PI)/M);
				
				for (double z = -this.h/2; z <= this.h/2; z+=dh) {
					
					double[] loc 		= cylindrical2cartesian(ro, phi, z);
					
					boolean isInImage 	= 
							(loc[0]>=0 && loc[0]<=image_calc.getImgHeight()	-1) 	&&
							(loc[1]>=0 && loc[1]<=image_calc.getImgWidth()	-1) 	&&
							(loc[2]>=0 && loc[2]<=image_calc.getImgLength()	-1);
					
					if(isInImage){
						
						extracted_samples++;
						
						float I 	= image_calc.interpolateAt_new((float)loc[0], (float)loc[1], (float)loc[2]); 
						
						if(ro<=in_out_ratio*this.r){
							
							avg_in += I;
							cnt_in ++;
							
						}
						else{
							
							avg_out += I;
							cnt_out++;
							
						}
						
					}
					
				}
				
			}
			
		}
		
		if(cnt_in>0){
			avg_in /= cnt_in;
		}
		
		if(cnt_out>0){
			avg_out /= cnt_out;
		}
		
		return new double[]{avg_in, avg_out, extracted_samples};
		
	}

	// TODO: should be deprecated in future in favor of extractAvgInOut_new 
	public double[] 	extractAvgInOut(ImagePlus input_img, double in_out_ratio, double dr, double dh){ 
		
		int 	cnt_in	= 0;
		int 	cnt_out	= 0;
		double 	avg_in 	= 0;
		double  avg_out	= 0;
		
		int extracted_samples = 0;
		
		IntensityCalc calc = new IntensityCalc(input_img.getStack());
		
		for (double ro = 0.0; ro <= this.r; ro+=dr) {
			
			int M = 1;
			if(ro>0){
				M = (int)Math.ceil((2*Math.PI)*(ro/dr));
			}
			
			for (int m = 1; m <= M; m++) {
				
				double phi = m*((2*Math.PI)/M);
				
				for (double z = -this.h/2; z <= this.h/2; z+=dh) {
					
					double[] loc 		= cylindrical2cartesian(ro, phi, z);
					boolean isInImage 	= 
							(loc[0]>=0 && loc[0]<=input_img.getStack().getHeight()-1) 	&&
							(loc[1]>=0 && loc[1]<=input_img.getStack().getWidth()-1) 	&&
							(loc[2]>=0 && loc[2]<=input_img.getStack().getSize()-1);
					
					if(isInImage){
						
						extracted_samples++;
						
						float I 		= calc.interpolateAt_new((float)loc[0], (float)loc[1], (float)loc[2]); 
						
						if(ro<=in_out_ratio*this.r){
							
							avg_in += I;
							cnt_in ++;
							
						}
						else{
							
							avg_out += I;
							cnt_out++;
							
						}
						
					}
					
				}
				
			}
			
		}
		
		if(cnt_in>0){
			avg_in /= cnt_in;
		}
		
		if(cnt_out>0){
			avg_out /= cnt_out;
		}
		
		return new double[]{avg_in, avg_out, extracted_samples};
		
	}
	
	public int numberOfCylinderPoints(double dr, double dh){
		int count_in_cyl = 0;
		for (double ro = 0.0; ro <= this.r; ro+=dr) {
			for (double z = -this.h/2; z <= this.h/2; z+=dh) {
				int M = 1;
				if(ro>0){
					M = (int)Math.ceil((2*Math.PI)*(ro/dr));
				}
				count_in_cyl += M;
			}
		}
		return count_in_cyl;
	}	

	/*	##################################
	 *  GRAPHICS
	 * 	##################################
	 */	
	
	// draws cylinder with specified grey-level value
	// if it is gray8  - valueToDraw 0-255
	// if it is gray16 - valueToDraw 0-65,535
	// use ImageConversions.setJet256ColorValue to draw over color images
	public ImagePlus drawOverGrayImage(ImagePlus template_image, int valueToDraw){
		
		int wt 	= template_image.getStack().getWidth();
		
		// check if it's grey image
		if(template_image.getType()==ImagePlus.GRAY8){
			valueToDraw = (valueToDraw>255) ? 255 : valueToDraw ;
			valueToDraw = (valueToDraw<  0) ?   0 : valueToDraw ;
		}
		else if(template_image.getType()==ImagePlus.GRAY16){
			valueToDraw = (valueToDraw>65535) ? 65535 : valueToDraw ;
			valueToDraw = (valueToDraw<  0  ) ?   0   : valueToDraw ;
		}
		else{
			System.err.println("Sphere:drawOverImage(): this image format is not suported - possible to use gray8, gray16 or rgb!");
			System.exit(-1);
		}
		
		int[][] template_image_array = ImageConversions.GraytoIntArray(template_image); 
		
		// extract the points
		int[][] cyl_coordinates 	= 
				new int[3][Sphere.numberOfVoxInSphere((int) Math.ceil(r*r+h*h))];
		
		int cnt = this.extractCoords(template_image, cyl_coordinates);
		
		if(cnt>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt; i++) {
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				template_image_array[layer][col+row*wt] = valueToDraw;
				
			}
			
		}else{
			
			System.out.println("Cylinder:drawOverImage()\nThere was not enough points to draw anything!");
		
		}

		return ImageConversions.toGray8(template_image_array, wt);     //toByteImage(template_image_array, wt);
	
	}
	
	public void 	drawOverColorImage(ImagePlus template_image, int valueToDraw){
		
		int w = template_image.getWidth();
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}

		byte[] rgb_values = ColourTransf.Jet256(valueToDraw);
		
		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
		// extract the points
		int[][] cyl_coordinates = new int[3][Sphere.numberOfVoxInSphere((int) Math.ceil(r*r+h*h))];
		
		int cnt = this.extractCoords(template_image, cyl_coordinates);
		
		if(cnt>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt; i++) {
				
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				
				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[0];
				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[1];
				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[2];
				
				//template_image_array[layer][col+row*wt] = valueToDraw;
				
			}
			
		}else{
			
			System.out.println("Cylinder:drawOverImage()\nThere was not enough points to draw anything!");
		
		}
		
		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}

	}

	public void 	drawOverColorImage(ImagePlus template_image, int r, int g, int b){
		
		int w = template_image.getWidth();
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}

		byte[] rgb_values = new byte[]{(byte)r, (byte)g, (byte)b};//ColourTransf.Jet256(valueToDraw);
		
		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
		// extract the points
		int[][] cyl_coordinates = new int[3][Sphere.numberOfVoxInSphere((int) Math.ceil(r*r+h*h))];
		
		int cnt = this.extractCoords(template_image, cyl_coordinates);
		
		if(cnt>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt; i++) {
				
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				
				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[0];
				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[1];
				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[2];
				
			}
			
		}else{
			
			System.out.println("Cylinder:drawOverImage()\nThere was not enough points to draw anything!");
		
		}
		
		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}

	}
	
	public void 	drawRedOverColorImage(ImagePlus template_image){
		
		int w = template_image.getWidth();
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}

		byte[] rgb_value_red = ColourTransf.Red();   //Jet256(valueToDraw);
		
		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
		// extract the points
		int[][] cyl_coordinates = new int[3][Sphere.numberOfVoxInSphere((int) Math.ceil(Math.pow(r, 2)+Math.pow(h/2, 2)))];
		
		int cnt = this.extractCoords(template_image, cyl_coordinates);
		
		if(cnt>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cnt; i++) {
				
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				
				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_value_red[0];
				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_value_red[1];
				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_value_red[2];
				
			}
			
		}else{
			
			System.out.println("Cylinder:drawOverImage()\nThere was not enough points to draw anything!");
		
		}
		
		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}

	}
	
	public static void drawOverColorImage(Cylinder[] cylinders_to_draw, ImagePlus template_image, int valueToDraw){
		for (int i = 0; i < cylinders_to_draw.length; i++) {
			cylinders_to_draw[i].drawOverColorImage(template_image, valueToDraw);
		}
	}
}

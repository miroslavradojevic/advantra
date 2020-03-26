package advantra.shapes;

import ij.ImageStack;
import advantra.general.Transf;


public class Cone extends RegionOfInterest {
	// 3 dimensional cone extracted from image	
	private double r;//cone basis radius
	// cone orientation
	private double ax, ay, az;
	// cone height
	private double h;
	
	public Cone(double x, double y, double z, double r, double ax, double ay, double az, double h){
		super(x, y, z);
		this.roi_type = RoiType.CONE;
		this.r = r;
		this.ax = ax;
		this.ay = ay;
		this.az = az;
		this.h  = h;
	}

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
	
	public void setOrientation(double ax, double ay, double az){
		this.ax = ax;
		this.ay = ay;
		this.az = az;
	}
	
	public double getOrientationX(){
		return this.ax;
	}
	public double getOrientationY(){
		return this.ay;
	}
	public double getOrientationZ(){
		return this.az;
	}	
	
	
	public int[][] extractVoxels(ImageStack stack){ 
		// extracts Sphere region of interest from particular image stack
		// sphere defines center and radius
		// output is nx4 integer array where columns 1-3 contain voxel indexes and column 4 voxel intensity
		int stack_width  = stack.getWidth();
		int stack_height = stack.getHeight();
		int stack_length = stack.getSize();

		// method will extract cube from the stack
		// around the cube center, with cube's length len
		
		int cone_radius 	= (int)Math.round(this.r);
		int cone_height 	= (int)Math.round(this.h);
		int cone_limit;
		if(cone_radius>cone_height){
			cone_limit = cone_radius;
		}
		else{
			cone_limit = cone_height;
		}
		// center in pixel coordinates
		int[] centerPix = {(int)Math.round(this.x), (int)Math.round(this.y), (int)Math.round(this.z)};
		
		// define the range
		int startX = 1; if (centerPix[0]-cone_limit > 1) startX = centerPix[0]-cone_limit;
		int startY = 1; if (centerPix[1]-cone_limit > 1) startY = centerPix[1]-cone_limit;
		int startZ = 1; if (centerPix[2]-cone_limit > 1) startZ = centerPix[2]-cone_limit;
			
		int endX   = stack_width;   if (centerPix[0]+cone_limit < stack_width)  endX = centerPix[0]+cone_limit;
		int endY   = stack_height;  if (centerPix[1]+cone_limit < stack_height) endY = centerPix[1]+cone_limit;
		int endZ   = stack_length;  if (centerPix[2]+cone_limit < stack_length) endZ = centerPix[2]+cone_limit;

		int count=0;
		// this consumes memory! - maybe not necessary to allocate here
		byte[] pixels = new byte[stack_width*stack_height]; //storage for pixels from one layer
		int[][] roi   = new int [2*cone_limit*2*cone_limit*2*cone_limit][4];//roi big enough to take cone at any orientation
		System.out.println("Cone height: "+cone_height);
		System.out.println("Cone radius: "+cone_radius);
		System.out.println("Cone limit: "+cone_limit);
		
		double[] b ={0, 0, 0};  
		double[] c ={0, 0, 0};

		
		Transf.cartesian(ax, ay, az, b, c);
				
		for(int z_layer=startZ; z_layer<endZ; z_layer++){
			pixels = (byte[])stack.getPixels(z_layer);
			for(int x_row=startX; x_row<endX; x_row++){
				for(int y_col=startY; y_col<endY; y_col++){
					
					double point_h, point_r;
					double vx, vy, vz;
					
					vx = x_row-x; 
					vy = y_col-y; 
					vz = z_layer-z;
					
					point_h = vx*ax+vy*ay+vz*az;  
					point_r = Math.sqrt(Math.pow(vx*b[0]+vy*b[1]+vz*b[2], 2)+Math.pow(vx*c[0]+vy*c[1]+vz*c[2], 2));
					
					boolean isInCone = false;
					isInCone = (point_h<=cone_height) && (point_h>=0) && (point_r<=-(cone_radius/cone_height)*point_h+cone_radius);
					
					if(isInCone){
						roi[count][3] = pixels[x_row + stack_width * y_col] & 0xff; //int
						roi[count][0] = x_row; 
						roi[count][1] = y_col; 
						roi[count][2] = z_layer;
						count ++;	
					}
				}
			}
		}	
		
		// quite inefficient for memory to extract twice... but left it this way by now
		int[][] extractedVoxels = new int [count][4];
		for(int i=0; i<count; i++){
			for(int j=0; j<4; j++){
				extractedVoxels[i][j]=roi[i][j];
			}
		}
		return extractedVoxels;
	}
}

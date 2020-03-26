package advantra.shapes;

import advantra.general.Transf;
import ij.ImageStack;

public class ConeCutoff extends RegionOfInterest {

	// 3 dimensional cutoff cone extracted from image - two 	
	private double r;//cone basis  radius
	private double r1;//cone radius at the cutoff
	// center of the cutoff
	private double x1, y1, z1;
	// cone height vector (connects two centers, hence defines orientation as well)
	private double[] h;
	
	public ConeCutoff(double x, double y, double z,  double r,     double x1, double y1, double z1, double r1){
		
		super(x, y, z);
		
		this.roi_type = RoiType.CONE_CUT;
		
		this.r = r;
		
		this.x1 = x1;
		this.y1 = y1;
		this.z1 = z1;
		this.r1 = r1;
		
		this.h = new double[3];
		this.h[0] = x1-x;
		this.h[1] = y1-y;
		this.h[2] = z1-z;
		
	}

	public void setR(double r, double r1){
		this.r = r;
		this.r1 = r1;
	}
	
	public double getR(){
		return this.r;
	}
	
	public double getR1(){
		return this.r1;
	}	
	
	public double getH(){
		return Math.sqrt(Math.pow(h[0], 2)+Math.pow(h[1], 2)+Math.pow(h[2], 2));
	}	
	
	public void setCenter(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
		// update h as well
		this.h[0] = x1-x;
		this.h[1] = y1-y;
		this.h[2] = z1-z;
		
	}
	
	public void setCenter1(double x1, double y1, double z1){
		this.x1 = x1;
		this.y1 = y1;
		this.z1 = z1;
		// update h as well
		this.h[0] = x1-x;
		this.h[1] = y1-y;
		this.h[2] = z1-z;
		
	}
	
	public double getOrientationX(){
		return h[0];
	}
	public double getOrientationY(){
		return h[1];
	}
	public double getOrientationZ(){
		return h[2];
	}	
	
	
	public int[][] extractVoxels(ImageStack stack){ 
		// extracts Sphere region of interest from particular image stack
		// sphere defines center and radius
		// output is nx4 integer array where columns 1-3 contain voxel indexes and column 4 voxel intensity
		int stack_width  = stack.getWidth();
		int stack_height = stack.getHeight();
		int stack_length = stack.getSize();

		// method will extract cone thath's been cutoff from the stack
		
		int cone_radius 	= (int)Math.round(this.r);
		//int cone_radius1 	= (int)Math.round(this.r1);
		int cone_height		= (int)Math.round(this.getH());
		int cone_limit		= (int)Math.ceil(Math.sqrt(Math.pow(r1, 2)+Math.pow(this.getH(), 2)));
		
		if(cone_radius>cone_limit){
			cone_limit = cone_radius;
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
		
		double[] b ={0, 0, 0};  
		double[] c ={0, 0, 0};
		double[] a ={0, 0, 0}; //this one will serve as unit h vector
		a[0] = h[0]/this.getH(); a[1] = h[1]/this.getH(); a[2] = h[2]/this.getH();
		
		Transf.cartesian(a[0], a[1], a[2], b, c);
		
		for(int z_layer=startZ; z_layer<endZ; z_layer++){
			pixels = (byte[])stack.getPixels(z_layer);
			for(int x_row=startX; x_row<endX; x_row++){
				for(int y_col=startY; y_col<endY; y_col++){
					
					double point_h, point_r;
					double vx, vy, vz;
					
					vx = x_row-x; 
					vy = y_col-y; 
					vz = z_layer-z;
					
					point_h = vx*a[0]+vy*a[1]+vz*a[2];  
					point_r = Math.sqrt(Math.pow(vx*b[0]+vy*b[1]+vz*b[2], 2)+Math.pow(vx*c[0]+vy*c[1]+vz*c[2], 2));
					
					boolean isInCone = true;
					isInCone = (point_h<=cone_height) && (point_h>=0) && (point_r<= ((r1-r)/this.getH())*point_h+r);//
					
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

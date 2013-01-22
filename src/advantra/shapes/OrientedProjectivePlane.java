package advantra.shapes;

import advantra.general.ArrayHandling;
import advantra.general.ImageConversions;
import advantra.general.Sort;
import advantra.general.Transf;
import advantra.processing.IntensityCalc;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.process.ColorProcessor;

public class OrientedProjectivePlane extends RegionOfInterest {

	/*
	 * 
	 * this is a square shaped plane that is defined with 
	 *  - reference centerpoint (from RegionOfInterest class)			:(cx, xy, cz)
	 * 	- orientation in 3d												:(vx, vy, vz)
	 * 	- and the distance from centerpoint following the orientation 	:(d)
	 * plane center is at 		:(cx, cy, cz)+d*(vx, vy, vz) -> in local (0, 0, z),
	 * local coordinate system is attached to center & the direction
	 * plane size is (2*d)X(2*d) in metric, and is covered with (numberOfPts x numberOfPts) resolution
	 * 
	 */
	
	private	double   	d;				// 	distance
	private double 		d_plane;		// thickness of the plane in metric
	private double 		vx, vy, vz; 	//	unit orientation vec.
	//private int			numberOfPts; 	// 	plane resolution will be (numberOfPts x numberOfPts)
	
	public 		OrientedProjectivePlane(){
		super(0, 0, 0); 			// roi mother class default initialize
		this.roi_type 		= RoiType.ORIENTED_PROJECTIVE_PLANE;
		this.d				= 1;
		this.d_plane		= 1;		
		this.vx 			= 1; 		// default orientation
		this.vy 			= 0;
		this.vz 			= 0;
		//this.numberOfPts 	= (int)Math.ceil(2*d);	 // 2x2
	}
	
	public		OrientedProjectivePlane(Sphere s, double d_plane, double[] v){
		super(s.getCenterX(), s.getCenterY(), s.getCenterZ());
		this.roi_type 		= RoiType.ORIENTED_PROJECTIVE_PLANE;
		this.d				= s.getR();
		this.d_plane		= d_plane;
		this.vx 			= v[0];
		this.vy 			= v[1];
		this.vz 			= v[2];
	}
	
	public		OrientedProjectivePlane(double[] c, double d, double d_plane, double[] v){
		super(c[0], c[1], c[2]);
		this.roi_type 		= RoiType.ORIENTED_PROJECTIVE_PLANE;
		this.d				= d;
		this.d_plane		= d_plane;
		this.vx 			= v[0];
		this.vy 			= v[1];
		this.vz 			= v[2];
	}
	
	public		OrientedProjectivePlane(double cx, double cy, double cz, double d, double d_plane, double vx, double vy, double vz){
		super(cx, cy, cz);
		this.roi_type 		= RoiType.ORIENTED_PROJECTIVE_PLANE;
		this.d				= d;
		this.d_plane		= d_plane;
		this.vx 			= vx;
		this.vy 			= vy;
		this.vz 			= vz;
	}
	
	public ImagePlus	extract(ImagePlus input_image, int resolution){
		
		ImagePlus output 		= NewImage.createByteImage(
				"exported_oriented_plane", resolution, resolution, 1, NewImage.FILL_BLACK);
		
		//ImageStack input_stack = input_image.getStack();
		int h = input_image.getHeight();
		int w = input_image.getWidth();
		int s = input_image.getStack().getSize();
		
		// allocate cartesian coordinates 
		double[][] plane_coords = new double[resolution*resolution][3];
		
		int idx = 0;
		IntensityCalc calc = new IntensityCalc(input_image.getStack());
		// 'resolution' corresponds to image height&width
		
		for (int r = 0; r < resolution; r++) {
			for (int c = 0; c < resolution; c++) {
				
				// take the distance towards center
				plane_coords[idx][0] = (r*(2f/(resolution-1)) - 1) * (2*d);
				plane_coords[idx][1] = (c*(2f/(resolution-1)) - 1) * (2*d);
				plane_coords[idx][2] = d;
				
				
				// transfer local coordinates to global for sampling from the image
				Transf.rotate(vx, vy, vz, plane_coords[idx]);
				Transf.shift(x, y, z, plane_coords[idx]);

//				System.out.format("%f, %f, %f \n", 
//						plane_coords[idx][0],
//						plane_coords[idx][1],
//						plane_coords[idx][2]);
				
				if(
						plane_coords[idx][0]>0 && plane_coords[idx][0]<=h-1 && 
								plane_coords[idx][1]>0 && plane_coords[idx][1]<=w-1	&& 
										plane_coords[idx][2]>0 && plane_coords[idx][2]<=s-1){
					
					int value = (int)Math.round(
							calc.interpolateAt(plane_coords[idx][0], plane_coords[idx][1], plane_coords[idx][2]));
					
					output.getProcessor().set(c, r, value);
					
				}
				
				/*
				int row_get = (int)Math.round(plane_coords[idx][0]);
				int col_get = (int)Math.round(plane_coords[idx][1]);
				int lay_get = (int)Math.round(plane_coords[idx][2]);
				*/
				
				idx++;
				
			}
		}
		
		return output;
		
	}
	
	public int[][] extractCoords(ImagePlus input_image){
		
		// coordinates of the 8 border points
		double[] u = new double[3];
		double[] w = new double[3];
		Transf.cartesian(vx, vy, vz, u, w); // u and w are unit vectors
		
		double[] plane_center = new double[3];
		plane_center[0] = x+vx*d;	plane_center[1] = y+vy*d;	plane_center[0] = z+vz*d;
		
		double[][] borders = new double[3][8];
		
		int idx = 0;
		borders[idx][0] = plane_center[idx]-2*d*u[idx]-2*d*w[idx]-(d_plane/2)*vx;
		borders[idx][1] = plane_center[idx]-2*d*u[idx]-2*d*w[idx]+(d_plane/2)*vx;
		borders[idx][2] = plane_center[idx]-2*d*u[idx]+2*d*w[idx]-(d_plane/2)*vx;
		borders[idx][3] = plane_center[idx]-2*d*u[idx]+2*d*w[idx]+(d_plane/2)*vx;
		borders[idx][4] = plane_center[idx]+2*d*u[idx]-2*d*w[idx]-(d_plane/2)*vx;
		borders[idx][5] = plane_center[idx]+2*d*u[idx]-2*d*w[idx]+(d_plane/2)*vx;
		borders[idx][6] = plane_center[idx]+2*d*u[idx]+2*d*w[idx]-(d_plane/2)*vx;
		borders[idx][7] = plane_center[idx]+2*d*u[idx]+2*d*w[idx]+(d_plane/2)*vx;
		idx = 1;
		borders[idx][0] = plane_center[idx]-2*d*u[idx]-2*d*w[idx]-(d_plane/2)*vy;
		borders[idx][1] = plane_center[idx]-2*d*u[idx]-2*d*w[idx]+(d_plane/2)*vy;
		borders[idx][2] = plane_center[idx]-2*d*u[idx]+2*d*w[idx]-(d_plane/2)*vy;
		borders[idx][3] = plane_center[idx]-2*d*u[idx]+2*d*w[idx]+(d_plane/2)*vy;
		borders[idx][4] = plane_center[idx]+2*d*u[idx]-2*d*w[idx]-(d_plane/2)*vy;
		borders[idx][5] = plane_center[idx]+2*d*u[idx]-2*d*w[idx]+(d_plane/2)*vy;
		borders[idx][6] = plane_center[idx]+2*d*u[idx]+2*d*w[idx]-(d_plane/2)*vy;
		borders[idx][7] = plane_center[idx]+2*d*u[idx]+2*d*w[idx]+(d_plane/2)*vy;
		idx = 2;
		borders[idx][0] = plane_center[idx]-2*d*u[idx]-2*d*w[idx]-(d_plane/2)*vz;
		borders[idx][1] = plane_center[idx]-2*d*u[idx]-2*d*w[idx]+(d_plane/2)*vz;
		borders[idx][2] = plane_center[idx]-2*d*u[idx]+2*d*w[idx]-(d_plane/2)*vz;
		borders[idx][3] = plane_center[idx]-2*d*u[idx]+2*d*w[idx]+(d_plane/2)*vz;
		borders[idx][4] = plane_center[idx]+2*d*u[idx]-2*d*w[idx]-(d_plane/2)*vz;
		borders[idx][5] = plane_center[idx]+2*d*u[idx]-2*d*w[idx]+(d_plane/2)*vz;
		borders[idx][6] = plane_center[idx]+2*d*u[idx]+2*d*w[idx]-(d_plane/2)*vz;
		borders[idx][7] = plane_center[idx]+2*d*u[idx]+2*d*w[idx]+(d_plane/2)*vz;
		
		// define the range
		
		int height = input_image.getStack().getHeight();
		int width  = input_image.getStack().getWidth();
		int len    = input_image.getStack().getSize();
		
		int startX = (int)Math.round(Sort.findMin(borders[0])); 
		startX = (startX < 0)? 0: startX;
		int startY = (int)Math.round(Sort.findMin(borders[1])); 
		startY = (startY < 0)? 0: startY;
		int startZ = (int)Math.round(Sort.findMin(borders[2])); 
		startZ = (startZ < 0)? 0: startZ;
		
		int endX = (int)Math.round(Sort.findMax(borders[0])); 
		endX = (endX > height-1)? height-1: endX;
		int endY = (int)Math.round(Sort.findMax(borders[1])); 
		endY = (endY > width-1)? width-1: endY;
		int endZ = (int)Math.round(Sort.findMax(borders[2])); 
		endZ = (endZ > len-1)? len-1: endZ;

		int count=0;
		int to_allocate = 16*(int)Math.ceil(d)*(int)Math.ceil(d)*(int)Math.ceil(d_plane);
		int[][] roi_coord = new int[3][to_allocate];
		for(int x=startX; x<=endX; x++){
			for(int y=startY; y<=endY; y++){
				for(int z=startZ; z<=endZ; z++){
					
					boolean isInPlane = 
							Math.abs(Transf.dotProd(new double[]{x,y,z}, new double[]{this.x, this.y, this.z}))<=2*d;
					
					if(isInPlane){
						
						roi_coord[0][count] = x; 
						roi_coord[1][count] = y; 
						roi_coord[2][count] = z;
						
						count ++;
						
						if(count>to_allocate){
							break;
						}
						
					}
				}
			}
		}
		
		if(count>to_allocate){
			return roi_coord;
		}
		else{
			// resize the output array - count will be the new size
			int[][] roi_coord_output = new int[3][count];
			for (int i = 0; i < count; i++) {
				roi_coord_output[0][i] = roi_coord[0][i];
				roi_coord_output[1][i] = roi_coord[1][i];
				roi_coord_output[2][i] = roi_coord[2][i];
			}
			
			return roi_coord_output;
		}
		
		
	}
	
	public void drawOverColorImage(ImagePlus template_image, int r, int g, int b){
		
		int w = template_image.getWidth();
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}

		byte[] rgb_values = new byte[]{(byte)r, (byte)g, (byte)b}; 
		
		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
		// extract the points
		int[][] cyl_coordinates = extractCoords(template_image);
		
		if(cyl_coordinates[0].length>0){
			
			// change the spots, assign them with valueToDraw
			for (int i = 0; i < cyl_coordinates[0].length; i++) {
				
				int row 	= cyl_coordinates[0][i];
				int col 	= cyl_coordinates[1][i];
				int layer 	= cyl_coordinates[2][i];
				
				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[0];
				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[1];
				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[2];
			}
			
		}else{
			
			System.out.println("OrientedProjecitvePlane:drawOverColorImage()\n" +
					"There was not enough points to draw anything!");
		
		}
		
		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}

	}
	
}

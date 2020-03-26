package advantra.shapes;

import advantra.general.ArrayHandling;
import advantra.general.ImageConversions;
import ij.process.ColorProcessor;
import ij.ImagePlus;

public class Point  extends RegionOfInterest  {

	public Point(){ 	// dummy construction
		super(0, 0, 0);
		this.roi_type = RoiType.POINT;
	}
	
	public Point(double x, double y, double z){
		super(x, y, z);
		this.roi_type = RoiType.POINT;
	}
	
	public Point(double[] p){
		super(p[0], p[1], p[2]);
		this.roi_type = RoiType.POINT;
	}
	
	public void setPoint(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/*	##################################
	 *  EXTRACTION METHODS
	 * 	##################################
	 */	


	
	/*	##################################
	 *  GRAPHICS
	 * 	##################################
	 */	
	public void drawOverColorImageStack(ImagePlus template_image, int r, int g, int b){
		
		if(template_image.getStack().getSize()==1){
			System.err.println("Point:drawOverColorImageStack() takes image stack as argument!");
			System.exit(1);
		}
		
		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
			template_image = ImageConversions.ImagePlusToRGB(template_image);
			System.out.println("converting ImagePlus to rgb...");
		}
		
		int pt_x = (int) Math.round(x);
		int pt_y = (int) Math.round(y);
		int pt_z = (int) Math.round(z);
		
		boolean isIn = 
				pt_x>=0 && pt_x <=template_image.getHeight()-1 &&
				pt_y>=0 && pt_y <=template_image.getWidth()-1 &&
				pt_z>=0 && pt_z <=template_image.getStack().getSize()-1;
		
		if(isIn){
			int w = template_image.getWidth();

			byte[] rgb_values = new byte[]{(byte)r, (byte)g, (byte)b}; 
		
			byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
		
			img_array[0][pt_z][ArrayHandling.sub2index_2d(pt_x, pt_y, w)] = rgb_values[0];
			img_array[1][pt_z][ArrayHandling.sub2index_2d(pt_x, pt_y, w)] = rgb_values[1];
			img_array[2][pt_z][ArrayHandling.sub2index_2d(pt_x, pt_y, w)] = rgb_values[2];
		
			for (int i = 1; i <= template_image.getStack().getSize(); i++) {
			
				((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
			}
		
		}else{
			
			System.out.println("Point:drawOverColorImageStack()\n" +
					"Point was out of the image to be drawn!");
		
		}

	}
	
	public static void drawOverColorImageStack(
			Point[] pts, 
			ImagePlus template_image, 
			int r, int g, int b){
		for (int i = 0; i < pts.length; i++) {
			pts[i].drawOverColorImageStack(template_image, r, g, b);
		}
	}
	
	public ImagePlus drawOverGrayImageStack(ImagePlus template_image, int valueToDraw){
		
		int wt = template_image.getStack().getWidth();
		
		// check if it's grey image
		if(template_image.getType()==ImagePlus.GRAY8){
			// expects 0-255 8-bit scale
			valueToDraw = (valueToDraw>255) ? 255 : valueToDraw ;
			valueToDraw = (valueToDraw<  0) ?   0 : valueToDraw ;
		}
		else if(template_image.getType()==ImagePlus.GRAY16){
			// expects 0-65,535 16-bit scale
			valueToDraw = (valueToDraw>65535) ? 65535 : valueToDraw ;
			valueToDraw = (valueToDraw<  0  ) ?   0   : valueToDraw ;
		}
		else{
			System.err.println("Point:drawOverImage()\nthis image format is not suported - possible to use gray8, gray16 or rgb!");
			System.exit(-1);
		}
		
		int[][] template_image_array = ImageConversions.GraytoIntArray(template_image);
		
		int pt_x = (int) Math.round(x);
		int pt_y = (int) Math.round(y);
		int pt_z = (int) Math.round(z);
		
		boolean isIn = 
				pt_x>=0 && pt_x <=template_image.getHeight()-1 &&
				pt_y>=0 && pt_y <=template_image.getWidth()-1 &&
				pt_z>=0 && pt_z <=template_image.getStack().getSize()-1;
		
		if(isIn){
			template_image_array[pt_z][pt_y+pt_x*wt] = valueToDraw;
		}else{
			System.out.println("Point:drawOverImage()\nThere was not enough points to draw anything!");
		}

		// to ImagePlus
		if(template_image.getType()==ImagePlus.GRAY8){
			return ImageConversions.toGray8(template_image_array, wt);
		}
		else if(template_image.getType()==ImagePlus.GRAY16){
			return ImageConversions.toGray16(template_image_array, wt);
			
		}
		else{
			
			System.err.println("Sphere:drawOverImage(): this image format is not suported - possible to use gray8, gray16 or rgb!");
			System.exit(-1);
			return ImageConversions.toGray8(template_image_array, wt); // just a dummy line to avoid error 
			
		}
		
	}



}

	
//	public void drawOverColorImage(ImagePlus template_image, int valueToDraw){
//		
//		int w = template_image.getWidth();
//		
//		if((template_image.getType()!=ImagePlus.COLOR_RGB)){
//			// set it to rgb
//			template_image = ImageConversions.ImagePlusToRGB(template_image);
//			System.out.println("converting ImagePlus to rgb...");
//		}
//
//		byte[] rgb_values = ColourTransf.Jet256(valueToDraw);
//		
//		byte[][][] img_array = ImageConversions.RgbToByteArray(template_image);
//		
//		// extract the points
//		int[][] cyl_coordinates = new int[3][Sphere.numberOfVoxInSphere((int)Math.ceil(r))]; 
//		
//		int cnt = this.extractCoords(template_image, cyl_coordinates);
//		
//		if(cnt>0){
//			
//			// change the spots, assign them with valueToDraw
//			for (int i = 0; i < cnt; i++) {
//				
//				int row 	= cyl_coordinates[0][i];
//				int col 	= cyl_coordinates[1][i];
//				int layer 	= cyl_coordinates[2][i];
//				
//				img_array[0][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[0];
//				img_array[1][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[1];
//				img_array[2][layer][ArrayHandling.sub2index_2d(row, col, w)] = rgb_values[2];
//			}
//			
//		}else{
//			
//			System.out.println("Sphere:drawOverImage()\nThere was not enough points to draw anything!");
//		
//		}
//		
//		for (int i = 1; i <= template_image.getStack().getSize(); i++) {
//			
//			((ColorProcessor)template_image.getStack().getProcessor(i)).setRGB(
//					img_array[0][i-1], 
//					img_array[1][i-1], 
//					img_array[2][i-1]
//							);
//		}
//
//	}
	

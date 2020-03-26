package advantra.tools;

import ij.ImagePlus;
import ij.gui.NewImage;
import advantra.general.ArrayHandling;
import advantra.general.ArrayHandling.IdxMode;

public class MeanShiftMasks {

	double 	angular_span;
	int		height;
	
	boolean[][][][] masks;
	
	public MeanShiftMasks(int height, double angular_span){ // template width will be 2xheight automatically
		
		this.angular_span 	= angular_span;
		this.height			= height;
		this.masks 			= new boolean[height][2*height][height][2*height];
		
	}
	
	public void createMasks(){
		
		System.out.print("Creating masks");
		
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < 2*height; j++) {
				
				double current_phi		= ArrayHandling.index2value(i, IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	height);
				double current_theta	= ArrayHandling.index2value(j, IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 	2*height);
				
				double current_x 		= Math.sin(current_phi)*Math.cos(current_theta);
				double current_y 		= Math.sin(current_phi)*Math.sin(current_theta);
				double current_z 		= Math.cos(current_phi);
				
				// mask at (i,j)... check for each direction from the mask whether it belongs to predefined angular span
				for (int k = 0; k < height; k++) {
					for (int l = 0; l < 2*height; l++) {
						
						double mask_phi		= ArrayHandling.index2value(k, IdxMode.LAST_INCLUDED, 0.0, Math.PI, 	height);
						double mask_theta	= ArrayHandling.index2value(l, IdxMode.LAST_EXCLUDED, 0.0, 2*Math.PI, 	2*height);
						
						double mask_x 		= Math.sin(mask_phi)*Math.cos(mask_theta);
						double mask_y 		= Math.sin(mask_phi)*Math.sin(mask_theta);
						double mask_z 		= Math.cos(mask_phi);
						
						// calculate angle
						double cos_mask_current = current_x*mask_x + current_y*mask_y + current_z*mask_z;
						//double cos_solid_angle = 1 - (angular_span/(2*Math.PI));
						
						if(cos_mask_current>Math.cos(angular_span)){
							masks[i][j][k][l] = true;
						}
					}
				}
				
				System.out.print(".");
			}
		}
		
		System.out.println("done.");
		
	}
	
	public ImagePlus getMask(int atRow, int atCol){
		
		if(atRow<0 || atRow>=height || atCol<0 || atCol>=(2*height)){
			System.err.println("index out of the range...");
			System.exit(1);
		}
		
		ImagePlus im = NewImage.createByteImage("extracted_mask", (2*height), height, 1, NewImage.FILL_BLACK);
		
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < 2*height; j++) {
				if(masks[atRow][atCol][i][j]){
					im.getStack().setVoxel(j, i, 0, 255);
				}
				
			}
		}
		
		return im;
		
	}
	
	public boolean[][][][] getMasks(){
		return masks;
	}
	
	public boolean isInNeighbourhood(int atRowPos, int atColPos, int referenceRowPos, int referenceColPos){
		return masks[referenceRowPos][referenceColPos][atRowPos][atColPos];
	}
	
	public boolean isInNeighbourhood(double[] atPos, double[] referencePos){
		
		int refX = (int)Math.round(referencePos[0]);
		int refY = (int)Math.round(referencePos[1]);
		int atX  = (int)Math.round(atPos[0]);
		int atY  = (int)Math.round(atPos[1]);
		
		return masks[refX][refY][atX][atY];
		
	}
	
}

package advantra.bifurcations;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import advantra.general.CreateDirectory;
import advantra.processing.IntensityCalc;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.NewImage;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class CreateTrainSet implements PlugInFilter, MouseListener, KeyListener {
	
	ImagePlus 		img;
	ImageCanvas 	canvas;
	
	ImagePlus[] 	img_patches;
	ImageCanvas[]	canvas_patches;
	
	int 			start_scale = 3;
	int 			end_scale 	= 3;
	int 			patch_size 	= 30; //()^3
	
	String 			out_path;
	IntensityCalc 	calc;
	
	int 			count_pts;
	
	public int setup(String arg, ImagePlus img) {
		this.img = img;
		return DOES_8G+NO_CHANGES;
	}
	
	public void run(ImageProcessor ip) {
		
		ImageWindow win = img.getWindow();
		canvas = win.getCanvas();
		canvas.addMouseListener(this);
		
		// dialog to enter input values
		GenericDialog gd = new GenericDialog("Create Train Set", IJ.getInstance());
		gd.addMessage("Extract cube patches :");
		gd.addNumericField( "start S :", start_scale, 0, 5, "2Sx2Sx2S");
		gd.addNumericField( "end   S :", end_scale  , 0, 5, "2Sx2Sx2S");
		gd.addMessage("Extracted patch dimensions :");
		gd.addNumericField( "size :", patch_size, 0, 5, "size x size x size");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		start_scale = (int)gd.getNextNumber();
		end_scale 	= (int)gd.getNextNumber();
		patch_size 	= (int)gd.getNextNumber();
		
		if(end_scale<start_scale){
			System.out.println("'end scale' has to be higher than the 'start_scale'");
			return;
		}
		
		out_path	= System.getProperty("user.home")+File.separator+"train_set_patches";
		CreateDirectory.createOneDir(out_path);
		
		calc = new IntensityCalc(img.getStack()); // for the interpolations
		count_pts = 0;
		
		img_patches = new ImagePlus[end_scale-start_scale+1];
		canvas_patches = new ImageCanvas[end_scale-start_scale+1];
		
		// images of the patches and their canvases to save
		for (int sc = 0; sc <= end_scale-start_scale; sc++) {
			String patch_name = String.format("patch_scale_%d", (2*(sc+start_scale)));
			img_patches[sc] = NewImage.createByteImage(patch_name, patch_size, patch_size, patch_size, NewImage.FILL_BLACK);
			canvas_patches[sc] 	= new ImageCanvas(img_patches[sc]);
			canvas_patches[sc].addMouseListener(this);
		}
		
	}
	
	public void mouseClicked(MouseEvent e) {
		
		int mouseY = canvas.offScreenX(e.getX());
		int mouseX = canvas.offScreenY(e.getY());
		int mouseZ =  	img.getCurrentSlice()-1;
		
		for (int s = start_scale; s <= end_scale; s++) {
			
			if( (mouseX-s)>=0 && (mouseX+s)<img.getHeight() && (mouseY-s)>=0 && (mouseY+s)<img.getWidth() && (mouseZ-s)>=0 && (mouseZ+s)<img.getStack().getSize()){
				// it's within the image
				byte[][] img_patch_array = new byte[patch_size][patch_size*patch_size];
				
				float step = (2*s)/((float)patch_size-1);
				
					for (int row = 0; row < patch_size; row++) {
						for (int col = 0; col < patch_size; col++) {
							for (int lay = 0; lay < patch_size; lay++) {
							
								float x = mouseX-s+row*step;
								float y = mouseY-s+col*step;
								float z = mouseZ-s+lay*step;
							
								img_patch_array[lay][row*patch_size+col] = (byte)( (int) (Math.round(calc.interpolateAt_new(x, y, z))) );
							
							}
						}
					}
				
					for (int i = 0; i < patch_size; i++) {
						img_patches[s-start_scale].getStack().setPixels(img_patch_array[i], (i+1));
					}
					
					String file_name = String.format("%s%s_(%d,%d,%d,%d).tif", 
							out_path, File.separator, (2*s), mouseX, mouseY, mouseZ);
					img_patches[s-start_scale].setTitle(file_name);
					img_patches[s-start_scale].show();
				
			}
		}
			
	}
		
//		else{
//			
//			// regular click sensed only on array of patch images
//			// dialog to save the patch
//			GenericDialog gd = new GenericDialog("Save patch", IJ.getInstance());
//			gd.addChoice("Save as:", new String[] {"background", "body", "bifurcation", "end-point", "cross-over", "cancel-save"}, "");
//			gd.showDialog();
//			if (gd.wasCanceled()) return;
//			
//			String save_name;
//			switch(gd.getNextChoiceIndex()){
//			
//			case 0:
//				System.out.println("background was chosen");
//				count_pts++;
//				save_name = String.format(
//						"%s%s%d_background_size_(%dx%dx%d)_pos_(%d,%d,%d).tif", 
//						out_path, File.separator, count_pts, (2*s), (2*s), (2*s), mouseX, mouseY, mouseZ);
//				(new FileSaver(img_patch)).saveAsTiffStack(save_name);
//				System.out.println(save_name+" is saved...");
//				break;
//			case 1:
//				System.out.println("body was chosen");
//				count_pts++;
//				save_name = String.format(
//						"%s%s%d_body_size_(%dx%dx%d)_pos_(%d,%d,%d).tif", 
//						out_path, File.separator, count_pts, (2*s), (2*s), (2*s), mouseX, mouseY, mouseZ);
//				(new FileSaver(img_patch)).saveAsTiffStack(save_name);
//				System.out.println(save_name+" is saved...");
//				break;
//			case 2:
//				System.out.println("bifurcation was chosen");
//				count_pts++;
//				save_name = String.format(
//						"%s%s%d_bifurcation_size_(%dx%dx%d)_pos_(%d,%d,%d).tif", 
//						out_path, File.separator, count_pts, (2*s), (2*s), (2*s), mouseX, mouseY, mouseZ);
//				(new FileSaver(img_patch)).saveAsTiffStack(save_name);
//				System.out.println(save_name+" is saved...");
//				break;
//			case 3:
//				System.out.println("end-point was chosen");
//				count_pts++;
//				save_name = String.format(
//						"%s%s%d_end-point_size_(%dx%dx%d)_pos_(%d,%d,%d).tif", 
//						out_path, File.separator, count_pts, (2*s), (2*s), (2*s), mouseX, mouseY, mouseZ);
//				(new FileSaver(img_patch)).saveAsTiffStack(save_name);
//				System.out.println(save_name+" is saved...");
//				break;
//			case 4:
//				System.out.println("cross-over was chosen");
//				count_pts++;
//				save_name = String.format(
//						"%s%s%d_cross-over_size_(%dx%dx%d)_pos_(%d,%d,%d).tif", 
//						out_path, File.separator, count_pts, (2*s), (2*s), (2*s), mouseX, mouseY, mouseZ);
//				(new FileSaver(img_patch)).saveAsTiffStack(save_name);
//				System.out.println(save_name+" is saved...");
//				break;
//			case 5:
//				System.out.println("cancelling the save...");
//				break;
//			default:
//				break;
//			}
//		
//		}
		
//	}		

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void keyTyped(KeyEvent e) {
		System.out.println(e.getKeyCode()+" was typed...");

		
	}
	
	public void keyPressed(KeyEvent e) {
		System.out.println(e.getKeyCode()+" was pressed...");

	}

	public void keyReleased(KeyEvent e) {
		
		
	}

}

package advantra.critpoint;

import java.io.File;
import java.io.FilenameFilter;

import advantra.feature.DifferentialStructure;
import advantra.file.AnalyzeCSV;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class ExtractFeatures implements PlugIn {

	String folder_name = System.getProperty("user.dir");
	ImagePlus 	train_img;
	double[][] 	train_loc;
	int[]		train_cls;
	
	// extraction parameters (standard deviations of the Gaussan used for smoothing/scaling)
	double 	sigma_1 = 2.0;
	double 	sigma_2 = 2.0;
	int		nr 		= 1;
	
	
	public void run(String arg0) {
		GenericDialog gd = new GenericDialog("Critical Point features");
		
		gd.addMessage("Locate folder...");
		gd.addStringField("folder with annotations", folder_name, 100);
		
		gd.addNumericField( "sigma start:", sigma_1, 	0, 5, "pix" );
		gd.addNumericField( "sigma end  :", sigma_2, 	0, 5, "pix" );	
		gd.addNumericField( "number of scales : ", nr,  0, 5, "");
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		folder_name 	= gd.getNextString();
		sigma_1 		= (double)gd.getNextNumber();
		sigma_2 		= (double)gd.getNextNumber();
		nr 				= (int)gd.getNextNumber();
		
		File dir = new File(folder_name);
		folder_name = dir.getAbsolutePath();
		
		System.out.println("loading files from : "+folder_name);
		
		if(!dir.isDirectory()){
			IJ.error(folder_name+" is not a directory!");
			return;
		}
		
		File[] csv_files = dir.listFiles(new FilenameFilter() {
		    	public boolean accept(File dir, String name) {
		    		return name.toLowerCase().endsWith(".csv");
		    	}
			}
		);
		
		File[] tif_files = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
				}
			}
		);
		
		for (int i = 0; i < csv_files.length; i++) {
			
			String 	current_csv_name 	= csv_files[i].getName();
			System.out.print("extracting "+current_csv_name+" ... ");
			boolean found = false;
			int idx_found = 0;
			// is there a .tif pair, remember it's index
			for (int j = 0; j < tif_files.length; j++) {
				String current_tif_name = tif_files[j].getName();
				if(removeExt(current_csv_name).equals(removeExt(current_tif_name))){
					found = true;
					idx_found = j;
					break;
				}
			}
			
			if(!found){
				System.out.println("FAILED");
				continue; // try other .csv-s
			}
			
			System.out.println();
			System.out.println("processing "+tif_files[idx_found].getName()+"...");
			
			
			// match was found, extract image, locations and class value
			train_img = new ImagePlus(tif_files[idx_found].getAbsolutePath());
			AnalyzeCSV reader_csv = new AnalyzeCSV(csv_files[i].getAbsolutePath());
			
			train_loc = reader_csv.readLn(2);
			train_cls = reader_csv.readLastCol();
			DifferentialStructure df = new DifferentialStructure(train_img, sigma_1, sigma_2, nr);
			df.getHessianDeterminant().show();
			df.getDoH().show();
			df.getBallness().show();
			double[][] feats  = df.exportFeatures(train_loc);
			
			
		}
		
		
		
	}
	
	private String removeExt(String in_name){
		int name_len = in_name.length();
		if(name_len>4){
			name_len -= 4;
			return in_name.substring(0, name_len);
		}
		return "";
	}
	
	private double[][] concatenateCols(double[][] in11, double[][] in12){
		 
		double[][] out = new double[in11.length][in11[0].length+in12[0].length];
		 
		 for (int i = 0; i < in11.length; i++) {
			for (int j = 0; j < in11[0].length; j++) {
				out[i][j] = in11[i][j];
			}
		 }
		 
		 for (int i = 0; i < in12.length; i++) {
			 for (int j = 0; j < in12[0].length; j++) {
				 out[i][j+in11[0].length] = in12[i][j];
			 }
		 }
		 
		 return out;
	}
	
	private double[][] concatenateRows(double[][] in11, double[][] in21){
		
		double[][] out = new double[in11.length+in21.length][in11[0].length];
		
		 for (int i = 0; i < in11.length; i++) {
			for (int j = 0; j < in11[0].length; j++) {
				out[i][j] = in11[i][j];
			}
		 }
		 
		 for (int i = 0; i < in21.length; i++) {
			 for (int j = 0; j < in21[0].length; j++) {
				 out[i+in11.length][j] = in21[i][j];
			 }
		 }
		
		return out;
		
	}
	
}	
//	ImagePlus 		img;
//	ImageCanvas 	canvas;
//	
//	ImagePlus[] 	img_patches;
//	ImageCanvas[]	canvas_patches;
//	
//	int 			start_scale = 3;
//	int 			end_scale 	= 3;
//	int 			patch_size 	= 30; //()^3
//	
//	String 			out_path;
//	IntensityCalc 	calc;
//	
//	int 			count_pts;
	
//	public int setup(String arg, ImagePlus img) {
//		this.img = img;
//		return DOES_8G+NO_CHANGES;
//	}
	
//	public void run(ImageProcessor ip) {
//		
//		ImageWindow win = img.getWindow();
//		canvas = win.getCanvas();
//		canvas.addMouseListener(this);
//		
//		// dialog to enter input values
//		GenericDialog gd = new GenericDialog("Create Train Set", IJ.getInstance());
//		gd.addMessage("Extract cube patches :");
//		gd.addNumericField( "start S :", start_scale, 0, 5, "2Sx2Sx2S");
//		gd.addNumericField( "end   S :", end_scale  , 0, 5, "2Sx2Sx2S");
//		gd.addMessage("Extracted patch dimensions :");
//		gd.addNumericField( "size :", patch_size, 0, 5, "size x size x size");
//		gd.showDialog();
//		if (gd.wasCanceled()) return;
//		
//		start_scale = (int)gd.getNextNumber();
//		end_scale 	= (int)gd.getNextNumber();
//		patch_size 	= (int)gd.getNextNumber();
//		
//		if(end_scale<start_scale){
//			System.out.println("'end scale' has to be higher than the 'start_scale'");
//			return;
//		}
//		
//		out_path	= System.getProperty("user.home")+File.separator+"train_set_patches";
//		CreateDirectory.createOneDir(out_path);
//		
//		calc = new IntensityCalc(img.getStack()); // for the interpolations
//		count_pts = 0;
//		
//		img_patches = new ImagePlus[end_scale-start_scale+1];
//		canvas_patches = new ImageCanvas[end_scale-start_scale+1];
//		
//		// images of the patches and their canvases to save
//		for (int sc = 0; sc <= end_scale-start_scale; sc++) {
//			String patch_name = String.format("patch_scale_%d", (2*(sc+start_scale)));
//			img_patches[sc] = NewImage.createByteImage(patch_name, patch_size, patch_size, patch_size, NewImage.FILL_BLACK);
//			canvas_patches[sc] 	= new ImageCanvas(img_patches[sc]);
//			canvas_patches[sc].addMouseListener(this);
//		}
//		
//	}
	
//	public void mouseClicked(MouseEvent e) {
//		
//		int mouseY = canvas.offScreenX(e.getX());
//		int mouseX = canvas.offScreenY(e.getY());
//		int mouseZ =  	img.getCurrentSlice()-1;
//		
//		for (int s = start_scale; s <= end_scale; s++) {
//			
//			if( (mouseX-s)>=0 && (mouseX+s)<img.getHeight() && (mouseY-s)>=0 && (mouseY+s)<img.getWidth() && (mouseZ-s)>=0 && (mouseZ+s)<img.getStack().getSize()){
//				// it's within the image
//				byte[][] img_patch_array = new byte[patch_size][patch_size*patch_size];
//				
//				float step = (2*s)/((float)patch_size-1);
//				
//					for (int row = 0; row < patch_size; row++) {
//						for (int col = 0; col < patch_size; col++) {
//							for (int lay = 0; lay < patch_size; lay++) {
//							
//								float x = mouseX-s+row*step;
//								float y = mouseY-s+col*step;
//								float z = mouseZ-s+lay*step;
//							
//								img_patch_array[lay][row*patch_size+col] = (byte)( (int) (Math.round(calc.interpolateAt_new(x, y, z))) );
//							
//							}
//						}
//					}
//				
//					for (int i = 0; i < patch_size; i++) {
//						img_patches[s-start_scale].getStack().setPixels(img_patch_array[i], (i+1));
//					}
//					
//					String file_name = String.format("%s%s_(%d,%d,%d,%d).tif", 
//							out_path, File.separator, (2*s), mouseX, mouseY, mouseZ);
//					img_patches[s-start_scale].setTitle(file_name);
//					img_patches[s-start_scale].show();
//				
//			}
//		}
//			
//	}
		
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

//	public void mousePressed(MouseEvent e) {}
//
//	public void mouseReleased(MouseEvent e) {}
//
//	public void mouseEntered(MouseEvent e) {}
//
//	public void mouseExited(MouseEvent e) {}
//
//	public void keyTyped(KeyEvent e) {
//		System.out.println(e.getKeyCode()+" was typed...");
//
//		
//	}
//	
//	public void keyPressed(KeyEvent e) {
//		System.out.println(e.getKeyCode()+" was pressed...");
//
//	}
//
//	public void keyReleased(KeyEvent e) {
//		
//		
//	}



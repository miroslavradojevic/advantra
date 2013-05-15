package advantra.plugins;

import java.io.File;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import advantra.tools.BranchModel;

public class GenerateBranch implements PlugIn {
	
	/*
	 * generates the branch and saves it in a output .tif file
	 * can be used as terminal command 
	 * java -cp advantra_.jar:ij-1.47g.jar:flanagan.jar advantra.plugins.GenerateBranch 64 64 32 2.0
	 * or from fiji-imagej as a plugin
	 */
	
	public static void main(String[] args){
		
		// terminal command
		// parameters
		int 	stack_height 	= 64;
		int 	stack_width 	= 64;
		int 	stack_size 		= 32;
		double 	radius_std		= 1.0;
		
		String out_path_h = 
				System.getProperty("user.home")+File.separator+"gen_branch_hor.tif";
		String out_path_v = 
				System.getProperty("user.home")+File.separator+"gen_branch_ver.tif";
		
		if(args.length==4){
			stack_height 	= (int) 	Integer.parseInt(	args[0]);
			stack_width  	= (int) 	Integer.parseInt(	args[1]);
			stack_size		= (int) 	Integer.parseInt(	args[2]);
			radius_std		= (double)	Double.parseDouble( args[3]);
		}
		else{
			System.out.println("Enter again, wrong arguments...\n" +
					"1 - height \n" +
					"2 - width  \n" +
					"3 - size \n" +
					"4 - radius std. \n");
			return;
		}
		
		BranchModel bm = new BranchModel(stack_height, stack_width, stack_size, radius_std, 40, 150);
		
		bm.drawHorizontalModel();
		ImagePlus im_out = bm.getModelAsImage();
		(new FileSaver(im_out)).saveAsTiffStack(out_path_h);
		
		bm.drawVerticalModel();
		im_out = bm.getModelAsImage();
		(new FileSaver(im_out)).saveAsTiffStack(out_path_v);
		
		System.out.println("files exported in :\n" +
				""+out_path_h+" and \n" +
						""+out_path_v);
		
	}

	public void run(String arg0) {

        System.out.println("generating branch...");
		
		// parameters
		int 	stack_height 	= 64;
		int 	stack_width 	= 64;
		int 	stack_size 		= 32;
		double 	radius_std		= 1.0;
		
		String out_path_h = 
				System.getProperty("user.home")+File.separator+"gen_branch_hor.tif";
		String out_path_v = 
				System.getProperty("user.home")+File.separator+"gen_branch_ver.tif";
		
		// dialog to enter input values
		GenericDialog gd = new GenericDialog("Generate Simple Bifurcation", IJ.getInstance());
		gd.addNumericField("image stack height :", stack_height, 0, 5, "pix");
		gd.addNumericField("image stack width  :", stack_width,  0, 5, "pix");
		gd.addNumericField("image stack size   :", stack_size, 	 0, 5, "pix");  
		gd.addNumericField("radius std 		   :", radius_std,   2, 5, "pix" );
		gd.showDialog();
		if (gd.wasCanceled()) return;
		stack_height = (int)gd.getNextNumber();
		stack_width  = (int)gd.getNextNumber();
		stack_size   = (int)gd.getNextNumber();
		radius_std   = (double)gd.getNextNumber();
		
		// imagej plugin
		BranchModel bm = new BranchModel(stack_height, stack_width, stack_size, radius_std, 40, 150);
		bm.drawHorizontalModel();
        IJ.showMessage("Done!");
		ImagePlus im_out = bm.getModelAsImage();
		im_out.setTitle("horizontal_model");
		im_out.show();
		(new FileSaver(im_out)).saveAsTiffStack(out_path_h);
		
		bm.drawVerticalModel();
		im_out = bm.getModelAsImage();
		im_out.setTitle("vertical_model");
		im_out.show();
		(new FileSaver(im_out)).saveAsTiffStack(out_path_v);
		
		System.out.println("files exported in :\n" +
				""+out_path_h+" and \n" +
						""+out_path_v);

		
	}
	
}

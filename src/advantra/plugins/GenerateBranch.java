package advantra.plugins;

import java.io.File;

import advantra.tools.BranchModel2D;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

public class GenerateBranch implements PlugIn {
	
	/*
	 * generates the branch, shows it, and saves it in a output .tif file
	 * terminal command commented out, but was possible:
	 * java -cp advantra_.jar:ij-1.47g.jar:flanagan.jar advantra.plugins.GenerateBranch 64 64 32 2.0
	 * or from fiji-imagej as a plugin
	 */

    String outDirPath;

	public void run(String arg0) {

        int H = 257, W = 257, N = 5;

        ImageStack isOut = new ImageStack(W, H);

        BranchModel2D bm2d = new BranchModel2D(W, H);

        for (int i = 0; i < N; i ++){
            ImageProcessor ip_to_add = new ByteProcessor(W, H);

            ImagePlus imShow = new ImagePlus(""+N, bm2d.generateRandomBranch());

            Overlay ovly = new Overlay();
            ovly.add(new PointRoi(bm2d.pc[0]+0.5, bm2d.pc[1]+0.5));
            ovly.add(new PointRoi(bm2d.p1[0]+0.5, bm2d.p1[1]+0.5));
            ovly.add(new PointRoi(bm2d.p2[0]+0.5, bm2d.p2[1]+0.5));
            ovly.add(new PointRoi(bm2d.p3[0]+0.5, bm2d.p3[1]+0.5));

            imShow.setOverlay(ovly);
            imShow.show();

            isOut.addSlice("", ip_to_add);
        }

        new ImagePlus("", isOut).show();

        System.out.println("done");

/*		// parameters
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
		BranchModel3D bm = new BranchModel3D(stack_height, stack_width, stack_size, radius_std, 40, 150);
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
						""+out_path_v);*/

		
	}

    // terminal command option (empty)
    public static void main(
            String[] args
    )
    {
        // terminal command
//		// parameters
//		int 	stack_height 	= 64;
//		int 	stack_width 	= 64;
//		int 	stack_size 		= 32;
//		double 	radius_std		= 1.0;
//
//		String out_path_h =
//				System.getProperty("user.home")+File.separator+"gen_branch_hor.tif";
//		String out_path_v =
//				System.getProperty("user.home")+File.separator+"gen_branch_ver.tif";
//
//		if(args.length==4){
//			stack_height 	= (int) 	Integer.parseInt(	args[0]);
//			stack_width  	= (int) 	Integer.parseInt(	args[1]);
//			stack_size		= (int) 	Integer.parseInt(	args[2]);
//			radius_std		= (double)	Double.parseDouble( args[3]);
//		}
//		else{
//			System.out.println("Enter again, wrong arguments...\n" +
//					"1 - height \n" +
//					"2 - width  \n" +
//					"3 - size \n" +
//					"4 - radius std. \n");
//			return;
//		}
//
//		BranchModel3D bm = new BranchModel3D(stack_height, stack_width, stack_size, radius_std, 40, 150);
//
//		bm.drawHorizontalModel();
//		ImagePlus im_out = bm.getModelAsImage();
//		(new FileSaver(im_out)).saveAsTiffStack(out_path_h);
//
//		bm.drawVerticalModel();
//		im_out = bm.getModelAsImage();
//		(new FileSaver(im_out)).saveAsTiffStack(out_path_v);
//
//		System.out.println("files exported in :\n" +
//				""+out_path_h+" and \n" +
//						""+out_path_v);
    }

}

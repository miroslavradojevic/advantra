package profile;

import java.io.File;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class ConnectedRegions  implements PlugInFilter {

	/*
	 * uses Find_Connected_Regions class
	 */

	ImagePlus img;

	public static void main(String[] args){

		if(args.length!=1) {
			System.out.println(
									  "/* finding connected regions */ \n" +
											  "/* with the same gray-level value */ \n" +
											  "/* recommendable not more than several gray-levels (computation time and usage reasons) */ \n" +
											  "usage: set image path or name as argument\n"
			);
			return;
		}

		if(!(new File(args[0])).exists() || (new File(args[0])).isDirectory()){
			System.out.println("there is no file named ' "+args[0]+" '");
			return;
		}

		ImagePlus img = new ImagePlus((new File(args[0]).getAbsolutePath()));

		Find_Connected_Regions conn_reg = new Find_Connected_Regions(img, true);
		conn_reg.run("");

		System.out.println(conn_reg.getNrConnectedRegions()+" connected regions extracted...");
		ImagePlus regs = conn_reg.showLabels();

		String out_name = "labeled_regions.tif";

		if(img.getStack().getSize()>1){
			new FileSaver(regs).saveAsTiffStack(out_name);
		}
		else{
			new FileSaver(regs).saveAsTiff(out_name);
		}

		System.out.println((new File(out_name)).getAbsolutePath()+" saved...");

	}

	public void run(ImageProcessor arg0) {

		// dialog to enter input values
		GenericDialog gd = new GenericDialog("Connected Regions", IJ.getInstance());

		gd.addMessage(
							 "Plugin will extract connected regions that contain the same GRAY8 value\n" +
									 "it is recommendable that the image does not contain too many intensities."
		);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		long t1 = System.currentTimeMillis();
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(img, true);
		conn_reg.run("");
		long t2 = System.currentTimeMillis();

		int nr_regions = conn_reg.getNrConnectedRegions();
		IJ.log(nr_regions+" connected regions extracted.\n" +
					   "elapsed: "+((t2-t1)/1000f)+ " seconds.");

		conn_reg.showLabels().show();

	}

	public int setup(String arg0, ImagePlus img) {

		this.img = img;// ImageConversions.ImagePlusToGray8(img);
		return DOES_8G+NO_CHANGES; //+DOES_RGB+

	}




}

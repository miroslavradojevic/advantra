package advantra.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
//import ij.gui.NewImage;
//import ij.io.FileSaver;
import ij.plugin.PlugIn;
//import ij.process.ColorProcessor;
//import ij.process.ImageProcessor;

public class StackToPng implements PlugIn {
	
	public void run(String arg){
		
		System.out.println("converting to png...");
		
		IJ.run("Close All", "");
		
		//IJ.run("Image Sequence...", "convert_to_rgb");//, "open=/home/miroslav/Desktop/stack/01.tif number=71 starting=1 increment=1 scale=100 file=[] or=[] sort");
		//ImagePlus 		stack_layer_img = NewImage.createRGBImage("BlankImage", img_stack.getWidth(), img_stack.getHeight(), stack_len, NewImage.FILL_BLACK);
		//byte[] stack_layer = new byte[img_stack.getHeight()*img_stack.getWidth()];
				
		IJ.run("Image Sequence...", "open=[] number=[] starting=1 increment=1 scale=[] file=[] or=[] convert_to_rgb sort");
		
		ImagePlus img_stack 	= new ImagePlus("InputStack", IJ.getImage().getStack());
		int stack_len 			= img_stack.getStackSize();
		
		String out_folder  = 	"/home/miroslav/";
		String file_prefix = 	"stack";
		GenericDialog gd = new GenericDialog("Conversion");
		gd.addStringField("output folder: ", out_folder);
		gd.addStringField("file prefix  : ", file_prefix);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		out_folder = gd.getNextString();
		file_prefix = gd.getNextString();
	    
		System.out.println("stack size   : "+stack_len);
		System.out.println("output folder: "+out_folder);
		System.out.println("prefix       : "+file_prefix);
		System.out.println("starting...");
		
//		ColorProcessor stack_layer_proc;
//		for (int i = 0; i < stack_len; i++) {
//			System.out.println("converting stack "+(i+1)+"/"+(stack_len)+"...");
//			//extract layer
//			stack_layer_proc 	= (ColorProcessor) img_stack.getStack().getProcessor(i+1).convertToRGB();
//			stack_layer_img.setProcessor(stack_layer_proc);
//		}
		
		img_stack.show();
		
		/*
		//save stack
		String output_destination = out_folder+file_prefix+".tif";
		System.out.println("saving... "+output_destination);
		
		FileSaver fs = new FileSaver(img_stack);
		
		if( fs.saveAsTiff(output_destination) ){
			System.out.println("saved...");
		}
		else{
			System.out.println("problem...");
		}
		*/
	}

}

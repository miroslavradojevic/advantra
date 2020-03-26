package advantra.plugins;

	import ij.IJ;
	import ij.ImagePlus;
	import ij.gui.Roi;
	import ij.plugin.filter.PlugInFilter;
	import ij.plugin.frame.RoiManager;
	import ij.process.ByteProcessor;
	import ij.process.ImageProcessor;

	public class CreateMaskFromROI implements PlugInFilter {

	private ImagePlus imp;

	public int setup(String arg, ImagePlus imp) {
	    this.imp = imp;
	    if (imp == null) {
	        IJ.error("Input image required!");
	        return DONE;
	    }
	    
	    if (RoiManager.getInstance() == null) {
	        IJ.error("Roi Manager is not opened!");
	        return DONE;
	    }

	    return DOES_8G + DOES_16 + DOES_32 + NO_CHANGES;
	}

	public void run(ImageProcessor ip) {

//		String filename = imp.getTitle();
//	    String filenamebase = filename.substring(0, filename.length() - 4);
//	    String dirname = imp.getOriginalFileInfo().directory;
//	    String newpath = dirname + filenamebase;

	    RoiManager rm = RoiManager.getInstance();
	    Roi roi [] = rm.getRoisAsArray();
	    
	    ByteProcessor ipb = new ByteProcessor(imp.getWidth(), imp.getHeight());
	    for (int i = 0; i < imp.getWidth(); i++) {
			for (int j = 0; j < imp.getHeight(); j++) {
			    for (int ri = 0; ri < roi.length; ri++) {
					Roi r = roi[ri];
					if (r.getType() < 4) {
						if (r.contains(i, j)) {
							ipb.set(i, j, 255);
						}
					}
				}
				
			}
	    }

	    ImagePlus impb = new ImagePlus("mask", ipb);
	    impb.show();

	    IJ.setTool("hand");
	}
}
	

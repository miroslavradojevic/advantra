package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/7/13
 * Time: 8:13 PM
 */
public class ProfilerDemo implements PlugInFilter {

	ImagePlus 	inimg;
	String 		inimgPath;
	ImagePlus   inmask;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) {
			IJ.showMessage("needs image opened"); return DONE; }
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;

	}

	public void run(ImageProcessor imageProcessor) {


		String inmaskPath = Tools.removeExtension(inimgPath)+".mask";

		if (new File(inmaskPath).exists()) {
			inmask = new ImagePlus(inmaskPath);
			inmask.setTitle("inmask");
		}
		else {
			byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
			for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
			inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
		}

		inmask.show();

		Profiler pf = new Profiler(0, 5);
		Profiler.load(inimg.getProcessor(), (ByteProcessor) inmask.getProcessor(), 4, 2);
//		IJ.log("angles: "+Profiler.offsets.size()+" , angularStep="+Profiler.resolDeg);



	}
}

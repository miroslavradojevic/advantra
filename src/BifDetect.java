import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/22/13
 * Time: 12:21 PM
 */
public class BifDetect implements PlugInFilter {

	ImagePlus 	inimg;
	String      inimgPath;
	ImagePlus   inmask;
	Feat f;

	public void run(ImageProcessor imageProcessor) {

//		inmaskPath = MyOpener.open("Open mask file", false);
//		if (inmaskPath==null) {
//			IJ.log(inmaskPath);
//			return;
//		}

		int     t       		= 4;
		double  scale   		= 2.0;
		double  D       		= 10;
		double  E       		= 20;
		String 	inmaskPath		= Tools.removeExtension(inimgPath)+".mask";
		boolean useMask 		= true;

		t     	= (int) Prefs.get("advantra.critpoint.neuron_diam", t);
		scale   = Prefs.get("advantra.critpoint.scale", scale);

		GenericDialog gd = new GenericDialog("Fit Features");
		gd.addMessage("feature parameters");
		gd.addNumericField("neuron diameter min", t, 0, 5, "pix");
		gd.addNumericField("n'hood", scale, 1, 5, "x diameter");

        gd.addMessage("detection parameters");
		gd.addNumericField("D", D, 1, 5, "");
		gd.addNumericField("E", E, 1, 5, "");

        gd.addMessage("detection parameters - configuration file");

        gd.addMessage("mask (avoid processing all)");
		gd.addStringField("mask path", inmaskPath, 50);
		gd.addCheckbox("", useMask);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		t 		=  	(int)gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuron_diam", 	t);
		scale   =   gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale", 	scale);
		D   	=   gd.getNextNumber();
		E   	=   gd.getNextNumber();
		inmaskPath = gd.getNextString();
		useMask = gd.getNextBoolean();

		f= new Feat(t, scale);

		inimg.show();

		if (new File(inmaskPath).exists() && useMask) {
			inmask = new ImagePlus(inmaskPath);
			inmask.setTitle("inmask");
		}
		else {
			byte[] array = new byte[inimg.getWidth() * inimg.getHeight()];
			for (int a = 0; a < array.length; a++) array[a] = (byte) 255;
			inmask = new ImagePlus("inmask", new ByteProcessor(inimg.getWidth(), inimg.getHeight(), array));
		}

//		inmask.show();
//		IJ.selectWindow("inimg");
//		IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
//		inmask.close();

		IJ.log("filtering...");
		long t1 = System.currentTimeMillis();
		new ImagePlus("detection,nd="+f.diam+","+"r="+f.r+",D="+D+",E="+E, filterWithMask(inmask.getProcessor(), (FloatProcessor) inimg.getProcessor(), D, E)).show();
		long t2 = System.currentTimeMillis();
		IJ.log("done. "+((t2-t1)/1000f)+" sec.");

	}

	public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) {
			IJ.showMessage("needs opened image"); return DONE;
		}
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;
	}

	public FloatProcessor filterWithMask(ImageProcessor msk, FloatProcessor input, double D, double E)
	{

		FloatProcessor ipOut = new FloatProcessor(input.getWidth(), input.getHeight());

		for (int x=0; x<input.getWidth(); x++) {
			for (int y=0; y<input.getHeight(); y++) {
				if (msk.getf(x, y)==255) {
					double sc = f.bifurcationess1(x, y, input, D, E);
					ipOut.setf(x, y, (float)sc );
				}
			}
		}
		return ipOut;

	}

}

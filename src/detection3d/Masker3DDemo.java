package detection3d;

import detection.Masker;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.Arrays;

/**
 * Input image is 3D ImageStack
 *
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/7/13
 * Time: 2:46 PM
 */
public class Masker3DDemo implements PlugInFilter {

	ImagePlus inimg;    // with float processor image stack
	ImagePlus inmask;   // with byte processor image stack

	/*
    parameters
     */
	float       nhoodRadius, iDiff, zDist;
	int         CPU_NR;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		if(imagePlus.getStackSize()<=1) return DONE;

		//inimg = toFloatImage(imagePlus);
		inimg = imagePlus;

		return DOES_8G+DOES_32+NO_CHANGES;
	}

//	public void run (ImageProcessor imageProcessor) {
//		for (int i=0; i<100; i++) {
//			IJ.log("\\Update:"+"line "+IJ.d2s(i, 0));
//			try {
//				Thread.sleep(1000);
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//		}
//	}

	public void run(ImageProcessor imageProcessor) {

		/*
        	Generic Dialog
         */
		nhoodRadius             = (float)   Prefs.get("critpoint.mask.nhoodRadius", 5);
		iDiff 					= (float)   Prefs.get("critpoint.mask.iDiff", 5);
		zDist					= (float)	Prefs.get("critpoint.mask.zDist", 1);

		GenericDialog gd = new GenericDialog("MASK EXTRACTOR (3D)");
		gd.addNumericField("radius ", 	nhoodRadius, 	1, 10, "spatial neighbourhood");
		gd.addNumericField("iDiff ", 	iDiff, 			0, 10, "intensity margin");
		gd.addNumericField("zDist ", 	zDist, 			2, 10, "z distance [pix]");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		nhoodRadius      = (float) gd.getNextNumber();
		Prefs.set("critpoint.mask.nhoodRadius",	nhoodRadius);

		iDiff       	= (float) gd.getNextNumber();
		Prefs.set("critpoint.mask.iDiff", 	    	iDiff);

		zDist       	= (float) gd.getNextNumber();
		Prefs.set("critpoint.mask.zDist", 	    	zDist);

		CPU_NR = Runtime.getRuntime().availableProcessors();

        /*
        	main
         */

		IJ.log("extracting background...    ");

		long t1 = System.currentTimeMillis();

		Masker3D.loadTemplate(inimg.getStack(), zDist, nhoodRadius, iDiff);

		int totalJobs = inimg.getHeight()*inimg.getWidth()*inimg.getStackSize();

		Masker3D ms_jobs[] = new Masker3D[CPU_NR];
		for (int i = 0; i < ms_jobs.length; i++) {
			ms_jobs[i] = new Masker3D(i*totalJobs/CPU_NR,  (i+1)*totalJobs/CPU_NR);
			ms_jobs[i].start();
		}
		for (int i = 0; i < ms_jobs.length; i++) {
			try {
				ms_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		long t2 = System.currentTimeMillis();
		int fgLocations = Masker3D.getNrLocations();
		IJ.log("done "+fgLocations+" locations ("+((float)fgLocations/(Masker3D.image_width*Masker3D.image_height*Masker3D.image_length))+" %) found in "+((t2-t1)/1000f)+" sec.");

		inmask = new ImagePlus("inmask", Masker3D.getMask());
		inmask.show();


//		int size = Masker3D.sizeCircularNbhood(nhoodRadius);

//		float[] storage = new float[size];

		//IJ.log("before" + Arrays.toString(storage));

//		System.out.println();

//		// per point
//		for (int z=0; z<inimg.getStackSize(); z++) {
//
//			System.out.print("\rlayer "+IJ.d2s(z,0)+" / "+(inimg.getStackSize()-1));
//
//			long t1, t2;
//
//			t1 = System.currentTimeMillis();
//
//			for (int x=0; x<inimg.getWidth(); x++) {
//				for (int y=0; y<inimg.getHeight(); y++) {
//
//					Masker3D.extractCircularNbhood(x, y, z, nhoodRadius, storage);
//					Masker3D.median(storage);
//
//				}
//			}
//
//			t2 = System.currentTimeMillis();
//
//			System.out.print("\t\t "+((t2-t1)/1000f)+" sec.");
//
//		}
//		System.out.println();

		//IJ.log("after" + Arrays.toString(storage));

		//new ImagePlus("loaded image", Masker3D.inis).show();

	}

	private ImagePlus toFloatImage(ImagePlus inputImage){

		int W = inputImage.getWidth();
		int H = inputImage.getHeight();
		int L = inputImage.getStackSize();

		ImageStack outIs = new ImageStack(W, H);

		// convert every element of the input to float
		for (int l=1; l<L; l++) {

			// convertToFloat will convert any ImageProcessor to FloatProcessor
			float[] floatArray = (float[]) inputImage.getStack().getProcessor(l).convertToFloat().getPixels();
			float[] floatArrayNew = new float[floatArray.length];

			for (int i=0; i<floatArray.length; i++) floatArrayNew[i] = floatArray[i];

			FloatProcessor outFp = new FloatProcessor(W, H, floatArrayNew);
			outIs.addSlice(outFp);

		}

		ImagePlus outImg = new ImagePlus("inputFloatImage", outIs);
		return outImg;

	}

}

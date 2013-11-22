package detection3d;

import ij.ImagePlus;
import ij.ImageStack;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/22/13
 * Time: 5:18 AM
 */
public class Detector3D {

	// inpput 3d image


	Sphere3D sph3D; // will be tool of thread classes - used to make geometric computations


	public Detector3D() {

		// store parameters


	}

	public void run(ImagePlus inimg, double neuronDiameter, float iDiff1) { // ArrayList<ArrayList<int[]>>



	}

	/*
		convert image to the array form (float[][][]) that will be used by all thread classes (Masker3D, Profiler3D, PeakExtractor, Analyzer3D) in run()
	 */
	private static float[][][] stackToZxyArray(ImageStack inis) {

		int W = inis.getWidth();
		int H = inis.getHeight();
		int L = inis.getSize();

		float[][][] img3d_zxy = new float[L][][];

		for (int l=0; l<L; l++) {

			img3d_zxy[l] = new float[W][H];
			float[] readSlice = (float[]) inis.getProcessor(l+1).convertToFloat().getPixels();

			for (int ww=0; ww<W; ww++) {
				for (int hh=0; hh<H; hh++) {
					img3d_zxy[l][ww][hh] = readSlice[hh*W+ww];
				}
			}

		}

		return img3d_zxy;
	}

}

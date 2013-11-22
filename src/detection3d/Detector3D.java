package detection3d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/22/13
 * Time: 5:18 AM
 */
public class Detector3D {

    float[][][] img3d_zxy;  // inpput 3d image in convenient form (z,x,y)
    float zDist;

	Sphere3D sph3D;         // will be tool of thread classes - used to make geometric computations

    byte [][][] foregroundLocations_zxy;    // = new byte[1][1][1];         // runMasker() output
    float[][][] backgroundEst_zxy;          // = new float[1][1][1];        // runMasker() output

    double 	D;
    float 	iDiff;
    float 	r;
    float 	scale = 1.5f;                   // heuristics: 1.5 is resonable

    public static float 	MIN_COS_ANG = .6f;
    public static float 	MIN_FUZZY_SCORE = .7f;
    public static float 	SCATTER_D2 	= 5;
    public static int		W_STD_RATIO_WRT_TO_D = 3;

    int CPU_NR;

    public Detector3D() {
        System.out.println("loading default parameters...");
    }

	public Detector3D( // constructor initializes the parameter, run() calls different images with predefined parameters
            float   minCosAng,
            float   minFuzzyScore,
            float   scatterDistance2,
            int     wStdRatioToD
    ) {

		// store parameters
        MIN_COS_ANG = minCosAng;
        MIN_FUZZY_SCORE = minFuzzyScore;
        SCATTER_D2 = scatterDistance2;
        W_STD_RATIO_WRT_TO_D = wStdRatioToD;

    }

	public void run(ImagePlus inimg, float zDist1, double neuronDiameter, float iDiff1) { // ArrayList<ArrayList<int[]>>

        // load image into float[][][] structure
        if(inimg==null || inimg.getStackSize()<=1) return ;

        this.img3d_zxy = stackToZxyArray(inimg.getStack());
        this.zDist = zDist1;
        this.D = neuronDiameter;
        this.iDiff = iDiff1;
        this.r = (float) (scale * D);

        // some important parameters
        this.CPU_NR = Runtime.getRuntime().availableProcessors();


        /*
            main processing
         */
//        int 	totalJobs; // paralelization variable used by Thread modules to spread the tasks

        /*
            general pipeline is:

            1 find mask containing foreground voxels    (Masker3D)

            2 extract profiles in 3d for each location  (Profiler3D)

            3 find peaks of the profiles and correlate them to subpixel locations in 3d ()

            4 analyze recursive lists starting from each location
         */

        // 1. MASK
        System.out.print("separating foreground...");
        int L = img3d_zxy.length;
        int W = img3d_zxy[0].length;
        int H = img3d_zxy[0][0].length;
        foregroundLocations_zxy     = new byte[L][W][H];
        backgroundEst_zxy           = new float[L][W][H];

        runMasker(  img3d_zxy, zDist, 1.5f*r, iDiff,
                    foregroundLocations_zxy, backgroundEst_zxy); // Masker3D threading within method call
        System.out.println("done.");

        //System.out.println("new "+backgroundEst_zxy.length+" , "+backgroundEst_zxy[0].length+" , "+backgroundEst_zxy[0][0].length);

    }

    /*
        run threaded masker within the method
     */
    private void runMasker(
            float[][][] img3d_zxy, float zDist, float radiusMask, float iDiff,
            byte [][][] foregroundLocations_zxy, float[][][] backgroundEst_zxy)
    {

//        IJ.log("radius="+radiusMask);

        long t1 = System.currentTimeMillis();

        Masker3D.loadTemplate(img3d_zxy, zDist, radiusMask, iDiff);

        int totalJobs = img3d_zxy.length * img3d_zxy[0].length * img3d_zxy[0][0].length; // inimg.getHeight()*inimg.getWidth()*inimg.getStackSize();
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
        System.out.println("done " + fgLocations + " locations (" + (((float) fgLocations / (Masker3D.image_width * Masker3D.image_height * Masker3D.image_length))  * 100) + " %) found in " + ((t2 - t1) / 1000f) + " sec.");

        new ImagePlus("", Masker3D.getMask()).show();


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

    // TODO add method that converts opposite way float[][][] to ImageStack

}

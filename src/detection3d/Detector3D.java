package detection3d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;

import java.util.Arrays;

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

//    byte [][][] foregroundLocations_zxy;    // = new byte[1][1][1];         // runMasker() output
//    float[][][] backgroundEst_zxy;          // = new float[1][1][1];        // runMasker() output

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
            general pipeline is:

            1 find mask containing foreground voxels    (Masker3D)

            2 extract profiles in 3d for each location  (Profiler3D)

            3 find peaks of the profiles and correlate them to subpixel locations in 3d ()

            4 analyze recursive lists starting from each location
         */

        // 1. MASK
        System.out.print("separating foreground... ");
		long t1 = System.currentTimeMillis();

		MaskerOutput mo = new MaskerOutput(); // output is stored here

		runMasker(img3d_zxy, zDist, 1.5f * r, iDiff, mo);

		long t2 = System.currentTimeMillis();

		System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");

		// 2. PROFILE


		// 3. PROFILE PEAK(S)



		// 4. ASSOCIATE PEAKS



	}

    /*
        run threaded masker within the method
     */
    private void runMasker(float[][][] img3d_zxy, float zDist, float radiusMask, float iDiff, MaskerOutput mo)
    {

		/*
			following bit will execute run() in parallel
		 */

        Masker3D.loadTemplate(img3d_zxy, zDist, radiusMask, iDiff);

		int Zdim = img3d_zxy.length;
		int Xdim = img3d_zxy[0].length;
		int Ydim = img3d_zxy[0][0].length;

        int totalJobs = Zdim * Xdim * Ydim;
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

		/*
		 	skipped locations are filled,
		  */
		Masker3D.fill();

		/*
		 	store locations, estimated background and binary mask them in output class (take from Masker3D into MaskerOutput instance
		  */
		mo.isForeground = new boolean[Zdim][Xdim][Ydim];
		mo.locIndex = new int[Zdim][Xdim][Ydim];
		int cnt = 0;
		for (int zz=0; zz<Zdim; zz++) {
			for (int xx=0; xx<Xdim; xx++) {
				for (int yy=0; yy<Ydim; yy++) {

					mo.isForeground[zz][xx][yy] = Masker3D.mask3[zz][xx][yy];

					if (Masker3D.mask3[zz][xx][yy]) {
						mo.locIndex[zz][xx][yy] = cnt;
						cnt++;
					}
					else {
						mo.locIndex[zz][xx][yy] = -1;
					}

				}
			}
		}
		// align
		mo.foregroundLocsZXY = new int[cnt][3];
		mo.backgroundEst = new float[cnt];

		cnt =0;
		for (int zz=0; zz<Zdim; zz++) {
			for (int xx=0; xx<Xdim; xx++) {
				for (int yy=0; yy<Ydim; yy++) {

					if (Masker3D.mask3[zz][xx][yy]) {

						mo.foregroundLocsZXY[cnt][0] = zz;
						mo.foregroundLocsZXY[cnt][1] = xx;
						mo.foregroundLocsZXY[cnt][2] = yy;

						mo.backgroundEst[cnt] = Masker3D.back3[zz][xx][yy];

						cnt++;

					}

				}
			}
		}


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

class MaskerOutput { // passed as an argument to void method

	int[][] 		foregroundLocsZXY 		= null; // nr loc x 3  (list of foreground locations, for paralellization later)
	float[] 		backgroundEst	= null; 		// nr. loc x 1, list of background estimates at foreground locations
	boolean[][][] 	isForeground 		= null; 	// for easier access when checking whether the location belongs to foreground
	int[][][]       locIndex = null; 				// for easier access when finding an index of some location or background value

}

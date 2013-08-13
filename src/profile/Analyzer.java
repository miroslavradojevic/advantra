package profile;

import ij.IJ;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/8/13
 * Time: 6:33 PM
 */
public class Analyzer extends Thread  {

    private int begN, endN;

	/*
		shared variables are static
	 */
	// input
    public static ArrayList<ArrayList<float[]>> profiles;

    // outputs
    public static float[][][] peakIdx; // array index

    // this is functionally not necessary - just to visualize how it converged (it's calculated anyway)
    public static ArrayList<ArrayList<float[]>> convIdx;

    public static int       nrPoints = 100;
    double[]  start 		= new double[nrPoints];
    double[]  msFinish 		= new double[nrPoints];

    public static int       maxIter = 100;
    public static double    epsilon = 0.00001;//Double.MIN_VALUE;//
    public static int       h = 2;              // in indexes
    public static double    minD = 0.5;
    public static int       M = 2;

    public Analyzer (int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadProfiles(ArrayList<ArrayList<float[]>> profiles1){

        /*
        profiles
         */
		profiles = new ArrayList<ArrayList<float[]>>(profiles1.size());
		for (int i=0; i<profiles1.size(); i++) {
            ArrayList<float[]> temp = new ArrayList<float[]>(profiles1.get(i).size());
            for (int j=0; j<profiles1.get(i).size(); j++) {

                float[] toAdd = new float[profiles1.get(i).get(j).length];

                for (int k=0; k<toAdd.length; k++ ) {
                    toAdd[k] = profiles1.get(i).get(j)[k];
                }

                temp.add(toAdd);

            }
            profiles.add(temp);
        }

		/*
		peakIdx
		 */
        peakIdx = new float[profiles.size()][profiles.get(0).size()][]; // keeps three corresp indexes of profile values

		/*
		convIdx
		 */
		convIdx = new ArrayList<ArrayList<float[]>>(); // this is where msFinish will be stored
        for (int i=0; i<profiles1.size(); i++) {
            ArrayList<float[]> temp = new ArrayList<float[]>(profiles1.get(i).size());
            for (int j=0; j<profiles1.get(i).size(); j++) {

                float[] toAdd = new float[nrPoints];
                temp.add(toAdd);

            }
            convIdx.add(temp);
        }

    }

	public static float[] extractPeakIdxs(float[] profile1)
	{

        double[] start11 	= new double[nrPoints];
		double[] finish11 	= new double[nrPoints];

	    int profileLength = profile1.length;

		for (int k=0; k<nrPoints; k++) {
		    start11[k] = ((float) k / nrPoints) * profileLength;
		}

		Tools.runMS(start11,
				    profile1,
					maxIter,
					epsilon,
					h,
					finish11);

	    Vector<float[]> cls = Tools.extractClusters(finish11, minD, M); // TODO: use extractClusters1() here

        if (cls.size()==0) return null;

        else if (cls.size()==1) return new float[]{cls.get(0)[0]};

        else if (cls.size()==2) return new float[]{cls.get(0)[0], cls.get(1)[0]};

        else if (cls.size()==3) return new float[]{cls.get(0)[0], cls.get(1)[0], cls.get(2)[0]};

        else if (cls.size()==4) return new float[]{cls.get(0)[0], cls.get(1)[0], cls.get(2)[0], cls.get(3)[0]};

        else return bestN(cls, 4);

	}

    private static float[] bestN(Vector<float[]> cls1, int N){

        // cls1.size should be more or equal than N
        boolean[] checked = new boolean[cls1.size()];
        float[] out = new float[N];

        for (int k = 0; k<N; k++) {
            // reset max search
            double  currMax = Double.MIN_VALUE;
            int     currMaxIdx = -1;

            for (int i=0; i<cls1.size(); i++) {

                // find max in this round
                if (!checked[i]) {
                    if (cls1.get(i)[1]>currMax) {

                        currMax = cls1.get(i)[1];
                        currMaxIdx = i;

                    }
                }
            }

            checked[currMaxIdx] = true;
            // set the output value
            out[k] = cls1.get(currMaxIdx)[0];

        }

        return out;
    }

    public void run(){ // considers begN and endN

        for (int locIdx = begN; locIdx < endN; locIdx++) {

            // do mean-shift for profiles at these locations

            for (int profileIdx=0; profileIdx<profiles.get(locIdx).size(); profileIdx++) {

                // access the profile
                int profileLength = profiles.get(locIdx).get(profileIdx).length;

				for (int k=0; k<nrPoints; k++) {
                    start[k] = ((float) k / nrPoints) * profileLength;
                }

				Tools.runMS(  	start,
                              	profiles.get(locIdx).get(profileIdx),
                              	maxIter,
                              	epsilon,
                              	h,
								msFinish);

                for (int i1=0; i1<nrPoints; i1++) {
                    convIdx.get(locIdx).get(profileIdx)[i1] = (float) msFinish[i1];
                }

                int inputProfileLength = profiles.get(locIdx).get(profileIdx).length;
                Vector<float[]> cls = Tools.extractClusters1(msFinish, minD, M, inputProfileLength);
				extractPeakIdx(cls, locIdx, profileIdx); // to extract 3 major ones (if there are three)

            }

        }

    }

	private static void extractPeakIdx(Vector<float[]> cls, int locIdx, int profileIdx)  // will update peakIdx static variable
	{
		// store the values
		if (cls.size()<=2) {
			//peakIdx[locIdx][profileIdx] = null;
			for (int a1=profileIdx; a1<profiles.get(locIdx).size(); a1++) {
				peakIdx[locIdx][profileIdx] = null;
			}

			//break; // break the loop for the rest of configurations

		}
		else if (cls.size()==3) {
			peakIdx[locIdx][profileIdx] = new float[3];
			peakIdx[locIdx][profileIdx][0]  = cls.get(0)[0];
			peakIdx[locIdx][profileIdx][1]  = cls.get(1)[0];
			peakIdx[locIdx][profileIdx][2]  = cls.get(2)[0];
		}
		else { // >3

			boolean[] checked = new boolean[cls.size()]; // all to false

			// extract 3 angles with most convergence points
			// cls.get(i)[0] - convergence point
			// cls.get(i)[1] - nr. points
			peakIdx[locIdx][profileIdx] = new float[3];

			// find top 3
			for (int k = 0; k<3; k++) {
				// reset max search
				double  currMax = Double.MIN_VALUE;
				int     currMaxIdx = -1;

				for (int i=0; i<cls.size(); i++) {

					// find max in this round
					if (!checked[i]) {
						if (cls.get(i)[1]>currMax) {

							currMax = cls.get(i)[1];
							currMaxIdx = i;

						}
					}
				}

				checked[currMaxIdx] = true;
				// set the output value
				peakIdx[locIdx][profileIdx][k] = cls.get(currMaxIdx)[0];

			}

		}

	}


}
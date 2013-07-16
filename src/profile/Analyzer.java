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

    public static int       maxIter = 200;
    public static double    epsilon = 0.0001;
    public static int       h = 4;
    public static double    minD = 0.5;
    public static int       M = 1;

    public Analyzer (int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadProfiles(ArrayList<ArrayList<float[]>> profiles1){

        // read & store profiles

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

        peakIdx = new float[profiles.size()][profiles.get(0).size()][]; // keeps three corresp indexes of profile values

		convIdx = new ArrayList<ArrayList<float[]>>(); // this is where msFinish will be stored

        for (int i=0; i<profiles1.size(); i++) {
            ArrayList<float[]> temp = new ArrayList<float[]>(profiles1.get(i).size());
            for (int j=0; j<profiles1.get(i).size(); j++) {

                float[] toAdd = new float[100];

//                for (int k=0; k<toAdd.length; k++ ) {
//                    toAdd[k] = profiles1.get(i).get(j)[k];
//                }

                temp.add(toAdd);

            }
            convIdx.add(temp);
        }

    }

    public void run(){ // considers begN and endN

        for (int locIdx = begN; locIdx < endN; locIdx++) {

            // do mean-shift for profiles at these locations

            for (int profileIdx=0; profileIdx<profiles.get(locIdx).size(); profileIdx++) {

                // access the profile
                int profileLength = profiles.get(locIdx).get(profileIdx).length;

				for (int k=0; k<nrPoints; k++) { // profileLength
                    start[k] = ((float) k / nrPoints) * profileLength;
                }

				Tools.runMS(  	start,
                              	profiles.get(locIdx).get(profileIdx),
                              	maxIter,
                              	epsilon,
                              	h,
								msFinish);

				// loop to store it in convIdx

                for (int i1=0; i1<100; i1++) {
                    convIdx.get(locIdx).get(profileIdx)[i1] = (float) msFinish[i1];
                }

                Vector<float[]> cls = Tools.extractClusters(msFinish, minD, M);

                // store the values
                if (cls.size()<=2) {
//                    peakVal[locIdx][profileIdx] = null;
                    peakIdx[locIdx][profileIdx] = null;

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

    }


}

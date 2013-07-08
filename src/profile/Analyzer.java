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

    public static ArrayList<ArrayList<float[]>> profiles; // input
    // outputs
    public static float[][][] peakVal; //
    public static float[][][] peakIdx; // array index

    public static int       maxIter = 150;
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

        peakVal = new float[profiles.size()][profiles.get(0).size()][]; // keeps three peak values
        peakIdx = new float[profiles.size()][profiles.get(0).size()][]; // keeps three corresp indexes

    }

    public void run(){ // considers begN and endN

        for (int locIdx = begN; locIdx < endN; locIdx++) {

            // do mean-shift for profiles at these locations
            // will fill
            // peakVal[locIdx][][3]
            // peakIdx[locIdx][][3]

            for (int profileIdx=0; profileIdx<profiles.get(locIdx).size(); profileIdx++) {

                // access the profile
                int profileLength = profiles.get(locIdx).get(profileIdx).length;
                double[] start = new double[profileLength];
                for (int k=0; k<profileLength; k++) {
                    start[k] = k;
                }

                double[] finish = Tools.runMS(  start,
                                                profiles.get(locIdx).get(profileIdx),
                                                maxIter,
                                                epsilon,
                                                h);

                Vector<float[]> cls = Tools.extractClusters(finish, minD, M);

                // store the values
                if (cls.size()<=2) {
                    peakVal[locIdx][profileIdx] = null;
                    peakIdx[locIdx][profileIdx] = null;

                }
                else if (cls.size()==3) {
                    peakVal[locIdx][profileIdx] = new float[3];
                    peakVal[locIdx][profileIdx][0]  = cls.get(0)[0];
                    peakVal[locIdx][profileIdx][1]  = cls.get(1)[0];
                    peakVal[locIdx][profileIdx][2]  = cls.get(2)[0];

                    peakIdx[locIdx][profileIdx] = new float[3];
                    peakIdx[locIdx][profileIdx][0]  = (float) Tools.interp1Darray(peakVal[locIdx][profileIdx][0], profiles.get(locIdx).get(profileIdx));
                    peakIdx[locIdx][profileIdx][1]  = (float) Tools.interp1Darray(peakVal[locIdx][profileIdx][1], profiles.get(locIdx).get(profileIdx));
                    peakIdx[locIdx][profileIdx][2]  = (float) Tools.interp1Darray(peakVal[locIdx][profileIdx][2], profiles.get(locIdx).get(profileIdx));

                }
                else { // >3

                    boolean[] checked = new boolean[cls.size()]; // all to false
                    // extract 3 angles with most convergence points

                    peakVal[locIdx][profileIdx] = new float[3];
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
                        peakVal[locIdx][profileIdx][k] = (float) Tools.interp1Darray(peakIdx[locIdx][profileIdx][k], profiles.get(locIdx).get(profileIdx));

                    }

                }

            }

        }

    }


}

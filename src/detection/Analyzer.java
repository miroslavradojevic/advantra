package detection;

import aux.Tools;

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

	// input
    public static ArrayList<ArrayList<float[]>> 	profiles;  	// nr. loc x nr.scales
	// aux
	public static ArrayList<double[]> 				startIdx; 	// nr. scales
	public static ArrayList<ArrayList<double[]>> 	finishIdx;  // nr. loc x nr.
    // output
    public static ArrayList<ArrayList<ArrayList<Float>>> peakIdx; // one detection will give ArrayList<Float>, different number of peaks can be detected

    public static int       nrPoints 	= 200;
    public static int       maxIter 	= 20;
    public static double    epsilon 	= 1e-8;
    public static int       h 			= 2;              						// in indexes
    public static double    minD 		= 0.5;     								// separation
    public static int       M = (int) Math.round(0.05*nrPoints);        		// 0.05 of the nrPoints

    public Analyzer (int n0, int n1)
	{
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadProfiles(ArrayList<ArrayList<float[]>> profiles1)
	{

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
        peakIdx = new ArrayList<ArrayList<ArrayList<Float>>>(profiles.size());
        for (int i=0; i<profiles.size(); i++) {

            ArrayList<ArrayList<Float>> profilesAtLoc = new ArrayList<ArrayList<Float>>(2); // ring A & B, two profiles

            for (int j=0; j<profiles.get(0).size(); j++) {

                ArrayList<Float> peaksForOneProfile = new ArrayList<Float>(4); // max 4
                profilesAtLoc.add(peaksForOneProfile);

            }

            peakIdx.add(profilesAtLoc);

        }

		/*
		finishIdx
		*/
		finishIdx = new ArrayList<ArrayList<double[]>>(profiles.size());
        for (int i=0; i<profiles1.size(); i++) {
            ArrayList<double[]> temp = new ArrayList<double[]>(profiles1.get(i).size());
            for (int j=0; j<profiles1.get(i).size(); j++) {

				double[] toAdd = new double[nrPoints];
                temp.add(toAdd);

            }
			finishIdx.add(temp);
        }

		/*
		startIdx
		 */
		startIdx = new ArrayList<double[]>(profiles.get(0).size());
		for (int j=0; j<profiles.get(0).size(); j++) {

			int profileLength = profiles.get(0).get(j).length;
			double[] startToAdd = new double[nrPoints];

			for (int k=0; k<nrPoints; k++) {
				startToAdd[k] = ((float) k / nrPoints) * profileLength;
			}

			startIdx.add(startToAdd);

		}

    }

	public static float[] extractPeakIdxs(float[] profile1, double[] startPts, double[] finishPts)
	{

		int convPoints = startPts.length;
	    int profileLength = profile1.length;

		for (int k=0; k<convPoints; k++) {
			startPts[k] = ((float) k / convPoints) * profileLength;
		}

		//Tools.runMeanShift(startPts, profile1, maxIter, epsilon, h,	finishPts);
		Tools.runMaxShift(startPts, profile1, maxIter, epsilon, h, finishPts);

	    Vector<float[]> cls = Tools.extractClusters1(finishPts, minD, M, profileLength);

        if (cls.size()==0) return null;

        else if (cls.size()==1) return new float[]{cls.get(0)[0]};

        else if (cls.size()==2) return new float[]{cls.get(0)[0], cls.get(1)[0]};

        else if (cls.size()==3) return new float[]{cls.get(0)[0], cls.get(1)[0], cls.get(2)[0]};

        else if (cls.size()==4) return new float[]{cls.get(0)[0], cls.get(1)[0], cls.get(2)[0], cls.get(3)[0]};

        else return bestN(cls, 4);

	}

	public static ArrayList<Float> extractPeakIdxsList(float[] profile1, double[] start1, double[] finish1)
	{

//		, double[] start11, double[] finish11
//		double[] start11 	= new double[nrPoints];
//		double[] finish11 	= new double[nrPoints];

		int profileLength = profile1.length;

//		for (int k=0; k<nrPoints; k++) {
//			start11[k] = ((float) k / nrPoints) * profileLength;
//		}

		Tools.runMaxShift(start1, profile1, maxIter, epsilon, h, finish1);

		Vector<float[]> cls = Tools.extractClusters1(finish1, minD, M, profileLength);

		if (cls.size()==0) return new ArrayList<Float>();

		else if (cls.size()==1) {
			ArrayList<Float> out = new ArrayList<Float>(1);
			out.add(0, cls.get(0)[0]);
			return out;
		}

		else if (cls.size()==2) {
            ArrayList<Float> out = new ArrayList<Float>(2);
            out.add(0, cls.get(0)[0]);
            out.add(1, cls.get(1)[0]);
            return out;
        }

		else if (cls.size()==3) {
            ArrayList<Float> out = new ArrayList<Float>(3);
            out.add(0, cls.get(0)[0]);
            out.add(1, cls.get(1)[0]);
            out.add(2, cls.get(2)[0]);
            return out;
        }

		else if (cls.size()==4) {
            ArrayList<Float> out = new ArrayList<Float>(4);
            out.add(0, cls.get(0)[0]);
            out.add(1, cls.get(1)[0]);
            out.add(2, cls.get(2)[0]);
            out.add(3, cls.get(3)[0]);
            return out;
        }

		else {
            return bestNList(cls, 4);
        }

	}

    private static float[] bestN(Vector<float[]> cls1, int N)
	{

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

	private static ArrayList<Float> bestNList(Vector<float[]> cls1, int N)
	{

		// cls1.size should be more or equal than N
		boolean[] checked = new boolean[cls1.size()];
		ArrayList<Float> out = new ArrayList<Float>(N);

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
			out.add(k, cls1.get(currMaxIdx)[0]);

		}

		return out;

	}

    public void run(){ // considers begN and endN

        for (int locIdx = begN; locIdx < endN; locIdx++) {

            // do mean-shift for profiles at these locations

            for (int profileIdx=0; profileIdx<profiles.get(locIdx).size(); profileIdx++) {

                // access the detection
//                int profileLength = profiles.get(locIdx).get(profileIdx).length;

                // calculate peaks for the ring 'profileIdx'
                ArrayList<Float> currPeaks = extractPeakIdxsList(profiles.get(locIdx).get(profileIdx), startIdx.get(profileIdx), finishIdx.get(locIdx).get(profileIdx));

                //if (currPeaks.size()<3) {
                //    // it is not a bifurcation according to MS for this ring, don't calculate further, leave empty fields of peakIdx at this location
                //    break;
                //}
                //else {
                    // add those points
                    for (int pp=0; pp<currPeaks.size(); pp++){
                        peakIdx.get(locIdx).get(profileIdx).add(pp, currPeaks.get(pp));
                    }
                //}

/*
				for (int k=0; k<nrPoints; k++) {
                    start[k] = ((float) k / nrPoints) * profileLength;
                }

				Tools.runMS(  	start,
                              	profiles.get(locIdx).get(profileIdx),
                              	maxIter,
                              	epsilon,
                              	h,
								msFinish);
*/

                /*
                for (int i1=0; i1<nrPoints; i1++) {
                    convIdx.get(locIdx).get(profileIdx)[i1] = (float) msFinish[i1];
                }
*/
/*
                int inputProfileLength = profiles.get(locIdx).get(profileIdx).length;
                Vector<float[]> cls = Tools.extractClusters1(msFinish, minD, M, inputProfileLength);
				extractPeakIdx(cls, locIdx, profileIdx); // to extract 3 major ones (if there are three)
*/

            }

        }

    }

	public static ArrayList<ArrayList<ArrayList<Float>>> exportPeakIdx()
	{

		ArrayList<ArrayList<ArrayList<Float>>> out = new ArrayList<ArrayList<ArrayList<Float>>>(peakIdx.size());
		for (int i=0; i<peakIdx.size(); i++) {

			ArrayList<ArrayList<Float>> out1 = new ArrayList<ArrayList<Float>>(peakIdx.get(i).size());

			for (int j=0; j<peakIdx.get(i).size(); j++) {

				ArrayList<Float> out2 = new ArrayList<Float>(peakIdx.get(i).get(j).size());

				for (int k=0; k<peakIdx.get(i).get(j).size(); k++) {

					out2.add(k, peakIdx.get(i).get(j).get(k));

				}

				out1.add(j, out2);

			}

			out.add(i, out1);
		}

		return out;

	}

}
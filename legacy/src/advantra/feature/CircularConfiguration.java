package advantra.feature;

import advantra.critpoint.Calc;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class CircularConfiguration {

    public int			nrPeaks;

	public int[] 		angResDeg;
	public double[]		angResRad;

	public int[] 		angBtwPeakDeg;
	public double[] 	angBtwPeakRad;
	
	public int 			nrRot;
	public double 		rotStepRad;
	
	public double 		minAngResRad;   // only used to define how many rotations there are
    public int          minAngResDeg;   // only used to define how many rotations there are

	public double[][] 	peaksRad;
	public double[]		ringR;          // array that defines the rings
                                        // ringR[0]<ON<ringR[1]
                                        // length 2 is enough for borders of the ring where
                                        // directions are filtered

    public double       innerRing;      // defines the inner ring

    public float[]      sumON;          // will be used for weighted score calculation
    public float[]      sumOFF;
    public int[]        nrON;
    public int[]        nrOFF;

    public static float TwoPi = (float) (2*Math.PI);

//	double[] 			filt;           // this filter will make a score
//                                        // it will have a different version with each rotation
//                                        // important to call initFilter() before using it
//                                        // before calling score() method
//                                        // to allocate enough values to cover the patch
	
	public CircularConfiguration(
            int[]       angular_resolution_deg,
            int[]       angles_between_peaks_deg,
			double[]    radius,
            double      innerRingRadius
    )
	{

        innerRing = innerRingRadius;

        ringR = new double[2];
		ringR[0] = radius[0];
        ringR[1] = radius[1];
		// this configuration has one ring (2 thresholds)
        // inner ring for (ringR[0]=0)
        // outer ring for (ringR[1]=1.0)

		// filter used when calculating score(), only to initialize
//		filt = new double[1];    // length will depend on th profile, same

		angResDeg = new int[angular_resolution_deg.length];
        for (int i = 0; i < angular_resolution_deg.length; i++){
            angResDeg[i] = angular_resolution_deg[i];
        }

		angResRad = new double[angular_resolution_deg.length];
        for (int i = 0; i < angular_resolution_deg.length; i++){
            angResRad[i] = (angResDeg[i]/360f)*TwoPi;
        }

		angBtwPeakDeg = new int[angles_between_peaks_deg.length];
		for (int i = 0; i < angles_between_peaks_deg.length; i++) {
			angBtwPeakDeg[i] = angles_between_peaks_deg[i];
		}
		
		angBtwPeakRad = new double[angles_between_peaks_deg.length];
		for (int i = 0; i < angles_between_peaks_deg.length; i++) {
			angBtwPeakRad[i] = (angles_between_peaks_deg[i]/360.0)*TwoPi;
		}

		nrPeaks = angles_between_peaks_deg.length;

        // initialize arrays
        sumON   = new float[nrPeaks+1];
        sumOFF  = new float[nrPeaks];
        nrON    = new int[nrPeaks+1];
        nrOFF   = new int[nrPeaks];

        // find the smallest angular resolution to carry out scores
        // for different rotations
        minAngResDeg = 360;
        for (int i = 0; i < angResDeg.length; i++){
            if (angResDeg[i]<minAngResDeg){
                minAngResDeg = angResDeg[i];
            }
        }
        minAngResRad = (minAngResDeg/360f)*TwoPi;

        // use the current minimal resolution to step when scoring rotations
        // bigger minimal angular scale => bigger steps

		nrRot = (int)Math.round( TwoPi / (minAngResRad/2) );    // number of rotations
		rotStepRad = minAngResRad/2;                            // step when rotating

		peaksRad = new double[nrRot][nrPeaks];
		// each configuration will be determined with
        // - rotations,
        // - peaks (and their angular width),
        // - angles between peaks
        // they will influence filter design
        // for this configuration
        // this is the core information of the config
        // (where the peaks are)
		
		for (int cnt_rots = 0; cnt_rots < nrRot; cnt_rots++) {
			
			double start_pos = cnt_rots*rotStepRad;

			for (int cnt_pks = 0; cnt_pks < nrPeaks; cnt_pks++) {
				
				if(cnt_pks==0)
					peaksRad[cnt_rots][cnt_pks] 	= start_pos;
				else
					peaksRad[cnt_rots][cnt_pks] 	=
                            peaksRad[cnt_rots][cnt_pks-1] +
                            angBtwPeakRad[cnt_pks-1];
				
			}
			
		}

	}

    public void print()
    {

        System.out.print("circular configuration (base rotation):  ");
        for (int j = 0; j < nrPeaks; j++){
            System.out.print(
                            IJ.d2s((peaksRad[0][j] / TwoPi) * 360, 2)+
                            " ( w "+
                            IJ.d2s((angResRad[j]/TwoPi)*360 ,2)+
                            " ) ; ring_range: "+
                            IJ.d2s(ringR[0],2)+" <-> "+
                            IJ.d2s(ringR[1],2)+" <-> inner ring "+innerRing
            );
        }

    }

    public float score(
            float[] val,
            float[] ang,
			float[] rad
    )
	{

        /*
        calculate score on this CircularConfiguration as
        highest score on all rotated versions of
        that configuration to be rotationally invariant
         */

        float 	score = Float.NEGATIVE_INFINITY;  // lowest possible score

        for (int r = 0; r < nrRot; r++) {

//            // set to zero
//            for (int k = 0; k < nrPeaks; k++) {
//
//                sumOFF[k] = 0;
//                nrOFF[k]  = 0;
//
//                sumON[k]  = 0;
//                nrON[k]   = 0;
//
//            }
//
//            //sumON[nrPeaks]= 0;
//            sumON[nrPeaks]= 0;
//            nrON[nrPeaks] = 0;

			float gmOFF = 0;
			int cntOFF = 0;
			float gmON  = 0;
			int cntON = 0;

            for (int i = 0; i < val.length; i++) {

				if (rad[i]<=ringR[1] && rad[i]>=ringR[0]) {

                    int regId = regionId(ang[i], r); // will give out index (+1, +nrPeaks), and (-1, -nrPeaks), and (nrPeaks+1)

					if(regId>0){
						gmON += Math.log(val[i]);
						cntON ++;
//						nrON[regId-1]++;
//                      sumON[regId-1]+= Math.log(val[i]);
					}
					else{
						gmOFF += Math.log(val[i]);
						cntOFF++;
//						nrOFF[-regId-1]++;
//						sumOFF[-regId-1]+=Math.log(val[i]);
					}
				}
//                else if (rad[i] <= innerRing) {
//					gmON += Math.log(val[i]);
//					cntON ++;
////                    nrON[nrPeaks]++;
////                    sumON[nrPeaks]+=Math.log(val[i]);
//                }

            }

//            float mulONAll = 1;
//            int nrONAll = 0;
//
//            float mulOFFAll = 1;
//            int nrOFFAll = 0;

//            for (int k = 0; k < nrPeaks; k++){
//                if (nrON[k]>0)  {mulONAll  *= Math.exp(sumON[k]/nrON[k]);  nrONAll++;}    //  sumONAll  += sumON[k]/nrON[k];
//                if (nrOFF[k]>0) {mulOFFAll *= Math.exp(sumOFF[k]/nrOFF[k]);  nrOFFAll++;} //  sumOFFAll += sumOFF[k]/nrOFF[k];
//            }

//            mulONAll  *= Math.exp(sumON[nrPeaks]/nrON[nrPeaks]); nrONAll++; // sumONAll += nrON[nrPeaks]/sumON[nrPeaks];

//            System.out.println("A"+nrPeaks + " : " + (sumON[nrPeaks]/nrON[nrPeaks]) +","+ nrON[nrPeaks] + " samples");

            //float sc =  ((sumONAll/nrONAll)-(sumOFFAll/nrOFFAll));
            //System.out.println("SCORE : "+sc+" (max "+score+"), (+)"+(sumONAll/nrONAll)+", (-)"+(sumOFFAll/nrOFFAll));

            double gmOn = Math.exp(gmON/cntON);
            double gmOff = Math.exp(gmOFF/cntOFF);
            //System.out.println("SCORE (geometric mean) : "+(gmOn-gmOff)+"), (+gm)"+gmOn+", (-gm)"+gmOff);

//            float sc = (float) Math.sqrt(gmOn/gmOff);
			float sc = (gmOff>0)? (float) Math.sqrt(gmOn/gmOff) : (float) Math.sqrt(gmOn);
            if(sc>score || r==0){
                score = sc;
            }

        }

        return score;

    }

	public float score(
			float[] val,
			float[] ang,
			float[] rad,
			int 	rotIdx
	)
	{

        // score for particular rotation
		float sc;
		int r = rotIdx;

		if (r >=0 && r < nrRot){

		//for (int r = 0; r < nrRot; r++) {

			// calculate scores for each rotation

//			// set to zero
//			for (int k = 0; k < nrPeaks; k++) {
//				sumOFF[k] = 0;
//				nrOFF[k]  = 0;
//				sumON[k]  = 0;
//				nrON[k]   = 0;
//			}
//			nrON[nrPeaks] = 0;
//			sumON[nrPeaks]= 0;

			float gmOFF = 0;
			int cntOFF = 0;
			float gmON  = 0;
			int cntON = 0;

			for (int i = 0; i < val.length; i++) {

				if (rad[i]<=ringR[1] && rad[i]>=ringR[0]) {

					int regId = regionId(ang[i], r); // will give out index (+1, +nrPeaks), and (-1, -nrPeaks), and (nrPeaks+1)

					if(regId>0){
						gmON += Math.log(val[i]);
						cntON ++;
//						nrON[regId-1]++;
//						sumON[regId-1]+=Math.log(val[i]);
					}
					else{
						gmOFF += Math.log(val[i]);
						cntOFF++;
//						nrOFF[-regId-1]++;
//						sumOFF[-regId-1]+=Math.log(val[i]);
					}
				}
//				else if (rad[i] <= innerRing) {
//					nrON[nrPeaks]++;
//					sumON[nrPeaks]+=Math.log(val[i]);
//				}
			}

//            float mulONAll = 1;
//			int nrONAll = 0;
//
//			float mulOFFAll = 1;
//            int nrOFFAll = 0;

			double gmOn = Math.exp(gmON/cntON);
			double gmOff = Math.exp(gmOFF/cntOFF);

//			for (int k = 0; k < nrPeaks; k++){
//				if (nrON[k]>0)  {mulONAll *= Math.exp(sumON[k]/nrON[k]);   nrONAll++;}  // sumONAll += sumON[k]/nrON[k];
//				if (nrOFF[k]>0) {mulOFFAll*= Math.exp(sumOFF[k]/nrOFF[k]); nrOFFAll++;} // sumOFFAll += sumOFF[k]/nrOFF[k];
//			}
//            mulONAll *= Math.exp(sumON[nrPeaks]/nrON[nrPeaks]); nrONAll++;   // sumONAll += nrON[nrPeaks]/sumON[nrPeaks];

//            System.out.println("A"+nrPeaks + " : " + (sumON[nrPeaks]/nrON[nrPeaks]) );

//            double gmOn = Math.pow(mulONAll, 1f/nrONAll);
//            double gmOff = Math.pow(mulOFFAll, 1f/nrOFFAll);

			//sc =  ((sumONAll/nrONAll)-(sumOFFAll/nrOFFAll));
            sc = (gmOff>0)? (float) Math.sqrt(gmOn/gmOff) : (float) Math.sqrt(gmOn);

		}
		else {
			System.out.println("illegal rotation index!");
			sc = Float.NaN;
		}

		return sc;

	}


//	public float[] scoreAllRot(
//            float[] val,
//            float[] ang,
//            float[] rad
//    )
//    {
//
//        float[] 	score = new float[nrRot];
//
//        for (int r = 0; r < nrRot; r++) {
//
//            // calculate scores for each rotation
//
//            // set to zero
//            for (int k = 0; k < nrPeaks; k++) {
//                sumOFF[k] = 0;
//                nrOFF[k]  = 0;
//                sumON[k]  = 0;
//                nrON[k]   = 0;
//            }
//            nrON[nrPeaks] = 0;
//            sumON[nrPeaks]= 0;
//
//            for (int i = 0; i < val.length; i++){
//
//                if (rad[i]<=ringR[1] && rad[i]>=ringR[0]){      // rad[] is radius actually, normalized 0-1.0, angle is calculated as (float)(atan2(r,c)+pi)
//
//                    int regId = regionId(ang[i], r); // will give out index (+1, +nrPeaks), and (-1, -nrPeaks)
//
//                    if(regId>0){
//                        nrON[regId-1]++;
//                        sumON[regId-1]+=Math.log(val[i]);
//                    }
//                    else{
//                        nrOFF[-regId-1]++;
//                        sumOFF[-regId-1]+=Math.log(val[i]);
//                    }
//                }
//                else if (rad[i] <= innerRing){
//                    nrON[nrPeaks]++;
//                    sumON[nrPeaks]+=Math.log(val[i]);
//                }
//
//            }
//
////            float sumONAll = 0;
//            float mulONAll = 1;
//            int nrONAll = 0;
//
////            float sumOFFAll = 0;
//            float mulOFFAll = 1;
//            int nrOFFAll = 0;
//
//            for (int k = 0; k < nrPeaks; k++){
//                if (nrON[k]>0)  {mulONAll  *= Math.exp(sumON[k]/nrON[k]);       nrONAll++;}  // sumONAll += sumON[k]/nrON[k];
//                if (nrOFF[k]>0) {mulOFFAll *= Math.exp(sumOFF[k]/nrOFF[k]);    nrOFFAll++;} // sumOFFAll += sumOFF[k]/nrOFF[k];
//            }
//
//            mulONAll *= Math.exp(nrON[nrPeaks]/sumON[nrPeaks]); nrONAll++;  // sumONAll += nrON[nrPeaks]/sumON[nrPeaks];
//
//            double gmOn = Math.pow(mulONAll, 1f/nrONAll);
//            double gmOff = Math.pow(mulOFFAll, 1f/nrOFFAll);
//
//            score[r] = (float) (gmOn-gmOff);//((sumONAll/nrONAll)-(sumOFFAll/nrOFFAll));
//
//        }
//
//        return score;
//
//    }


    private boolean isOn(
            float angle,
            int rotIdx
    )
	{

		for (int p = 0; p < nrPeaks; p++) {
			if(Math.abs(wrap_PI(angle-peaksRad[rotIdx][p])) <= angResRad[p]/2){
                // only peaks will change with rotations, angular resolutions stay the same
				return true;
			}
		}

		return false;

	}

    private int regionId(
            float angle,
            int rotIdx
    )
    {
        // give out positive index for region with ON
        // negative index for region with OFF
        // indexes range from 1 to nrPeaks and from -1 to -nrPeaks

        for (int p = 0; p < nrPeaks; p++) {
            if(Math.abs(wrap_PI(angle-peaksRad[rotIdx][p])) <= angResRad[p]/2){
                // only peaks will change with rotations, angular resolutions stay the same
                return p+1;
            }
        }

        // it is not ON, assign index now, check where it fits
        for (int p = 0; p < nrPeaks-1; p++) {
            // check if it falls in the (p, p+1)
            if(wrap_PI(angle-peaksRad[rotIdx][p]) > 0 && wrap_PI(angle-peaksRad[rotIdx][p+1])<0){
                // only peaks will change with rotations, angular resolutions stay the same
                return -(p+1);
            }
        }

        return -(nrPeaks);

    }
	
	public ImageProcessor plot(
            int N
    )
	{

        // N will represent the resolution, can be arbitrary value,
        // it is useful to set the real resolution of the profile (usually 2*radius+1)
        // where radius was used for profile extraction

        int centerX = (N-1)/2;
        int centerY = (N-1)/2;
		int Rout		= (N-1)/2;
        double R1      = ringR[0]*Rout;
        double R2      = ringR[1]*Rout;
        double Rinn     = innerRing*Rout;

	    ImageProcessor fp = new FloatProcessor(N, N);

            // fill the values
            for (int c = 0; c < fp.getWidth(); c++){
                for (int r = 0; r < fp.getHeight(); r++){

                    double ro = Math.pow((c-centerX), 2)+Math.pow((r-centerY), 2);

                    if (ro <= R2*R2 && ro >= R1*R1 ){

                        float theta = (float) (Math.atan2(r-centerY, c-centerX) + Math.PI);
                        int regId = regionId(theta, 0);

                        fp.setf(c, r, regId);

                    }
                    else if (ro<=Rinn*Rinn) {
                        fp.setf(c, r, nrPeaks+1);
                    }

                }
            }

		return fp;

	}

    public ImageStack plotAllRotations(
            int N
    )
    {

        // N will represent the resolution, can be arbitrary value,
        // it is useful to set the real resolution of the profile (usually 2*radius+1)
        // where radius was used for profile extraction,
        // to see how it really works

        ImageStack viz_is = new ImageStack(N, N);

        int centerX = (N-1)/2;
        int centerY = (N-1)/2;
        int Rout		= (N-1)/2;
        double R1      = ringR[0]*Rout;
        double R2      = ringR[1]*Rout;
        double Rinn     = innerRing*Rout;

        for (int rotIdx = 0; rotIdx < nrRot; rotIdx++) {

            ImageProcessor fp = new FloatProcessor(N, N);

            // fill the values
            for (int c = 0; c < fp.getWidth(); c++){
                for (int r = 0; r < fp.getHeight(); r++){

                    double ro = Math.pow((c-centerX), 2)+Math.pow((r-centerY), 2);

                    if (ro <= R2*R2 && ro >= R1*R1 ){

                        float theta = (float) (Math.atan2(r-centerY, c-centerX) + Math.PI);
                        int regId = regionId(theta, rotIdx);

                        fp.setf(c, r, regId);

                    }
                    else if(ro<=Rinn*Rinn){
                        fp.setf(c, r, nrPeaks+1);
                    }

                }
            }

            viz_is.addSlice("rot."+rotIdx+"/"+nrRot, fp);

        }

        return viz_is;

    }

	private static double wrap_PI(
            double in
    )
	{
		
		double out = in;
		
		while(out<=-Math.PI){
			out += 2*Math.PI;
		}
		while(out>Math.PI){
			out -= 2*Math.PI;
		}
		
		return out;
	}
	
}
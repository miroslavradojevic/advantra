package advantra.feature;

import ij.IJ;
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

    public static float TwoPi = (float) (2*Math.PI);

//	double[] 			filt;           // this filter will make a score
//                                        // it will have a different version with each rotation
//                                        // important to call initFilter() before using it
//                                        // before calling score() method
//                                        // to allocate enough values to cover the patch
	
	public CircularConfiguration(
            int[] angular_resolution_deg,
            int[] angles_between_peaks_deg,
			double[] radius,
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

//	public void initFilter(int length)
//	{
//		filt = new double[length];  // length will correspond to the number of locations in the patch
//	}

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

        double  sumPos, sumNeg;
        int     nrPos, nrNeg;

        for (int r = 0; r < nrRot; r++) {

            // calculate scores for each rotation

            sumPos  = sumNeg    = 0;
            nrPos   = nrNeg     = 0;

            for (int i = 0; i < val.length; i++){
				if (rad[i]<=ringR[1] && rad[i]>=ringR[0]){      // rad[] is radius actually, normalized 0-1.0
					if(isOn(ang[i], r)){                        // considering that the angle is calculated as (float)(atan2(r,c)+pi)
						nrPos++;
						sumPos+=val[i];
					}
					else{
						nrNeg++;
						sumNeg+=val[i];
					}
				}
                else if (rad[i] <= innerRing) {
                    nrPos++;
                    sumPos += val[i];
                }
            }

            float sc = (float) ((sumPos/nrPos)-(sumNeg/nrNeg));

            if(sc>score || r==0){
                score = sc;
            }

        }

        return score;

    }

    public float[] scoreAllRot(
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

        float[] 	score = new float[nrRot];
        double      sumPos, sumNeg;
        int         nrPos, nrNeg;

        for (int r = 0; r < nrRot; r++) {

            // calculate scores for each rotation

            sumPos  = sumNeg = 0;
            nrPos   = nrNeg  = 0;

            for (int i = 0; i < val.length; i++){

                if (rad[i]<=ringR[1] && rad[i]>=ringR[0]){      // rad[] is radius actually, normalized 0-1.0
                    if(isOn(ang[i], r)){                        // considering that the angle is calculated as (float)(atan2(r,c)+pi)
                        nrPos++;
                        sumPos+=val[i];
                    }
                    else{
                        nrNeg++;
                        sumNeg+=val[i];
                    }
                }
                else if (rad[i]<=innerRing){
                        nrPos++;
                        sumPos+=val[i];
                }


            }

//            System.out.println("sumPos "+sumPos+" sumNeg "+sumNeg+" nrPos "+nrPos+" nrNeg "+nrNeg);

            score[r] = (float) ((sumPos/nrPos)-(sumNeg-nrNeg));

        }

        return score;

    }


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

						if(isOn(theta, 0)){

							fp.setf(c, r, +1);

						}
						else {

							fp.setf(c, r, -1);

						}

                    }
                    else if (ro<=Rinn) {
                        fp.setf(c, r, +1);
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

                        if(isOn(theta, rotIdx)){

                            fp.setf(c, r, +1);

                        }
                        else {

                            fp.setf(c, r, -1);

                        }

                    }
                    else if(ro<=Rinn){
                        fp.setf(c, r, +1);
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
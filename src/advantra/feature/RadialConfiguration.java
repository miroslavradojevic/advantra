package advantra.feature;

import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/10/13
 * Time: 12:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class RadialConfiguration {

    public int          nrRings;

    public double[]     ringRes;      // the width of the ring

    public double[]     ringStart;    // radial distance at which the ring starts

	double[] 			filt;

    public RadialConfiguration(
            double[] ring_resol,
            double[] ring_btw_resol
    )
	{

		// filter default
		filt = new double[1];

        nrRings = ring_resol.length;

        ringRes = new double[nrRings];
        for (int i = 0; i < nrRings; i++){
            ringRes[i] = ring_resol[i];
        }

        ringStart = new double[nrRings];

        /*
        assign initial radial distances wrt. given resolution
        and predefined minimal radial distance between
         */

        ringStart[0] = 0;

        for (int i = 1; i < nrRings; i++){
            ringStart[i] = ringStart[i-1] + ringRes[i-1] + ring_btw_resol[i-1];
        }

    }

	public float score(
		float[] val,
		float[] rad
	)
	{

		int sumPos, sumNeg;
		sumPos = sumNeg = 0;

		for (int i = 0; i < val.length; i++){

			if(isOn(rad[i])){
				filt[i] = +1;
				sumPos++;
			}
			else{
				filt[i] = -1;
				sumNeg++;
			}

		}

		// normalize filt[]
		for (int i1 = 0; i1 < filt.length; i1++) {
			if(filt[i1]>0){
				filt[i1] = filt[i1] / sumPos;// * (sumNeg/(float)sumPos);
			}
			else{
				filt[i1] = filt[i1] / sumNeg;
			}
		}

//		System.out.println(""+val.length+" , "+filt.length);
//		System.out.println("(+): "+sumPos+","+(1f/sumPos)+" (-): "+sumNeg+", "+(1f/sumNeg));

		// score for this filt[]
		float sc = 0;
		for (int i = 0; i < filt.length; i++) {
			sc += filt[i]*val[i];
		}

		return sc;

	}

	public void initFilter(int length)
	{

		filt = new double[length];

	}

    public void print()
	{

         for (int i = 0; i < ringStart.length; i++){
             System.out.println("ring  :  "+ringStart[i]+" < "+ringRes[i]+" > ");
         }

    }

	private boolean isOn(float radius)
	{

		for (int i = 0; i < nrRings; i++){
			if(radius>=ringStart[i] && radius<=ringStart[i]+ringRes[i]){
				return true;
			}
		}

		return false;

	}

	public ImageStack plot()
	{
		int N = 101;
		ImageStack viz_is = new ImageStack(N, N);

		int centerX = (N-1)/2;
		int centerY = (N-1)/2;
		int R2		= (N-1)/2;

		ImageProcessor fp = new FloatProcessor(N, N);

		// fill the values
		for (int c = 0; c < fp.getWidth(); c++){
			for (int r = 0; r < fp.getHeight(); r++){

				float ro;
				ro = (float) (Math.sqrt(Math.pow((c-centerX), 2)+Math.pow((r-centerY), 2))/R2);

				if (ro<=1.0){

					fp.setf(c, r, -1);

					if(isOn(ro)){

						fp.setf(c, r, +1);

					}

				}

			}
		}

		viz_is.addSlice(fp);

		return viz_is;
	}

}

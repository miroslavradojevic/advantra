package advantra.filter;

import java.lang.IllegalArgumentException;

public class Lowpass {
	
	public static float[] createRow(float[] radius, float cutOff, int n){
		
		// allocates low pass butterworth filt.
		// cutOff       - cutoff frequency [0, 0.5]
		// n 			- order - defines the shape
		
		if(cutOff<0 || cutOff>0.5)
			throw new IllegalArgumentException("ERROR: cutoff frequency out of the [0,0.5] range");
		
		if(n<1)
			throw new IllegalArgumentException("ERROR: n must be an integer >= 1");
		
		float[] filt = new float[radius.length];
		
		//float radius;
		for (int i = 0; i < radius.length; i++) {

			filt[i] = (float) ( 1f / (1 + Math.pow(radius[i]/cutOff, 2*n)));
			
		}
		
		return filt;
	}

}
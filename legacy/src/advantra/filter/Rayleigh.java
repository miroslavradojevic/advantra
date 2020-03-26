package advantra.filter;

import flanagan.analysis.Stat;

public class Rayleigh {
	public static float mode(double[] data, int nbins){
		
		Stat st = new Stat(data);
		
		double mx 	= st.getMaximum_as_double();
		
		double[][] hist = Stat.histogramBins(data, mx/nbins, 0, mx);

		int idx = (new Stat(hist[1])).maximumIndex();
		
		float out = 0;
		out = (float) hist[0][idx];
		return out;
		
	}
}

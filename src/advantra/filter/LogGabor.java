package advantra.filter;

public class LogGabor {
	
	public static void createRow(float[] radius, float minWvl, float mult, float sigma, int scale, float[] output){
		
		// take the radius values from the predefined template
		
		for (int i = 0; i < radius.length; i++) {
			
			double f0 		= 1 / (minWvl * Math.pow(mult,scale));
			output[i]		= (float) Math.exp( -Math.pow(Math.log(radius[i]/f0), 2) / (2 * Math.pow(Math.log(sigma), 2)) );
			output[0]		= 0.f;
		
		}
	}
}
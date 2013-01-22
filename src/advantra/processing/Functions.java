package advantra.processing;

public class Functions {
	
	public static double gaussianFn(double input, double mean, double std){
		input = 1.0 * Math.exp((-1/(2*std*std)) * (input-mean) * (input-mean)); // (1/(Math.sqrt(2*Math.PI)*std))
		return input;
	}
	
	public static double nr_most_probable(double[] weights){
		// finds the diversity in weights
		double Ne = 0;
		
		for (int i = 0; i < weights.length; i++) {
			Ne += Math.pow(weights[i], 2);
		}
		
		return 1/Ne;
	}

	public static double nr_most_probable(double[] weights, int N){
		// finds the diversity in weights (first N weights)
		if(N>=weights.length){
			N = weights.length;
		}
		
		double Ne = 0;
		
		for (int i = 0; i < N; i++) {
			Ne += Math.pow(weights[i], 2);
		}
		
		return 1/Ne;
	}
	
}

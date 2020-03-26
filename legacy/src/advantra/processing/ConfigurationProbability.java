package advantra.processing;

public class ConfigurationProbability {
	double AngleLimit; // angle limit is set to 90 degrees
	int N; // discretize space
	double normalizationSum;
	
	public ConfigurationProbability(int N){
		// angle ranges from 0 to 90 degrees, discretized to N levels
		this.N = N;
		AngleLimit = Math.PI/2;
		this.normalizationSum = 0;
		for (int i = 1; i <= N; i++) {
			double a_i = (Math.PI/2)/(2*N) * (2*i-1);
			//System.out.println("a[i] :  "+a_i);
			normalizationSum += Math.cos(a_i);
		}
		//System.out.println("sum :  "+normalizationSum);
		//System.out.println("class created!");
	}
	
	private double discretizeAngle(double angle){
		
		double i = Math.floor(angle/(AngleLimit/N))+1;
		angle = (AngleLimit/(2*N)) * (2*i-1);
		angle = (angle > AngleLimit) ? AngleLimit : angle ;
		angle = (angle < 0) 		 ? 0 		  : angle ;
		return angle;
		
	}
	
	public double configurationPty(double angle){
		
		// discretize it...
		angle = discretizeAngle(angle);
		//return angle;
		return Math.cos(angle)/normalizationSum;
		
	}

}

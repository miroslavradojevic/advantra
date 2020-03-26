package advantra.feature;

public class PC_Params {

	public int nscale;
	public int norient;
	public int minWvl;
	public float mult;
	public float sigmaOnf;
	public float k;
	public float cutOff;
	public float g;
	public int noiseMethod;
	
	public PC_Params(){
		nscale = 4;
		norient = 12;
		minWvl = 3;
		mult = 2.1f;
		sigmaOnf = 0.55f;
		k = 10f;
		cutOff = 0.5f;
		g = 10;
		noiseMethod = -1;
	}

	public void printAll(){
		System.out.println("Parameters......");
		System.out.println("nscale="+nscale);
		System.out.println("norient="+norient);
		System.out.println("minimal wavelength="+minWvl);
		System.out.println("filter scaling(mult)="+mult);
		System.out.println("sigmaOnf="+sigmaOnf);
		System.out.println("k="+k);
		System.out.println("cutOff="+cutOff);
		System.out.println("g="+g);
		System.out.println("noiseMethod="+noiseMethod);
	}
}
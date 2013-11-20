package detection3d;

import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/13/13
 * Time: 11:40 AM
 */
public class Profiler3D extends Thread {

	private int begN, endN;

	public static double 	neuronDiam;
	public static double 	scale;
	public static int 		outerRange;
	public static int 		weightStdRatioToD;

	private static int 		resolDegTheta;
	private static int 		resolDegPhi;

	//public static ArrayList<ArrayList<double[]>> 	offsets;
	//public static ArrayList<ArrayList<Double>> 		weights;

	// profile resolution (4N+1) x (2N+1)  4N+1 cover 360 deg, 2N+1 cover 180 deg.
	public int 				N;
	public int[]			indexMap;
	public float[] 			theta; 			// azimuth - horizontal   -PI, +PI
	public float[] 			phi;    		// inclination - vertical -PI/2, +PI/2

	public Profiler3D(int n0, int n1){

		this.begN = n0; // split oriented filters
		this.endN = n1;

	}

	public static void loadTemplate(){

	}

	public static float[] sequence(float beginValue, float endValue, int Nelements, boolean lastIncluded)
	{

		float[] seqArray = new float[Nelements];

		if (lastIncluded) {
			float step = (endValue-beginValue)/(Nelements-1);
			for (int ii=0; ii<Nelements; ii++) seqArray[ii] = beginValue + ii * step;
		}
		else {
			float step = (endValue-beginValue)/Nelements;
			for (int ii=0; ii<Nelements; ii++) seqArray[ii] = beginValue + ii * step;
		}

		return seqArray;

	}

}

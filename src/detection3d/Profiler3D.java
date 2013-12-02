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

	public static Sphere3D 	sph3;
	public static int[][] 	listLocs3D;
	public static short[][]	prof3;

	public Profiler3D(int n0, int n1)
	{
		this.begN = n0; // split oriented filters
		this.endN = n1;
	}

	public static void loadTemplate(Sphere3D sph3_init, int[][] listLocs3D_init){

		sph3 = sph3_init; // just assign link, no allocation necessary since there will be one Sphere3D instance for all
		listLocs3D = listLocs3D_init;

		/*
			set prof3
		 */
		System.out.println("allocated for profiles: "+listLocs3D.length+" locations x "+sph3.getProfileLength()+" (profile length)");

		prof3 = new short[listLocs3D.length][sph3.getProfileLength()];

		System.out.println("input "+listLocs3D.length+" locations x "+listLocs3D[0].length+" coordinates");


	}

	public void run()
	{

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

package detection3d;

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/13/13
 * Time: 11:40 AM
 */

public class Profiler3DDemo implements PlugInFilter {

	private static float TwoPi = (float) (2*Math.PI);
	private static float OnePi = (float) (Math.PI);
	private static float HlfPi = (float) (Math.PI/2);

	public int setup(String s, ImagePlus imagePlus) {

		return 0;
	}

	public void run(ImageProcessor imageProcessor) {


	}

	/*
		terminal call
	 */

	public static void main(String[] args) {

		System.out.println("testing...");
		boolean lastIncluded = true;
		int N = 10;
		System.out.println(Arrays.toString(Profiler3D.sequence(-HlfPi, HlfPi, N, lastIncluded)));
		System.out.println(Arrays.toString(Profiler3D.sequence(0, TwoPi, N, !lastIncluded)));



	}
}

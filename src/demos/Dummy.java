package demos;

import aux.Stat;
import detection2d.FuzzyBranch2D;
import ij.plugin.PlugIn;

import java.util.Arrays;

/**
 * Created by miroslav on 3/25/14.
 * java -cp "$HOME/critpoint/*:$HOME/jarlib/*" demos.Dummy
 */
public class Dummy implements PlugIn {



	public static void main(String[] args) {

		demo();

	}

	public void run(String s) {

		demo();

	}

	public static void demo()
	{
		System.out.println("Dummy working...");
		float[] ar = new float[]{1, 2, 3, 4, 5};
		System.out.println(Arrays.toString(ar));
		Stat.normalize(ar);
		System.out.println(Arrays.toString(ar));

		FuzzyBranch2D test = new FuzzyBranch2D(			40,
													  // ncc
													  1, 0.2f,      // high
													  0.4f, 0.2f,   // low
													  // lhood
													  0.85f, 0.25f, // high
													  0.0f, 0.25f,  // low
													  .25f  // std output membership functions - defines separation
		);
		test.showFuzzification();
		test.showDefuzzification();


		float[] temp = new float[2];
		System.out.print("(1, 1) -> " + test.branchStrengthDefuzzified(1f, 1f));
		test.branchStrengthFuzzified(1f, 1f, temp);
		System.out.println(" -> " + Arrays.toString(temp));
		test.showAgg();

		System.out.print("(.9, .9) -> " + test.branchStrengthDefuzzified(.9f, .9f));
		test.branchStrengthFuzzified(.9f, .9f, temp);
		System.out.println(" -> " + Arrays.toString(temp));
		test.showAgg();

		test.showDefuzzificationSurface();
	}

}

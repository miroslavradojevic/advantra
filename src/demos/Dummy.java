package demos;

import detection2d.Fuzzy2D;
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

		Fuzzy2D test = new Fuzzy2D(			40,
													  // ncc
													  1, 0.2f,      // high
													  0.4f, 0.2f,   // low
													  // lhood
													  0.85f, 0.25f, // high
													  0.0f, 0.25f,  // low
													  .4f  // std output membership functions - defines separation
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

        System.out.print("endpoint on membership (.9, .9, .5, .3) -> " + test.endpointFuzzified(.9f, .9f, .5f, .3f));

	}

}

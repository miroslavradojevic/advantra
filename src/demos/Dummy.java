package demos;

import aux.Stat;

import java.util.Arrays;

/**
 * Created by miroslav on 3/25/14.
 * java -cp "$HOME/critpoint/*:$HOME/jarlib/*" demos.Dummy
 */
public class Dummy {
	public static void main(String[] args) {
		System.out.println("Dummy working...");
		float[] ar = new float[]{1, 2, 3, 4, 5};
		System.out.println(Arrays.toString(ar));
		Stat.normalize(ar);
		System.out.println(Arrays.toString(ar));
	}
}

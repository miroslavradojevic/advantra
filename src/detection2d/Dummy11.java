package detection2d;

import java.util.ArrayList;

/**
 * Created by miroslav on 25-6-14.
 */
public class Dummy11 {

	public static void main(String [] args) {
		System.out.println("testing..." + Float.NaN);

		ArrayList<CritpointRegion> regs = new ArrayList<CritpointRegion>(10);
		System.out.println(regs.size() + " elements");
		for (int i = 0; i < 10; i++)
			regs.add(new CritpointRegion(CritpointRegion.RegionType.END, 1f, 1f, 1f, 1f, new float[2][2], 1));

		System.out.println(regs.size() + " elements");
		for (int i = 0; i < regs.size(); i++) System.out.println(i + " : " + regs.get(i).type);
		System.out.println("---");


		regs.set(2, null);

		System.out.println(regs.size() + " elements");
		System.out.println(regs.size() + " elements");
		for (int i = 0; i < regs.size(); i++)
			if (regs.get(i)!=null)
				System.out.println(i + " : " + regs.get(i).type);
			else
				System.out.println("NULL");
		System.out.println("---");



	}

}

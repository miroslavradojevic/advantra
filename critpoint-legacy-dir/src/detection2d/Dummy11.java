package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 25-6-14.
 */
public class Dummy11 {

	public static void updateIt(FloatProcessor fp_to_update)
	{



		for (int i = 0; i < fp_to_update.getHeight(); i++) {
			fp_to_update.setf(i, i, (float)i);
		}


	}

	public static void main(String [] args) {

		ImagePlus ip = new ImagePlus();
		FloatProcessor fp = new FloatProcessor(new float[255][255]);
		ByteProcessor bp = new ByteProcessor(255, 255);

		updateIt(fp);

		float[] kernel = new float[9];
		Arrays.fill(kernel, 1/9f);


		ip.setProcessor(fp);

		System.out.println("old max: " + fp.getMax());
		IJ.saveAs(ip, "Tiff", "/home/miroslav/Desktop/orig.tif");

//		IJ.run(ip, "Smooth", "");
//		IJ.run(ip, "Smooth", "");
//		IJ.run(ip, "Smooth", "");
		fp.convolve(kernel, 3, 3);

		float[] pix = (float[]) fp.getPixels();

		float max_val = Float.NEGATIVE_INFINITY;
		for (int i = 0; i <pix.length; i++) {
			if (pix[i]>max_val)
				max_val = pix[i];
		}

		System.out.println("new max: " + max_val);


		for (int i = 0; i <pix.length; i++) {
			if (pix[i]>=0.5*max_val)
				bp.set(i, (int)255);
		}

		ip.setProcessor(bp);


		IJ.saveAs(ip, "Tiff", "/home/miroslav/Desktop/orig_th.tif");

		if (true) return;

		ArrayList[][] table = new ArrayList[10][10];

		table[0][0] = new ArrayList(); // add another ArrayList object to [0,0]
		table[0][0].add(new float[2]); // add object to that ArrayList


		for (int i = 0; i < table.length; i++) {
			for (int j = 0; j < table[0].length; j++) {
				//System.out.println(i+","+j+" -> "+table[i][j].size());
				table[i][j] = null;
			}
		}


		System.out.println("is it null? " + table[0][0]==null);
		table[1][2] = new ArrayList(1);
		table[1][2].add(new float[2]);
		table[3][4] = new ArrayList(3);
		table[3][4].add(new float[2]);
		table[1][2].add(new float[]{3,3});
		table[1][2].add(new float[]{3,4,4,56});

		for (int i = 0; i < table.length; i++) {
			for (int j = 0; j < table[0].length; j++) {
				if (table[i][j] !=null) System.out.println("at "+i+","+j+" -> "+table[i][j].size());
//				table[i][j] = null;
			}
		}

		System.out.println(Arrays.toString((float[])table[1][2].get(2)));

		if (true) return;

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

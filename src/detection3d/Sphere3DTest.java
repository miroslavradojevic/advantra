package detection3d;

import ij.ImagePlus;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/16/13
 * Time: 4:15 PM
 */
public class Sphere3DTest {

	public static void main(String[] args) {

		/*
			demo for the Sphere3D class
		 */

		float neuronDiameter = 6.0f;
		float scale = 1.5f;
		float nbhoodDeg = 20f;

		String dir = System.getProperty("user.home")+File.separator+"Sphere3D_Demo";
		File theDir = new File(dir);

		// if the directory does not exist, create it
		if (!theDir.exists()) {
			System.out.println("creating directory: " + dir);
			boolean result = theDir.mkdir();

			if(result) {
				System.out.println("dir created");
			}
		}

		// this is the demo on Sphere3D class
		Sphere3D s3 = new Sphere3D(neuronDiameter, scale);
		System.out.println("Sphere3D with neuron diameter "+neuronDiameter+" and scale "+scale+ " formed.");

		/*
			export sphere sampling
		 */
		String fileName = dir+File.separator+"S3D_direction_sampling.swc";
		s3.printFilterDirections(fileName);
		System.out.println(fileName);

		/*
			export sampling used to filtering the profiles
		 */
		System.out.println("\nexport filter offsets...\n");
		s3.printFilterOffsets(dir);

		/*
			export masks
		 */
		System.out.println("\nexport filter neighbourhood(s)...\n");
		s3.printMasks(dir);

		/*
			extract & show profile on sample image at some location
		 */


		/*
			detect peaks on selected profile
		 */


		/*
			extract profile element on sample image at some location (necessary for parallelisation)
		 */


		//new ImagePlus("test", s3.getProfile(new float[])).show();

	}

}

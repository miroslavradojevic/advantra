package detection3d;

import conn.Find_Connected_Regions;
import ij.ImagePlus;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/14/13
 * Time: 12:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class Trial {

	public static void main(String[] args){
		ImagePlus testSt = new ImagePlus("/home/miroslav/Desktop/testSt.tif");
		testSt.show();
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(testSt, true);  // true means save locations
		conn_reg.run("");
		System.out.println(" done.");
		ArrayList<ArrayList<int[]>> a =  conn_reg.getConnectedRegions3D_XYZ();
		System.out.println(a.size()+" regions");

		for (int i=0; i<a.size(); i++) {

			System.out.println("\n"+i + " : " + a.get(i).size() + " elements: ");
			for (int j=0; j<a.get(i).size(); j++) {
				System.out.print(Arrays.toString(a.get(i).get(j))+"\t");
			}

		}

		conn_reg.showLabels().show();

	}

	}

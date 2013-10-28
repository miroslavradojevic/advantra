package generate;

import aux.ReadSWC;
import aux.ReadSWC;
import ij.ImagePlus;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:46 PM
 */
public class Generator3D {
    // will be used to generate 3D image stacks

    private static int margin = 10;

    public static ImagePlus fromSWC(String swcPath, float k, float snr) {

        // read swc
        ReadSWC readerSWC = new ReadSWC(swcPath);   // read file into list

		//System.out.println(readerSWC.nodes.size()+" ELEMENTS");

		// determine boundaries
		int xMin, xMax, yMin, yMax, zMin, zMax, rMax;


		// allocate 8bit image

		// loop through the list to form cones,


		return new ImagePlus();

    }

	private void drawCone(float x, float y, float z, float r, float xPrev, float yPrev, float zPrev, float rPrev, byte[] image3d) {



	}

}

package generate;

import aux.ReadSWC;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import imagescience.random.PoissonGenerator;

import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:46 PM
 */
public class GeneratorSwc {
    // will be used to generate 2D/3D image stacks

    private static int bg = 20; // background level  (used when defining snr)
	private static int MARGIN_MIN = 50; // pixels
	private static float R_MIN = 1;

    /*
    	generate ImageStack from swc & export corresponding swc
     */
    public static ImageStack swc2image(String inSwc, boolean is2D, float k, float snr, String outRec, String outBif, String outEnd)
	{

		/*
		define the dimensions of the output stack based on readouts from inSwc, loop once to find boundaries
		 */
		ReadSWC readerSWC = new ReadSWC(inSwc); // read swc TODO: take this out, pass ReadSWC as argument

//		IJ.log("radius range: "+readerSWC.minR+" - "+readerSWC.maxR);

		float margin = MARGIN_MIN + ((readerSWC.maxR<1f)? 1f : readerSWC.maxR);
		// offsets when filling the values up
		float dX = -readerSWC.minX + margin;
		float dY = -readerSWC.minY + margin;
		float dZ = -readerSWC.minZ + margin;

//		IJ.log("margin: "+margin);
//		IJ.log("X : "+readerSWC.minX+" , "+readerSWC.maxX);
//		IJ.log("Y : "+readerSWC.minY+" , "+readerSWC.maxY);
//		IJ.log("Z : "+readerSWC.minZ+" , "+readerSWC.maxZ);

		// dimensions
		int W = (int)Math.ceil(readerSWC.maxX + dX + margin);
		int H = (int)Math.ceil(readerSWC.maxY + dY + margin);
		int L = (int)Math.ceil(readerSWC.maxZ + dZ + margin);
//		IJ.log("allocated image: "+W+" x "+H+" x "+L);

		// allocate 8 bit image
		byte[] outIm;
		if (is2D) outIm = new byte[W * H * 1]; // allocate 2d image
		else outIm = new byte[W * H * L]; // allocate 3d image

		float fg = foregroundLevel(bg, snr); // this is where snr calculations are embedded

        // initialize writers for reconstruction, bifurcations, and endpoints
        PrintWriter logRecWriter=null, logBifWriter=null, logEndWriter=null;
        // empty the files
        try {
            logRecWriter = new PrintWriter(outRec); logRecWriter.print(""); logRecWriter.close();
            logBifWriter = new PrintWriter(outBif); logBifWriter.print(""); logBifWriter.close();
            logEndWriter = new PrintWriter(outEnd); logEndWriter.print(""); logEndWriter.close();
        } catch (FileNotFoundException ex) {}

        // initialize output files
        try {
            logRecWriter = new PrintWriter(new BufferedWriter(new FileWriter(outRec, true))); //logRecWriter.println("# reconstruction\n# source swc file "+ inSwc);
            logBifWriter = new PrintWriter(new BufferedWriter(new FileWriter(outBif, true))); //logBifWriter.println("# bifurcations\n# source swc file "+ inSwc);
            logEndWriter = new PrintWriter(new BufferedWriter(new FileWriter(outEnd, true))); //logEndWriter.println("# endpoints\n# source swc file "+ inSwc);
        } catch (IOException e) {}


		// fill up the values and output swc files
		for (int ii=0; ii<readerSWC.nodes.size(); ii++) {  // loop through the list to draw cones and export bifs and ends

			int currId, currMotherId, laterMotherId, count;
			float currX, currY, currZ, currR;
			boolean isEnd, isBif;

			currId = Math.round(readerSWC.nodes.get(ii)[readerSWC.ID]);
			currMotherId = Math.round(readerSWC.nodes.get(ii)[readerSWC.MOTHER]);
			currX = readerSWC.nodes.get(ii)[readerSWC.XCOORD] + dX;
			currY = readerSWC.nodes.get(ii)[readerSWC.YCOORD] + dY;
			currZ = readerSWC.nodes.get(ii)[readerSWC.ZCOORD] + dZ;
			currR = readerSWC.nodes.get(ii)[readerSWC.RADIUS];

			//// critpoint check
			count = 0;
			isEnd = true;
			isBif = false;

			if (currMotherId!=-1) { // it is not root, assign it as endpoint by default if it is root, no need to loop
				for (int jj=ii+1; jj<readerSWC.nodes.size(); jj++) {

					laterMotherId = Math.round(readerSWC.nodes.get(jj)[readerSWC.MOTHER]);
					if (laterMotherId==currId) {
						isEnd = false;
						count++;
						if (count==2) {
							isBif = true;
							break;
						}
					}

				}
			}

			if (isEnd) {
				if(is2D) {
					assert logEndWriter != null;
					logEndWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", currId, 6, currX, currY, 0f, k*currR));
				}
				else {
					assert logEndWriter != null;
					logEndWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", currId, 6, currX, currY, currZ, k*currR));
				}
			}

			if (isBif) {
				if (is2D) {
					assert logBifWriter != null;
					logBifWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", currId, 5, currX, currY, 0f, k*currR));
				}
				else {
					assert logBifWriter != null;
					logBifWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", currId,	5, currX, currY, currZ, k*currR));
				}
			}

			if (is2D) {
				assert logRecWriter != null;
				logRecWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", currId, Math.round(readerSWC.nodes.get(ii)[readerSWC.TYPE]), currX, currY, 0f, k*currR, currMotherId));
			}
			else {
				assert logRecWriter != null;
				logRecWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", currId, Math.round(readerSWC.nodes.get(ii)[readerSWC.TYPE]), currX, currY, currZ, k*currR, currMotherId));
			}
			// fill the cones' values to array if it is not the root node
			if (currMotherId!=-1) {

				// loop backwards till the mother node is found
				int jj;
				for (jj=ii-1; jj>=0; jj--) {
					if (readerSWC.nodes.get(jj)[readerSWC.ID]==currMotherId) {
						break;
					}
				}

				float prevX = readerSWC.nodes.get(jj)[readerSWC.XCOORD] + dX;
				float prevY = readerSWC.nodes.get(jj)[readerSWC.YCOORD] + dY;
				float prevZ = readerSWC.nodes.get(jj)[readerSWC.ZCOORD] + dZ;
				float prevR = readerSWC.nodes.get(jj)[readerSWC.RADIUS];

				// draw
				currR = R_MIN + (currR - readerSWC.minR);
				prevR = R_MIN + (prevR - readerSWC.minR);

				if (is2D) {
					int countElements = drawCone(currX, currY, currR, prevX, prevY, prevR, outIm, W, k, fg);
//					System.out.println(" done(2d). "+countElements+" pixels added");
				}
				else {
					int countElements = drawCone(currX, currY, currZ, currR, prevX, prevY, prevZ, prevR, outIm, W, H, k, fg);
//					System.out.println(" done(3d). "+countElements+" pixels added");
				}

			}

		}

		/*
			add poisson noise
		*/

		System.out.print("adding poisson noise... ");

		PoissonGenerator poiss 	= new PoissonGenerator();

		int total_size = (is2D)? W*H : W*H*L;
		for (int j=0; j<total_size; j++) {

			// add background to be able to emulate poisson noise everywhere
			int currVal = (int)(outIm[j] & 0xff);
			currVal += bg;
			// poisson
			outIm[j] = (byte) ((int) Math.round(poiss.next(currVal)));

		}

		System.out.println("done");

		logRecWriter.close();  System.out.println(outRec+" exported.");
		logBifWriter.close();  System.out.println(outBif+" exported.");
		logEndWriter.close();  System.out.println(outEnd+" exported.");

		ImageStack  isOut       = toImageStack(outIm, W, H);
		return isOut;

    }

	private static int drawCone(float x, float y, float r, float xPrev, float yPrev, float rPrev, byte[] image2d, int W, float k, float fg)
	{

		// x,y are expected to be valid 2d coordinates of the image stack layer
		// cone defined with (x1,y1,r1),(x2,y2,r2)

		int H = image2d.length / W;

		// count how many values were added
		int count = 0;

		// range in x
		int xMin = (int)Math.floor(Math.min(x-r, xPrev-rPrev));
		int xMax = (int)Math.ceil(Math.max(x+r, xPrev+rPrev));

		// range in y
		int yMin = (int)Math.floor(Math.min(y-r, yPrev-rPrev));
		int yMax = (int)Math.ceil(Math.max(y+r, yPrev+rPrev));

//		// range in z
//		int zMin = (int)Math.floor(Math.min(z-r, zPrev-rPrev));
//		int zMax = (int)Math.ceil(Math.max(z+r, zPrev+rPrev));

		// loop
		for (int xLoop=xMin; xLoop<=xMax; xLoop++) {
			for (int yLoop=yMin; yLoop<=yMax;yLoop++) {
//				for (int zLoop=zMin; zLoop<=zMax; zLoop++) {
					if ((xLoop>=0 && xLoop<W) && (yLoop>=0 && yLoop<H)) { // is inside the image && (zLoop>=0 && zLoop<L)

						int currentIndex        = xy2ind(xLoop, yLoop, W);
						int currentValue        = (image2d[currentIndex] & 0xff);           // byte to int
						int calculatedValue     = //255;
								coneIntensity(
													 xLoop, yLoop,// zLoop,
													 x, y, r, //z,
													 xPrev, yPrev, rPrev, //, zPrev
													 k, fg
								);

						if (calculatedValue>currentValue) {
							image2d[currentIndex] = (byte)calculatedValue;
							count++;
						}

					}
//				}
			}
		}

		return count;

	}


	private static int drawCone(float x, float y, float z, float r, float xPrev, float yPrev, float zPrev, float rPrev, byte[] image3d, int W, int H, float k, float fg)
    {

        // x,y,z are expected to be valid 3d coordinates of the image stack
        // cone defined with (x1,y1,z1,r1),(x2,y2,z2,r2)

		int L = image3d.length / (W*H);

        // count how many values were added
        int count = 0;

        // range in x
        int xMin = (int)Math.floor(Math.min(x-r, xPrev-rPrev));
        int xMax = (int)Math.ceil(Math.max(x+r, xPrev+rPrev));

        // range in y
        int yMin = (int)Math.floor(Math.min(y-r, yPrev-rPrev));
        int yMax = (int)Math.ceil(Math.max(y+r, yPrev+rPrev));

        // range in z
        int zMin = (int)Math.floor(Math.min(z-r, zPrev-rPrev));
        int zMax = (int)Math.ceil(Math.max(z+r, zPrev+rPrev));

        // loop
        for (int xLoop=xMin; xLoop<=xMax; xLoop++) {
            for (int yLoop=yMin; yLoop<=yMax;yLoop++) {
                for (int zLoop=zMin; zLoop<=zMax; zLoop++) {
					if ((xLoop>=0 && xLoop<W) && (yLoop>=0 && yLoop<H) && (zLoop>=0 && zLoop<L)) { // is inside the image

                        int currentIndex        = xyz2ind(xLoop, yLoop, zLoop, W, H);
                        int currentValue        = (image3d[currentIndex] & 0xff);           // byte to int
                        int calculatedValue     = //255;
                                coneIntensity(
                                xLoop, yLoop, zLoop,
                                x, y, z, r,
                                xPrev, yPrev, zPrev, rPrev,
                                k, fg
                                );

                        if (calculatedValue>currentValue) {
                            image3d[currentIndex] = (byte)calculatedValue;
                            count++;
                        }

					}
                }
            }
        }

        return count;

	}

	private static int coneIntensity( // will output value at xyLoc based on cone geometry and gaussian cross-profile outputs -1 in case it is out of the cone
									  int xLoc, int yLoc, //int zLoc,
									  float x1, float y1, float r1, // float z1,
									  float x2, float y2, float r2, //  float z2,
									  float k, // scales the gaussian sigma
									  float fg // foreground scaled wrt snr and chosen background
	)
	{
		// returns 0-255 value for gaussian cross-profile
		// of the cone based on the normal distance from the axis
		// if it's out of the cone or the cone is really short, returns -1
		float p21_x = x2 - x1;
		float p21_y = y2 - y1;
//		float p21_z = z2 - z1;

		double P12 = Math.sqrt(p21_x*p21_x + p21_y*p21_y); //  + p21_z*p21_z

		if (P12>0.0001) {// don't check below some small distance
			p21_x /= P12;   // normalize
			p21_y /= P12;
//			p21_z /= P12;
		}
		else
			return 0;

		float p12_x = -p21_x;
		float p12_y = -p21_y;
//		float p12_z = -p21_z;

		float c1_x = xLoc - x1;
		float c1_y = yLoc - y1;
//		float c1_z = zLoc - z1;

		float c2_x = xLoc - x2;
		float c2_y = yLoc - y2;
//		float c2_z = zLoc - z2;

		if ( p21_x * c1_x + p21_y * c1_y < 0 ) { // + p21_z * c1_z
			// take the spherical distance from (x1,y1,z1)

			double d1 		= Math.sqrt( Math.pow(c1_x, 2) + Math.pow(c1_y, 2) ); //+ Math.pow(c1_z, 2)
			double sigma1 	= r1 * k;

			int val = 0;
			if (d1<=3*sigma1)
				val =  (int) (fg * Math.exp(-(d1 * d1) / (2 * sigma1 * sigma1)));

			return val;
		}
		if ( p12_x * c2_x + p12_y * c2_y < 0 ) { //  + p12_z * c2_z
			// take the spherical distance from (x2,y2,z2)

			double d2 		= Math.sqrt( Math.pow(c2_x, 2) + Math.pow(c2_y, 2) ); // + Math.pow(c2_z, 2)
			double sigma2 	= r2 * k;

			int val = 0;
			if (d2<=3*sigma2)
				val =  (int) (fg * Math.exp(-(d2 * d2) / (2 * sigma2 * sigma2)));

			return val;
		}


		float d2 = (c2_x * p12_x) + (c2_y * p12_y) ; // projection c2 on p12      + (c2_z * p12_z)

		double d = Math.sqrt( Math.pow(-p12_x * d2 + c2_x, 2) + Math.pow(-p12_y * d2 + c2_y, 2) ); // + Math.pow(-p12_z * d2 + c2_z, 2)

		// interpolate radius based on r1 and r2
		double r = (d2*r1 + (P12-d2)*r2)/P12;
		double sigma =  r * k;
		// todo: add poisson here too, so that r or sigma change

		// use d, sigma to calculate the value and scale it with fg
		int val = 0;
		if (d<=3*sigma)
			//val = (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
			//val =  (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
			val =  (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
		return val;

	}

    private static int coneIntensity( // will output value at xyzLoc based on cone geometry and gaussian cross-profile outputs -1 in case it is out of the cone
        int xLoc, int yLoc, int zLoc,
        float x1, float y1, float z1, float r1,
        float x2, float y2, float z2, float r2,
        float k, // scales the gaussian sigma
        float fg // foreground scaled wrt snr and chosen background
    )
    {
        // returns 0-255 value for gaussian cross-profile
        // of the cone based on the normal distance from the axis
        // if it's out of the cone or the cone is really short, returns -1
        float p21_x = x2 - x1;
        float p21_y = y2 - y1;
        float p21_z = z2 - z1;

        double P12 = Math.sqrt(p21_x*p21_x + p21_y*p21_y + p21_z*p21_z);

        if (P12>0.0001) {// don't check below some small distance
            p21_x /= P12;   // normalize
            p21_y /= P12;
            p21_z /= P12;
        }
        else
            return 0;

        float p12_x = -p21_x;
        float p12_y = -p21_y;
        float p12_z = -p21_z;

        float c1_x = xLoc - x1;
        float c1_y = yLoc - y1;
        float c1_z = zLoc - z1;

        float c2_x = xLoc - x2;
        float c2_y = yLoc - y2;
        float c2_z = zLoc - z2;

        if ( p21_x * c1_x + p21_y * c1_y + p21_z * c1_z < 0 ) {
			// take the spherical distance from (x1,y1,z1)

			double d1 		= Math.sqrt( Math.pow(c1_x, 2) + Math.pow(c1_y, 2) + Math.pow(c1_z, 2) );
			double sigma1 	= r1 * k;

			int val = 0;
			if (d1<=3*sigma1)
				val =  (int) (fg * Math.exp(-(d1 * d1) / (2 * sigma1 * sigma1)));

			return val;
		}
        if ( p12_x * c2_x + p12_y * c2_y + p12_z * c2_z < 0 ) {
			// take the spherical distance from (x2,y2,z2)

			double d2 		= Math.sqrt( Math.pow(c2_x, 2) + Math.pow(c2_y, 2) + Math.pow(c2_z, 2) );
			double sigma2 	= r2 * k;

			int val = 0;
			if (d2<=3*sigma2)
				val =  (int) (fg * Math.exp(-(d2 * d2) / (2 * sigma2 * sigma2)));

			return val;
		}


        float d2 = (c2_x * p12_x) + (c2_y * p12_y) + (c2_z * p12_z); // projection c2 on p12

        double d = Math.sqrt( Math.pow(-p12_x * d2 + c2_x, 2) + Math.pow(-p12_y * d2 + c2_y, 2) + Math.pow(-p12_z * d2 + c2_z, 2) );

        // interpolate radius based on r1 and r2
        double r = (d2*r1 + (P12-d2)*r2)/P12;
        double sigma =  r * k;
		// todo: add poisson here too, so that r or sigma change

        // use d, sigma to calculate the value and scale it with fg
        int val = 0;
        if (d<=3*sigma)
            //val = (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
            //val =  (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
            val =  (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
        return val;

    }

    private static float foregroundLevel(float bg, float snr)
    {
        // explanation b=background, a=added on top of b
        // a        signal
        // a+b      background
        // both poisson distribution (mean = lambda, variance = lambda, sigma = sqrt(lambda))
        // snr = mean(signal)/sigma(signal) = a / sqrt(a+b)
        // solve quadratic equation for a
        // quadratic equation is:
        // a^2  -snr^2  *a  -snr^2*b   = 0
        // a  = (1/2) (snr^2 + snr * sqrt(snr^2+4b))
        return (float) (0.5 * (snr * snr + snr * Math.sqrt(4 * bg + snr * snr)));

    }

    private static boolean isCone(
            int xLoc,
            int yLoc,
            int zLoc,

            float x1,
            float y1,
            float z1,
            float r1,

            float x2,
            float y2,
            float z2,
            float r2
    )
	{

        float p21_x = x2 - x1;
        float p21_y = y2 - y1;
        float p21_z = z2 - z1;

        double P12 = Math.sqrt(p21_x*p21_x + p21_y*p21_y + p21_z*p21_z);

        if (P12>0.1) {// don't check below some small distance
            p21_x /= P12;   // normalize
            p21_y /= P12;
            p21_z /= P12;
        }
        else {
            return false;
        }

        float p12_x = -p21_x; //x1 - x2; // -p1_x
        float p12_y = -p21_y; //y1 - y2; // -p2_y
        float p12_z = -p21_z; //z1 - z2; // -p2_z

        float c1_x = xLoc - x1;
        float c1_y = yLoc - y1;
        float c1_z = zLoc - z1;

        float c2_x = xLoc - x2;
        float c2_y = yLoc - y2;
        float c2_z = zLoc - z2;

        // dot prod p21*c1
        if ( p21_x * c1_x + p21_y * c1_y + p21_z * c1_z >= 0 ) {
            // dot prod p12*c2
            if ( p12_x * c2_x + p12_y * c2_y + p12_z * c2_z >= 0 ) {

                float d2 = (c2_x * p12_x) + (c2_y * p12_y) + (c2_z * p12_z); // projection c2 on p12

                double dist = Math.sqrt( Math.pow(-p12_x * d2 + c2_x, 2) + Math.pow(-p12_y * d2 + c2_y, 2) + Math.pow(-p12_z * d2 + c2_z, 2) );

                // use d2 to check dist threshold
                double distLimit = (d2*r1 + (P12-d2)*r2)/P12;

                if (dist<=distLimit) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }

    }

    private static int[] idx2xyz(int idx, int W, int H){  // will contain the allocation inside
        int[] xyz = new int[3];

        // z = layer
        xyz[2] = idx / (W*H);
        // x ~ width
        xyz[0] = (idx - xyz[2] * H * W) % W;
        // y ~ height
        xyz[1] = (idx - xyz[2] * H * W) / W;

        return xyz;
    }

    private static int xyz2ind(int x, int y, int z, int W, int H){
        return H*W*z + y*W + x;
    }

	private static int xy2ind(int x, int y, int W){
		return y*W + x;
	}

    private static ImageStack toImageStack(byte[] image3d, int W, int H) {

        int L = image3d.length / (W*H);

        ImageStack isOut = new ImageStack(W, H);

        for (int ii=0; ii<L; ii++) {

            ByteProcessor ipOut = new ByteProcessor(W, H);

            for (int jj=0; jj<H*W; jj++) {
                ipOut.set(jj, (int)image3d[ii*H*W+jj]&0xff);
            }

            isOut.addSlice(ipOut);

        }
        return isOut;

    }

}

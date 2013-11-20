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
public class Generator3D {
    // will be used to generate 3D image stacks

    private static int bg = 20; // background level  (used when defining snr)

    /*
    	generate ImageStack from swc & export corresponding swc
     */
    public static ImageStack swc2Stack(String inSwc, float k, float snr, String outSwc, String outBif, String outEnd, String outTif)
	{

		// TODO fix cases when there are several root locations (with -1) - it considers only one now

        /*
            initialize writer
         */
        PrintWriter logWriter = null;

        /*
            empty the file
         */
        try {
            logWriter = new PrintWriter(outSwc);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        /*
            initialize detection log file
         */
        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outSwc, true)));
            logWriter.println("# source "+ inSwc);
        } catch (IOException e) {}

        // read swc
        ReadSWC readerSWC = new ReadSWC(inSwc);

		System.out.println(readerSWC.nodes.size()+" nodes found in " + inSwc);

		// determine boundaries based on the swc node coordinates
		float
				xMin=Float.MAX_VALUE, xMax=Float.MIN_VALUE,
                yMin=Float.MAX_VALUE, yMax=Float.MIN_VALUE,
                zMin=Float.MAX_VALUE, zMax=Float.MIN_VALUE,
                rMax=Float.MIN_VALUE;

        for (int ii=0; ii<readerSWC.nodes.size(); ii++) {

            double readVal; //, readY, readZ, readR;

            // 2 -> x
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[2]));
            if (readVal>xMax) xMax = (int)readVal;
			readVal = Math.round(Math.floor(readerSWC.nodes.get(ii)[2]));
            if (readVal<xMin) xMin = (int)readVal;
            // 3 -> y
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[3]));
            if (readVal>yMax) yMax = (int)readVal;
			readVal = Math.round(Math.floor(readerSWC.nodes.get(ii)[3]));
            if (readVal<yMin) yMin = (int)readVal;
            // 4 -> z
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[4]));
            if (readVal>zMax) zMax = (int)readVal;
			readVal = Math.round(Math.floor(readerSWC.nodes.get(ii)[4]));
            if (readVal<zMin) zMin = (int)readVal;
            // 5 -> r
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[5]));
            if (readVal>rMax) rMax = (int)readVal;

        }

		// xMin + margin will correspond to start (index 0)
		float margin = rMax + 5;

        float dX = -xMin + margin;
        float dY = -yMin + margin;
        float dZ = -zMin + margin;

        int W = (int)Math.ceil(xMax + dX + margin);
        int H = (int)Math.ceil(yMax + dY + margin);
        int L = (int)Math.ceil(zMax + dZ + margin);

        byte[] outIm = new byte[W*H*L]; // allocate 8 bit image

        System.out.println("allocated image: "+W+" x "+H+" x "+L);

        float fg = foregroundLevel(bg, snr);

        // first line
		logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1",
											   Math.round(readerSWC.nodes.get(0)[0]+1),
											   Math.round(readerSWC.nodes.get(0)[1]),
											   readerSWC.nodes.get(0)[2]+dX,
											   readerSWC.nodes.get(0)[3]+dY,
											   readerSWC.nodes.get(0)[4]+dZ,
											   readerSWC.nodes.get(0)[5],
											   -1)
		);

		// loop through the list to form cones
        for (int coneId=1; coneId<readerSWC.nodes.size(); coneId++) {

            float x = readerSWC.nodes.get(coneId)[2] + dX; //- rMax + shiftX;
            float y = readerSWC.nodes.get(coneId)[3] + dY; //- rMax + shiftY;
            float z = readerSWC.nodes.get(coneId)[4] + dZ; //- rMax + shiftZ;
            float r = readerSWC.nodes.get(coneId)[5]; 	   //- rMax + shiftZ;

            int indexPrev = Math.round(readerSWC.nodes.get(coneId)[6]);

            float xPrev = readerSWC.nodes.get(indexPrev)[2] + dX;   //- rMax + shiftX;
            float yPrev = readerSWC.nodes.get(indexPrev)[3] + dY;   //- rMax + shiftY;
            float zPrev = readerSWC.nodes.get(indexPrev)[4] + dZ;   //- rMax + shiftZ;
            float rPrev = readerSWC.nodes.get(indexPrev)[5];        //- rMax + shiftZ;

			// clamp th radius values to be minimum 1 when generating and storing the output
			r = 1 + (r-readerSWC.minR);
			rPrev = 1 + (rPrev-readerSWC.minR);

            int countElements = drawCone(x,y,z,r,   xPrev,yPrev,zPrev,rPrev, outIm, W, H, k, fg);

            //System.out.println(" done. "+countElements+" pixels added");

			/*
            	add corresponding swc
			 */
            logWriter.println( // TODO String.format
                        IJ.d2s(readerSWC.nodes.get(coneId)[0]+1, 0) + " " +
                        IJ.d2s(readerSWC.nodes.get(coneId)[1], 0) + " " +
                        IJ.d2s(x, 3)                              + " " +
                        IJ.d2s(y, 3)                              + " " +
                        IJ.d2s(z, 3)                              + " " +
                        IJ.d2s(r, 3)                              + " " +
                        IJ.d2s(indexPrev+1, 0)

            );



        }

		/*
			add poisson noise
		*/

		System.out.print("adding poisson noise... ");

		PoissonGenerator poiss 	= new PoissonGenerator();

		for (int j=0; j<W*H*L; j++) {

			// add background to be able to emulate poisson noise everywhere
			int currVal = (int)(outIm[j] & 0xff);
			currVal += bg;
			// poisson
			outIm[j] = (byte) ((int) Math.round(poiss.next(currVal)));

		}

		System.out.println("done");

        logWriter.close();

        System.out.println(outSwc+" exported");

		/*
			critical points
		 */
		ReadSWC readNewSwc = new ReadSWC(outSwc);
		readNewSwc.exportBifurcations(outBif);   	// extract bifurcations
		readNewSwc.exportEndpoints(outEnd);         // extract endpoints

        ImageStack  isOut       = toImageStack(outIm, W, H);

        FileSaver fs = new FileSaver(new ImagePlus("", isOut));
        fs.saveAsTiffStack(outTif);

        System.out.println(outTif+" exported");

		return isOut;

    }

	public static ImageStack yjunctions2Stack(float snr, float diam1, float diam2, float diam3, String outSwc, String outBif, String outEnd, String outTif)
	{

		// generate y juctions with branches of different scales in 3d
		return new ImageStack();

	}




	private static int drawCone(float x, float y, float z, float r, float xPrev, float yPrev, float zPrev, float rPrev, byte[] image3d, int W, int H, float k, float fg)
    {

        // define the range to loop for plotting
        // x,y,z are expected to be valid 3d coordinates of the image stack
        // cone defined with (x1,y1,z1,r1),(x2,y2,z2,r2) should not be out of image (but that's taken care of when defining the size)

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

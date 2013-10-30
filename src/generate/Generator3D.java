package generate;

import aux.ReadSWC;
import aux.ReadSWC;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;

import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:46 PM
 */
public class Generator3D {
    // will be used to generate 3D image stacks

    private static int bg = 20; // background level

    public static double Swc2MinRadius(
            String swcPath
    )
    {
        // read swc
        ReadSWC readerSWC = new ReadSWC(swcPath);   // read file into list

        double rMin = Double.MAX_VALUE;

        for (int ii=0; ii<readerSWC.nodes.size(); ii++) {
            double readR = Math.round(Math.ceil(readerSWC.nodes.get(ii)[3])); // 3rd is the radius
            if (readR<rMin)
                rMin = readR;
        }

        return rMin;

    }

    /*
    generate ImageStack from swc & export corresponding swc
     */
    public static ImageStack Swc2Stack(String inSwc, float k, float snr, String outSwc) {

        // outSwc will be recorded

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
            logWriter.println("# "+ inSwc + " reconstruction");
        } catch (IOException e) {}


        // read swc
        ReadSWC readerSWC = new ReadSWC(inSwc);

		System.out.println(readerSWC.nodes.size()+" nodes in SWC");

		// determine boundaries based on the swc node coordinates
		float
				xMin=Float.MAX_VALUE, xMax=Float.MIN_VALUE,
                yMin=Float.MAX_VALUE, yMax=Float.MIN_VALUE,
                zMin=Float.MAX_VALUE, zMax=Float.MIN_VALUE,
                rMax=Float.MIN_VALUE;


        for (int ii=0; ii<readerSWC.nodes.size(); ii++) {

            double readVal;//, readY, readZ, readR;

            // 2 -> x
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[0]));
            if (readVal>xMax) xMax = (int)readVal;
			readVal = Math.round(Math.floor(readerSWC.nodes.get(ii)[0]));
            if (readVal<xMin) xMin = (int)readVal;
            // 3 -> y
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[1]));
            if (readVal>yMax) yMax = (int)readVal;
			readVal = Math.round(Math.floor(readerSWC.nodes.get(ii)[1]));
            if (readVal<yMin) yMin = (int)readVal;
            // 4 -> z
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[2]));
            if (readVal>zMax) zMax = (int)readVal;
			readVal = Math.round(Math.floor(readerSWC.nodes.get(ii)[2]));
            if (readVal<zMin) zMin = (int)readVal;
            // 5 -> r
			readVal = Math.round(Math.ceil(readerSWC.nodes.get(ii)[3]));
            if (readVal>rMax) rMax = (int)readVal;

        }

		// xMin + margin will correspond to start (index 0)
		float margin = rMax + 10;

        System.out.println("margin: "+margin);

        float dX = -xMin + margin;
        float dY = -yMin + margin;
        float dZ = -zMin + margin;

        int W = (int)Math.ceil(xMax + dX + margin);
        int H = (int)Math.ceil(yMax + dY + margin);
        int L = (int)Math.ceil(zMax + dZ + margin);

        byte[] outIm = new byte[W*H*L]; // allocate 8 bit image

        System.out.println("allocated image: "+W+" x "+H+" x "+L);

        float fg = foregroundLevel(bg, snr);

        System.out.println("foreground set to: "+fg);

        logWriter.println("first line");

		// loop through the list to form cones
        for (int coneId=1; coneId<readerSWC.nodes.size(); coneId++) {

            System.out.print("drawing cone "+coneId+"/"+(readerSWC.nodes.size()-1));

            // set value at the node center location to 255
//            int locX = (int)Math.round(readerSWC.nodes.get(coneId)[0] - rMax + shiftX);
//            int locY = (int)Math.round(readerSWC.nodes.get(coneId)[1] - rMax + shiftY);
//            int locZ = (int)Math.round(readerSWC.nodes.get(coneId)[2] - rMax + shiftZ);
//            outIm[xyz2ind(locX, locY, locZ, W, H)] = (byte)255;

            float x = readerSWC.nodes.get(coneId)[0] + dX; //- rMax + shiftX;
            float y = readerSWC.nodes.get(coneId)[1] + dY; //- rMax + shiftY;
            float z = readerSWC.nodes.get(coneId)[2] + dZ; //- rMax + shiftZ;
            float r = readerSWC.nodes.get(coneId)[3]; 	   //- rMax + shiftZ;

            int indexPrev = Math.round(readerSWC.nodes.get(coneId)[4]);

            float xPrev = readerSWC.nodes.get(indexPrev)[0] + dX; //- rMax + shiftX;
            float yPrev = readerSWC.nodes.get(indexPrev)[1] + dY; //- rMax + shiftY;
            float zPrev = readerSWC.nodes.get(indexPrev)[2] + dZ; //- rMax + shiftZ;
            float rPrev = readerSWC.nodes.get(indexPrev)[3];// - rMax + shiftZ;

            int countElements = drawCone(x,y,z,r,   xPrev,yPrev,zPrev,rPrev, outIm, W, H, k, fg);

            System.out.println(" done. "+countElements+" pixels added");

            logWriter.println("modified line...");

        }

        logWriter.close();

        System.out.println(outSwc+" exported");


        // draw the cone
		return toImageStack(outIm, W, H);//new ImagePlus("out", );

    }

//    private static int drawCone(float x, float y, float z, float r, float xPrev, float yPrev, float zPrev, float rPrev, byte[] image3d, int W, int H, float k, float fg) {
//        image3d[10] = (byte) 18;
//        return 1;
//    }

	private static int drawCone(float x, float y, float z, float r, float xPrev, float yPrev, float zPrev, float rPrev, byte[] image3d, int W, int H, float k, float fg) {


        // aux
        int currentIndex        = xyz2ind((int)Math.round(x), (int)Math.round(y), (int)Math.round(z), W, H);
        int currentValue        = (image3d[currentIndex] & 0xff);           // byte to int
        int calculatedValue     = 255;
        if (calculatedValue>currentValue)
            image3d[currentIndex] = (byte)calculatedValue;

        System.out.println("at "+x+" , "+y+" , "+z);
        System.out.println("at "+xPrev+" , "+yPrev+" , "+zPrev);


        currentIndex        = xyz2ind((int)Math.round(xPrev), (int)Math.round(yPrev), (int)Math.round(zPrev), W, H);
        currentValue        = (image3d[currentIndex] & 0xff);           // byte to int
        calculatedValue     = 255;
        if (calculatedValue>currentValue)
            image3d[currentIndex] = (byte)calculatedValue;

        return 1;


        }
//        {
//
//
//        // define the range to loop for plotting
//        // x,y,z are expected to be valid 3d coordinates of the image stack
//        // cone defined with (x1,y1,z1,r1),(x2,y2,z2,r2) should not be out of image (but that's taken care of when defining the size)
//
//		int L = image3d.length / (W*H);
//
//        // count how many values were added
//        int count = 0;
//
//        // range in x
//        int xMin = (int)Math.floor(Math.min(x-r, xPrev-rPrev));
//        int xMax = (int)Math.ceil(Math.max(x+r, xPrev+rPrev));
//
//        // range in y
//        int yMin = (int)Math.floor(Math.min(y-r, yPrev-rPrev));
//        int yMax = (int)Math.ceil(Math.max(y+r, yPrev+rPrev));
//
//        // range in z
//        int zMin = (int)Math.floor(Math.min(z-r, zPrev-rPrev));
//        int zMax = (int)Math.ceil(Math.max(z+r, zPrev+rPrev));
//
//        System.out.println(xPrev+" , "+x+ "    x range: "+xMin+" -- "+xMax);
//
//        // loop
//        for (int xLoop=xMin; xLoop<=xMax; xLoop++) {
//            for (int yLoop=yMin; yLoop<=yMax;yLoop++) {
//                for (int zLoop=zMin; zLoop<=zMax; zLoop++) {
//					if ((xLoop>=0 && xLoop<W) && (yLoop>=0 && yLoop<H) && (zLoop>=0 && zLoop<L)) { // is inside the image
//
//
//
//
//
//                        int currentIndex        = xyz2ind(xLoop, yLoop, zLoop, W, H);
//                        int currentValue        = (image3d[currentIndex] & 0xff);           // byte to int
//                        int calculatedValue     = 255;
////                                coneIntensity(
////                                xLoop, yLoop, zLoop,
////                                x, y, z, r,
////                                xPrev, yPrev, zPrev, rPrev,
////                                k, fg
////                                );
//
//                        if (calculatedValue>currentValue) {
//                            image3d[currentIndex] = (byte)calculatedValue;
//                            //System.out.print(" im["+currentIndex+"] "+calculatedValue+" byte: "+image3d[currentIndex]+"");
//                            count++;
//                        }
//
//
//
//
////						if (isCone(
////										  xLoop,  yLoop,  zLoop,
////										  x,      y,      z,      r,
////										  xPrev,  yPrev,  zPrev,  rPrev
////
////						)) {
////							//System.out.println("x: "+xLoop+" y: "+yLoop+" z: "+zLoop);
////							image3d[xyz2ind(xLoop, yLoop, zLoop, W, H)] = (byte)255;
////						}
//
//
//					}
//                }
//            }
//        }
//
//        return count;
//
//	}

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

        if (P12>0.1) {// don't check below some small distance
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

        if ( p21_x * c1_x + p21_y * c1_y + p21_z * c1_z < 0 )   return 0;
        if ( p12_x * c2_x + p12_y * c2_y + p12_z * c2_z < 0 )  return 0; // no need to calculate the value

        return 255;

//        float d2 = (c2_x * p12_x) + (c2_y * p12_y) + (c2_z * p12_z); // projection c2 on p12
//
//        double d = Math.sqrt( Math.pow(-p12_x * d2 + c2_x, 2) + Math.pow(-p12_y * d2 + c2_y, 2) + Math.pow(-p12_z * d2 + c2_z, 2) );
//
//        // interpolate radius based on r1 and r2
//        double r = (d2*r1 + (P12-d2)*r2)/P12;
//        double sigma =  r * k;
//
//        // use d, sigma to calculate the value and scale it with fg
//        int val = 0;
//        if (d<=1*sigma)
//            //val = (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
//            val = 255;// (int) (fg * Math.exp(-(d * d) / (2 * sigma * sigma)));
//        return val;

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

    ) {

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

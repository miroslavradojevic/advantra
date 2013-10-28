package generate;

import aux.ReadSWC;
import aux.ReadSWC;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:46 PM
 */
public class Generator3D {
    // will be used to generate 3D image stacks

//    private static int margin = 10;

    public static ImagePlus fromSWC(String swcPath, float k, float snr) {

        // read swc
        ReadSWC readerSWC = new ReadSWC(swcPath);   // read file into list

		System.out.println(readerSWC.nodes.size()+" ELEMENTS");

		// determine boundaries based on the swc locations
		int xMin=Integer.MAX_VALUE, xMax=Integer.MIN_VALUE,
                yMin=Integer.MAX_VALUE, yMax=Integer.MIN_VALUE,
                zMin=Integer.MAX_VALUE, zMax=Integer.MIN_VALUE,
                rMax=Integer.MIN_VALUE;


        for (int ii=0; ii<readerSWC.nodes.size(); ii++) {

            double readX, readY, readZ, readR;

            readX = Math.round(Math.ceil(readerSWC.nodes.get(ii)[0]));
            if (readX>xMax) xMax = (int)readX;
            readX = Math.round(Math.floor(readerSWC.nodes.get(ii)[0]));
            if (readX<xMin) xMin = (int)readX;

            readY = Math.round(Math.ceil(readerSWC.nodes.get(ii)[1]));
            if (readY>yMax) yMax = (int)readY;
            readY = Math.round(Math.floor(readerSWC.nodes.get(ii)[1]));
            if (readY<yMin) yMin = (int)readY;

            readZ = Math.round(Math.ceil(readerSWC.nodes.get(ii)[2]));
            if (readZ>zMax) zMax = (int)readZ;
            readZ = Math.round(Math.floor(readerSWC.nodes.get(ii)[2]));
            if (readZ<zMin) zMin = (int)readZ;

            readR = Math.round(Math.ceil(readerSWC.nodes.get(ii)[3]));
            if (readR>rMax) rMax = (int)readR;

        }

        System.out.println("X "+xMin+" -- "+xMax);
        System.out.println("Y "+yMin+" -- "+yMax);
        System.out.println("Z "+zMin+" -- "+zMax);
        System.out.println("R "+rMax);

        int shiftX = -xMin + rMax;
        int shiftY = -yMin + rMax;
        int shiftZ = -zMin + rMax;

        int W = xMax + rMax + shiftX + 1;
        int H = yMax + rMax + shiftY + 1;
        int L = zMax + rMax + shiftZ + 1;


//        System.out.println(""+sizeX+" , "+sizeY+" , "+sizeZ);
//        for (int z=0; z<3; z++) {
//            for (int y=0; y<3; y++) {
//                for (int x=0; x<4; x++) {
//                }
//            }
//        }

        byte[] outIm = new byte[W*H*L]; // allocate 8 bit image

		// loop through the list to form cones
        for (int coneId=1; coneId<readerSWC.nodes.size(); coneId++) {

            // set value at the node center location to 255
//            int locX = (int)Math.round(readerSWC.nodes.get(coneId)[0] - rMax + shiftX);
//            int locY = (int)Math.round(readerSWC.nodes.get(coneId)[1] - rMax + shiftY);
//            int locZ = (int)Math.round(readerSWC.nodes.get(coneId)[2] - rMax + shiftZ);

//            System.out.println(locX+", "+locY+", "+locZ);
//            outIm[xyz2ind(locX, locY, locZ, W, H)] = (byte)255;

            float x = readerSWC.nodes.get(coneId)[0] - rMax + shiftX;
            float y = readerSWC.nodes.get(coneId)[1] - rMax + shiftY;
            float z = readerSWC.nodes.get(coneId)[2] - rMax + shiftZ;
            float r = readerSWC.nodes.get(coneId)[3];// - rMax + shiftZ;

            int indexPrev = (int)Math.round(readerSWC.nodes.get(coneId)[4]);

            float xPrev = readerSWC.nodes.get(indexPrev)[0] - rMax + shiftX;
            float yPrev = readerSWC.nodes.get(indexPrev)[1] - rMax + shiftY;
            float zPrev = readerSWC.nodes.get(indexPrev)[2] - rMax + shiftZ;
            float rPrev = readerSWC.nodes.get(indexPrev)[3];// - rMax + shiftZ;

            drawCone(x,y,z,r,   xPrev,yPrev,zPrev,rPrev, outIm, W, H);

        }

        // draw the cone
		return new ImagePlus("out", toImageStack(outIm, W, H));

    }

    private void setValueAt() {

    }

	private static void drawCone(float x, float y, float z, float r, float xPrev, float yPrev, float zPrev, float rPrev, byte[] image3d, int W, int H) {

        // define the range to loop for plotting
        // x,y,z are expected to be valid 3d coordinates of the image stack
        // cone defined with (x1,y1,z1,r1),(x2,y2,z2,r2) should not be out of image (but that's taken care of when defining the size)

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
                for (int zLoop=0; zLoop<=zMax; zLoop++) {
                    if (isCone(
                            xLoop,  yLoop,  zLoop,
                            x,      y,      z,      r,
                            xPrev,  yPrev,  zPrev,  rPrev

                    )) {
                        image3d[xyz2ind(xLoop, yLoop, zLoop, W, H)] = (byte)255;
                    }
                }
            }
        }

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

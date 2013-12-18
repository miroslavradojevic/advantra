package detection2d;

import detection.Interpolator;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 10:37 AM
 * Computation unit for 2D processing
 */
public class Sphere2D {

    private static float    arcRes 	        = 0.7f;
    private static float 	arcNbhood       = 2*arcRes;
    private static float    samplingStep    = 0.7f;
    public static float     R_FULL 	        = 1.00f;
    public static float     T_HALF 	        = 0.50f;
    public static int 		weightStdRatioToD = 4;                  // could be a parameter

    private float 	radius;
    private float   scale;
    private float   neuronDiameter;
    private int     N;

    private int 	limR, limT;

    ArrayList<Float>        theta = new ArrayList<Float>(); 	    // list of elements (theta) covering the circle

    ArrayList<int[]> 		masks = new ArrayList<int[]>(); 	    // list of list indexes of the neighbours for each list element

    ArrayList<float[][]> 	offstXY = new ArrayList<float[][]>(); 	// list of filter offsets for each direction

    float[] weights;

    /*
    *********************************************************************
     */

    public Sphere2D(float neuronDiam, float scale) {

        this.radius 	= scale * neuronDiam;
        this.scale      = scale;
        this.neuronDiameter = neuronDiam;
        this.N 	= (int) Math.ceil( ( (2 * Math.PI * radius) / arcRes) );
        this.limR = (int) Math.ceil(R_FULL*neuronDiameter/samplingStep);    // how many to take radially with given sampling step
        this.limT = (int) Math.ceil(T_HALF*neuronDiameter/samplingStep);    //


        theta.clear();
        for (int i=0; i<N; i++) {
            theta.add( i * ( (float)(2*Math.PI) / N ) );
        }

        masks.clear();
        for (int ii=0; ii<theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            // loop the rest of the elements to see if they are in the angNeighbourhood range
            ArrayList<Integer> nbrs = new ArrayList<Integer>();
            for (int jj = 0; jj<theta.size(); jj++) {

                if (jj!=ii) {

                    // check
                    float theta_test = theta.get(jj);

                    float arc_btw = arcBetweenDirections(curr_theta, theta_test);
                    if (arc_btw<=arcNbhood) {
                        nbrs.add(jj);
                    }

                }

            }

            // convert list to regular array and add
            int[] nbrsArray = new int[nbrs.size()];
            for (int iii=0; iii<nbrs.size(); iii++) {
                nbrsArray[iii] = nbrs.get(iii);
            }

            masks.add(nbrsArray);

        }

        offstXY.clear();
        weights = new float[(2*limT+1)*(limR+1)];
        float sumWgt = 0;

        for (int ii = 0; ii<theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            /*
				form sampling (offset) points
			 */

            float[][] offsetsPerDirection = new float[(2*limT+1)*(limR+1)][2];

            int cnt = 0;

            for (int k=-limR; k<=0; k++) {

                for (int i=-limT; i<=limT; i++) {

//                    for (int j = -limT; j<=limT; j++) {

                        float px = i * samplingStep;
//                        float py = j * samplingStep;
                        float py = k * samplingStep;

                        offsetsPerDirection[cnt][0] = px;
                        offsetsPerDirection[cnt][1] = py;
//                        offsetsPerDirection[cnt][2] = pz;

                        //*** HERE IT DEFINES THE FILTER PROFILE WEIGHTS (only in first iteration, the rest are the same)
                        if (ii==0) {
                                float dstAxis = point2line(0, 0,        0, 1,       px, py);
                            weights[cnt] = (float) Math.exp(-(dstAxis*dstAxis)/(2*(neuronDiam/weightStdRatioToD)*(neuronDiam/weightStdRatioToD)));
                            sumWgt += weights[cnt];
                        }

                        cnt++;

//                    }

                }

            }

            /*
				transformation for offsets before adding
			 */
            transY(radius, offsetsPerDirection);
            //rotY(-phi+HalfPI, offsetsPerDirection);
            rotZ(curr_theta, offsetsPerDirection);
            offstXY.add(offsetsPerDirection); //store

        }

        /*
				normalize weights
	    */
        for (int iii=0; iii<weights.length; iii++) {
            weights[iii] /= sumWgt;
        }

    }

    public ImagePlus showSampling(){

        int DIM = 2 * (int) Math.ceil(Math.sqrt(Math.pow(neuronDiameter*T_HALF,2)+Math.pow(radius, 2))) + 1;
        int CX = DIM/2;
        int CY = CX;

        ImageStack isOut = new ImageStack(DIM, DIM);
        Overlay ov = new Overlay();

        for (int i=0; i<offstXY.size(); i++) {
            isOut.addSlice(new ByteProcessor(DIM, DIM));

            // center
            OvalRoi p = new OvalRoi(CX+0.5 - .5, CY+0.5 -.5, 1, 1);
            p.setPosition(i+1);
            p.setFillColor(Color.RED);
            p.setFillColor(Color.RED);
            ov.add(p);

            // sampling
            for (int i1=0; i1<offstXY.get(i).length; i1++) {

                float offX = offstXY.get(i)[i1][0];
                float offY = offstXY.get(i)[i1][1];

                PointRoi p1 = new PointRoi(CX+offX+.5, CY+offY+.5);
                p1.setPosition(i+1);
                ov.add(p1);

            }

        }

        ImagePlus outIm = new ImagePlus("offsets", isOut);
        outIm.setOverlay(ov);
        return outIm;

    }

    public ImagePlus showWeights(){

        float sum = 0;
        for (int k=0; k<weights.length; k++) sum += weights[k];
        return new ImagePlus("weights", new FloatProcessor((2*limT+1), (limR+1), weights));

    }

    public int getProfileLength()
    {
        return offstXY.size();
    }

    public int getOuterRadius(){
        return (int) Math.ceil(Math.sqrt(Math.pow(neuronDiameter*T_HALF, 2)+Math.pow(radius, 2)));
    }

    public short extractProfile(int profileIdx, float atX, float atY, float[][] _inimg_xy) {
        // one element filter output (indexed with profileIdx)

        float value = 0;

        for (int offsetIdx=0; offsetIdx<offstXY.get(profileIdx).length; offsetIdx++) {

            float x_offs_pix = offstXY.get(profileIdx)[offsetIdx][0];
            float y_offs_pix = offstXY.get(profileIdx)[offsetIdx][1];
//            float z_offs_lay = offstXY.get(profileIdx)[offsetIdx][2] / zDist; // convert to layers

            float imgVal = Interpolator.interpolateAt(atX + x_offs_pix, atY + y_offs_pix, _inimg_xy); // , atZ + z_offs_lay
            value += weights[offsetIdx] * imgVal;

        }

        return  (short) ((int) ((value/255f)*65535f));//&  0xffff

    }

	public void peakCoords_4xXY(short[] _profile, ){

	}

    /*
    *********************************************************************
     */

    private void rotZ(float ang, float[][] coords) {
        for (int i=0; i<coords.length; i++) {
            float x_temp = coords[i][0];
            float y_temp = coords[i][1];
            coords[i][0] = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
            coords[i][1] = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
        }
    }


    private void transY(float ty, float[][] coords){
        for (int i=0; i<coords.length; i++){
            coords[i][1] += ty;
        }
    }


    private float point2line(float n1x, float n1y,  // float n1z,
                             float n2x, float n2y,  // float n2z,
                             float px,  float py    //, float pz
    ) {

        float d = 0;

        double[] p_b = new double[2];

        //double[] n21 = new double[3];
        float n21Len = (float) Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2)); // +Math.pow(n2z-n1z,2)
        float n21x = (n2x-n1x)/n21Len;
        float n21y = (n2y-n1y)/n21Len;
//        float n21z = (n2z-n1z)/n21Len;

        float proj = (px - n1x) * n21x + (py - n1y) * n21y; // + (pz - n1z) * n21z; // dot prod

        p_b[0] = -(px - n1x) + proj * n21x;
        p_b[1] = -(py - n1y) + proj * n21y;
//        p_b[2] = -(pz - n1z) + proj * n21z;

        return (float) Math.sqrt(p_b[0]*p_b[0] + p_b[1]*p_b[1]); // + p_b[2]*p_b[2]

    }

    private float arcBetweenDirections(float theta1, float theta2){

        float x1 = getX(radius, theta1);
        float y1 = getY(radius, theta1);

        float x2 = getX(radius, theta2);
        float y2 = getY(radius, theta2);

        return radius * (float) Math.acos( (x1*x2+y1*y2)/(radius * radius) );

    }

    private float getX(float r, float theta){return r * (float) Math.cos(theta);}

    private float getY(float r, float theta){return r * (float) Math.sin(theta);}

}

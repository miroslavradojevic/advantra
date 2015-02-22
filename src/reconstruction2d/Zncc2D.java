package reconstruction2d;

import aux.Interpolator;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 15-2-15.
 */
public class Zncc2D {

    private static float arcRes = 1.5f;
    private static float samplingStep = .7f;

    private float   radius;
    private float   diameter;
    private int     N;
    private int 	limR, limT;

    private float[] sigmas;
    private float sigma_step = .5f;
    private float sg_min = 1.0f;

    public static ArrayList<Float>          theta = new ArrayList<Float>(); 	        // list of elements (theta) covering the circle
    private static ArrayList<float[][]> 	offstXY = new ArrayList<float[][]>(); 	    // list of filter offsets for each direction

    private float[][]   tplt;           // template
    private float[]     tplt_avg;       // template average
    private float[][]   tplt_hat;       // template - template average
    private float[]     tplt_hat_sum_2; // sum(template- template average)^2

    private static float TWO_PI = (float) (2 * Math.PI);

    public Zncc2D(float _radius)
    {

        radius = _radius;
        diameter = 2*radius;
        N 	= (int) Math.ceil(((TWO_PI * radius)/arcRes));
        limT = (int) Math.ceil(radius/samplingStep);    // transversal sampling limits
        limR = 2 * limT + 1; // how many to take radially with given sampling step

        int rr = (int) Math.ceil(_radius);
        rr = (rr<1)? 1 : rr;

        // sigmas define
        int cnt = 0;
        for (float sg = sg_min; sg <= rr; sg+=sigma_step) cnt++;
        sigmas = new float[cnt];
        cnt = 0;
        for (float sg = sg_min; sg <= rr; sg+=sigma_step) sigmas[cnt++] = sg;

        theta.clear();
        for (int i=0; i<N; i++) theta.add(i * (TWO_PI/N));

        offstXY.clear();
//        float sumWgt = 0;
        for (int ii = 0; ii<theta.size(); ii++) {

            // define offsetsPerDirection, weights (at one direction only)
            float[][] offsetsPerDirection = new float[limR*limR][2]; // (limR+1)

            cnt = 0;
            for (int k=-limT; k<=limT; k++) {
                for (int i=-limT; i<=limT; i++) {
//                    for (int j = -limT; j<=limT; j++) {

                    float px = i * samplingStep;
//                        float py = j * samplingStep;
                    float py = k * samplingStep;

                    offsetsPerDirection[cnt][0] = px;
                    offsetsPerDirection[cnt][1] = py;
//                        offsetsPerDirection[cnt][2] = pz;

                    cnt++;

//                    }
                }
            }

            // transformation for offsets before adding
//            transY(radius, offsetsPerDirection);
            //rotY(-phi+HalfPI, offsetsPerDirection);
            rotZ(theta.get(ii), offsetsPerDirection);
            offstXY.add(offsetsPerDirection); //store

        }

        tplt            = new float[sigmas.length][limR*limR]; // weights at different sigmas only for the first theta, the rest are the same
        tplt_avg        = new float[sigmas.length];
        tplt_hat        = new float[sigmas.length][limR*limR];
        tplt_hat_sum_2  = new float[sigmas.length];

        for (int sig_idx = 0; sig_idx < sigmas.length; sig_idx++) {

            float minWgt = Float.POSITIVE_INFINITY;
            float maxWgt = Float.NEGATIVE_INFINITY;



            cnt = 0;
            for (int k=-limT; k<=limT; k++) {
                for (int i = -limT; i <= limT; i++) {

                    float px = i * samplingStep;
                    float py = k * samplingStep;
                    float dstAxis = point2line(0,0,        0,1,       px,py);

                    tplt[sig_idx][cnt] = (float) Math.exp(-(dstAxis*dstAxis)/(2*Math.pow(sigmas[sig_idx],2)));


                    if (tplt[sig_idx][cnt]>maxWgt) maxWgt = tplt[sig_idx][cnt];
                    if (tplt[sig_idx][cnt]<minWgt) minWgt = tplt[sig_idx][cnt];

                    cnt++;
                }
            }

            // normalize template values (min/max) & calculate average of normalized template values
            tplt_avg[sig_idx] = 0;
            for (int iii=0; iii<tplt[sig_idx].length; iii++) {
                tplt[sig_idx][iii] = (tplt[sig_idx][iii] - minWgt) / (maxWgt - minWgt);
                tplt_avg[sig_idx] += tplt[sig_idx][iii];
            }
            tplt_avg[sig_idx] /= (float) tplt[sig_idx].length;

            // calculate mean subtracted template values and their sum of squares
            tplt_hat_sum_2[sig_idx] = 0;
            for (int iii=0; iii<tplt[sig_idx].length; iii++) {
                tplt_hat[sig_idx][iii] = tplt[sig_idx][iii] - tplt_avg[sig_idx];
                tplt_hat_sum_2[sig_idx] += Math.pow(tplt_hat[sig_idx][iii], 2);
            }

        }

    }

    private float point2line(float n1x, float n1y,  // float n1z,
                             float n2x, float n2y,  // float n2z,
                             float px,  float py    //, float pz
    )
    {

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

    private void rotZ(
            float ang,
            float[][] coords
    )
    {
        for (int i=0; i<coords.length; i++) {
            float x_temp = coords[i][0];
            float y_temp = coords[i][1];
            coords[i][0] = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
            coords[i][1] = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
        }
    }

    public ImagePlus getSampling()
    {

        int DIM = 2 * (int) Math.ceil(Math.sqrt(Math.pow(radius,2)+Math.pow(radius, 2))) + 1;
        float CX = DIM/2f;
        float CY = CX;

        ImageStack isOut = new ImageStack(DIM, DIM);
        Overlay ov = new Overlay();

        for (int i=0; i<offstXY.size(); i++) {

            isOut.addSlice(new ByteProcessor(DIM, DIM));

            // center
            OvalRoi p = new OvalRoi(CX+.5-.25, CY+.5-.25, .5, .5);
            p.setPosition(i + 1);
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

            // arrow direction
            Line l = new Line(CX+.5, CY+.5, CX+.5+radius*(-Math.sin(theta.get(i))), CY+.5+radius*Math.cos(theta.get(i)));
            l.setFillColor(new Color(1,0,0,.5f));
            l.setStrokeColor(new Color(1,0,0,.5f));
            l.setStrokeWidth(0.25);
            l.setPosition(i+1);
            ov.add(l);
        }

        ImagePlus outIm = new ImagePlus("offsets", isOut);
        outIm.setOverlay(ov);
        return outIm;

    }

    public ImagePlus getTemplates()
    {

        ImageStack is_weights = new ImageStack(limR, limR);

        for (int sig_idx = 0; sig_idx < sigmas.length; sig_idx++)
            is_weights.addSlice("sigma="+ IJ.d2s(sigmas[sig_idx],2), new FloatProcessor(limR, limR, tplt[sig_idx]));
        
        return new ImagePlus("weights", is_weights);

    }

    private void transY(
            float ty,
            float[][] coords
    )
    {
        for (int i=0; i<coords.length; i++){
            coords[i][1] += ty;
        }
    }

    // do zncc calculations at one location index
    public void extract(
            int         _idx,           // index from the foreground list
            int[][]     _i2xy,          // foreground map using Masker2D
            float[][]   _inimg_xy,      // image in array form
            float[]     _i2zncc,        // output - list of correlations    (from Profiler2D)
            float[]     _i2sigma,       // output - index of the sigma      (from Profiler2D)
            float[][]   _i2vxy          // output - vector, 2D direction    (from Profiler2D)
    )
    {

        int at_x = _i2xy[_idx][0];
        int at_y = _i2xy[_idx][1];

        // outputs...
        _i2zncc[_idx]   = Float.NEGATIVE_INFINITY; // because we look for the highest
        _i2sigma[_idx]  = Float.NaN;
        _i2vxy[_idx][0] = Float.NaN;
        _i2vxy[_idx][1] = Float.NaN;

        float curr_zncc;
        float[] vals = new float[limR*limR]; // can be allocated outside BUT independently for each thread
        float vals_avg = 0;

        for (int dir_idx = 0; dir_idx < offstXY.size(); dir_idx++) {

//            Arrays.fill(vals, 0);
            // extract values at this direction by interpolating at real 2d locations
            // min/max normalized values - maybe not necessary!
            float vals_min = Float.POSITIVE_INFINITY;
            float vals_max = Float.NEGATIVE_INFINITY;
            for (int cnt = 0; cnt < offstXY.get(dir_idx).length; cnt++) {
                float dx = offstXY.get(dir_idx)[cnt][0];
                float dy = offstXY.get(dir_idx)[cnt][1];
                vals[cnt] = Interpolator.interpolateAt(at_x+dx, at_y+dy, _inimg_xy);
                if (vals[cnt]<vals_min) vals_min = vals[cnt];
                if (vals[cnt]>vals_max) vals_max = vals[cnt];
            }

            vals_avg = 0;
            for (int cnt = 0; cnt < offstXY.get(dir_idx).length; cnt++) {
                vals[cnt] = (vals[cnt]-vals_min)/(vals_max-vals_min);
                vals_avg += vals[cnt];
            }
            vals_avg /= offstXY.get(dir_idx).length;

            for (int sigma_idx = 0; sigma_idx <sigmas.length; sigma_idx++) {

                // calculate zncc
                curr_zncc = zncc(vals, vals_avg, tplt_hat[sigma_idx], tplt_hat_sum_2[sigma_idx]);

                if (curr_zncc>_i2zncc[_idx]) {
                    _i2zncc[_idx] = curr_zncc;
                    _i2sigma[_idx] = sigmas[sigma_idx];
                    _i2vxy[_idx][0] = (float) - Math.sin(theta.get(dir_idx));
                    _i2vxy[_idx][1] = (float)   Math.cos(theta.get(dir_idx));
                }

            }
        }

    }

    private float zncc(
            float[] v,
            float   v_avg,
            float[] tmplt_hat,
            float   tmplt_hat_sum_sqr
    )
    {
        // TODO: place the reference here

        float num = 0;
        float den = 0;

        for (int i = 0; i < v.length; i++) {
            num += (v[i] - v_avg) * tmplt_hat[i];
            den += Math.pow(v[i] - v_avg, 2);
        }

        return (float) (num / Math.sqrt(den * tmplt_hat_sum_sqr));

    }

}

package reconstruction;

import aux.Interpolator;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.ArrayList;
import java.util.Random;

/**
 * module that calculates correlation in different directions
 * this time sampoing of different directions in 3d is done using the
 *
 * Created by miroslav on 1-3-15.
 */
public class Zncc3D {

//    private static int      N;
private static float arcRes = 1.5f;
    private static float    samplingStep = 1f;
    private float           radius;
    private float           zDist;
    public int 	        limR, limT;

    private float[] sigmas;
    private float   sigma_step = .5f;
    private float   sg_min = 0.5f;


    public static ArrayList<float[]> u                  = new ArrayList<float[]>(); // list of unit directions ux, uy, uz randomly distributed on the sphere
    private static ArrayList<float[][]> 	offstXYZ    = new ArrayList<float[][]>(); 	    // list of filter offsets for each direction u

    private float[][]       tplt;           // template
    private float[]         tplt_avg;       // template average
    private float[][]       tplt_hat;       // template - template average
    private float[]         tplt_hat_sum_2; // sum(template- template average)^2

    public ImagePlus       ipviz;

    public Zncc3D(float _radius, float _zDist) {

        // _show_weights will do additional allocation to visualize the stack with weights

        radius = _radius;
//        N = _N;
        zDist = _zDist;

        limT = (int) Math.ceil(radius/samplingStep);
        limR = 2 * limT + 1;

        int rr = (int) Math.ceil(_radius);
        rr = (rr<1)? 1 : rr;

        // sigmas define
        int cnt = 0;
        for (float sg = sg_min; sg <= rr/2; sg+=sigma_step) cnt++;
        sigmas = new float[cnt];
        cnt = 0;
        for (float sg = sg_min; sg <= rr/2; sg+=sigma_step) sigmas[cnt++] = sg;

        u.clear();
        // form N vxyz (N directions), sphere point picking
        // saaf kuijlars method to uniformly sample on the sphere
//        u = fullSpherePts(N);

        // http://mathworld.wolfram.com/SpherePointPicking.html   Marsaglia (1972)  method
//        Random rg = new Random();
//        cnt = 0;
//        while (cnt<N) {
//            float x1 = 2 * rg.nextFloat() - 1;
//            float x2 = 2 * rg.nextFloat() - 1;
//            if (Math.pow(x1,2)+Math.pow(x2,2)<1) {
//                float x = (float) (2 * x1 * Math.sqrt(1 - Math.pow(x1,2) - Math.pow(x2,2)));
//                float y = (float) (2 * x2 * Math.sqrt(1 - Math.pow(x1,2) - Math.pow(x2,2)));
//                float z = (float) (1 - 2 * (Math.pow(x1,2)+Math.pow(x2,2)));
////                System.out.println(IJ.d2s(x,1) + ", " + IJ.d2s(y,1) + ", " + IJ.d2s(z,1));
//                u.add(new float[]{x, y, z});
//                cnt++;
//            }
//        }

        // none of those two was optimal - need to focus on 2d layers and not oversample as it is problem for the profiler

        int zlim = (int) Math.floor(radius/zDist);
        System.out.println("+/-" + zlim + " neighbouring layers considered in Zncc3D");

        float x_temp = 1;
        float y_temp = 0;

        for (int dz = -zlim; dz <= zlim; dz++) {

            float rds = (float) (radius * Math.cos(  ((dz * zDist) / radius) * (Math.PI/2)   ));

            for (float darc = 0; darc < 2*rds*Math.PI; darc+=arcRes) {

                float ang = darc/rds;

                float vx = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
                float vy = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
                float vz = dz;

                // normalize to make it direction
                float vnorm = (float) Math.sqrt(vx*vx+vy*vy+vz*vz);
                vx /= vnorm;
                vy /= vnorm;
                vz /= vnorm;

                u.add(new float[]{vx, vy, vz});


            }

        }

        System.out.println(u.size() + " directions uzed by Zncc3D");

        // form offsets for each direction

        offstXYZ.clear();

        float[][] R  = new float[3][3]; // transformation
        float[] Rvec = new float[3]; // vector after transformation

        for (int ii = 0; ii<u.size(); ii++) {

            float ux = u.get(ii)[0];
            float uy = u.get(ii)[1];
            float uz = u.get(ii)[2];

            rotation_matrix(1,0,0, ux,uy,uz, R); // R contains transformation from 1,0,0 to generated random direction
            rotation_apply(R, 0, 1, 0, Rvec);

            float vx = Rvec[0];
            float vy = Rvec[1];
            float vz = Rvec[2];

            rotation_apply(R, 0, 0, 1, Rvec);

            float wx = Rvec[0];
            float wy = Rvec[1];
            float wz = Rvec[2];

            // define offsetsPerDirection, weights (at one direction only)
            float[][] offsetsPerDirection = new float[limR*limR*limR][3]; // (limR+1)

            cnt = 0;

            for (int i=-limT; i<=limT; i++) {
                for (int j = -limT; j<=limT; j++) {
                    for (int k=-limT; k<=limT; k++) {

                        float A = i * samplingStep;
                        float B = j * samplingStep;
                        float C = k * samplingStep;

                        // Au+Bv+Cw is the coordinate in the new coordinate system
                        // (ux,uy,uz), (vx,vy,vz), (wx,xy,xz) are unit vectors of the new coordinate system
                        // apply rotation here immediately
                        offsetsPerDirection[cnt][0] = A*ux + B*vx + C*wx; // x
                        offsetsPerDirection[cnt][1] = A*uy + B*vy + C*wy; // y
                        offsetsPerDirection[cnt][2] = A*uz + B*vz + C*wz; // z

                        // scale z component if necessary
                        offsetsPerDirection[cnt][2] /= zDist;

                        cnt++;

                    }
                }
            }

            uz /= zDist;
            float aa = (float) Math.sqrt(ux*ux+uy*uy+uz*uz);
            u.get(ii)[0] = (float) (ux/aa);
            u.get(ii)[1] = (float) (uy/aa);
            u.get(ii)[2] = (float) (uz/aa);

            offstXYZ.add(offsetsPerDirection); //store

        }

        // fill in the templates
        tplt            = new float[sigmas.length][limR*limR*limR]; // weights at different sigmas only for the first theta, the rest are the same
        tplt_avg        = new float[sigmas.length];
        tplt_hat        = new float[sigmas.length][limR*limR*limR];
        tplt_hat_sum_2  = new float[sigmas.length];

        ImageStack isviz = new ImageStack(limR, limR);
        float[][][]  viz = new float[limR][limR][limR];

        for (int sig_idx = 0; sig_idx < sigmas.length; sig_idx++) {

            float minWgt = Float.POSITIVE_INFINITY;
            float maxWgt = Float.NEGATIVE_INFINITY;

            cnt = 0;

            for (int i=-limT; i<=limT; i++) {
                for (int j = -limT; j<=limT; j++) {
                    for (int k=-limT; k<=limT; k++) {

                        float A = i * samplingStep; // ~ slice
                        float B = j * samplingStep; // ~ row
                        float C = k * samplingStep; // ~ column

                        // distance from the axis in pix - using B,C
                        float dstAxis = (float) Math.sqrt(B*B+C*C); // distance in pixels
                        tplt[sig_idx][cnt] = (float) Math.exp(-(dstAxis*dstAxis)/(2*Math.pow(sigmas[sig_idx],2)));
                        if (tplt[sig_idx][cnt]>maxWgt) maxWgt = tplt[sig_idx][cnt];
                        if (tplt[sig_idx][cnt]<minWgt) minWgt = tplt[sig_idx][cnt];

                        viz[i+limT][j+limT][k+limT] = tplt[sig_idx][cnt];

                        cnt++;

                    }
                }
            }

            // add viz by layers
            for (int i = 0; i < limR; i++) isviz.addSlice(new FloatProcessor(viz[i].clone()));

//            System.out.println(isviz.getSize() + " is size " + limR + " , " + sigmas.length);

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

        ipviz = new ImagePlus("weights", isviz);
        ipviz.setDimensions(1, limR, sigmas.length);
//        ipviz.show();

    }

    private void rotation_matrix(
            float a1, float a2, float a3,
            float b1, float b2, float b3,
            float[][] R
    )
    {

        // v is cross product of (a1, a2, a3) and (b1, b2, b3)
        float v1 = a2*b3 - b2*a3;
        float v2 = -(a1*b3-b1*a3);
        float v3 = a1*b2-b1*a2;

        float tt = (1-(a1*b1+a2*b2+a3*b3))/(v1*v1+v2*v2+v3*v3);

        // from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        R[0][0] = 1 + 0     + tt * (-v3*v3-v2*v2);
        R[0][1] = 0 + (-v3) + tt * (v1*v2);
        R[0][2] = 0 + (v2)  + tt * (-v1*v3);

        R[1][0] = 0 + (v3)  + tt * (v1*v2);
        R[1][1] = 1 + 0     + tt * (-v3*v3-v1*v1);
        R[1][2] = 0 + (-v1) + tt * (v2*v3);

        R[2][0] = 0 + (-v2) + tt * (v1*v3);
        R[2][1] = 0 + (v1)  + tt * (v2*v3);
        R[2][2] = 1 + 0     + tt * (-v2*v2-v1*v1);

    }

    private void rotation_apply(
            float[][] R,
            float v1, float v2, float v3,
            float[] Rvec
    ){

        Rvec[0] = R[0][0]*v1 + R[0][1]*v2 + R[0][2]*v3;
        Rvec[1] = R[1][0]*v1 + R[1][1]*v2 + R[1][2]*v3;
        Rvec[2] = R[2][0]*v1 + R[2][1]*v2 + R[2][2]*v3;

    }

    public void exportTemplates(String out_dir_path){
        IJ.saveAs(ipviz, "Tiff", out_dir_path + "templates.tif");
    }

    public void exportFilterSampling(String out_dir_path){

        File f = new File(out_dir_path);
        if (!f.exists()) f.mkdirs();

        for(File file: f.listFiles()) file.delete();

        for (int i = 0; i < offstXYZ.size(); i++) {

            String det_path = f.getAbsolutePath() + File.separator + "" + IJ.d2s(u.get(i)[0],1) + "," + IJ.d2s(u.get(i)[1],1) + "," + IJ.d2s(u.get(i)[2],1) + ".swc";

            PrintWriter logWriter = null;
            try {
                logWriter = new PrintWriter(det_path);
                logWriter.print("");
                logWriter.close();
            } catch (FileNotFoundException ex) {}

            try {
                logWriter = new PrintWriter(new BufferedWriter(new FileWriter(det_path, true)));
                logWriter.println("# ");
            } catch (IOException e) {}

            for (int j = 0; j < offstXYZ.get(i).length; j++)
                logWriter.println(String.format(
                        "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                        (j + 1), 2, offstXYZ.get(i)[j][0], offstXYZ.get(i)[j][1], offstXYZ.get(i)[j][2], .1f, -1));

            // center
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    (offstXYZ.get(i).length + 1), 1, 0f, 0f, 0f, .5f, -1));


            // main direction (line)
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    (offstXYZ.get(i).length + 2), 3, 0f, 0f, 0f, .1f, -1));
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    (offstXYZ.get(i).length + 3), 3, radius*u.get(i)[0], radius*u.get(i)[1], radius*u.get(i)[2], .05f, (offstXYZ.get(i).length + 2)));

            logWriter.close();

            System.out.println("done. "+det_path);
        }

        String det_path = f.getAbsolutePath() + File.separator + "" + "directions" + ".swc";

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(det_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(det_path, true)));
            logWriter.println("# ");
        } catch (IOException e) {}


        int cnt = 1;
        for (int i = 0; i < u.size(); i++) {

            // main direction (line)
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt, 3, 0f, 0f, 0f, .1f, -1));
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+1, 3, radius*u.get(i)[0], radius*u.get(i)[1], radius*u.get(i)[2], .05f, cnt));

            cnt+=2;

        }

        logWriter.close();

    }

    public ArrayList<float[]> fullSpherePts(int N){

        ArrayList<float[]> points = new ArrayList<float[]>(N);// float[N][3];

        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < N; k++) {

            h_k = -1 + 2 * (double)k/(N-1); // -1 : 1

            theta_k = Math.acos(h_k);

            if(k==0 || k==(N-1)){

                phi_k   = 0;
                phi_k_1 = 0;

            }
            else{

                phi_k = phi_k_1 + 3.6 / ( Math.sqrt(N) * Math.sqrt(1-h_k*h_k));
                phi_k_1 = phi_k;

            }


            float px = (float) (Math.sin(theta_k) * Math.cos(phi_k));
            float py = (float) (Math.sin(theta_k) * Math.sin(phi_k));
            float pz = (float) Math.cos(theta_k);

            points.add(new float[]{px, py, pz});

        }

        return points;

    }





    // do zncc calculations at one location index
    public void extract(
            int         _idx,           // index from the foreground list
            int[][]     _i2xyz,         // foreground map
            float[][][] _inimg_xyz,     // image in array form

            float[]     _i2zncc,        // output - list of correlations    (from Profiler3D)
            float[]     _i2sigma,       // output - index of the sigma      (from Profiler3D)
            float[][]   _i2vxyz,        // output - vector, 3D direction    (from Profiler3D)

            // auxilliary - storage for the values, needs to be properly allocated length=limR*limR*limR
            float[] vals
    )
    {

        int at_x = _i2xyz[_idx][0];
        int at_y = _i2xyz[_idx][1];
        int at_z = _i2xyz[_idx][2];

        // outputs...
        _i2zncc[_idx]   = Float.NEGATIVE_INFINITY; // because we look for the highest
        _i2sigma[_idx]  = Float.NaN;
        _i2vxyz[_idx][0] = Float.NaN;
        _i2vxyz[_idx][1] = Float.NaN;
        _i2vxyz[_idx][2] = Float.NaN;

        float curr_zncc;
        float vals_avg = 0;

        for (int dir_idx = 0; dir_idx < offstXYZ.size(); dir_idx++) {

//            Arrays.fill(vals, 0);
            // extract values at this direction by interpolating at real 2d locations
            float vals_min = Float.POSITIVE_INFINITY;
            float vals_max = Float.NEGATIVE_INFINITY;
            for (int cnt = 0; cnt < offstXYZ.get(dir_idx).length; cnt++) {

                float dx = offstXYZ.get(dir_idx)[cnt][0];
                float dy = offstXYZ.get(dir_idx)[cnt][1];
                float dz = offstXYZ.get(dir_idx)[cnt][2];

                vals[cnt] = Interpolator.interpolateAt(at_x + dx, at_y + dy, at_z + dz, _inimg_xyz);
                if (vals[cnt]<vals_min) vals_min = vals[cnt];
                if (vals[cnt]>vals_max) vals_max = vals[cnt];
            }

            // todo perhaps it is not necessary to normalize here
            vals_avg = 0;
            for (int cnt = 0; cnt < offstXYZ.get(dir_idx).length; cnt++) {
                vals[cnt] = (vals[cnt]-vals_min)/(vals_max-vals_min);
                vals_avg += vals[cnt];
            }
            vals_avg /= offstXYZ.get(dir_idx).length;
            // todo perhaps it is not necessary to normalize here

            for (int sigma_idx = 0; sigma_idx <sigmas.length; sigma_idx++) {

                // calculate zncc
                curr_zncc = zncc(vals, vals_avg, tplt_hat[sigma_idx], tplt_hat_sum_2[sigma_idx]);

                if (curr_zncc>_i2zncc[_idx]) {
                    _i2zncc[_idx] = curr_zncc;
                    _i2sigma[_idx] = sigmas[sigma_idx];
                    _i2vxyz[_idx][0] = u.get(dir_idx)[0];
                    _i2vxyz[_idx][1] = u.get(dir_idx)[1];
                    _i2vxyz[_idx][2] = u.get(dir_idx)[2];
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

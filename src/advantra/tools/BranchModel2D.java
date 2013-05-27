package advantra.tools;

import advantra.general.ArrayHandling;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/16/13
 * Time: 10:12 AM
 * To change this template use File | Settings | File Templates.
 */

public class BranchModel2D {

    /*
    creates artificial 2D branch, as 8bit image
     */
    public byte[] values; // to store the values
    public byte[] mask;   // object versus background mask

    public int h, w;

    // 4 points define the branch point: pc (midpoint), p1, p2, p3
    public int[] p1, p2, p3, pc;

    public double[] v1, v2, v3;

    // standard deviation of the radius
    public double rstd;

    // 3 angles define the branch point
    public int  alfa12, alfa23, alfa31;

    // 3 lengths define branches taking different directions
    public double l1, l2, l3;

    // background, foreground
    int bg, fg;

    public static float TwoPi = (float) (2 * Math.PI);

    int minAngle = 30;
	int maxAngle = 150;
    double toFix = 0.5;

    public BranchModel2D(int height, int width)
    {

        h = height;
        w = width;

        pc = new int[2];

        pc[0] = w/2;
        pc[1] = h/2;

        p1 = new int[2];
        p2 = new int[2];
        p3 = new int[2];

        v1 = new double[2];
        v2 = new double[2];
        v3 = new double[2];

        values  = new byte[h*w];
        mask    = new byte[h*w];

    }

    public byte[] generateRandomBranch(int bgBias, int bgRange, int fgBias, int fgRange)
    {

        Random gen = new Random();

        // random parameters
        bg = gen.nextInt(bgRange) + bgBias;
        for (int i = 0; i < values.length; i++) values[i] = (byte) bg;
        for (int i = 0; i < mask.length; i++) mask[i] = (byte) 0;

        fg = gen.nextInt(fgRange) + fgBias;

        boolean angCorrect = false;
        while (!angCorrect){
            alfa12 = gen.nextInt(maxAngle); // avoid angles higher than maxAngle degrees - they don't make sense
			alfa23 = gen.nextInt(maxAngle);

			angCorrect = alfa12>=minAngle && alfa23>=minAngle && (360-alfa12-alfa23)>=minAngle && (360-alfa12-alfa23)<maxAngle;

        }

		//alfa12 = minAngle + gen.nextInt(maxAngle-minAngle);
		//alfa23 = alfa12 + minAngle + gen.nextInt(360 - (alfa12 + minAngle));


        double maxLength = Math.min(0.8*w/2, 0.8*h/2);

        l1 = Math.random() * maxLength * (1-toFix) + maxLength * toFix;
        l2 = Math.random() * maxLength * (1-toFix) + maxLength * toFix;
        l3 = Math.random() * maxLength * (1-toFix) + maxLength * toFix;

		int startAngle = gen.nextInt(360);

		float startAngleRad = (startAngle/360f)*TwoPi;
        p1[0] = pc[0] + (int) (Math.cos(startAngleRad) * l1);
        p1[1] = pc[1] + (int) (Math.sin(startAngleRad) * l1);

		float startAngleRad1 = ((startAngle+alfa12)/360f)*TwoPi;
        p2[0] = pc[0] + (int) (Math.cos(startAngleRad1) * l2);
        p2[1] = pc[1] + (int) (Math.sin(startAngleRad1) * l2);

		// third will be between
		// choose random angle between startAngle+180 and startAngle+alfa12+180
		//int endAngle = startAngle + 180 + gen.nextInt(alfa12);
		float endAngleRad = ((startAngle+alfa12+alfa23)/360f)*TwoPi;

        p3[0] = pc[0] + (int) (Math.cos(endAngleRad) * l3); // row
        p3[1] = pc[1] + (int) (Math.sin(endAngleRad) * l3); // col

        /*
		 * memorize the unit directions as well
		 */

        double 	v_norm;

        // v1
        v1[0] 		= p1[0] - pc[0];
        v1[1] 		= p1[1] - pc[1];
        v_norm 	= Math.sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
        v1[0]        = v1[0] / v_norm;
        v1[1]        = v1[1] / v_norm;

        // v2
        v2[0] 		= p2[0] - pc[0];
        v2[1] 		= p2[1] - pc[1];
        v_norm 	= Math.sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
        v2[0]        = v2[0] / v_norm;
        v2[1]        = v2[1] / v_norm;

        // v3
        v3[0] 		= p3[0] - pc[0];
        v3[1] 		= p3[1] - pc[1];
        v_norm 	= Math.sqrt(v3[0]*v3[0] + v3[1]*v3[1]);
        v3[0]        = v3[0] / v_norm;
        v3[1]        = v3[1] / v_norm;

        rstd = gen.nextInt(2) + 1;

        writeLineBetween(pc, p1);
        writeLineBetween(pc, p2);
        writeLineBetween(pc, p3);

/*        // bor1
        int bor_row = (int) ((p1[0]+pc[0])/2 + 2*rstd*(-v1[1]));
        int bor_col = (int) ((p1[1]+pc[1])/2 + 2*rstd*v1[0]);
		values[bor_col+bor_row*w] = (byte) 255;
        bor_row = (int) ((p1[0]+pc[0])/2 - 2*rstd*(-v1[1]));
        bor_col = (int) ((p1[1]+pc[1])/2 - 2*rstd*v1[0]);
        values[bor_col+bor_row*w] = (byte) 255;

        // bor2
        bor_row = (int) ((p2[0]+pc[0])/2 + 2*rstd*(-v2[1]));
        bor_col = (int) ((p2[1]+pc[1])/2 + 2*rstd*v2[0]);
        values[bor_col+bor_row*w] = (byte) 255;
        bor_row = (int) ((p2[0]+pc[0])/2 - 2*rstd*(-v2[1]));
        bor_col = (int) ((p2[1]+pc[1])/2 - 2*rstd*v2[0]);
        values[bor_col+bor_row*w] = (byte) 255;

        // bor3
        bor_row = (int) ((p3[0]+pc[0])/2 + 2*rstd*(-v3[1]));
        bor_col = (int) ((p3[1]+pc[1])/2 + 2*rstd*v3[0]);
        values[bor_col+bor_row*w] = (byte) 255;
        bor_row = (int) ((p3[0]+pc[0])/2 - 2*rstd*(-v3[1]));
        bor_col = (int) ((p3[1]+pc[1])/2 - 2*rstd*v3[0]);
        values[bor_col+bor_row*w] = (byte) 255;*/

        return values;
    }

    private void writeLineBetween(
            int[] point_1,
            int[] point_2
    )
    {

        // correct row (0-index) if they're out
        point_1[0] = (point_1[0]<0)?0:point_1[0];	point_1[0] = (point_1[0]>=h)?(h-1):point_1[0];
        point_2[0] = (point_2[0]<0)?0:point_2[0];	point_2[0] = (point_2[0]>=h)?(h-1):point_2[0];

        // correct col (1-index) if they're out
        point_1[1] = (point_1[1]<0)?0:point_1[1];	point_1[1] = (point_1[1]>=w )?(w -1):point_1[1];
        point_2[1] = (point_2[1]<0)?0:point_2[1];	point_2[1] = (point_2[1]>=w )?(w -1):point_2[1];

		/*
		 * unit orientation
		 */

        double[] n 	= new double[2];
        n[0] 		= point_2[0] - point_1[0];
        n[1] 		= point_2[1] - point_1[1];

        double 	n_len 	= Math.sqrt(n[0]*n[0] + n[1]*n[1]);

        if(n_len<=0.001){
            System.err.println("Error: distance p1-p2 was too small tocalculate direction: "+n_len);
            return;
        }
        // normalize it
        n[0]        = n[0] / n_len;
        n[1]        = n[1] / n_len;

        // assign it to the class member

        //go through every point
//        for (int layer = 0; layer < values.length; layer++) {
            for (int coord = 0; coord < values.length; coord++) {

                int[] p = new int[2];
                p[0] = ArrayHandling.index2row_2d(coord, w);
                p[1] = ArrayHandling.index2col_2d(coord, w);

				/*
				 * check where within the tube it is
				 */

                // calculate (p1-p) dotprod (p2-p1)
                double projection_1 =
                        (point_1[0]-p[0])*(point_2[0]-point_1[0]) +
                                (point_1[1]-p[1])*(point_2[1]-point_1[1]);
//                +               (point_1[2]-p[2])*(point_2[2]-point_1[2]);

                // calculate p2-p
                double projection_2 =
                        (point_2[0]-p[0])*(point_1[0]-point_2[0]) +
                                (point_2[1]-p[1])*(point_1[1]-point_2[1]);
//                                +                                (point_2[2]-p[2])*(point_1[2]-point_2[2]);

                double distance = 0;

                if(projection_1<=0 && projection_2<=0){

                    //image_array[layer][ArrayHandling.sub2index_2d(p[0], p[1], image_width)] = (byte)255;
                    // inside the line distance from line with a=p1 and n=unitVec(p1,p2)
                    distance = distance_point_to_line_2d(point_1, n, p);

                }
                else if(projection_1>0 && projection_2<=0){

                    // outside p1 - measure the distance to p1
                    distance = dist(p, point_1);

                }
                else if(projection_1<=0 && projection_2>0){

                    // outside p2 - measure the distance to p2
                    distance = dist(p, point_2);

                }
                else{
                    System.err.println("this case is not possible when printing line between p1 and p2");
                    return;
                }

                if(distance<3*rstd) { // some reasonable distance from where we consider intensity is ~0

                    byte value_to_write = (byte)Math.round(fg * Math.exp(-Math.pow(distance,2)/(2*Math.pow(rstd, 2))));

                    values[ArrayHandling.sub2index_2d(p[0], p[1], w)] =
                            (byte) Math.max(
                                    (int)(value_to_write & 0xff),
                                    (int)(values[ArrayHandling.sub2index_2d(p[0], p[1], w)] & 0xff)
                            );

                }
                if (distance<3.5*rstd) {

                    mask[ArrayHandling.sub2index_2d(p[0], p[1], w)] = (byte) 255;

                }
            }
    }

    private double 	distance_point_to_line_2d(int[] line_base_point, double[] line_unit_direction, int[] point){

        // line is in vector form (point+unit_direction)

        double distance_from_line = 0;

        if		(
                line_base_point.length              != 2 	||
                        line_unit_direction.length  != 2    ||
                        point.length                != 2
                ){

            System.err.println(
                    "vectors are not of length 2"
            );
            return Double.NaN;

        }

//        if(vectorNorm(line_unit_direction)< 0.99 || vectorNorm(line_unit_direction)> 1.01){
//            System.err.println(
//                    "Transf:distance_point_to_line():Line's unit orientation vector must have unit length! \nCurrent length is "+
//                            vectorNorm(line_unit_direction));
//            return Double.NaN;
//        }

        double[] distance_2d = new double[2];

        // line_base_point - point
        distance_2d[0] = line_base_point[0] - point[0];
        distance_2d[1] = line_base_point[1] - point[1];
//        distance_2d[2] = line_base_point[2] - point[2]; 	// it contains difference a-p

        double proj =
                        distance_2d[0] * line_unit_direction[0] +
                        distance_2d[1] * line_unit_direction[1];// +
//                        distance_2d[2] * line_unit_direction[2]; 	// dotproduct((a-p), n)

        distance_2d[0] = distance_2d[0] - proj * line_unit_direction[0];
        distance_2d[1] = distance_2d[1] - proj * line_unit_direction[1];
//        distance_2d[2] = distance_2d[2] - proj * line_unit_direction[2]; // || (a-p) - ((a-p)*n) * n ||

        distance_from_line = vectorNorm(distance_2d);

        return distance_from_line;
    }

    private double 	vectorNorm(
            double[] vector
    )
    {
        double norm = 0;

        for (int i = 0; i < vector.length; i++) {
            norm += vector[i] * vector[i];
        }

        return Math.sqrt(norm);
    }

    private double  dist(
            int[]       p1,
            int[]       p2
    )
    {
        return Math.sqrt( Math.pow(p1[0]-p2[0], 2) + Math.pow(p1[1]-p2[1], 2) );
    }

	int wrap_360(
						int in
	)
	{
		while(in>=360) in-=360;
		while (in<0) in+=360;
		return in;
	}

}

import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/3/13
 * Time: 12:32 PM
 *
 *  reference:
 *  Mean shift, mode seeking, and clustering, Yizong Cheng
 *  doi: 10.1109/34.400568
 */

public class MeanShift2D {

    int 		image_width;
    int 		image_height;

    double[][]	S;					// finite set, data, sample
    double[][]	T;					// cluster centers

    double[] 	w; 					// weights (image intensities)

    int 		h_spatial;
    int 		h_spatial2; 		// squared values

    FloatProcessor inip;

    static double Tdist 	= 0.001;

    public MeanShift2D(FloatProcessor imgp, int h_spatial){ // add an option for ByteProcessor

        image_height 	= imgp.getHeight();
        image_width 	= imgp.getWidth();

        float[] image_values = (float[])imgp.getPixels();

//        w				= new double[image_height*image_width];
//        for (int i = 0; i < image_height*image_width; i++) {
//            w[i] 		= (double)(image_values[i] & 0xff);
//        }

        inip = new FloatProcessor(image_width, image_height, image_values);

        int count = 0;
        S				= new double[image_height*image_width][2];
        T				= new double[image_height*image_width][2];
        for (int i = 0; i < image_height; i++) {
            for (int j = 0; j < image_width; j++) {
                S[count][0] = (double)i;
                S[count][1] = (double)j;
                count++;
            }
        }

        this.h_spatial 		= h_spatial;
        h_spatial2 			= h_spatial*h_spatial;

    }

//    public void set(ByteProcessor imgp, int h_spatial) {
//
//        image_height 	= imgp.getHeight();
//        image_width 	= imgp.getWidth();
//
//        byte[] image_values = (byte[])imgp.getPixels();
//
//        w				= new double[image_height*image_width];
//        for (int i = 0; i < image_height*image_width; i++) {
//            w[i] 		= (double)(image_values[i] & 0xff);
//        }
//
//        inip = new FloatProcessor(image_width, image_height, w);
//
//        int count = 0;
//        S				= new double[image_height*image_width][2];
//        T				= new double[image_height*image_width][2];
//        for (int i = 0; i < image_height; i++) {
//            for (int j = 0; j < image_width; j++) {
//                S[count][0] = (double)i;
//                S[count][1] = (double)j;
//                count++;
//            }
//        }
//
//        this.h_spatial 		= h_spatial;
//        h_spatial2 			= h_spatial*h_spatial;
//
//    }

    private double[] 	runOneIterAtPos(double[] curr_pos){
        double[] 	new_pos		= new double[2];
        double 		sum 		= 0;
        for (int x = -h_spatial; x <= h_spatial; x++) {
            if ((  curr_pos[0] + x >=0) && ( curr_pos[0] + x <= image_height-1) ){
                for (int y = -h_spatial  ; y <= h_spatial; y++) {
                    if(( curr_pos[1] + y >=0) && (curr_pos[1] + y <=image_width-1)) {
                        if (x * x + y * y <= h_spatial2) {
                            //double value_read 		= w[x0 * image_width + y0];
                            double value_read = Interpolator.interpolateAt(curr_pos[1]+y, curr_pos[0]+x, inip);//interpolateIntensity(curr_pos[0]+x, curr_pos[1]+y);
                            // possible kernel
                            // / (1 + x * x) / (1 + y * y);
                            // / (1 + x * x) / (1 + y * y);
                            // / (1 + x * x) / (1 + y * y);
                            //switch(knType){
//                                case FLAT:
                                    new_pos[0] 	+= value_read * x;
                                    new_pos[1] 	+= value_read * y;
                                    sum 		+= value_read 	 ;

//                                    break;
//                                case TRUNCATED_GAUSSIAN_1:
//                                    new_pos[0] 	+= value_read * x * Math.exp(-(x*x+y*y));
//                                    new_pos[1] 	+= value_read * y * Math.exp(-(x*x+y*y));
//                                    sum 		+= value_read 	 ;
//                                    break;
//                                default:
//                                    new_pos[0] 	+= value_read * x;
//                                    new_pos[1] 	+= value_read * y;
//                                    sum 		+= value_read 	 ;
//                                    break;
                            //}
                        }
                    }
                }
            }
        }
        if(sum>0){
            new_pos[0] = new_pos[0]/sum + curr_pos[0];
            new_pos[1] = new_pos[1]/sum + curr_pos[1];
        }
        else {
            new_pos[0] = curr_pos[0];
            new_pos[1] = curr_pos[1];
        }
        return new_pos;
    }

    public double[][] 	run(int max_iter, double epsilon){

        T				= new double[image_height*image_width][2];
        for (int i = 0; i < T.length; i++) {
            T[i][0] = S[i][0];	T[i][1] = S[i][1];
        }

        for (int i = 0; i < image_height*image_width; i++) {
            int iter = 0;
            double d = Double.MAX_VALUE;
            do{
                double[] new_pos = runOneIterAtPos(T[i]);
                d = d2(new_pos, T[i]);
                T[i][0] = new_pos[0];
                T[i][1] = new_pos[1];
                iter++;
            }
            while(iter < max_iter && d > epsilon*epsilon);
        }

        return T;
    }

    public double[][] extractConvPoints(double range, int M){

        // 'range' describes neighborhood range size,
        // 'M' number of samples within the cluster

        // T will be checked and extracted clusters stored in it, starting from the
        // element with first index

        boolean[] 	clustered 		= new boolean[T.length]; // initialized with false
        int[] 		cluster_size 	= new int[T.length];
        int nr_clusters = 0;
        int nr_clusters_nigher_than_M = 0;

        for (int i = 0; i < T.length; i++) {

            if(!clustered[i]){

                clustered[i] = true;
                cluster_size[nr_clusters]++;
                T[nr_clusters][0] = T[i][0];
                T[nr_clusters][1] = T[i][1];

                for (int j = i+1; j < T.length; j++) {

                    if(!clustered[j] && dist(T[i], T[j])<=range){

                        clustered[j] = true;
                        cluster_size[nr_clusters]++;

                    }

                }

                if(cluster_size[nr_clusters]>M){
                    nr_clusters_nigher_than_M++;
                }

                nr_clusters++;

            }

        }

        // take them out
        double[][] conv_pts = null;
        int cnt = 0;
        conv_pts = new double[nr_clusters_nigher_than_M][2];

        for (int i = 0; i < nr_clusters; i++) {
            if(cluster_size[i]>M){

                conv_pts[cnt][0] = T[i][0]; conv_pts[cnt][1] = T[i][1];
                cnt++;
            }
        }

        return conv_pts; // nr_pointsx2

    }

    public int[][] extractClust(int M){

        // create list
        ArrayList<int[]> T_;
        T_			= new ArrayList<int[]>();

        // initialize the list with T values
        for (int i = 0; i < T.length; i++) {
            T_.add(new int[]{ (int)Math.round(T[i][0]) , (int)Math.round(T[i][1]) });
        }


        // final list - after cropping those with small number of elements
        ArrayList<int[]>	out_list;
        out_list	= new ArrayList<int[]>();

        while(!T_.isEmpty()){

            // take the first element
            int[] 	check_element 		= T_.get(0);
            int		check_element_count	= 1;
            T_.remove(0);

            int i = 0;
            while(i<T_.size()){
                int[] take_element = T_.get(i);

                if(equal(take_element, check_element)){
                    T_.remove(i);
                    check_element_count++;
                }
                else{
                    i++;
                }

            }

            if(check_element_count>M){
                out_list.add(new int[]{ check_element[0], check_element[1], check_element_count });
            }

        }
        // out_list contains the selection of the coordinates, each corersponding to a direction

        int[][] out = new int[out_list.size()][3];
        for (int i = 0; i < out_list.size(); i++) {
            int[] take_dir_spherical = out_list.get(i);

            out[i][0] = take_dir_spherical[0];
            out[i][1] = take_dir_spherical[1];
            out[i][2] = take_dir_spherical[2];

        }

        return out;

    }

    public void printRoundedT(){
        for (int i = 0; i < T.length; i++) {
            System.out.format("%d, %d \n", (int)Math.round(T[i][0]), (int)Math.round(T[i][1]));
        }
    }

    private double d2(double[] a, double[] b){
        return Math.pow(a[0]-b[0], 2)+Math.pow(a[1]-b[1], 2);
    }

    private boolean equal(int[] a, int[] b){
        return a[0]==b[0] && a[1]==b[1];
    }

    private double dist(double[] a, double[] b){
        return Math.sqrt(Math.pow((a[0]-b[0]), 2)+Math.pow((a[1]-b[1]), 2));
    }

}

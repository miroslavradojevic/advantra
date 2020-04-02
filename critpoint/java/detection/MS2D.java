package detection;

import aux.Interpolator;
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

public class MS2D extends Thread {

    private int begN, endN;

    public static int 		image_width;
    public static int 		image_height;
    public static int 		toProcess;


    public static double[][]	S;					// x,y
    public static double[][]	T;					// after convergence

    public static int 		h_spatial = 4;
    public static int 		h_spatial2 = 4*4; 		    // squared

    public static FloatProcessor    inip;
    public static ByteProcessor     maskip;

    static double Tdist 	= 0.001;
    public static int       max_iter    = 150;
    public static double    epsilon     = 0.0001;

    public MS2D (int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(ImageProcessor inip1, ByteProcessor maskip1, int h1, int max_iter1, double epsilon1)
    {

        /*
        set inip
         */
        inip = new FloatProcessor(inip1.getWidth(), inip1.getHeight());
        for (int i=0; i<inip1.getWidth()*inip1.getHeight(); i++) {
            inip.setf(i, inip1.getf(i));
        }

        /*
        set maskip
         */
        int cnt = 0;
        maskip = new ByteProcessor(maskip1.getWidth(), maskip1.getHeight());
        byte[] arr = (byte[]) maskip1.getPixels();
        for (int i=0; i<arr.length; i++) {
            if (arr[i]==(byte)255) {
                maskip.set(i, 255);
                cnt++;
            }
            else {
                maskip.set(i, 0);
            }
        }

        toProcess = cnt;

        image_height 	= inip.getHeight();
        image_width 	= inip.getWidth();

        /*
        set locations
         */
        S = new double[cnt][2];
        T = new double[cnt][2];

        cnt = 0;
        for (int i=0; i<arr.length; i++) {
            if (arr[i]==(byte)255) {
                S[cnt][1] = i%maskip.getWidth();// column
                S[cnt][0] = i/maskip.getWidth();// row
                cnt++;
            }
        }

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

        h_spatial 		    = h1;
        h_spatial2 			= h_spatial*h_spatial;

        max_iter = max_iter1;
        epsilon = epsilon1;

    }

    private double[] 	runOneIterAtPos(double[] curr_pos){
        double[] 	new_pos		= new double[2];
        double 		sum 		= 0;
        for (int x = -h_spatial; x <= h_spatial; x++) {
            if ((  curr_pos[0] + x >=0) && ( curr_pos[0] + x <= image_height-1) ){
                for (int y = -h_spatial  ; y <= h_spatial; y++) {
                    if(( curr_pos[1] + y >=0) && (curr_pos[1] + y <=image_width-1)) {
                        if (x * x + y * y <= h_spatial2) {
                            //double value_read 		= w[x0 * image_width + y0];
                            double value_read = Interpolator.interpolateAt(curr_pos[1] + y, curr_pos[0] + x, inip);//interpolateIntensity(curr_pos[0]+x, curr_pos[1]+y);
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

    public void 	    run(){

        // will work out range of values from T

        //T				= new double[S.length][2];
        for (int i = begN; i < endN; i++) {
            T[i][0] = S[i][0];	T[i][1] = S[i][1];
        }

        for (int i = begN; i < endN; i++) {
            int iter = 0;
            double d;// = Double.MAX_VALUE;
            do{
                double[] new_pos = runOneIterAtPos(T[i]);
                d = d2(new_pos, T[i]);
                T[i][0] = new_pos[0];
                T[i][1] = new_pos[1];
                iter++;
            }
            while(iter < max_iter && d > epsilon*epsilon);
        }

        //return T;
    }

    public static ArrayList<ArrayList<double[]>> extractClusters(double minDist){

        // 'range' describes neighborhood range size,
        // 'M' number of samples within the cluster

        //boolean[] 	clustered 		= new boolean[T.length]; // initialized with false
//		int[]		clusterIdx		= new int[T.length];
//		for (int i1=0; i1<T.length; i1++) clusterIdx[i1] = -1;

		ArrayList<ArrayList<double[]>> clusters = new ArrayList<ArrayList<double[]>>(T.length); // allocate T.length

        int nr_clusters = 0;

        for (int i = 0; i < T.length; i++) {

//            if(clusterIdx[i]==-1){

				// it was not clustered: check if it belongs to some formed cluster
				// or make a new one if not

				boolean assignedToClstr = false;
				for (int a1=0; a1<clusters.size(); a1++) {

					for (int b1=0; b1<clusters.get(a1).size(); b1++) {

						if (d2(T[i], clusters.get(a1).get(b1))<=minDist*minDist) {

							clusters.get(a1).add(T[i]);
//							clusterIdx[i] = a1;
							assignedToClstr = true;
							break;

						}

					}

					if (assignedToClstr) break;

				}

				if (!assignedToClstr) {

					//initialize one cluster
					ArrayList<double[]> newClstr = new ArrayList<double[]>();
					newClstr.add(T[i]);
//					clusterIdx[i] = clusters.size();
					clusters.add(newClstr);

				}

//				clusterIdx[i] = nr_clusters;
                //cluster_size[nr_clusters]++;
//                T[nr_clusters][0] = T[i][0];
//                T[nr_clusters][1] = T[i][1];

//                for (int j = i+1; j < T.length; j++) {
//
//                    if(clusterIdx[j]==-1 && d2(T[i], T[j])<=minDist*minDist){
//
//						clusterIdx[j] = nr_clusters;
////                        cluster_size[nr_clusters]++;
//
//                    }
//
//                }

//                if(cluster_size[nr_clusters]>M){
//                    nr_clusters_nigher_than_M++;
//                }

//                nr_clusters++;

//            }

        }

//        // take them out
//        double[][] conv_pts = null;
//        int cnt = 0;
//        conv_pts = new double[nr_clusters_nigher_than_M][2];
//
//        for (int i = 0; i < nr_clusters; i++) {
//            if(cluster_size[i]>M){
//
//                conv_pts[cnt][0] = T[i][0]; conv_pts[cnt][1] = T[i][1];
//                cnt++;
//            }
//        }

        return clusters;

    }

    public static ArrayList<double[]> extractClusterCentroids(ArrayList<ArrayList<double[]>> clusters, int M){

        ArrayList<double[]> selection = new ArrayList<double[]>();

        // row, col, nr_points, peak value
        for (int clusterIdx=0; clusterIdx<clusters.size(); clusterIdx++) {

            if (clusters.get(clusterIdx).size()>M) {

                // find centroid
                double xc = 0;
                double yc = 0;
                double nr = clusters.get(clusterIdx).size();


                for (int i1=0; i1<clusters.get(clusterIdx).size(); i1++) {
                    xc += clusters.get(clusterIdx).get(i1)[1];
                    yc += clusters.get(clusterIdx).get(i1)[0];
                }

                xc/=nr;
                yc/=nr;

                double pk = inip.getf((int) Math.round(xc), (int) Math.round(yc));

                selection.add(new double[]{xc, yc, nr, pk});

            }

        }

        return selection;

    }

    private static double d2(double[] a, double[] b){
        return (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]);
    }

    private boolean equal(int[] a, int[] b){
        return a[0]==b[0] && a[1]==b[1];
    }

}

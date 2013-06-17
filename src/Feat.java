import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/13/13
 * Time: 10:30 AM
 */
public class Feat {

	public int			r;
	public int			d;
	public int          diam;

	public int          rInner;
	public int          rLower;

	int					resolDeg;
	int					h;

	ArrayList<ArrayList<int[]>> offsets = new ArrayList<ArrayList<int[]>>();

	float TwoPI = (float) (Math.PI*2);

	public Feat(int diam, double scale){

		this.diam = diam;
		r = (int) (diam*scale);
		r = (r<6)? 6 : r; // lower limit
		d = 2*r+1;

		rInner = diam/2;
		rLower = r/2;

		resolDeg = 10;

		double  angStep     = 2 * Math.asin((float)diam/(2*r));
		int     angStepDeg  = (int) Math.round( ((angStep/(2*Math.PI))*360) );
		angStepDeg = (angStepDeg<resolDeg)?resolDeg:angStepDeg;
		angStepDeg  = (angStepDeg/resolDeg)*resolDeg; // will be used for mean-shift, defines the neighbourhood to search
		h = angStepDeg/resolDeg;
		h = (h<1)?1:h;

		// form the kernels
		int xc = d/2;
		int yc = d/2;
		int[] cent = new int[]{0, 0};
		float[] n = new float[2];
		int[] p = new int[2];

		for (int aDeg = 0; aDeg<360; aDeg+=resolDeg) {

			float aRad = ((float)aDeg/360)*TwoPI;

			ArrayList<int[]> offsetsAngle = new ArrayList<int[]>();

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					float d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (Math.floor(d2) <= r*r && Math.floor(d2) >= rLower*rLower) {

						p[0] = x-xc;
						p[1] = -(y-yc);

						float ang = aRad;  		// read the angle

						n[0] = (float) Math.cos(ang);
						n[1] = (float) Math.sin(ang);

						float dst = point2dir(cent, n, p);

						if (Math.round(dst)<=diam/2) { // belongs to pIdx ON peak and not filled

							offsetsAngle.add(new int[]{p[0], p[1]});

						}

					}
				}
			}

			offsets.add(offsetsAngle);
//			System.out.println(offsets.size()+" -> adding "+offsetsAngle.size()+" elements");

		}

//		System.out.println("total "+offsets.size());

//		String[] s_all = new String[]{"'ro'", "'go'", "'bo'", "'yo'"};
//		for (int c = 0; c<offsets.size(); c+=4) {
//			System.out.println("a = [...");
//			for (int c1 = 0; c1<offsets.get(c).size(); c1++)
//				System.out.println(""+offsets.get(c).get(c1)[0]+" , "+offsets.get(c).get(c1)[1]+";...");
//			System.out.println("];");
//			String s =  s_all[new Random().nextInt(3)];
//			System.out.println("plot(a(:,1), a(:,2), "+s+"); axis equal; grid on; hold on;");
//		}

	}

    public ImageStack showOffsets(){

        ImageStack isOut = new ImageStack(d, d);

        for (int a = 0; a<offsets.size(); a++) {

            ImageProcessor ip = new ByteProcessor(d, d);

            for (int b = 0; b<offsets.get(a).size(); b++) {

                int offX = offsets.get(a).get(b)[0];
                int offY = offsets.get(a).get(b)[1];

                ip.setf(d/2+offX, d/2+offY, +255);

            }

            isOut.addSlice(ip);

        }

        return isOut;

    }

//	public Vector<Double> getPeakAnglesDeg360(int atX, int atY, FloatProcessor inip) {
//
//		float[] sumsPerOrt  = new float[offsets.size()];
//
//		for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {
//
//			for (int locIdx = 0; locIdx<offsets.get(dirIdx).size(); locIdx++) {
//
//				sumsPerOrt[dirIdx] += inip.getf(
//					atX+offsets.get(dirIdx).get(locIdx)[0], atY+offsets.get(dirIdx).get(locIdx)[1]
//				);
//
//			}
//
//		}
//
//		// find peaks - initialize start points for ms iterations
//		int Npoints = 200;
//		double eps = 0.0001;
//		int Niter = 200;
//		int M = 1;
//
//		double[] start = new double[Npoints];
//		double[] finish = new double[start.length];
//
//		for (int i=0; i<Npoints; i++)
//			start[i] = ((double)i/Npoints)*sumsPerOrt.length;
//
//		finish = runMS(start, sumsPerOrt, Niter, eps, h);
//
//		Vector<Double> cls = extractClusters(finish, 0.5, M);
//
//
//
//	}

    public ImageStack plotSums(int atX, int atY, FloatProcessor inip) {

        // calculate sums at each orientation
        float[] sumsPerOrt  = new float[offsets.size()];
        float[] angles      = new float[offsets.size()];

        for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {

            angles[dirIdx] = dirIdx;

            for (int locIdx = 0; locIdx<offsets.get(dirIdx).size(); locIdx++) {

                sumsPerOrt[dirIdx] += inip.getf(
                        atX+offsets.get(dirIdx).get(locIdx)[0], atY+offsets.get(dirIdx).get(locIdx)[1]
                );

            }

        }

		// plot angles vs. sums
		Plot p;
		ImageStack isOut;
		p = new Plot("", "", "", angles, sumsPerOrt);
		isOut = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
		isOut.addSlice(p.getProcessor());

		// find peaks - initialize start points for ms iterations
		int Npoints = 200;
		double[] start = new double[Npoints];
		for (int i=0; i<Npoints; i++)
			start[i] = ((double)i/Npoints)*sumsPerOrt.length;

		double[] plotStart = new double[Npoints];
		for (int i=0; i<Npoints; i++)
			plotStart[i]=interp1Darray(start[i], sumsPerOrt);//  sumsPerOrt[(int)Math.floor(start[i])];
		p = new Plot("", Npoints+" points start", "", angles, sumsPerOrt);
		p.addPoints(start, plotStart, Plot.CIRCLE);
		isOut.addSlice(p.getProcessor());

		double[] finish = new double[start.length];
		double[] plotFinish;//= new double[start.length];

		double eps = 0.0001;
		int maxIter = 200;
		int M = 5;
		double minD = 0.1;
		int h = 4;

		for (maxIter = 1; maxIter<=15; maxIter++) {

			//finish = runMS(start, sumsPerOrt, Niter, eps, h);
			finish = runBlurring(start, sumsPerOrt, maxIter, eps, h, M, minD);

			plotFinish = new double[finish.length];
			for (int i=0; i<finish.length; i++)
				plotFinish[i]=interp1Darray(finish[i], sumsPerOrt);

			Vector<Double> cls = extractClusters(finish, minD, M);
			p = new Plot("", maxIter+" iters,h="+h+","+cls.size()+"clusters", "", angles, sumsPerOrt);
			p.addPoints(finish, plotFinish, Plot.CROSS);
			if (cls.size()>0) {
				double[] plotClsY = new double[cls.size()];
				double[] plotClsX = new double[cls.size()];
				for (int i=0; i<plotClsY.length; i++) {
					plotClsY[i] = interp1Darray(cls.get(i), sumsPerOrt);
					plotClsX[i] = cls.get(i);
				}
				p.addPoints(plotClsX, plotClsY, Plot.BOX);
			}

			isOut.addSlice(p.getProcessor());

		}

        return isOut;

    }

	private float 	point2dir(
		int[] 	b,    	// direciton base point
		float[] n, 		// direction vector
		int[] 	p     	// point considered
	)
	{
		// line is in vector from b in n direction

		float d = 0;

		float[] p_b = new float[2];

		// p - b
		p_b[0] = p[0] - b[0];
		p_b[1] = p[1] - b[1];

		float proj = p_b[0] * n[0] + p_b[1] * n[1];

		if(proj<0){
			// "behind" the orientation
			return (float) Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
		}

		// || (p-b) - dot(p-b,n) * n ||
		p_b[0] = p_b[0] - proj * n[0];
		p_b[1] = p_b[1] - proj * n[1];

//		distance_from_line = vectorNorm(distance_2d);

		return (float) Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
	}


    public double[] runMS(double[] start, float[] inputProfile, int max_iter, double epsilon, int h){

        double[] T				= new double[start.length];

        for (int i = 0; i < T.length; i++) {
            T[i] = start[i];//	T[i][1] = S[i][1];
        }

        for (int i = 0; i < T.length; i++) {

			int iter = 0;
            double d = Double.MAX_VALUE;

            do{

                double new_pos = runOne(T[i], h, inputProfile);
                d = Math.abs(new_pos - T[i]);
                T[i] = new_pos;
                iter++;
            }
            while(iter < max_iter && d > epsilon);
//			if (d<=epsilon) {
//				System.out.println(iter+" iterations, "+max_iter+",eps="+epsilon);
//			}
        }

        return T;

    }

	public 	double[] 	runBlurring(
		double[] 	start,
		float[] 	inputProfile,
		int 		maxIter,
		double 		epsilon,
		int 		h,
		int 		M,
		double		intraDist
	){

		double[]	S = new double[start.length];
		double[] 	T = new double[start.length];
		double 		dmax;
		int 		iter = 0;

		for (int i=0; i<start.length; i++) {
			S[i] = T[i] = start[i];

		}

		do{

			dmax = Double.MIN_VALUE;

			for (int i = 0; i < S.length; i++) {      //

				double sum = 0;
				double weightedShift = 0;

				for (int j = 0; j < S.length; j++) {

					double d1 = S[j] - S[i];
					double d2 = (S[j]+inputProfile.length) - S[i];
					double d3 = (S[j]-inputProfile.length) - S[i];

					double d12 = d1*d1;//Math.pow(S[i]-S[j], 2);
					double d22 = d2*d2;//Math.pow(S[i]-(S[j]+inputProfile.length), 2);
					double d32 = d3*d3;//Math.pow(S[i]-(S[j]-inputProfile.length), 2);

					if(d12<=h*h){
						double w = interp1Darray(S[j], inputProfile);
						weightedShift 	+= w * d1;
						sum 	+= w;
					}
					else if (d22<=h*h) {
						double w = interp1Darray(S[j], inputProfile);
						weightedShift 	+= w * d2;
						sum 	+= w;
					}
					else if (d32<=h*h) {
						double w = interp1Darray(S[j], inputProfile);
						weightedShift 	+= w * d3;
						sum 	+= w;
					}


				}
				if(sum>0){
					T[i] = T[i] + weightedShift/sum;
					// wrapping
					T[i] = (T[i]<-0.5)?(T[i]+inputProfile.length):((T[i]>=inputProfile.length-0.5)?(T[i]-inputProfile.length):T[i]);
				}
				else{
					T[i] = S[i];
				}

				if(Math.abs(S[i]-T[i])>dmax){
					dmax = Math.abs(S[i]-T[i]);
				}

			}

			// check the amount of clusters
			//Vector<Double> cls = extractClusters(T, intraDist, M);

//			if (cls.size()<=3) {
//				IJ.log("found <= 3");
//				break;
//			}

			for (int i = 0; i < S.length; i++) {
				S[i] = T[i];
			}

			iter++;

		}
		while(iter<=maxIter && dmax>epsilon);

		return T;

	}

	public Vector<Double> extractClusters(double[] msConvergence, double range, int M){

		// 'range' describes neighborhood range size,
		// 'M' number of samples within the cluster

		boolean[] 	seed_clustered 		= new boolean[msConvergence.length];
		int[] 		seed_cluster_idx 	= new int[msConvergence.length];
		Vector<Integer> cluster_sizes	= new Vector<Integer>();
		int nr_clusters 				= 0;

		for (int i = 0; i < msConvergence.length; i++) {
			for (int j = 0; j < msConvergence.length; j++) {

				if(!seed_clustered[j] && j!=i){

					if(Math.abs(msConvergence[i]-msConvergence[j])<range){

						if(seed_clustered[i]){
							// i was clustered
							seed_clustered[j] = true;
							seed_cluster_idx[j] = seed_cluster_idx[i];

							int tmp = cluster_sizes.get(seed_cluster_idx[i]);
							tmp++;
							cluster_sizes.setElementAt(tmp, seed_cluster_idx[i]);

						}
						else{
							// i was not clustered
							seed_clustered[i] = seed_clustered[j] = true;
							seed_cluster_idx[i] = seed_cluster_idx[j] = nr_clusters;
							cluster_sizes.add(2);
							nr_clusters++;
						}
					}
				}
			}

			if(!seed_clustered[i]){
				seed_clustered[i] = true;
				seed_cluster_idx[i] = nr_clusters;
				cluster_sizes.add(1);//set(nr_clusters, 1);
				nr_clusters++;
			}

		}

		//System.out.println("\n---\ncluster sizes for "+msConvergence.length+" conv. pts\n---\n");
		Vector<Double> clust = new Vector<Double>();
		for (int k = 0; k < cluster_sizes.size(); k++) {
			if (cluster_sizes.get(k)>M) {

				double sum = 0;
				int cnt = 0;

				for (int l=0; l<seed_cluster_idx.length; l++) {

					if (seed_cluster_idx[l]==k) {

						sum+=msConvergence[l];
						cnt++;

					}

				}

				if (cnt>0) {
					clust.add(sum/cnt);
//					System.out.println("\n\ncluster[" + clust.size()+"]="+clust.get(clust.size()-1));
				}


			}
		}

		return clust;


	}

    private double interp1Darray(double x, float[] inarray) {

        int wth = inarray.length;

        // x is a real position of a pixel in an image
        // x = [-0.5, 	W-0.5)
        // width is array width - indexes will be wrapped - first corresponds to the last one

        // find four surrounding points
        int p11 = (int) Math.floor(x);
        int p12 = (int) Math.ceil(x);
        // bilinear coeffs a,b
        double a = ((p12 -p11)>0)?(x-p11)/(p12-p11) : 0.5;
        // wrap cols of surrounding pts
        p11 = (p11<-0.5)? p11+wth : ((p11>=wth-0.5)?(p11-wth) : p11);
        p12 = (p12<-0.5)? p12+wth : ((p12>=wth-0.5)?(p12-wth) : p12);

        double I11 = inarray[p11];
        double I12 = inarray[p12];

        return (1-a)*I11 + a*I12;
    }

	private double 	runOne(double curr_pos, int h, float[] inputProfile){
		double 	    new_pos     = 0;
		double 		sum 		= 0;

		for (int x = -h; x <= h; x++) {
			if (true) {// ((  curr_pos + x >=0) && ( curr_pos + x <= inputProfile.length-1) ){
				if (true) {
					double value_read = interp1Darray(curr_pos+x, inputProfile);
					new_pos 	+= value_read * x;
					sum 		+= value_read 	 ;
				}
			}
		}
		if(sum>0){
			new_pos = new_pos/sum + curr_pos;
//            new_pos[1] = new_pos[1]/sum + curr_pos[1];
			// wrap it again, this time as a position
			new_pos = (new_pos<-0.5)?(new_pos+inputProfile.length):((new_pos>=inputProfile.length-0.5)?(new_pos-inputProfile.length):new_pos);
		}
		else {
			new_pos = curr_pos;
		}

		return new_pos;
	}

//	public float scoreAt(){
//
//	}
//
//	public ImageProcessor template(
//										   float[] 	directionAngles,
//										   int 		Rpix,
//										   int 		Rlower,
//										   int 		Rinner,
//										   float 	Tpix,
//										   // outputs
//										   int[][][]     regionMap,
//										   //int[][][]     region_size,
//										   float[][][]   kernel
//	)
//	{
//
//		int[][]     idxMapLocal = new int[nRot][];
//		int[][]     regSzsLocal = new int[nRot][];
//		float[][]   kernelLocal = new float[nRot][];
//
//		float[][] peaksRad = new float[nRot][nPeaks];
//
//		for (int cnt_rots = 0; cnt_rots < nRot; cnt_rots++) {
//
//			float start_pos = cnt_rots*rotStp;
//
//			for (int cnt_pks = 0; cnt_pks < nPeaks; cnt_pks++) {
//
//				if(cnt_pks==0)
//					peaksRad[cnt_rots][cnt_pks] 	= start_pos;
//				else
//					peaksRad[cnt_rots][cnt_pks] 	=
//							peaksRad[cnt_rots][cnt_pks-1] +
//									angles[cnt_pks-1];
//
//			}
//		}
//
//		System.out.println("forming with angles...");
//		for (int cnt_rots = 0; cnt_rots < nRot; cnt_rots++) {
//			for (int cnt_pks = 0; cnt_pks < nPeaks; cnt_pks++) {
//				System.out.print(peaksRad[cnt_rots][cnt_pks]+", ");
//			}
//			System.out.println();
//		}
//		System.out.println("--\n");
//
//		// form the kernels
//		int xc = d/2;
//		int yc = d/2;
//		int[] cent = new int[]{0, 0};
//		float[] n = new float[2];
//		int[] p = new int[2];
//		float[] ap = new float[nPeaks];
//
//		for (int rIdx = 0; rIdx < nRot; rIdx++) {
//
//			idxMapLocal[rIdx] = new int[d*d];
//			regSzsLocal[rIdx] = new int[2*nPeaks+1];
//			kernelLocal[rIdx] = new float[d*d];
//
//
//
//			for (int pIdx = 0; pIdx < nPeaks; pIdx++)
//				ap[pIdx] = ( peaksRad[rIdx][pIdx] );
//
//			for (int i = 0; i<3; i++) {System.out.print("ap in:"+ap[i]);}
//
//			Arrays.sort(ap);
//
//			for (int i = 0; i<3; i++) {System.out.print("ap srt:"+ap[i]);}
//
//			for (int pIdx = 0; pIdx < nPeaks; pIdx++)
//				ap[pIdx] = wrap_0_2PI( peaksRad[rIdx][pIdx] );
//
//			// wrapped and sorted angles for this rotation
//
//			for (int i = 0; i<3; i++) {System.out.print("ap wrapped:"+ap[i]);}
//
//			int sumON = 0;
//			int sumOFF = 0;// for kernel normalization
//
//			for (int x = 0; x < d; x++) {
//				for (int y = 0; y < d; y++) {
//
//					int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);
//
//					if (d2 <= Rpix*Rpix && d2 >= Rlower*Rlower) {
//
//						p[0] = x-xc;
//						p[1] = -(y-yc);
//
//						boolean isON = false;
//
//						for (int pIdx = 0; pIdx < nPeaks; pIdx++) {      // check first peak first always
//
//							float ang = peaksRad[rIdx][pIdx];  		// read the angle
//
//							n[0] = (float) Math.cos(ang);
//							n[1] = (float) Math.sin(ang);
//
//							float dst = point2dir(cent, n, p);
//
//							if (dst<=Tpix/2 && idxMapLocal[rIdx][x+d*y]==0) { // belongs to pIdx ON peak and not filled
//
//								idxMapLocal[rIdx][x+d*y] = pIdx+1;
//								regSzsLocal[rIdx][pIdx+1]++;
//								kernelLocal[rIdx][x+d*y] = +1;
//								sumON++;
//								isON = true;
//								//break; check all three peaks independently
//
//							}
//
//						}
//
//						if(!isON) {
//
//							float a = wrap_0_2PI((float) Math.atan2(p[1], p[0]));
//
//							boolean added = false;
//							for (int pIdx = 0; pIdx < nPeaks-1; pIdx++) {
//
//								if (wrap_PI(a-ap[pIdx])>0 && wrap_PI(a-ap[pIdx+1])<=0 && idxMapLocal[rIdx][x+d*y]==0) {
//
//									idxMapLocal[rIdx][x+d*y] = nPeaks+1+pIdx;
//									regSzsLocal[rIdx][nPeaks+1+pIdx]++;
//									kernelLocal[rIdx][x+d*y] = -1;
//									sumOFF++;
//									added = true;
//									//break;
//								}
//
//							}
//							if (!added) {
//								idxMapLocal[rIdx][x+d*y] = 2*nPeaks;
//								regSzsLocal[rIdx][2*nPeaks]++;
//								kernelLocal[rIdx][x+d*y] = -1;
//								sumOFF++;
//							}
//
//						}
//
//					}
//					else if ( d2 <= Rinner*Rinner ) {
//						idxMapLocal[rIdx][x+d*y] = 0;
//						regSzsLocal[rIdx][0]++;
//						kernelLocal[rIdx][x+d*y] = +1;
//						sumON++;
//					}
//					else {
//						idxMapLocal[rIdx][x+d*y] = -1; // nothing
//						kernelLocal[rIdx][x+d*y] = 0;
//					}
//				}
//			}
//
//			System.out.println("\nrot. "+rIdx);
//			for (int i = 0; i<2*nPeaks+1; i++) {
//				System.out.print("["+i+" -> "+regSzsLocal[rIdx][i]+"],");
//			}
//			System.out.println();
//
//			// if it happens that it's too small region in one of the rotations, don't add the whole rotation
//			for (int k = 0; k < 2*nPeaks+1; k++) {
//
//				if (regSzsLocal[rIdx][k] <= 3) { // allow second sum to be 0 elements (overlaping with first)  && k!=2 && k!=4
//
//					System.out.println("some reg was small "+k+"th sum   : "+regSzsLocal[rIdx][k]+" elements");
//					return false;
//				}
//
//			}
//
//			// normalize kernel for this rotation
//			for (int k = 0; k<d*d; k++) {
//				if (kernelLocal[rIdx][k]>0) {
//					kernelLocal[rIdx][k] /= sumON;
//				}
//				else if (kernelLocal[rIdx][k]<0) {
//					kernelLocal[rIdx][k] /= sumOFF;
//				}
//			}
//
//		}
//
//		index_map[0] = idxMapLocal;
//		region_size[0] = regSzsLocal;
//		kernels_rot[0] = kernelLocal;
//
//		return true;
//
//	}

}

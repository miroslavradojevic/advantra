import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Arrays;
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

	ArrayList<ArrayList<int[]>> offsets;
	ArrayList<ArrayList<int[]>> template;

	float[] sumsPerOrt;

	float TwoPI = (float) (Math.PI*2);

	// will be used when forming kernels
	int 		xc;
	int 		yc;
	float[] 	n;
	int[] 		p;
	double[] 	ap;

	public Feat(
					   int diam,
					   double scale
	)
	{

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

		// form the list of offsets
		offsets = new ArrayList<ArrayList<int[]>>();

		xc = d/2;
		yc = d/2;
		n = new float[2];
		p = new int[2];
		ap 	= new double[3];

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

						float dst = point2dir(n, p);

						if (Math.round(dst)<=diam/2) { // belongs to pIdx ON peak and not filled

							offsetsAngle.add(new int[]{p[0], p[1]});

						}

					}
				}
			}

			offsets.add(offsetsAngle);

		}

		template = new ArrayList<ArrayList<int[]>>(7); // nrPeaks+1

		sumsPerOrt  = new float[offsets.size()];

	}

    public ImageStack plotOffsets()
	{
		// offsets have to be defined
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

	public  double[] getAngles(
        int             atX,
        int             atY,
        FloatProcessor  inip,
        int             Npoints,
        int             Niter,
        double          eps,
        int             M,
    	double          minD
    )
    {

		// check if values can be taken for the profile
		int margin = d/2+1;
		if ( atX<=margin || atY>=inip.getWidth()-margin ) {
			IJ.log("it's out of margin");
			return null; // out
		}
		if ( atY<=margin || atY>=inip.getHeight()-margin ) {
			IJ.log("it's out of margin");
			return null; // out
		}

		for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {
			for (int locIdx = 0; locIdx<offsets.get(dirIdx).size(); locIdx++) {
				sumsPerOrt[dirIdx] += inip.getf(
					atX+offsets.get(dirIdx).get(locIdx)[0], atY+offsets.get(dirIdx).get(locIdx)[1]
				);
			}
		}

		double[] start = new double[Npoints];
		double[] finish;

		for (int i=0; i<Npoints; i++)
			start[i] = ((double)i/Npoints)*sumsPerOrt.length;

		finish = runMS(start, sumsPerOrt, Niter, eps, h);

        Vector<double[]> cls = extractClusters(finish, minD, M);

        double[] anglesDirections;

        if (cls.size()>=3) {

            // extract 3 angles with most convergence points
            boolean[] checked = new boolean[cls.size()]; // all to false
            anglesDirections = new double[3];

            // find top 3
            for (int k = 0; k<3; k++) {

                // reset max search
                double  currMax = Double.MIN_VALUE;
                int     currMaxIdx = -1;

                for (int i=0; i<cls.size(); i++) {

                    // find max in this round
                    if (!checked[i]) {
                        if (cls.get(i)[1]>currMax) {

                            currMax = cls.get(i)[1];
                            currMaxIdx = i;

                        }
                    }
                }

                checked[currMaxIdx] = true;
                anglesDirections[k] = wrap_0_2PI( cls.get(currMaxIdx)[0]* (((double)resolDeg/360)*Math.PI*2) );

            }

        }
        else {
			IJ.log("clustered out "+cls.size()+" elements...");
            return null;
        }

        return anglesDirections;

	}

    public ImageStack plotProfilesWithMSDetection(
		int 			atX,
		int 			atY,
		FloatProcessor 	inip,
		int				Npoints,
		int 			maxIterBlurring,
		int 			maxIterNonBlurring,
		double			eps,
		double			minD,
		int 			M
	)
	{

        // calculate sums at each orientation
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
		p = new Plot("", "orientation (degrees)", "integrated response", angles, sumsPerOrt);
		isOut = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
		isOut.addSlice(p.getProcessor());

		// find peaks - initialize start points for ms iterations
		double[] start = new double[Npoints];
		for (int i=0; i<Npoints; i++)
			start[i] = ((double)i/Npoints)*sumsPerOrt.length;

		double[] plotStart = new double[Npoints];
		for (int i=0; i<Npoints; i++)
			plotStart[i]=interp1Darray(start[i], sumsPerOrt);

		p = new Plot("", Npoints+" points start", "", angles, sumsPerOrt);
		p.addPoints(start, plotStart, Plot.CIRCLE);
		isOut.addSlice(p.getProcessor());

		double[] finish;
		double[] plotFinish;

		finish = runBlurring(start, sumsPerOrt, maxIterBlurring, eps, h);
		plotFinish = new double[finish.length];
		for (int i=0; i<finish.length; i++)
				plotFinish[i]=interp1Darray(finish[i], sumsPerOrt);

		p = new Plot("", "blurring,"+maxIterBlurring+" iters,h="+h, "", angles, sumsPerOrt);
		p.addPoints(finish, plotFinish, Plot.BOX);
        isOut.addSlice(p.getProcessor());

        finish = runMS(start, sumsPerOrt, maxIterNonBlurring, eps, h);
        plotFinish = new double[finish.length];
        for (int i=0; i<finish.length; i++)
            plotFinish[i]=interp1Darray(finish[i], sumsPerOrt);
        p = new Plot("", "regularMS,"+maxIterNonBlurring+" iters,h="+h, "", angles, sumsPerOrt);
        p.addPoints(finish, plotFinish, Plot.X);

        Vector<double[]> cls = extractClusters(finish, minD, M);

        double[] anglesDirections;

        if (cls.size()>=3) { // location -> cls.get(i)[0] elements -> cls.get(i)[1]

            // extract directions
            boolean[] checked = new boolean[cls.size()]; // all to false
            anglesDirections = new double[3];

            // find top 3
            for (int k = 0; k<3; k++) {

                // reset max search
                double  currMax = Double.MIN_VALUE;
                int     currMaxIdx = -1;

                for (int i=0; i<cls.size(); i++) {

                    // find max in this round
                    if (!checked[i]) {
                        if (cls.get(i)[1]>currMax) {

                            currMax = cls.get(i)[1];
                            currMaxIdx = i;

                        }
                    }
                }

                checked[currMaxIdx] = true;

                anglesDirections[k] = cls.get(currMaxIdx)[0];

            }

            double[] plotClsY = new double[3];
            double[] plotClsX = new double[3];

			for (int i=0; i<3; i++) {
                plotClsY[i] = interp1Darray(anglesDirections[i], sumsPerOrt);
                plotClsX[i] = anglesDirections[i];
            }

			p.addPoints(plotClsX, plotClsY, Plot.BOX);

		}
        else {

            anglesDirections = null;

        }

		isOut.addSlice(p.getProcessor());

        return isOut;

    }

	private float 	point2dir(
		float[] n, 		// direction vector
		int[] 	p     	// point considered
	)
	{
		// line is in vector from b in n direction

		float d = 0;

		float[] p_b = new float[2];

		// p - b
		p_b[0] = p[0];// - b[0];
		p_b[1] = p[1];// - b[1];

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

    public double[] runMS(
		double[] 	start,
		float[] 	inputProfile,
		int 		max_iter,
		double 		epsilon,
		int 		h
	){

        double[] T				= new double[start.length];

        for (int i = 0; i < T.length; i++) {
            T[i] = start[i];
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

        }

        return T;

    }

	private	double[] runBlurring(
		double[] 	start,
		float[] 	inputProfile,
		int 		maxIter,
		double 		epsilon,
		int 		h
	)
    {

		double[]	S = new double[start.length];
		double[] 	T = new double[start.length];
		double 		dmax;
		int 		iter = 0;

		for (int i=0; i<start.length; i++) {
			S[i] = T[i] = start[i];

		}

		do{

			dmax = Double.MIN_VALUE;

			for (int i = 0; i < S.length; i++) {

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

			for (int i = 0; i < S.length; i++) {
				S[i] = T[i];
			}

			iter++;

		}
		while(iter<=maxIter && dmax>epsilon);

		return T;

	}

	public Vector<double[]> extractClusters(
            double[] msConvergence,
            double range,
            int M
    )
    {
		// 'range' describes intra neighborhood range size,
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

		Vector<double[]> clust = new Vector<double[]>();
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
					clust.add(new double[] {sum/cnt, cnt});
				}

			}
		}

		return clust;

	}

    private double interp1Darray(
            double 	x,
            float[] inarray
    )
    {

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

	private double 	runOne(
            double  curr_pos,
            int     h,
            float[] inputProfile){
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

    private double wrap_0_2PI(
            double  in
    )
    {

        double out = in;

        while(out<0){
            out += 2*Math.PI;
        }
        while(out>=2*Math.PI){
            out -= 2*Math.PI;
        }

        return out;
    }

	public double[] scores(
		int             atX,
		int             atY,
		FloatProcessor  inip,
		int             Npoints,
		int             Niter,
		double          eps,
		int             M,
		double          minD
	)
	{

		long t1 = System.currentTimeMillis(); // calculation time

		double[] angs = getAngles( 	atX, atY,
		  							inip,
		             				Npoints,
		             				Niter,
		          					eps,
		             				M,
		          					minD
		); // will give 3 angles or nothing

		if (angs!=null) {



			int sumON, sumOFF, nrON, nrOFF;
			sumON = sumOFF = nrON = nrOFF = 0;

			int sumA0, nrA0, sumA1, nrA1, sumA2, nrA2, sumA3, nrA3;
			sumA0 = nrA0 = sumA1 = nrA1 = sumA2 = nrA2 = sumA3 = nrA3 = 0;

			int sumB1, nrB1, sumB2, nrB2, sumB3, nrB3;
			sumB1 = nrB1 = sumB2 = nrB2 = sumB3 = nrB3 = 0;

			// skeleton is similar to the one from  calculateTemplate( angles )
			// it is not necessary to calculate the template ArrayList

			// form the template elements
			for (int pIdx = 0; pIdx < 3; pIdx++)
				ap[pIdx] = angs[pIdx];

			Arrays.sort(ap);

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if ( d2 <= rInner*rInner ) {

						p[0] = x-xc;
						p[1] = -(y-yc);

						sumA0 += inip.getf(atX+p[0], atY+p[1]); nrA0++;
						sumON  += inip.getf(atX+p[0], atY+p[1]);nrON++;
						//templateCenter.add(new int[]{p[0], p[1]});

					}
				}
			}

			// template(1, 2, 3)
			for (int pIdx = 0; pIdx < 3; pIdx++) {

				double ang = angs[pIdx];

//				ArrayList<int[]> templateAngle = new ArrayList<int[]>();

				for (int x = 0; x < d; x++) {
					for (int y = 0; y < d; y++) {

						int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

						if (d2 <= r*r && d2 >= rLower*rLower) {

							p[0] = x-xc;
							p[1] = -(y-yc);

							n[0] = (float) Math.cos(ang);
							n[1] = (float) Math.sin(ang);

							float dst = point2dir(n, p);

							if (dst<=diam/2) { // belongs to pIdx ON peak and not filled  // && idxMapLocal[x+d*y]==0

								sumON  += inip.getf(atX+p[0], atY+p[1]);nrON++;
								if (pIdx==0) {
									sumA1 += inip.getf(atX+p[0], atY+p[1]); nrA1++;
								}
								else if (pIdx==1) {
									sumA2 += inip.getf(atX+p[0], atY+p[1]); nrA2++;
								}
								else if (pIdx==2) {
									sumA3 += inip.getf(atX+p[0], atY+p[1]); nrA3++;
								}
								//templateAngle.add(new int[]{p[0], p[1]});
							}
						}
					}
				}
				//template.add(templateAngle); // added template(1, 2, 3)
			}

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (d2 <= r*r && d2 >= rLower*rLower) {

						p[0] = x-xc;
						p[1] = -(y-yc);

						boolean isON = false;

						for (int pIdx = 0; pIdx < 3; pIdx++) {

							double ang = angs[pIdx];

							n[0] = (float) Math.cos(ang);
							n[1] = (float) Math.sin(ang);

							float dst = point2dir(n, p);

							if (dst<=diam/2) {
								isON = true;
								break;
							}

						}

						if (!isON) {

							// check where it is
							double a = wrap_0_2PI(Math.atan2(p[1], p[0]));

							boolean added = false;

							sumOFF  += inip.getf(atX+p[0], atY+p[1]);nrOFF++;

							if (wrap_PI(a-ap[0])>0 && wrap_PI(a-ap[1])<=0 ) {  // 0-1
								sumB1 += inip.getf(atX+p[0], atY+p[1]); nrB1++;
								//templateBtwAngle12.add(new int[]{p[0], p[1]});
							}
							else if (wrap_PI(a-ap[1])>0 && wrap_PI(a-ap[2])<=0 ) { // 1-2
								sumB2 += inip.getf(atX+p[0], atY+p[1]); nrB2++;
								//templateBtwAngle23.add(new int[]{p[0], p[1]});
							}
							else {
								sumB3 += inip.getf(atX+p[0], atY+p[1]); nrB3++;
								//templateBtwAngle31.add(new int[]{p[0], p[1]}); // 2-0
							}

						}

					}
				}
			}

			long t2 = System.currentTimeMillis(); // calculation time

			IJ.log(((t2-t1)/1000f)+" sec. calculated: "+sumON+" ("+nrON+"), "+sumOFF+" ("+nrON+")");

			return new double[]{};
		}
		else {
			return null;
		}

	}

	public ImageStack exportTemplate(
	    double[] 	directionAngles
	)
	{

		return plotTemplate(calculateTemplate(directionAngles));

	}

	private ArrayList<ArrayList<int[]>> calculateTemplate(
		double[] 	directionAngles
	)
	{
		// form the template elements
		for (int pIdx = 0; pIdx < 3; pIdx++)
			ap[pIdx] = directionAngles[pIdx];

		Arrays.sort(ap);

		template.clear();

		// template(0)
		ArrayList<int[]> templateCenter = new ArrayList<int[]>();
		for (int x = 0; x < d; x++) {
			for (int y = 0; y < d; y++) {

				int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

				if ( d2 <= rInner*rInner ) {

					p[0] = x-xc;
					p[1] = -(y-yc);

					templateCenter.add(new int[]{p[0], p[1]});

				}
			}
		}
		template.add(templateCenter);  // template(0) added

		// template(1, 2, 3)
		for (int pIdx = 0; pIdx < 3; pIdx++) {

			double ang = directionAngles[pIdx];

			ArrayList<int[]> templateAngle = new ArrayList<int[]>();

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (d2 <= r*r && d2 >= rLower*rLower) {

						p[0] = x-xc;
						p[1] = -(y-yc);

						n[0] = (float) Math.cos(ang);
						n[1] = (float) Math.sin(ang);

						float dst = point2dir(n, p);

						if (dst<=diam/2) { // belongs to pIdx ON peak and not filled  // && idxMapLocal[x+d*y]==0

							templateAngle.add(new int[]{p[0], p[1]});

						}


					}
				}
			}

			template.add(templateAngle); // added template(1, 2, 3)

		}

		ArrayList<int[]> templateBtwAngle12 = new ArrayList<int[]>();
		ArrayList<int[]> templateBtwAngle23 = new ArrayList<int[]>();
		ArrayList<int[]> templateBtwAngle31 = new ArrayList<int[]>();

		for (int x = 0; x < d; x++) {
			for (int y = 0; y < d; y++) {

				int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

				if (d2 <= r*r && d2 >= rLower*rLower) {

					p[0] = x-xc;
					p[1] = -(y-yc);

					boolean isON = false;

					for (int pIdx = 0; pIdx < 3; pIdx++) {

						double ang = directionAngles[pIdx];

						n[0] = (float) Math.cos(ang);
						n[1] = (float) Math.sin(ang);

						float dst = point2dir(n, p);

						if (dst<=diam/2) {
							isON = true;
							break;
						}

					}

					if (!isON) {

						// check where it is
						double a = wrap_0_2PI(Math.atan2(p[1], p[0]));

						boolean added = false;

						if (wrap_PI(a-ap[0])>0 && wrap_PI(a-ap[1])<=0 ) {  // 0-1
							templateBtwAngle12.add(new int[]{p[0], p[1]});
						}
						else if (wrap_PI(a-ap[1])>0 && wrap_PI(a-ap[2])<=0 ) { // 1-2
							templateBtwAngle23.add(new int[]{p[0], p[1]});
						}
						else {
							templateBtwAngle31.add(new int[]{p[0], p[1]}); // 2-0
						}

					}

				}
			}
		}

		template.add(templateBtwAngle12);
		template.add(templateBtwAngle23);
		template.add(templateBtwAngle31);  // added template(4, 5, 6)

		return template;
	}

	private ImageStack plotTemplate(ArrayList<ArrayList<int[]>> template)
	{

		ImageStack isOut = new ImageStack(d, d);
		ImageProcessor ip = new ByteProcessor(d, d);

		// inner template(0)
		if(template.size()>0) {
			for (int b = 0; b<template.get(0).size(); b++) {
				int offX = template.get(0).get(b)[0];
				int offY = template.get(0).get(b)[1];
				ip.setf(d/2+offX, d/2+offY, +255);
			}
		}

		// on
		//if(template.size()>=4) {
		for (int templateIdx=1; templateIdx<=3; templateIdx++) {
			for (int b = 0; b<template.get(templateIdx).size(); b++) {
				int offX = template.get(templateIdx).get(b)[0];
				int offY = template.get(templateIdx).get(b)[1];
				ip.setf(d/2+offX, d/2+offY, +255);
			}
		}
		//}

		// off
//		if(template.size()<=7) {
		for (int templateIdx=4; templateIdx<=6; templateIdx++) {
			for (int b = 0; b<template.get(templateIdx).size(); b++) {
				int offX = template.get(templateIdx).get(b)[0];
				int offY = template.get(templateIdx).get(b)[1];
				ip.setf(d/2+offX, d/2+offY, +0);
			}
		}
//		}

		isOut.addSlice(ip);

		return isOut;

	}

    private static double wrap_PI(
        double in
    )
    {
        double out = in;

        while(out<=-Math.PI){
            out += 2*Math.PI;
        }
        while(out>Math.PI){
            out -= 2*Math.PI;
        }

        return out;
    }

}
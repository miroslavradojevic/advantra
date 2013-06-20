import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.RefineryUtilities;

import java.awt.*;
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

	int					resolDeg = 10;
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
    double[]    start;

    // ms
    int Niter = 100;
    double eps = 0.001;
    double minD = 0.5;
    int M = 1;
    int Kpts = 1;

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

		//resolDeg = 10;

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

        // ms-variables fixed parameters
        start = new double[offsets.size()*Kpts];
        for (int i=0; i<start.length; i++)
            start[i] = ((double)i/start.length)*sumsPerOrt.length;

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
        boolean         print
    )
    {

		// check if values can be taken for the profile
		int margin = d/2+1;
		if ( atX<=margin || atX>=inip.getWidth()-margin ) {
			return null;
		}
		if ( atY<=margin || atY>=inip.getHeight()-margin ) {
			return null;
		}

        for (int dirIdx = 0; dirIdx<sumsPerOrt.length; dirIdx++) {
            sumsPerOrt[dirIdx] = 0;     // reset
        }

		for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {
			for (int locIdx = 0; locIdx<offsets.get(dirIdx).size(); locIdx++) {
				sumsPerOrt[dirIdx] += inip.getf(
					atX+offsets.get(dirIdx).get(locIdx)[0], atY+offsets.get(dirIdx).get(locIdx)[1]);
			}
		}

		double[] finish;

		finish = runMS(start, sumsPerOrt, Niter, eps, h);

        Vector<double[]> cls = extractClusters(finish, minD, M);

        double[] anglesDirections = null;
        boolean[] checked = new boolean[cls.size()]; // all to false

        if (cls.size()<=0) {
            anglesDirections = null;
        }

        if (cls.size()==1) {
            anglesDirections = new double[1];
            anglesDirections[0] = (cls.get(0)[0]+.0f) * (double)resolDeg;
            anglesDirections[0] = (anglesDirections[0] * (TwoPI/360));
            anglesDirections[0] = wrap_0_2PI(anglesDirections[0]);
        }

        if (cls.size()==2) {
            anglesDirections = new double[2];
            anglesDirections[0] = (cls.get(0)[0]+.0f) * (double)resolDeg;
            anglesDirections[0] = (anglesDirections[0] * (TwoPI/360));
            anglesDirections[0] = wrap_0_2PI(anglesDirections[0]);

            anglesDirections[1] = (cls.get(1)[0]+.0f) * (double)resolDeg;
            anglesDirections[1] = (anglesDirections[1] * (TwoPI/360));
            anglesDirections[1] = wrap_0_2PI(anglesDirections[1]);
        }

        if (cls.size()>=3) {

            // extract 3 angles with most convergence points
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
                // set the output value
                anglesDirections[k] = (cls.get(currMaxIdx)[0]) * (double)resolDeg;
                anglesDirections[k] = (anglesDirections[k] * (TwoPI/360));
                anglesDirections[k] = wrap_0_2PI(anglesDirections[k]);

            }

        }

        if (print) {

            // find peaks - initialize start points for ms iterations
            double[] plotStart = new double[start.length];
            double[] plotFinish = new double[finish.length];

            for (int i=0; i<plotStart.length; i++)
                plotStart[i]=interp1Darray(start[i], sumsPerOrt);

            for (int i=0; i<plotFinish.length; i++)
                plotFinish[i]=interp1Darray(finish[i], sumsPerOrt);

            float[] angles = new float[offsets.size()];
            for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {
                angles[dirIdx] = dirIdx;
            }

            Plot p;
            p = new Plot("", start.length+" points start", "", start, plotStart);
            p.draw();
            p.setColor(Color.BLUE);
            p.addPoints(finish, plotFinish, Plot.X);

            if (cls.size()==1) {
                double[] cluster = new double[1];
                double[] plotCluster = new double[1];
                cluster[0] = cls.get(0)[0];
                plotCluster[0] = interp1Darray(cluster[0], sumsPerOrt);
                p.setColor(Color.RED);
                p.addPoints(cluster, plotCluster, Plot.BOX);
            }

            if (cls.size()==2) {
                double[] cluster = new double[2];
                double[] plotCluster = new double[2];
                cluster[0] = cls.get(0)[0];
                cluster[1] = cls.get(1)[0];
                plotCluster[0] = interp1Darray(cluster[0], sumsPerOrt);
                plotCluster[1] = interp1Darray(cluster[1], sumsPerOrt);
                p.setColor(Color.RED);
                p.addPoints(cluster, plotCluster, Plot.BOX);
            }

            if (cls.size()>=3) { // take max
                double[] cluster = new double[3];
                double[] plotCluster = new double[3];

                int cnt = 0;

                for (int k=0; k<cls.size(); k++) {

                    if (checked[k]) {
                        cluster[cnt] = cls.get(k)[0];
                        plotCluster[cnt] = interp1Darray(cls.get(k)[0], sumsPerOrt);
                        cnt++;
                    }

                }
                p.setColor(Color.RED);
                p.addPoints(cluster, plotCluster, Plot.BOX);


            }

            new ImagePlus("print", p.getProcessor()).show();
            //p.show();

        }

        return anglesDirections;

	}

    public ImageStack plotIntegResponse(
            int             atX,
            int             atY,
            FloatProcessor  inip
    )
    {
        // calculate sums at each orientation
        float[] angles      = new float[offsets.size()];
        float[] intSum  = new float[offsets.size()]; // substitute class variable

        for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {
            angles[dirIdx] = dirIdx;
            for (int locIdx = 0; locIdx<offsets.get(dirIdx).size(); locIdx++) {
                intSum[dirIdx] += inip.getf(
                        atX+offsets.get(dirIdx).get(locIdx)[0], atY+offsets.get(dirIdx).get(locIdx)[1]);
            }
        }

        // plot angles vs. sums
        Plot p;
        ImageStack isOut;
        p = new Plot("", "orientation", "integrated response", angles, intSum);
        p.setColor(Color.BLACK);
        isOut = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
        isOut.addSlice(p.getProcessor());

        return isOut;

    }

    public void plotProfilesWithMSDetection(
		int 			atX,
		int 			atY,
		FloatProcessor 	inip
	)
	{
        /*
		finish = runBlurring(start, intSum, maxIterBlurring, eps, h);
		plotFinish = new double[finish.length];
		for (int i=0; i<finish.length; i++)
				plotFinish[i]=interp1Darray(finish[i], intSum);

		p = new Plot("", "blurring,"+maxIterBlurring+" iters,h="+h, "", angles, intSum);
        p.setColor(Color.BLACK); // so that the color can be added later
		p.addPoints(finish, plotFinish, Plot.BOX);
        isOut.addSlice("with_blurring", p.getProcessor());
        */
        double[] anglesDirections = getAngles(atX, atY, inip, true);

        if (anglesDirections!=null) {
            IJ.log("\nfinal angles:");
            for (int g = 0; g<anglesDirections.length; g++) {
                IJ.log("angle"+g+" > "+(anglesDirections[g]*(360f/TwoPI))+" degrees");
            }
        }

        score(atX, atY, inip, true);

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

    private double[] runMS(
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
				cluster_sizes.add(1);
				nr_clusters++;
			}

		}

		Vector<double[]> clust = new Vector<double[]>();
		for (int k = 0; k < cluster_sizes.size(); k++) {
			if (cluster_sizes.get(k)>=M) {

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

    /*
    //
     */
    public int nrDirections(
            int             atX,
            int             atY,
            FloatProcessor  inip
    ){

        double[] angs = getAngles(atX, atY, inip, false);
        return (angs!=null)? angs.length : 0 ;

    }

    public ImageProcessor nrDirections(
            FloatProcessor  inip
    ){

        FloatProcessor sc = new FloatProcessor(inip.getWidth(), inip.getHeight());

        // fill the values in
        for (int x = 0; x<inip.getWidth(); x++) {
            for (int y = 0; y<inip.getHeight(); y++) {

                //System.out.println(""+x+","+y);
                int nrd = nrDirections(x, y,inip);
                sc.setf(x, y, ((nrd>=3)? 255f : 0f) );

            }
        }

        return sc;

    }

	public double score(
		int             atX,
		int             atY,
		FloatProcessor  inip,
		boolean 		print
	)
	{

		//long t1 = System.currentTimeMillis(); // calculation time

		double[] angs = getAngles(atX, atY, inip, false);

		if (angs!=null && angs.length>=3) {   // return null otherwise

//			int sumON, sumOFF, nrON, nrOFF;
//			sumON = sumOFF = 1;
//			nrON = nrOFF = 0;


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
						//sumON  += inip.getf(atX+p[0], atY+p[1]);nrON++;
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

								//sumON  += inip.getf(atX+p[0], atY+p[1]);nrON++;

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

							//sumOFF  += inip.getf(atX+p[0], atY+p[1]);nrOFF++;

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

			//long t2 = System.currentTimeMillis(); // calculation time

            //double avgON    = (nrON>0)? (sumON/nrON) : (0);
            //double avgOFF   = (nrOFF>0)? (sumOFF/nrOFF) : (0);

            double avgA0   = (nrA0>0)? (sumA0/nrA0) : (0);

			double avgA1   = (nrA1>0)? (sumA1/nrA1) : (0);
            double avgA2   = (nrA2>0)? (sumA2/nrA2) : (0);
            double avgA3   = (nrA3>0)? (sumA3/nrA3) : (0);

            double avgB1   = (nrB1>3)? (sumB1/nrB1) : (Double.MAX_VALUE);
            double avgB2   = (nrB2>3)? (sumB2/nrB2) : (Double.MAX_VALUE);
            double avgB3   = (nrB3>3)? (sumB3/nrB3) : (Double.MAX_VALUE);

			double L0 	= (3*avgA0-avgB1-avgB2-avgB3>0)? ((3*avgA0-avgB1-avgB2-avgB3)/3) : 0;

			double L11 = (2*avgA1-avgB1-avgB3>0)? ((2*avgA1-avgB1-avgB3)/2) : 0;
//			double L12 = (avgA1-avgB3>0)? (avgA1-avgB3) : 0;

			double L21 = (2*avgA2-avgB2-avgB1>0)? ((2*avgA2-avgB2-avgB1)/2) : 0;
//			double L22 = (avgA2-avgB1>0)? (avgA2-avgB1) : 0;

			double L31 = (2*avgA3-avgB3-avgB2>0)? ((2*avgA3-avgB3-avgB2)/2) : 0;
//			double L32 = (avgA3-avgB2>0)? (avgA3-avgB2) : 0;

			if(print) {

//				final double[][] data = new double[][] {
//                    {avgA1, 	avgA2, 	avgA3, 	avgA0},
//                    {avgB1, 	avgB2, 	avgB3, 	0},
//                    {L11, 	    L21, 	L31, 	avgA0-0},
//				};

                /*
                ArrayList<double[]> data = new ArrayList<double[]>();
                                    //x,  y,  z
                data.add(new double[]{0, 	1, 	2, 	0,  0});
                data.add(new double[]{1, 	2, 	3, 	1,  0});
                data.add(new double[]{2, 	3, 	4, 	2,  0});
                */
                JFreeChartTools chart = new JFreeChartTools("MyPlot");
                chart.plotVec("profile", "",      "int.sum.", sumsPerOrt);
                //chart.plotVec("title", "start", "",         start);

				/*
				MyBarChart chart = new MyBarChart("", data);
				chart.pack();
				RefineryUtilities.centerFrameOnScreen(chart);
				chart.setVisible(true);
				*/

			}

            int D = 20;
			return ( 1-Math.exp(-(L11)/D) ) * ( 1-Math.exp(-(L21)/D) ) * ( 1-Math.exp(-(L31)/D) ) * ( 1-Math.exp(-(avgA0)/D) );

		}
		else {
            // if angs was null or <=2 for some reason
			return 0;
		}

	}

    public ImageStack score(
            FloatProcessor  inip
    )
    {
        FloatProcessor score0 = new FloatProcessor(inip.getWidth(), inip.getHeight());
        //FloatProcessor score1 = new FloatProcessor(inip.getWidth(), inip.getHeight());

        // fill the values in
        for (int x = 0; x<inip.getWidth(); x++) {
            for (int y = 0; y<inip.getHeight(); y++) {

				double sc = score(x, y, inip, false);

				score0.setf(x, y, (float)sc);

                /*if (sc!=null && sc.length==2) {
                    score0.setf(x, y, (float)sc[0]);
                    score1.setf(x, y, (float)sc[1]);
                }
                else {
                    score0.setf(x, y, Float.MIN_VALUE);
                    score1.setf(x, y, Float.MIN_VALUE);
                }*/

            }
        }

        ImageStack isOut = new ImageStack(inip.getWidth(), inip.getHeight());
        isOut.addSlice(score0);
//        isOut.addSlice(score1);
        return isOut;

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
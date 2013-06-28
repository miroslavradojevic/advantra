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

	int					resolDeg;// = 10;
	int					h;

	ArrayList<ArrayList<int[]>> template;

	ArrayList<ArrayList<double[]>> offsetsReal;

	float[] sumsPerOrt;

	float TwoPI = (float) (Math.PI*2);

	// will be used when forming kernels
	double 		step = 0.5;
	int 		xc;
	int 		yc;
	float[] 	n;
	int[] 		p;
    double[]    start;

    // mean-shift
    int Niter   = 100;
    double eps  = 0.001;
    double minD = 0.5;
    int M       = 1;
    float Kpts  = 1.5f;

    // features used
    double A0=0, A1=0, A2=0, A3=0, B1=0, B2=0, B3=0;
    int nA0=0, nA1=0, nA2=0, nA3=0, nB1=0, nB2=0, nB3=0;
	double[] 	ap;
	double[] feats;

	public Feat(
	    int diam,
		double scale
	)
	{

		this.diam = diam;
		r = (int) Math.round(diam*scale);
		r = (r<6)? 6 : r; // lower limit
		d = 2*r+1;

		resolDeg = (int) ( Math.round( ( 2*Math.asin(1f/(2 * scale))*(1f/4) / TwoPI) * 360 ) );
		resolDeg = (resolDeg>=1)? resolDeg : 1;

		rInner = diam/2;
		rLower = 3;//(int) Math.round(r/2.0); // (1.0*scale) // 2*scale at max

		double  angStep     = 2 * Math.asin((float)diam/(2*r));
		int     angStepDeg  = (int) Math.round( ((angStep/(2*Math.PI))*360) );
		angStepDeg = (angStepDeg<resolDeg)?resolDeg:angStepDeg;
		angStepDeg  = (angStepDeg/resolDeg)*resolDeg; // will be used for mean-shift, defines the neighbourhood to search
		h = angStepDeg/resolDeg;
		h = (h<1)?1:h;

		xc = d/2;
		yc = d/2;
		n = new float[2];
		p = new int[2];
		ap 	= new double[3];
		feats = new double[9];

		// form the list of offsets
		//offsets = new ArrayList<ArrayList<int[]>>();
		offsetsReal = new ArrayList<ArrayList<double[]>>();

		for (int aDeg = 0; aDeg<360; aDeg+=resolDeg) {

			float aRad = ((float)aDeg/360)*TwoPI;

			// offsets real
			ArrayList<double[]> offsetsAngleReal = new ArrayList<double[]>();
			for (double x = 0; x < d; x+=step) {
				for (double y = 0; y < d; y+=step) {

					double d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (d2 <= r*r && Math.floor(d2) >= rLower*rLower) {

						double px = x-xc;
						double py = -(y-yc);

						double nx = (float) Math.cos(aRad);
						double ny = (float) Math.sin(aRad);

						double dst = point2dir(nx, ny, px, py);

						if (dst<=diam/2) { // belongs to pIdx ON peak and not filled

							offsetsAngleReal.add(new double[]{px, py});

						}

					}
				}
			}

			offsetsReal.add(offsetsAngleReal);

		}

		template = new ArrayList<ArrayList<int[]>>(7); // nrPeaks+1
		//templateReal = new ArrayList<ArrayList<double[]>>(7);

		sumsPerOrt  = new float[offsetsReal.size()];

        // ms-variables fixed parameters
        int nrStartPoints = (int) Math.round(sumsPerOrt.length*Kpts);
        start = new double[nrStartPoints];
        for (int i=0; i<start.length; i++) {
            start[i] = ((double)i/nrStartPoints)*sumsPerOrt.length;
        }

	}

    public ImageStack plotOffsets()
	{
        ImageStack isOut = new ImageStack(d, d);

        for (int a = 0; a<offsetsReal.size(); a++) {

            ImageProcessor ip = new ByteProcessor(d, d);

            for (int b = 0; b<offsetsReal.get(a).size(); b++) {

                int offX = (int) Math.round(offsetsReal.get(a).get(b)[0]);
                int offY = (int) Math.round(offsetsReal.get(a).get(b)[1]);

                ip.setf(d/2+offX, d/2+offY, +255);

            }

            isOut.addSlice(ip);

        }

        return isOut;

    }

	public  ImageProcessor getAngles(
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

        for (int v = 0; v<sumsPerOrt.length; v++) sumsPerOrt[v] = 0;     // reset

		for (int d = 0; d<offsetsReal.size(); d++) {
			for (int l = 0; l<offsetsReal.get(d).size(); l++) {
				sumsPerOrt[d] += Interpolator.interpolateAt(
					atX+offsetsReal.get(d).get(l)[0],
					atY+offsetsReal.get(d).get(l)[1],
					inip
				);
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

        if (cls.size()==3) {
            anglesDirections = new double[3];
            anglesDirections[0] = (cls.get(0)[0]+.0f) * (double)resolDeg;
            anglesDirections[0] = (anglesDirections[0] * (TwoPI/360));
            anglesDirections[0] = wrap_0_2PI(anglesDirections[0]);

            anglesDirections[1] = (cls.get(1)[0]+.0f) * (double)resolDeg;
            anglesDirections[1] = (anglesDirections[1] * (TwoPI/360));
            anglesDirections[1] = wrap_0_2PI(anglesDirections[1]);

            anglesDirections[2] = (cls.get(2)[0]+.0f) * (double)resolDeg;
            anglesDirections[2] = (anglesDirections[2] * (TwoPI/360));
            anglesDirections[2] = wrap_0_2PI(anglesDirections[2]);
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

		if (anglesDirections==null) {
			ap = null;
		}
		else {
			for (int pIdx = 0; pIdx < anglesDirections.length; pIdx++)
				ap[pIdx] = anglesDirections[pIdx];
		}

		Arrays.sort(ap);

		Plot p = null;

        if (print) {

            // find peaks - initialize start points for ms iterations
            double[] plotStart = new double[start.length];
            double[] plotFinish = new double[finish.length];

            for (int i=0; i<plotStart.length; i++)
                plotStart[i]=interp1Darray(start[i], sumsPerOrt);

            for (int i=0; i<plotFinish.length; i++)
                plotFinish[i]=interp1Darray(finish[i], sumsPerOrt);

            float[] angles = new float[offsetsReal.size()];
            for (int dirIdx = 0; dirIdx<offsetsReal.size(); dirIdx++) {
                angles[dirIdx] = dirIdx;
            }

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
			//p.show();

        }

        if (p==null) {
            return null;
        }
        else {
            return p.getProcessor();
        }


	}

    public void plotIntegResponse(
            int             atX,
            int             atY,
            FloatProcessor  inip
    )
    {
        // calculate sums at each orientation
        float[] angles      = new float[offsetsReal.size()];
        float[] intSum  	= new float[offsetsReal.size()]; // substitute class variable

        for (int d = 0; d<offsetsReal.size(); d++) {
            angles[d] = d;
            for (int l = 0; l<offsetsReal.get(d).size(); l++) {
                intSum[d] += Interpolator.interpolateAt(atX+offsetsReal.get(d).get(l)[0], atY+offsetsReal.get(d).get(l)[1], inip);
            }
        }

        // plot angles vs. sums
        Plot p;
        ImageStack isOut;
        p = new Plot("", "orientation", "integrated response", angles, intSum);
        p.setColor(Color.BLACK);
        p.show();
//		isOut = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
//      isOut.addSlice(p.getProcessor());
//      return isOut;

    }

	public void regionScores(
		int             atX,
		int             atY,
		FloatProcessor  inip,
		double[] 		angs
	)
	{

        A0 = A1 = A2 = A3 =     nA0 = nA1 = nA2 = nA3   = 0;
        B1 = B2 = B3 =          nB1 = nB2 = nB3         = 0;

		if (angs!=null && angs.length>=3) {   // return 0 otherwise

			for (int pIdx = 0; pIdx < 3; pIdx++)
				ap[pIdx] = angs[pIdx];

			Arrays.sort(ap);

			// real version
			for (double x = 0; x < d; x+=step) {
				for (double y = 0; y < d; y+=step) {

					double d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					double px = x-xc;
					double py = -(y-yc);

					if (d2 <= r*r && d2 >= rLower*rLower) {

						boolean isON = false;

						for (int pIdx = 0; pIdx < 3; pIdx++) {

							//double ang = angs[pIdx];

							double nx = (float) Math.cos(ap[pIdx]);
							double ny = (float) Math.sin(ap[pIdx]);

							double dst = point2dir(nx, ny, px, py);

                            if (dst<=(double)diam) {

                                if (dst<=(double)diam/2) {

                                    if (pIdx==0) {
                                        A1 += Interpolator.interpolateAt(atX+px, atY+py, inip);
                                        nA1++;
                                    }
                                    else if (pIdx==1) {
                                        A2 += Interpolator.interpolateAt(atX+px, atY+py, inip);
                                        nA2++;
                                    }
                                    else if (pIdx==2) {
                                        A3 += Interpolator.interpolateAt(atX+px, atY+py, inip);
                                        nA3++;
                                    }

                                }
                                else {

                                    if (pIdx==0) {
                                        B1 += Interpolator.interpolateAt(atX+px, atY+py, inip);
                                        nB1++;
                                    }
                                    else if (pIdx==1) {
                                        B2 += Interpolator.interpolateAt(atX+px, atY+py, inip);
                                        nB2++;
                                    }
									else if (pIdx==2) {
										B3 += Interpolator.interpolateAt(atX+px, atY+py, inip);
										nB3++;
									}

                                }


                            }

//							if (dst<=(double)diam/2) {
//								isON = true;
//								break;
//							}

						}

//						if (!isON) {
//
//							// check where it is
//							double a = wrap_0_2PI(Math.atan2(py, px));
//
//							if (wrap_PI(a-ap[0])>0 && wrap_PI(a-ap[1])<=0 ) {  // 0-1
//								B1 += Interpolator.interpolateAt(atX+px, atY+py, inip);
//								nB1++;
//							}
//							else if (wrap_PI(a-ap[1])>0 && wrap_PI(a-ap[2])<=0 ) { // 1-2
//								B2 += Interpolator.interpolateAt(atX+px, atY+py, inip);
//								nB2++;
//							}
//							else {
//								B3 += Interpolator.interpolateAt(atX+px, atY+py, inip);
//								nB3++;
//							}
//
//						}

					}
					else if ( d2 <= (double)rInner*(double)rInner ) {

						A0 += Interpolator.interpolateAt(atX+px, atY+py, inip);//inip.getf(atX+p[0], atY+p[1]);
						nA0++;

					}
				}
			}

		}

	}

	public void extractFeatures(int atX, int atY, FloatProcessor  inip)
	{

		getAngles(atX, atY, inip, false);

		regionScores(atX, atY, inip, ap); // forms ap, A0, A1, nA0...

        double aA0 = (nA0>0)? (A0/nA0) : Double.NaN;
        double aA1 = (nA1>0)? (A1/nA1) : Double.NaN;
        double aA2 = (nA2>0)? (A2/nA2) : Double.NaN;
        double aA3 = (nA3>0)? (A3/nA3) : Double.NaN;
        double aB1 = (nB1>0)? (B1/nB1) : Double.NaN;
        double aB2 = (nB2>0)? (B2/nB2) : Double.NaN;
        double aB3 = (nB3>0)? (B3/nB3) : Double.NaN;

        // change here!
        feats[0] = aA0-aB1;
        feats[1] = aA0-aB2;
        feats[2] = aA0-aB3;
        feats[3] = aA1-aB1;
        feats[4] = aA1-aB3;
        feats[5] = aA2-aB2;
        feats[6] = aA2-aB1;
        feats[7] = aA3-aB3;
        feats[8] = aA3-aB2;

	}

    public double bifurcationess(
        int     atX,
        int     atY,
        FloatProcessor  inip,
        double D
    )
    {
		getAngles(atX, atY, inip, false);

        regionScores(atX, atY, inip, ap);

        double avgA0   = (nA0>3)? (A0/nA0) : (Double.MIN_VALUE);
		double avgA1   = (nA1>3)? (A1/nA1) : (Double.MIN_VALUE);
        double avgA2   = (nA2>3)? (A2/nA2) : (Double.MIN_VALUE);
        double avgA3   = (nA3>3)? (A3/nA3) : (Double.MIN_VALUE);
        double avgB1   = (nB1>3)? (B1/nB1) : (Double.MAX_VALUE);
        double avgB2   = (nB2>3)? (B2/nB2) : (Double.MAX_VALUE);
        double avgB3   = (nB3>3)? (B3/nB3) : (Double.MAX_VALUE);

		double L0 	= (3*avgA0-avgB1-avgB2-avgB3)/3;
        L0 = (L0>0)? L0 : 0;

		double L11 = (2*avgA1-avgB1-avgB3)/2;
        L11 = (L11>0)? L11 : 0;
		double L21 = (2*avgA2-avgB2-avgB1)/2;
        L21 = (L21>0)? L21 : 0;
        double L31 = (2*avgA3-avgB3-avgB1)/2;
        L31 = (L31>0)? L31 : 0;

		double U1 = (avgA1+avgA2-2*avgB1)/2;
		U1 = (U1>0)? U1 : 0;
		double U2 = (avgA2+avgA3-2*avgB2)/2;
		U2 = (U2>0)? U2 : 0;
		double U3 = (avgA3+avgA1-2*avgB3)/2;
		U3 = (U3>0)? U3 : 0;

            return 	(1-Math.exp(-L11/D)) *
					(1-Math.exp(-L21/D)) *
                    (1-Math.exp(-L31/D)) *
					(1-Math.exp(-U1/D)) *
					(1-Math.exp(-U2/D)) *
					(1-Math.exp(-U3/D)) *
                    (1-Math.exp(-L0/D));

    }

	public double bifurcationess1(
	int     atX,
	int     atY,
	FloatProcessor  inip,
	double  D,
	double 	E
	)
	{
		getAngles(atX, atY, inip, false);

		regionScores(atX, atY, inip, ap);

		double avgA0   = (nA0>3)? (A0/nA0) : (Double.MIN_VALUE);
		double avgA1   = (nA1>3)? (A1/nA1) : (Double.MIN_VALUE);
		double avgA2   = (nA2>3)? (A2/nA2) : (Double.MIN_VALUE);
		double avgA3   = (nA3>3)? (A3/nA3) : (Double.MIN_VALUE);
		double avgB1   = (nB1>3)? (B1/nB1) : (Double.MAX_VALUE);
		double avgB2   = (nB2>3)? (B2/nB2) : (Double.MAX_VALUE);
		double avgB3   = (nB3>3)? (B3/nB3) : (Double.MAX_VALUE);

		double C1 	= avgA0-avgB1; C1 = (C1>0)? C1 : 0;
		double C2 	= avgA0-avgB2; C2 = (C2>0)? C2 : 0;
		double C3 	= avgA0-avgB3; C3 = (C3>0)? C3 : 0;

		double G11 	= avgA1-avgB1; G11 = (G11>0)? G11 : 0;
		double G12 	= avgA1-avgB3; G12 = (G12>0)? G12 : 0;

		double G21 	= avgA2-avgB2; G21 = (G21>0)? G21 : 0;
		double G22 	= avgA2-avgB1; G22 = (G22>0)? G22 : 0;

		double G31 	= avgA3-avgB3; G31 = (G31>0)? G31 : 0;
		double G32 	= avgA3-avgB2; G32 = (G32>0)? G32 : 0;

		double E12	= avgA1-avgA2;
		double E23	= avgA2-avgA3;
		double E31	= avgA3-avgA1;

		return 	(1-Math.exp(-C1/D)) *
				(1-Math.exp(-C2/D)) *
				(1-Math.exp(-C3/D)) *
				(1-Math.exp(-G11/D)) *
				(1-Math.exp(-G12/D)) *
				(1-Math.exp(-G21/D)) *
				(1-Math.exp(-G22/D)) *
				(1-Math.exp(-G22/D)) *
				(1-Math.exp(-G31/D)) *
				(1-Math.exp(-G32/D)) *
				Math.exp(-(E12*E12)/(2*E*E)) *
				Math.exp(-(E23*E23)/(2*E*E)) *
				Math.exp(-(E31*E31)/(2*E*E))
				;

	}

    public ImageStack bifurcationess(
            FloatProcessor  inip,
            double          D
    )
    {
        FloatProcessor score0 = new FloatProcessor(inip.getWidth(), inip.getHeight());

        // fill the values in
        for (int x = 0; x<inip.getWidth(); x++) {
            for (int y = 0; y<inip.getHeight(); y++) {

				double sc = bifurcationess(x, y, inip, D);

				score0.setf(x, y, (float)sc);

            }
        }

        ImageStack isOut = new ImageStack(inip.getWidth(), inip.getHeight());
        isOut.addSlice(score0);
        return isOut;

    }

	public ImageStack exportTemplate(
	    double[] 	directionAngles
	)
	{

		return plotTemplate(calculateTemplate(directionAngles));

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

	private double 	point2dir(
									  double nx, 		// direction vector
									  double ny,
									  double px,        // point considered
									  double py
	)
	{
		// line is in vector from b in n direction

		float d = 0;

		double[] p_b = new double[2];

		// p - b
		p_b[0] = px;// - b[0];
		p_b[1] = py;// - b[1];

		double proj = p_b[0] * nx + p_b[1] * ny;

		if(proj<0){
			// "behind" the orientation
			return Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
		}

		// || (p-b) - dot(p-b,n) * n ||
		p_b[0] = p_b[0] - proj * nx;
		p_b[1] = p_b[1] - proj * ny;

//		distance_from_line = vectorNorm(distance_2d);

		return Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
	}

	private double[] runMS(
								  double[] 	start,
								  float[] 	inputProfile,
								  int 		max_iter,
								  double 	epsilon,
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

	private Vector<double[]> extractClusters(
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


	private ArrayList<ArrayList<int[]>> calculateTemplate( // exportTemplate() uses
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

				double d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

				if ( d2 <= (double)rInner*(double)rInner ) {

					p[0] = x-xc;
					p[1] = -(y-yc);

					templateCenter.add(new int[]{p[0], p[1]});

				}
			}
		}
		template.add(templateCenter);  // template(0) added

        ArrayList<int[]> templateBtwAngle12 = new ArrayList<int[]>();
        ArrayList<int[]> templateBtwAngle23 = new ArrayList<int[]>();
        ArrayList<int[]> templateBtwAngle31 = new ArrayList<int[]>();

        ArrayList<int[]> AtAngle1 = new ArrayList<int[]>();
        ArrayList<int[]> AtAngle2 = new ArrayList<int[]>();
        ArrayList<int[]> AtAngle3 = new ArrayList<int[]>();

		for (int pIdx = 0; pIdx < 3; pIdx++) {

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					double d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (d2 <= (double)r*(double)r && d2 >= (double)rLower*(double)rLower) {

						double px = x-xc;
						double py = -(y-yc);

						double nx = (float) Math.cos(ap[pIdx]);
						double ny = (float) Math.sin(ap[pIdx]);

						//float dst = point2dir(n, p);
                        double dst = point2dir(nx, ny, px, py);

                        if (dst<=(double)diam) {

                            if (dst<=(double)diam/2) { // belongs to pIdx ON peak and not filled  // && idxMapLocal[x+d*y]==0

                                if (pIdx==0) {
                                    AtAngle1.add(new int[]{x-xc, -(y-yc)});
                                }
                                else if (pIdx==1) {
                                    AtAngle2.add(new int[]{x-xc, -(y-yc)});
                                }
                                else if (pIdx==2) {
                                    AtAngle3.add(new int[]{x-xc, -(y-yc)});
                                }

                            }
                            else {

                                if (pIdx==0) {
                                    templateBtwAngle12.add(new int[]{x-xc, -(y-yc)});
                                }
                                else if (pIdx==1) {
                                    templateBtwAngle23.add(new int[]{x-xc, -(y-yc)});
                                }
                                else if (pIdx==2) {
                                    templateBtwAngle31.add(new int[]{x-xc, -(y-yc)});
                                }


                            }

                        }

					}
				}
			}

		}

        //IJ.log("---> "+AtAngle1.size()+" direction 1");

		template.add(AtAngle1);
		template.add(AtAngle2);
		template.add(AtAngle3);

        template.add(templateBtwAngle12);
        template.add(templateBtwAngle23);
        template.add(templateBtwAngle31);

		return template;
	}

	private ImageStack plotTemplate(ArrayList<ArrayList<int[]>> template) // exportTemplate() uses
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
				ip.setf(d/2+offX, d/2+offY, +128);
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

}
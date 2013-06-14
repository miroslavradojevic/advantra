import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

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

		}

		System.out.println("formed "+offsets.size()+"   , "+offsets.get(0).size());

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

    public ImageStack plotSums(int atX, int atY, FloatProcessor inip) {


        // calculate sums
        float[] sumsPerOrt  = new float[offsets.size()];
        float[] angles      = new float[offsets.size()];

        for (int dirIdx = 0; dirIdx<offsets.size(); dirIdx++) {

            angles[dirIdx] = dirIdx*resolDeg;

            for (int locIdx = 0; locIdx<offsets.get(dirIdx).size(); locIdx++) {

                sumsPerOrt[dirIdx] += inip.getf(
                        atX+offsets.get(dirIdx).get(locIdx)[0], atY+offsets.get(dirIdx).get(locIdx)[1]
                );

            }

        }

        Plot p;
        ImageStack isOut;

        p = new Plot("", "", "", angles, sumsPerOrt);
        isOut = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
        isOut.addSlice(p.getProcessor());

        int res = 1000;
        float[] angles2             = new float[res];
        float[] sumsPerOrt2         = new float[res];

        for (int i=0; i<res; i++) {
            angles2[i] = i*((float) sumsPerOrt.length/res)-0.5f;//((float)i/res)*360;
            sumsPerOrt2[i] = (float) interp1Darray(angles2[i], sumsPerOrt);
        }

        p = new Plot("", "higher res", "", angles2, sumsPerOrt2);
        //isOut.addSlice(p.getProcessor());

        // find peaks
        double[] start = new double[sumsPerOrt.length];
        for (int i=0; i<start.length; i++) start[i] = i;
        double[] finish = runMS(sumsPerOrt, 100, 0, 3);
        double[] plotY = new double[start.length];
        for (int i=0; i<plotY.length; i++) plotY[i]=sumsPerOrt[0]+50;
        //for (int i=0; i<plotY.length; i++) plotY[i]=sumsPerOrt[1];

        //p.addPoints(start, plotY, Plot.BOX);// = new Plot("", "start", "", start, plotY);
        //isOut.addSlice(p.getProcessor());

        p.addPoints(finish, plotY, Plot.BOX);// = new Plot("", "finish", "", finish, plotY);
        isOut.addSlice(p.getProcessor());

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

    private double 	runOne(double curr_pos, int h, float[] inputProfile){
        double 	    new_pos     = 0;
        double 		sum 		= 0;

        for (int x = -h; x <= h; x++) {
            if ((  curr_pos + x >=0) && ( curr_pos + x <= inputProfile.length-1) ){
                    if (x * x <= h * h) {
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

    public double[] 	runMS(float[] inputProfile, int max_iter, double epsilon, int h){

        double[] T				= new double[inputProfile.length];

        for (int i = 0; i < T.length; i++) {
            T[i] = i;//	T[i][1] = S[i][1];
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

    private double interp2Darray(double x, double y, float[] inarray, int wth) {// not tested

        // (x, y) is a real position in an 2d matrix
        // y = [0.0, 	H-1.0)
        // x = [-0.5, 	W-0.5)
        // wth is array width - columns will be wrapped - first corresponds to the last one

        // find four surrounding points
        int[] p11 = {(int)Math.floor(y),	(int)Math.floor(x)};
        int[] p12 = {(int)Math.floor(y), 	(int)Math.ceil(x)};
        int[] p21 = {(int)Math.ceil(y), 	(int)Math.floor(x)};
        int[] p22 = {(int)Math.ceil(y), 	(int)Math.ceil(x)};
        // bilinear coeffs a,b
        double a = ((p12[1]-p11[1])>0)?(x-p11[1])/(p12[1]-p11[1]) : 0.5;
        double b = ((p21[0]-p11[0])>0)?(y-p11[0])/(p21[0]-p11[0]) : 0.5;
        // wrap cols of surrounding pts
        p11[1] = (p11[1]<-0.5)? p11[1]+wth : ((p11[1]>=wth-0.5)?(p11[1]-wth) : p11[1]);
        p12[1] = (p12[1]<-0.5)? p12[1]+wth : ((p12[1]>=wth-0.5)?(p12[1]-wth) : p12[1]);
        p21[1] = (p21[1]<-0.5)? p21[1]+wth : ((p21[1]>=wth-0.5)?(p21[1]-wth) : p21[1]);
        p22[1] = (p22[1]<-0.5)? p22[1]+wth : ((p22[1]>=wth-0.5)?(p22[1]-wth) : p22[1]);

        double I11 = inarray[p11[0]*wth+p11[1]];
        double I12 = inarray[p12[0]*wth+p12[1]];
        double I21 = inarray[p21[0]*wth+p21[1]];
        double I22 = inarray[p22[0]*wth+p22[1]];

        return (1-a)*(1-b)*I11 + (1-a)*b*I21 + a*(1-b)*I12 + a*b*I22;
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

}

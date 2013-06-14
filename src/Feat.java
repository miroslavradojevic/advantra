import java.util.ArrayList;
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

					int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (d2 <= r*r && d2 >= rLower*rLower) {

						p[0] = x-xc;
						p[1] = -(y-yc);

						float ang = aRad;  		// read the angle

						n[0] = (float) Math.cos(ang);
						n[1] = (float) Math.sin(ang);

						float dst = point2dir(cent, n, p);

						if (dst<=diam/2) { // belongs to pIdx ON peak and not filled

							offsetsAngle.add(new int[]{p[0], p[1]});

						}

					}
				}
			}

			offsets.add(offsetsAngle);

		}

		System.out.println("formed "+offsets.size()+"   , "+offsets.get(0).size());

		String[] s_all = new String[]{"'ro'", "'go'", "'bo'", "'yo'"};
		for (int c = 0; c<offsets.size(); c+=4) {    //

			System.out.println("a = [...");

			for (int c1 = 0; c1<offsets.get(c).size(); c1++) {
				System.out.println(""+offsets.get(c).get(c1)[0]+" , "+offsets.get(c).get(c1)[1]+";...");
			}

			System.out.println("];");
			String s =  s_all[new Random().nextInt(3)];
			System.out.println("plot(a(:,1), a(:,2), "+s+"); axis equal; grid on; hold on;");
		}

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


}

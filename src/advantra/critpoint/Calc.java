package advantra.critpoint;

import advantra.feature.FilterSet;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/20/13
 * Time: 6:10 PM
 * To change this template use File | Settings | File Templates.
 */

public class Calc {

	public static int circularProfileSize(
												 int radius
	)
	{

		int cnt = 0;

		for (int x = -radius; x <= radius; x++){
			for (int y = -radius; y <= radius; y++){
				if (x*x+y*y<=radius*radius)	cnt++;
			}
		}

		return cnt;

	}

	public static void getProfileLocations(
												  int xin,
												  int yin,
												  int rin,
												  int[] xloc,
												  int[] yloc
	)
	{

		int xs = xin-rin;
		int xe = xin+rin;

		int ys = yin-rin;
		int ye = yin+rin;

		int cnt = 0;
		for (int x = xs; x <= xe; x++){
			for (int y = ys; y <= ye; y++){
				int d = (x-xin)*(x-xin)+(y-yin)*(y-yin);
				if(d <= rin*rin){
					xloc[cnt] = x;
					yloc[cnt] = y;
					cnt++;
				}
			}
		}

	}

	public static void getProfile(
								   ImagePlus inimg,
								   int xin,
								   int yin,
								   int rin,
								   float[] vals,
								   float[] angs,
								   float[] rads
	)
	{

		// will extract out vals[], rads[], and angs[]
		// those will be used to score() on filterSet

		int xs = xin-rin;
		int xe = xin+rin;

		int ys = yin-rin;
		int ye = yin+rin;

		int cnt = 0;
		for (int x = xs; x <= xe; x++){
			for (int y = ys; y <= ye; y++){

				int d = (x-xin)*(x-xin)+(y-yin)*(y-yin);

				if(d <= rin*rin){

					vals[cnt] = inimg.getProcessor().getPixelValue(x, y);
					rads[cnt] = (float) (Math.sqrt(d) / rin);
					angs[cnt] = (float) (Math.atan2(y-yin, x-xin) + Math.PI);
					cnt++;

				}

			}
		}
	}

	public static ImageStack plotProfile(
            float[] vals,
            float[] angs,
            float[] rads)

	{


		ImageStack viz = new ImageStack(400, 200);

		// find max for plotting
		float max_val = vals[0];
		float min_val = vals[0];

		float max_rad = rads[0]; // this is not necessary, it is 1
		float min_rad = rads[0]; // this is not necessary, it is 0

		float max_ang = angs[0];
		float min_ang = angs[0];

		for (int i = 1; i < vals.length; i++){



			if (vals[i] > max_val) max_val = vals[i];
			if (vals[i] < min_val) min_val = vals[i];

			if (rads[i] > max_rad) max_rad = rads[i];
			if (rads[i] < min_rad) min_rad = rads[i];

			if (angs[i] > max_ang) max_ang = angs[i];
			if (angs[i] < min_ang) min_ang = angs[i];

		}

//		System.out.print(min_val+" -> "+max_val+"\n");
//		System.out.print(min_rad+" -> "+max_rad+"\n");
//		System.out.print(min_ang+" -> "+max_ang+"\n");

		Plot p = new Plot("circular_profile", "angle[rad]", "value");
		p.setLimits(min_ang, max_ang, min_val, max_val);
		p.setSize(400, 200);
		p.addPoints(angs, vals, Plot.BOX);
		viz.addSlice("circular_profile", p.getProcessor());

		Plot p1 = new Plot("radial_profile", "radius", "value");
		p1.setLimits(min_rad, max_rad, min_val, max_val);
		p1.setSize(400, 200);
		p1.addPoints(rads, vals, Plot.CIRCLE);
		viz.addSlice("radial_profile", p1.getProcessor());

		return viz;
	}

	public static ImageProcessor getProfilePatch(
														ImagePlus inimg,
														int xin,
														int yin,
														int rin
	){
		ImageProcessor ip = new FloatProcessor(2*rin+1, 2*rin+1);

		for (int x = 0; x < 2*rin+1; x++){
			for (int y = 0; y < 2*rin+1; y++){
				ip.setf(x, y, inimg.getProcessor().getPixelValue(x+xin-rin, y+yin-rin));
			}
		}

		return ip;
	}

	public static ImageProcessor plotAllResponses(
														 FilterSet fs,
														 float[] vals,
														 float[] angs,
														 float[] rads
	)
	{

		float[] c_idx = new float[fs.circConfs.size()+fs.radlConfs.size()];
		float[] c_sco = new float[fs.circConfs.size()+fs.radlConfs.size()];
		for (int i = 0; i < fs.circConfs.size(); i++){
			c_idx[i] = i;
            System.out.println("<OK?!>, score at circ conf. "+i);
			c_sco[i] = fs.circConfs.get(i).score(vals, angs, rads);
		}
		for (int i = fs.circConfs.size(); i < fs.circConfs.size()+fs.radlConfs.size(); i++){
			c_idx[i] = i;
			c_sco[i] = fs.radlConfs.get(i-fs.circConfs.size()).score(vals, rads);
		}

		Plot p = new Plot("filtering_score", "filter_idx", "filter_score", c_idx, c_sco);
		p.setSize(400, 200);
		return p.getProcessor();

	}

	public static ImageProcessor filterResponse(
												  ImagePlus template,
												  FilterSet fs,
												  int fsIdx,
												  int margin
	)
	{
		int W = template.getWidth();
		int H = template.getHeight();

		// storage for the profile
		int N = circularProfileSize(margin);
		float[] vals = new float[N];
		float[] angs = new float[N];
		float[] rads = new float[N];

		ImageProcessor ip = new FloatProcessor(W, H);
		for (int x = 0; x < W; x++){
			for (int y = 0; y < H; y++){
				if (x>margin && x<W-margin && y>margin && y<H-margin){

					getProfile(template, x, y, margin, vals, angs, rads);

					if (fsIdx<fs.circConfs.size()) {
						ip.setf(x, y, fs.circConfs.get(fsIdx).score(vals, angs, rads));
					}
					else {
						ip.setf(x, y, fs.radlConfs.get(fsIdx-fs.circConfs.size()).score(vals, rads));
					}

				}
			}
		}
		return ip;
	}

}

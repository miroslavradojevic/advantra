package advantra.feature;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ByteProcessor;

import java.awt.Color;
import java.util.ArrayList;

public class CircHaarFeat {

	private 	int 				M; 		// number of angles in 2*pi
	private		ArrayList<float[]>  feats;	// collection of features, each M elements
	private 	float				featsMax;
	private 	float				featsMin;
	
	public CircHaarFeat(int M){
		
		// M describes the length of the angular signature
		this.M = M;
		feats = new ArrayList<float[]>();
		
	}
	
	public void 		createFeatures(){
		
		System.out.println("Creating features with "+M+" angular responses per 360 degrees...");
		
		int scale_limit = 1;
		
		int cnt;
		
		System.out.print("creating 1xON feature...");
		
		cnt = 0;
		
		// consider only part of the scales range
		for (int s = 1; s <= scale_limit; s++) { 
			
			float[] vec = new float[M];
			
			for (int i = 0; i < vec.length; i++) {
				if(i<s){
					vec[i] = (float)(M-s)/(float)(s); // ON
				}
				else{
					vec[i] = -1; // OFF
				}
			}
			
			feats.add(vec);
			cnt++;
			
		}
		
		System.out.println(cnt+" created");
		
		System.out.print("creating 2xON feature...");
		
		cnt = 0;
		
		for (int s1 = 1; s1 <= scale_limit; s1++) { 
			for (int sep = 1; sep <= M; sep++) {
				for (int s2 = 1; s2 <= scale_limit; s2++) {
					
					if(s1+s2+sep<=M){
						
						float[] vec = new float[M];
						for (int i = 0; i < vec.length; i++) {
							if(i<s1){
								vec[i] = (float)(M-(s1+s2))/(float)(s1+s2); // ON
							}
							else if(i>=s1 && i<s1+sep){
								vec[i] = -1; // OFF
							}
							else if(i>=s1+sep && i<s1+sep+s2){
								vec[i] = (float)(M-(s1+s2))/(float)(s1+s2); // ON
							}
							else{
								vec[i] = -1;///(float)(M-(s1+s2+s3)); // OFF
							}
						}
						
						feats.add(vec);
						cnt++;
					}
				}
			}
		}
		
		
		
		System.out.println(cnt+" created");
		
		System.out.print("creating 3xON feature...");
		
		cnt = 0;
		
		for (int s1 = 1; s1 <= scale_limit; s1++) { // M-5 is the maximum width
			for (int sep1 = 1; sep1 <= M; sep1++) {
				for (int s2 = 1; s2 <= scale_limit; s2++) {
					for (int sep2 = 1; sep2 <= M; sep2++) {
						for (int s3 = 1; s3 <= scale_limit; s3++) {
							
							if(s1+s2+s3+sep1+sep2<=M){ // there's enough space to add it
						
								float[] vec = new float[M];
								
								for (int i = 0; i < vec.length; i++) {
									if(i<s1){
										vec[i] = (float)(M-(s1+s2+s3))/(float)(s1+s2+s3); // ON
									}
									else if(i>=s1 && i<s1+sep1){
										vec[i] = -1; // OFF
									}
									else if(i>=s1+sep1 && i<s1+sep1+s2){
										vec[i] = (float)(M-(s1+s2+s3))/(float)(s1+s2+s3); // ON
									}
									else if(i>=s1+sep1+s2 && i<s1+sep1+s2+sep2){
										vec[i] = -1; // OFF
									}
									else if(i>=s1+sep1+s2+sep2 && i<s1+sep1+s2+sep2+s3){
										vec[i] = (float)(M-(s1+s2+s3))/(float)(s1+s2+s3); // ON
									}
									else{
										vec[i] = -1; // OFF
									}
								}
						
								feats.add(vec);
								cnt++;
							}
						}
					}
				}
			}
		}
		
		System.out.println(cnt+" created");
		
		System.out.print("creating 4xON feature...");
		cnt = 0;
		for (int s1 = 1; s1 <= scale_limit; s1++) {
			for (int sep1 = 1; sep1 <= M; sep1++) {
				for (int s2 = 1; s2 <= scale_limit; s2++) {
					for (int sep2 = 1; sep2 <= M; sep2++) {
						for (int s3 = 1; s3 <= scale_limit; s3++) {
							for (int sep3 = 1; sep3 <= M; sep3++) {
								for (int s4 = 1; s4 <= scale_limit; s4++) {

									if(s1+s2+s3+s4+sep1+sep2+sep3<=M){
						
										float[] vec = new float[M];
										for (int i = 0; i < vec.length; i++) {
											if(i<s1){
												vec[i] = (float)(M-(s1+s2+s3+s4))/(float)(s1+s2+s3+s4); // ON
											}
											else if(i>=s1 && i<s1+sep1){
												vec[i] = -1; // OFF
											}
											else if(i>=s1+sep1 && i<s1+sep1+s2){
												vec[i] = (float)(M-(s1+s2+s3+s4))/(float)(s1+s2+s3+s4); // ON
											}
											else if(i>=s1+sep1+s2 && i<s1+sep1+s2+sep2){
												vec[i] = -1; // OFF
											}
											else if(i>=s1+sep1+s2+sep2 && i<s1+sep1+s2+sep2+s3){
												vec[i] = (float)(M-(s1+s2+s3+s4))/(float)(s1+s2+s3+s4); // ON
											}
											else if(i>=s1+sep1+s2+sep2+s3 && i<s1+sep1+s2+sep2+s3+sep3){
												vec[i] = -1; // OFF
											}
											else if(i>=s1+sep1+s2+sep2+s3+sep3 && i<s1+sep1+s2+sep2+s3+sep3+s4){
												vec[i] = (float)(M-(s1+s2+s3+s4))/(float)(s1+s2+s3+s4); // ON
											}
											else{
												vec[i] = -1; // OFF
											}
										}
						
										feats.add(vec);
										cnt++;
									}
								}
							}
						}
					}
				}
			}
		}
		
		System.out.println(cnt+" created");
		
		System.out.println("total " + feats.size() + " features created.");
		
		// max
		featsMin = Float.MAX_VALUE;
		featsMax = Float.MIN_VALUE;
		for (int i = 0; i < feats.size(); i++) {
			for (int j = 0; j < feats.get(i).length; j++) {
				if(feats.get(i)[j]>featsMax) featsMax = feats.get(i)[j];
				if(feats.get(i)[j]<featsMin) featsMin = feats.get(i)[j];
			}
		}
		
		System.out.println("coeffs range from "+featsMin+" to "+featsMax);
		
	}
	
	public int 			nrFeats(){
		
		return feats.size();
	
	}
	
	public void 			showFeature(int index){
		
		byte[] fill 	= new byte[M*M];
		float[] shifted = new float[M];
		ImageStack out_stk =  new ImageStack(M, M);
		
		int cnt = 0;
		for (int f1 = 0; f1 < M; f1++) {
			for (int f2 = 0; f2 < M; f2++) {
				
				int shift_index = f1;
				
				if(shift_index<1){
					fill[cnt] = (feats.get(index)[f2]>0)?(byte)255:(byte)0; 
					cnt++;
				}
				else{
					shifted = shiftVals(feats.get(index), shift_index);
					fill[cnt] = (shifted[f2]>0)?(byte)255:(byte)0;
					cnt++;
				}
				
			}
		}
		
		ByteProcessor ip = new ByteProcessor(M, M, fill);
		out_stk.addSlice("feat." + index, ip);
		
		new ImagePlus("feature", out_stk).show();
		
		// plot
		float[] x_axis = new float[M];
		for (int i = 0; i < M; i++) {
			x_axis[i] = i * (2*(float)Math.PI/(float)M);
		}
				
		ImageStack plot_stack = new ImageStack(600, 300);
			Plot p = new Plot("feat"+index, "angle", "filter");
			p.setSize(600, 300);
			p.setLimits(0, 2*Math.PI, featsMin, featsMax);
			p.setColor(Color.RED);
			p.addPoints(x_axis, feats.get( index ), Plot.BOX);
			plot_stack.addSlice("feat"+index, p.getProcessor());
		new ImagePlus("features", plot_stack).show();
		
	}
	
	public void 			showFeatures(int[] indexes){
		
		ImageStack out_stk =  new ImageStack(M, M);
		float[] shifted = new float[M]; // just to visualize
		
		for (int i = 0; i < indexes.length; i++) {
			
			byte[] fill = new byte[M*M];
			
			int cnt = 0;
			
			for (int f1 = 0; f1 < M; f1++) {
				for (int f2 = 0; f2 < M; f2++) {
					
					int shift_index = f1;
					
					if(shift_index<1){
						fill[cnt] = (feats.get(  indexes[i]  )[f2]>0)?(byte)255:(byte)0; 
						cnt++;
					}
					else{
						shifted = shiftVals(feats.get(  indexes[i]  ), shift_index);
						fill[cnt] = (shifted[f2]>0)?(byte)255:(byte)0;
						cnt++;
					}
					
				}
			}
			
			ByteProcessor ip = new ByteProcessor(M, M, fill);
			out_stk.addSlice("feat" + i, ip);
			
		}
		
		new ImagePlus("selected_features", out_stk).show();
		
		// plot each feature
		float[] x_axis = new float[M];
		for (int i = 0; i < M; i++) {
			x_axis[i] = i * (2*(float)Math.PI/(float)M);
		}
				
		ImageStack plot_stack = new ImageStack(600, 300);
		for (int i = 0; i < indexes.length; i++) {
			Plot p = new Plot("feat"+indexes[i] , "angle", "filter");
			p.setSize(600, 300);
			p.setLimits(0, 2*Math.PI, featsMin, featsMax);
			p.setColor(Color.RED);
			p.addPoints(x_axis, feats.get( indexes[i] ), Plot.BOX);
			plot_stack.addSlice("feat"+indexes[i] , p.getProcessor());
		}
				
		new ImagePlus("features", plot_stack).show();
		
	}
	
	public void 			showFeatures(){
		
		ImageStack out_stk =  new ImageStack(M, M);
		float[] shifted = new float[M];
		
		for (int i = 0; i < feats.size(); i++) {
			
			byte[] fill = new byte[M*M];
			
			int cnt = 0;
			
			for (int f1 = 0; f1 < M; f1++) {
				for (int f2 = 0; f2 < M; f2++) {
					
					int shift_index = f1;
					
					if(shift_index<1){
						fill[cnt] = (feats.get(i)[f2]>0)?(byte)255:(byte)0; 
						cnt++;
					}
					else{
						shifted = shiftVals(feats.get(i), shift_index);
						fill[cnt] = (shifted[f2]>0)?(byte)255:(byte)0;
						cnt++;
					}
					
				}
			}
			
			ByteProcessor ip = new ByteProcessor(M, M, fill);
			out_stk.addSlice("feat" + i, ip);
			
		}
		
		new ImagePlus("angular_features", out_stk).show();
		
		// plot each feature
		float[] x_axis = new float[M];
		for (int i = 0; i < M; i++) {
			x_axis[i] = i * (2*(float)Math.PI/(float)M);
		}
		
		ImageStack plot_stack = new ImageStack(600, 300);
		for (int i = 0; i < feats.size(); i++) {
			Plot p = new Plot("feat"+i, "angle", "filter");
			p.setSize(600, 300);
			p.setLimits(0, 2*Math.PI, featsMin, featsMax);
			p.setColor(Color.RED);
			p.addPoints(x_axis, feats.get(i), Plot.BOX);
			
			plot_stack.addSlice("feat"+i, p.getProcessor());
		}
		
		new ImagePlus("features", plot_stack).show();

	}
	
	public float[] 			allFeatScore(double[] profile, int shiftStep){
		
		float[] out = new float[feats.size()];
		
		for (int i = 0; i < feats.size(); i++) {
			out[i] = featScore(profile, feats.get(i), shiftStep);
		}
		
		return out;
		
	}
	
	private float[] shiftVals(float[] in, int shift){
		
		if(shift==0){
			return in;
		}
		
		float[] out = new float[in.length];
		float[] tmp = new float[shift];
		
		// store last shift values, starting from the last
		for (int i = 0; i < shift; i++) {
			tmp[i] = in[in.length-1-i];
		}
		
		for (int i = in.length-1; i >= shift; i--) {
			out[i] = in[i-shift];
		}
		
		int cnt = 0;
		for (int i = shift-1; i >= 0; i--) {
			out[i] = tmp[cnt]; 
			cnt++;
		}
		
		return out;
		
	}
	
	private float conv(double[] a, float[] b){
		
		float c = 0;
		
		for (int i = 0; i < b.length; i++) {
			c += a[i]*b[i];
		}
		
		return c;
	}
	
	private float featScore(double[] profile, float[] feature, int shiftStep){
		
		float total_score 	= Float.MIN_VALUE;
		float curr_score 	= 0;
		
		for (int d = 0; d < feature.length; d+=shiftStep) {

			if(d==0){
				curr_score = conv(profile, feature);
			}
			else{
				curr_score = conv(profile, shiftVals(feature, d));
			}
			
			if(curr_score>total_score){
				total_score = curr_score;
			}
			
		}
		
		return total_score;
		
	}
	
}

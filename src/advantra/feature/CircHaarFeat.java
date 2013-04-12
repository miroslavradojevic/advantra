package advantra.feature;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;

import java.util.ArrayList;

public class CircHaarFeat {

	int M;
	ArrayList<float[]>  feats;
	
	public CircHaarFeat(int M){
		
		// M describes the length of the angular signature
		this.M = M;
		feats = new ArrayList<float[]>();
		
	}
	
	public void createFeatures(){
		
		System.out.println("Creating features with "+M+" angular responses...");
		
		int cnt = 0;
		
		System.out.print("ON x1 feature...");
		cnt = 0;
		for (int s = 1; s <= M-1; s++) { // M-1 is the maximum width
			
			float[] vec = new float[M];
			for (int i = 0; i < vec.length; i++) {
				if(i<s){
					vec[i] = 1/(float)s; // ON
				}
				else{
					vec[i] = -1/(float)(M-s); // OFF
				}
			}
			
			feats.add(vec);
			cnt++;
			
		}
		
		System.out.println(cnt+" created");
		
		System.out.print("ON x2 feature...");
		cnt = 0;
		for (int s1 = 1; s1 <= M-3; s1++) { // M-3 is the maximum width
			for (int s2 = 1; s2 <= M-3; s2++) {
				for (int s3 = 1; s3 <= M-3; s3++) {
					
					if(s1+s2+s3<=M-1){
						
						float[] vec = new float[M];
						for (int i = 0; i < vec.length; i++) {
							if(i<s1){
								vec[i] = 1/(float)s1; // ON
							}
							else if(i>=s1 && i<s1+s2){
								vec[i] = -1/(float)s2; // OFF
							}
							else if(i>=s1+s2 && i<s1+s2+s3){
								vec[i] = 1/(float)s3; // ON
							}
							else{
								vec[i] = -1/(float)(M-(s1+s2+s3)); // OFF
							}
						}
						
						feats.add(vec);
						cnt++;
					}
				}
			}
		}
		
		System.out.println(cnt+" created");
		
		System.out.print("ON x3 feature...");
		cnt = 0;
		for (int s1 = 1; s1 <= M-5; s1++) { // M-5 is the maximum width
			for (int s2 = 1; s2 <= M-5; s2++) {
				for (int s3 = 1; s3 <= M-5; s3++) {
					for (int s4 = 1; s4 <= M-5; s4++) {
						for (int s5 = 1; s5 <= M-5; s5++) {
							
							if(s1+s2+s3+s4+s5<=M-1){ // there's enough space to add it
						
								float[] vec = new float[M];
								
								for (int i = 0; i < vec.length; i++) {
									if(i<s1){
										vec[i] = 1/(float)s1; // ON
									}
									else if(i>=s1 && i<s1+s2){
										vec[i] = -1/(float)s2; // OFF
									}
									else if(i>=s1+s2 && i<s1+s2+s3){
										vec[i] = 1/(float)s3; // ON
									}
									else if(i>=s1+s2+s3 && i<s1+s2+s3+s4){
										vec[i] = -1/(float)s4; // OFF
									}
									else if(i>=s1+s2+s3+s4 && i<s1+s2+s3+s4+s5){
										vec[i] = 1/(float)s5; // ON
									}
									else{
										vec[i] = -1/(float)(M-(s1+s2+s3+s4+s5)); // OFF
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
		
		System.out.print("ON x4 feature...");
		cnt = 0;
		for (int s1 = 1; s1 <= M-7; s1++) {
			for (int s2 = 1; s2 <= M-7; s2++) {
				for (int s3 = 1; s3 <= M-7; s3++) {
					for (int s4 = 1; s3 <= M-7; s3++) {
						for (int s5 = 1; s3 <= M-7; s3++) {
							for (int s6 = 1; s6 <= M-7; s6++) {
								for (int s7 = 1; s7 <= M-7; s7++) {

									if(s1+s2+s3+s4+s5+s6+s7<=M-1){
						
										float[] vec = new float[M];
										for (int i = 0; i < vec.length; i++) {
											if(i<s1){
												vec[i] = 1/(float)s1; // ON
											}
											else if(i>=s1 && i<s1+s2){
												vec[i] = -1/(float)s2; // OFF
											}
											else if(i>=s1+s2 && i<s1+s2+s3){
												vec[i] = 1/(float)s3; // ON
											}
											else if(i>=s1+s2+s3 && i<s1+s2+s3+s4){
												vec[i] = -1/(float)s4; // OFF
											}
											else if(i>=s1+s2+s3+s4 && i<s1+s2+s3+s4+s5){
												vec[i] = 1/(float)s5; // ON
											}
											else if(i>=s1+s2+s3+s4+s5 && i<s1+s2+s3+s4+s5+s6){
												vec[i] = -1/(float)s6; // OFF
											}
											else if(i>=s1+s2+s3+s4+s5+s6 && i<s1+s2+s3+s4+s5+s6+s7){
												vec[i] = 1/(float)s7; // ON
											}
											else{
												vec[i] = -1/(float)(M-(s1+s2+s3+s4+s5+s6+s7)); // OFF
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
		
	}
	
	public ImagePlus showFeature(int index){
		
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
		out_stk.addSlice("feat.nr." + index, ip);
		
		return new ImagePlus("feature", out_stk);
		
	}
	
	public ImagePlus showFeatures(){
		
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
		
		return new ImagePlus("angular_features", out_stk);

	}
	
	private float[] shiftVals(float[] in, int shift){
		
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
	
}

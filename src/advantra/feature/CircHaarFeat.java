package advantra.feature;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Vector;

import advantra.general.ArrayHandling;

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
		
		System.out.println("ON/OFF x1 feature...");
		
		for (int s = 1; s <= M-1; s++) {
			
			float[] vec = new float[M];
			for (int i = 0; i < vec.length; i++) {
				if(i<s){
					vec[i] = 1/(float)s; // ON
				}
				else{
					vec[i] = -1/(float)(M-s); // OFF
				}
			}
			
//			System.out.println("add feature:");
//			System.out.print(""+s+" -> ");
//			for (int i = 0; i < vec.length; i++) {
//				System.out.format("%.2f , ", vec[i]);
//			}
			
			System.out.println();
			feats.add(vec);
			
		}
		
		System.out.println("ON/OFF x2 feature...");
		
		int count = 0;
		for (int s1 = 1; s1 <= M-1; s1++) {
			for (int s2 = 0; s2 < M-1; s2++) {
				for (int s3 = 0; s3 < M-1; s3++) {
					// s1 1st ON
					if(s1+s2+s3<=M-1){
						
						float[] vec = new float[M];
						for (int i = 0; i < vec.length; i++) {
							if(i<s){
								vec[i] = 1/(float)s; // ON
							}
							else{
								vec[i] = -1/(float)(M-s); // OFF
							}
						}
						
						feats.add(new float[M]);
						count++;
					}
				}
			}
		}
		
		System.out.println(feats.size()+" features created.");
		
	}
	
	public ImagePlus showFeatures(int plot_size){
		
		int block_size = plot_size/M;
		
		int H = block_size*M;
		int W = block_size*M;
		
		ImageStack out_stk =  new ImageStack(W, H);
		float[] shifted = new float[M];
		
		
		for (int i = 0; i < feats.size(); i++) {
			
			byte[] fill = new byte[H*W];
			
			int cnt = 0;
			for (int f1 = 0; f1 < H; f1++) {
				for (int f2 = 0; f2 < W; f2++) {
					
					// f1 is row f1/block_size is rotation index
					// f2 is col f2/block_size is feature vector index
					
					int shift_index = f1/block_size;
					
					if(shift_index<1){
						fill[cnt] = (feats.get(i)[f2/block_size]>0)?(byte)255:(byte)0; 
						cnt++;
					}
					else{
						shifted = shiftVals(feats.get(i), shift_index);
						fill[cnt] = (shifted[f2/block_size]>0)?(byte)255:(byte)0;
						cnt++;
					}
					
				}
			}
			
			ByteProcessor ip = new ByteProcessor(W, H, fill);
			out_stk.addSlice("feat" + i, ip);
			
		}
		
		ImagePlus out_img = new ImagePlus("angular_features", out_stk);
		return out_img;

	}
	
	
	private void shift(float[] in){
		
		float temp = in[in.length-1]; // last
		
		for (int i = 1; i < in.length; i++) {
			in[i] = in[i-1];
		}
		
		in[0] = temp;
		
	}
	
	private float[] shiftVals(float[] in, int shift){
		
		float[] out = new float[in.length];
		float[] tmp = new float[shift];
		
		// store last shift values, starting from the last
		for (int i = 0; i < shift; i++) {
			tmp[i] = in[in.length-1-i];
		}
		
		float temp = in[in.length-1]; // last
		
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

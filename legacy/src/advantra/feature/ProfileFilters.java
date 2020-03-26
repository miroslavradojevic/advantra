package advantra.feature;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;

import java.util.Vector;

public class ProfileFilters {
	
	/*
	 * feature is a vector corresponding to directional response
	 * the aim of the feature is to count the number of peaks, emulate the distribution
	 * of the directional profile (response)
	 */
	
	private Vector<float[][]> 	feat;
	private int					featLen;
	
	private int					angStepLen;
	private int					nrRots;
	
	private float[] 			anglesRad; // will be used to guide the construction of steps
	
	public ProfileFilters(int featLen, double angStep){
		
		this.feat 		= new Vector<float[][]>();
		this.featLen 	= featLen;
		
		this.angStepLen = (int)Math.ceil((angStep/(2*Math.PI))*featLen);
		this.nrRots 	= (featLen%angStepLen==0)? (featLen/angStepLen-1) : (featLen/angStepLen) ; 
		
		this.anglesRad = new float[featLen];
		for (int i = 0; i < featLen; i++) {
			this.anglesRad[i] = (float) (i*2*Math.PI/featLen);
		}
		
		System.out.println(""+featLen);
		System.out.println(""+angStepLen);
		System.out.println(""+nrRots);
	}
	
	public void 		create(){
		
		int 	sLim 	= 1;
		float 	sumNeg, sumPos;
		int 	cnt 	= 0;
		
		// 1x
		for (int s = 1; s <= sLim; s++) { 
			
			float[][] rotd = new float[nrRots+1][];
			float[] base = new float[featLen];
			rotd[0] = new float[featLen];
			
			sumNeg = sumPos = 0;
			
			for (int i = 0; i < featLen; i++){
				
				// inflection point
				int flex1 = s*angStepLen;
				
				if(i<flex1){
					base[i] = +1;
					sumPos ++;
				}
				else{
					base[i] = -1;
					sumNeg ++;
				}
				
			}
			
			//j==0
			for (int i = 0; i < featLen; i++){
				rotd[0][i] = (base[i]>0)? base[i]*((float)sumNeg/sumPos) : base[i] ;
			}
			
			// rotated (for orientation invariance)
			for (int j = 1; j <= nrRots; j++) {
				rotd[j] = shiftVals(rotd[0], j*angStepLen);
			}
			
			feat.add(rotd);
			cnt++;
			
		}
		
		//2x
		for (int s1 = 1; s1 <= sLim; s1++) { 
			for (int s2 = 1; s2 <= sLim; s2++) {
				for (int w1 = 1; w1 <= nrRots; w1++) {
					
					int flex1 = s1*angStepLen;
					int flex2 = flex1 + w1*angStepLen;
					int flex3 = flex2 + s2*angStepLen;
					
					if(flex3<featLen){
						
						float[][] rotd = new float[nrRots+1][];
						float[] base = new float[featLen];
						rotd[0] = new float[featLen];
						
						sumNeg = sumPos = 0;
						
						for (int i = 0; i < featLen; i++) {
							if(i<flex1){
								base[i] = +1; 
								sumPos ++;
							}
							else if(i>=flex1 && i<flex2){
								base[i] = -1; 
								sumNeg ++;
							}
							else if(i>=flex2 && i<flex3){
								base[i] = +1; 
								sumPos ++;
							}
							else{
								base[i] = -1; 
								sumNeg ++;
							}
						}
						
						//j==0
						for (int i = 0; i < featLen; i++){
							rotd[0][i] = (base[i]>0)? base[i]*((float)sumNeg/sumPos) : base[i] ;
						}
						
						// rotated (for orientation invariance)
						for (int j = 1; j <= nrRots; j++) {
							rotd[j] = shiftVals(rotd[0], j*angStepLen);
						}
						
						feat.add(rotd);
						cnt++;
						
					}
				}
			}
		}
		
		//3x
		for (int s1 = 1; s1 <= sLim; s1++) { 
			for (int s2 = 1; s2 <= sLim; s2++) {
				for (int s3 = 1; s3 <= sLim; s3++) {
					for (int w1 = 1; w1 <= nrRots; w1++) {
						for (int w2 = 1; w2 <= nrRots; w2++) {
							
							int flex1 = s1*angStepLen;
							int flex2 = flex1 + w1*angStepLen;
							int flex3 = flex2 + s2*angStepLen;
							int flex4 = flex3 + w2*angStepLen;
							int flex5 = flex4 + s3*angStepLen;
							
							if(flex5<featLen){
								
								float[][] rotd = new float[nrRots+1][];
								float[] base = new float[featLen];
								rotd[0] = new float[featLen];
								
								sumNeg = sumPos = 0;
								
								for (int i = 0; i < featLen; i++) {
									if(i<flex1){
										base[i] = +1; 
										sumPos ++;
									}
									else if(i>=flex1 && i<flex2){
										base[i] = -1; 
										sumNeg ++;
									}
									else if(i>=flex2 && i<flex3){
										base[i] = +1; 
										sumPos ++;
									}
									else if(i>=flex3 && i<flex4){
										base[i] = -1; 
										sumNeg ++;
									}
									else if(i>=flex4 && i<flex5){
										base[i] = +1; 
										sumPos ++;
									}
									else{
										base[i] = -1; 
										sumNeg ++;
									}
								}
								
								//j==0
								for (int i = 0; i < featLen; i++){
									rotd[0][i] = (base[i]>0)? base[i]*((float)sumNeg/sumPos) : base[i] ;
								}
								
								// rotated (for orientation invariance)
								for (int j = 1; j <= nrRots; j++) {
									rotd[j] = shiftVals(rotd[0], j*angStepLen);
								}
								
								feat.add(rotd);
								cnt++;
								
							}
						}
					}
				}
			}
		}
		
		//4x
		for (int s1 = 1; s1 <= sLim; s1++) { 
			for (int s2 = 1; s2 <= sLim; s2++) {
				for (int s3 = 1; s3 <= sLim; s3++) {
					for (int s4 = 1; s4 <= sLim; s4++) {
						for (int w1 = 1; w1 <= nrRots; w1++) {
							for (int w2 = 1; w2 <= nrRots; w2++) {
								for (int w3 = 1; w3 <= nrRots; w3++) {
									
									int flex1 = s1*angStepLen;
									int flex2 = flex1 + w1*angStepLen;
									int flex3 = flex2 + s2*angStepLen;
									int flex4 = flex3 + w2*angStepLen;
									int flex5 = flex4 + s3*angStepLen;
									int flex6 = flex5 + w3*angStepLen;
									int flex7 = flex6 + s4*angStepLen;
									
									if(flex7<featLen){
										
										float[][] rotd = new float[nrRots+1][];
										float[] base = new float[featLen];
										rotd[0] = new float[featLen];
										
										sumNeg = sumPos = 0;
										
										for (int i = 0; i < featLen; i++) {
											if(i<flex1){
												base[i] = +1; 
												sumPos ++;
											}
											else if(i>=flex1 && i<flex2){
												base[i] = -1; 
												sumNeg ++;
											}
											else if(i>=flex2 && i<flex3){
												base[i] = +1; 
												sumPos ++;
											}
											else if(i>=flex3 && i<flex4){
												base[i] = -1; 
												sumNeg ++;
											}
											else if(i>=flex4 && i<flex5){
												base[i] = +1; 
												sumPos ++;
											}
											else if(i>=flex5 && i<flex6){
												base[i] = -1; 
												sumNeg ++;
											}
											else if(i>=flex6 && i<flex7){
												base[i] = +1; 
												sumPos ++;
											}
											else{
												base[i] = -1; 
												sumNeg ++;
											}
										}
										
										//j==0
										for (int i = 0; i < featLen; i++){
											rotd[0][i] = (base[i]>0)? base[i]*((float)sumNeg/sumPos) : base[i] ;
										}
										
										// rotated (for orientation invariance)
										for (int j = 1; j <= nrRots; j++) {
											rotd[j] = shiftVals(rotd[0], j*angStepLen);
										}
										
										feat.add(rotd);
										cnt++;
										
									}
								}
							}
						}
					}
				}
			}
		}		
		
		System.out.println("created total " + cnt + " feats");
		
	}
	
	public int 			nrFilters(){
		
		return feat.size();
	
	}
	
	public float[]		getBaseFilter(int filterIdx){
		return feat.get(filterIdx)[0];
	}
	
	public void 		showFilters(){
		
		ImageStack viz_stk 		=  	new ImageStack(featLen, nrRots+1);
		for (int i = 0; i < feat.size(); i++) {
			
			FloatProcessor fp = new FloatProcessor(featLen, nrRots+1);
			
			float[][] current_filt = feat.get(i);
			for (int k = 0; k < current_filt.length; k++) {
				for (int l = 0; l < current_filt[0].length; l++) {
					fp.setf(l, k, current_filt[k][l]);
				}
			}
			
			viz_stk.addSlice(fp);
		}
		
		ImagePlus viz_img =  new ImagePlus("filters", viz_stk);
		
		viz_img.show();
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		
		ImageStack plot_stack 	= 	new ImageStack(800, 400);
		for (int i = 0; i < feat.size(); i++) {
			Plot p = new Plot("feat"+i, "angle", "weight", anglesRad, feat.get(i)[0]);
			p.setSize(800, 400);
			plot_stack.addSlice("feat"+i, p.getProcessor());
		}
		new ImagePlus("filters", plot_stack).show();
		
	}

	public void 		showFilters(int idx){
		
//		ImageStack viz_stk 		=  	new ImageStack(featLen, nrRots);
//		for (int i = 0; i < feat.size(); i++) {
			FloatProcessor fp = new FloatProcessor(featLen, nrRots);
			
			float[][] current_filt = feat.get(idx);
			for (int k = 0; k < current_filt.length; k++) {
				for (int l = 0; l < current_filt[0].length; l++) {
					fp.setf(l, k, current_filt[k][l]);
				}
			}
			
//			viz_stk.addSlice(fp);
//		}
		
		ImagePlus viz_img =  new ImagePlus("filter", fp);
		
		viz_img.show();
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		
//		ImageStack plot_stack 	= 	new ImageStack(800, 400);
//		for (int i = 0; i < feat.size(); i++) {
			Plot p = new Plot("feat"+idx, "angle", "weight", anglesRad, feat.get(idx)[0]);
			p.setSize(800, 400);
//			p.setLimits(0, 2*Math.PI, -1, +1);
//			p.setColor(Color.RED);
//			p.addPoints(angProfile, feat.get(i), Plot.LINE);
//			plot_stack.addSlice("feat"+idx, p.getProcessor());
//		}
		
		new ImagePlus("feat"+idx, p.getProcessor()).show();
		
		
	}
	
	public double[]		calculateProfileFeatures(double[] profile){
		
		double[] score = new double[feat.size()];
		
		for (int i = 0; i < feat.size(); i++) {
			
			double max_score = Double.MIN_VALUE;
			
			for (int j = 0; j < feat.get(i).length; j++) {
				
				double scr = 0;
				for (int k = 0; k < feat.get(i)[0].length; k++) {
					scr  += feat.get(i)[j][k] * profile[k];
				}	
				
				if(scr>max_score){
					max_score = scr;
				}
				
			}
			
			score[i] = max_score;
			
		}
		
		return score;
		
	}

//	public float[] 			allFeatScore(double[] profile){ //, int shiftStep
//		
//		float[] out = new float[feat.size()];
//		
//		for (int i = 0; i < 1; i++) {//feat.size()
//			out[i] = featScore(profile, feat.get(i)); //, shiftStep
//		}
//		
//		return out;
//		
//	}
	
//	private float featScore(double[] profile, float[] feature){//, int shiftStep
//		
//		float curr_score 	= conv(profile, feature);
//		float total_score 	= curr_score;//conv(profile, feature);//Float.MIN_VALUE;
//		
//		ImageStack feat_shifts = new ImageStack(800, 400);
//		feat_shifts.addSlice("shf"+0, VizFeatures.plotValues(feature, Plot.BOX));
////		
//		
//		for (int d = angStepLen; d < feature.length; d+=angStepLen) {
//			
//			float[] f_sh = shiftVals(feature, d);
//			
//			feat_shifts.addSlice("shf"+d, VizFeatures.plotValues(f_sh, Plot.BOX));
//			
//			curr_score = conv(profile, f_sh);
//			
//			if(curr_score>total_score){
//				total_score = curr_score;
//			}
//			
//		}
//		
//		new ImagePlus ("shf"+0, ).show();
//		
//		return total_score;
//		
//	}
	
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
	
}

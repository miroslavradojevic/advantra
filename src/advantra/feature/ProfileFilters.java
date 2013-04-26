package advantra.feature;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;

import java.awt.Color;
import java.util.Vector;

public class ProfileFilters {
	
	/*
	 * feature is a vector corresponding to directional response
	 * the aim of the feature is to count the number of peaks
	 * of the directional response, borders and filter levels
	 * example feature 
	 * ON OFF 
	 * ON OFF ON OFF
	 * ON OFF ON OFF ON OFF
	 */
	
	private Vector<float[]> 	feat;
	private int					featLen;
	
//	private double				angStep;
	private int					angStepLen;
	
	private float[] 			angProfile; // will be used to guide the construction of steps
	
	public ProfileFilters(int featLen, double angStep){
		
		this.feat 		= new Vector<float[]>();
		this.featLen 	= featLen;
		
//		this.angStep 	= angStep;
		this.angStepLen = (int)Math.ceil((angStep/(2*Math.PI))*featLen);
		
		this.angProfile = new float[featLen];
		for (int i = 0; i < featLen; i++) {
			this.angProfile[i] = (float) (i*2*Math.PI/featLen);
		}
		
	}
	
	public void 		create(){
		
		int 	sLim 	= 1;
		float 	sumNeg, sumPos;
		int 	cnt 	= 0;
		
		// 1x
		for (int s = 1; s <= sLim; s++) { 
			
			float[] vec = new float[featLen];
			sumNeg = sumPos = 0;
			
			for (int i = 0; i < featLen; i++){
				
				if(i<s*angStepLen){
					vec[i] = +1;
					sumPos ++;
				}
				else{
					vec[i] = -1;
					sumNeg ++;
				}
				
			}
			
//			System.out.println("poss: "+sumPos+"  negs "+sumNeg+ " total "+(sumPos+sumNeg));
			
			for (int i = 0; i < featLen; i++){
				vec[i] = (vec[i]>0)? vec[i]/sumPos : vec[i]/sumNeg ;
			}
			
			feat.add(vec);
			cnt++;
			
		}
		
		//2x
		for (int s1 = 1; s1 <= sLim; s1++) { 
			for (int s2 = 1; s2 <= sLim; s2++) {
				for (int w1 = 1; w1 <= featLen; w1++) {
				
					
					if(s1*angStepLen+w1*angStepLen+s2*angStepLen<featLen){
						
						float[] vec = new float[featLen];
						sumNeg = sumPos = 0;
						
						for (int i = 0; i < vec.length; i++) {
							if(i<s1*angStepLen){
								vec[i] = +1; sumPos ++;
							}
							else if(i>=s1*angStepLen && i<s1*angStepLen+w1*angStepLen){
								vec[i] = -1; sumNeg ++;
							}
							else if(i>=s1*angStepLen+w1*angStepLen && i<s1*angStepLen+w1*angStepLen+s2*angStepLen){
								vec[i] = +1; sumPos ++;
							}
							else{
								vec[i] = -1; sumNeg ++;
							}
						}
						
//						System.out.println("poss: "+sumPos+"  negs "+sumNeg+ " total "+(sumPos+sumNeg));
						
						for (int i = 0; i < featLen; i++){
							vec[i] = (vec[i]>0)? vec[i]/sumPos : vec[i]/sumNeg ;
						}
						
						feat.add(vec);
						cnt++;
						
					}
				}
			}
		}
		
		//3x
		for (int s1 = 1; s1 <= sLim; s1++) { 
			for (int s2 = 1; s2 <= sLim; s2++) {
				for (int s3 = 1; s3 <= sLim; s3++) {
					for (int w1 = 1; w1 <= featLen; w1++) {
						for (int w2 = 1; w2 <= featLen; w2++) {
							if(s1*angStepLen+w1*angStepLen+s2*angStepLen+w2*angStepLen+s3*angStepLen<featLen){
								
								float[] vec = new float[featLen];
								sumNeg = sumPos = 0;
								
								for (int i = 0; i < vec.length; i++) {
									if(i<s1*angStepLen){
										vec[i] = +1; sumPos ++;
									}
									else if(i>=s1*angStepLen && i<s1*angStepLen+w1*angStepLen){
										vec[i] = -1; sumNeg ++;
									}
									else if(i>=s1*angStepLen+w1*angStepLen && i<s1*angStepLen+w1*angStepLen+s2*angStepLen){
										vec[i] = +1; sumPos ++;
									}
									else if(i>=s1*angStepLen+w1*angStepLen+s2*angStepLen && i<s1*angStepLen+w1*angStepLen+s2*angStepLen+w2*angStepLen){
										vec[i] = -1; sumNeg ++;
									}
									else if(i>=s1*angStepLen+w1*angStepLen+s2*angStepLen+w2*angStepLen && i<s1*angStepLen+w1*angStepLen+s2*angStepLen+w2*angStepLen+s3*angStepLen){
										vec[i] = +1; sumPos ++;
									}
									else{
										vec[i] = -1; sumNeg ++;
									}
								}
								
//								System.out.println("poss: "+sumPos+"  negs "+sumNeg+ " total "+(sumPos+sumNeg));
								
								for (int i = 0; i < featLen; i++){
									vec[i] = (vec[i]>0)? vec[i]/sumPos : vec[i]/sumNeg ;
								}
								
								feat.add(vec);
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
						for (int w1 = 1; w1 <= featLen; w1++) {
							for (int w2 = 1; w2 <= featLen; w2++) {
								for (int w3 = 1; w3 <= featLen; w3++) {
									
									int Y1 = s1*angStepLen;
									int N1 = w1*angStepLen;
									int Y2 = s2*angStepLen;
									int N2 = w2*angStepLen;
									int Y3 = s3*angStepLen;
									int N3 = w3*angStepLen;
									int Y4 = s4*angStepLen;
									
									if(Y1+N1+Y2+N2+Y3+N3+Y4<featLen){
										
										float[] vec = new float[featLen];
										sumNeg = sumPos = 0;
										
										for (int i = 0; i < vec.length; i++) {
											if(i<Y1){
												vec[i] = +1; sumPos ++;
											}
											else if(i>=Y1 && i<Y1+N1){
												vec[i] = -1; sumNeg ++;
											}
											else if(i>=Y1+N1 && i<Y1+N1+Y2){
												vec[i] = +1; sumPos ++;
											}
											else if(i>=Y1+N1+Y2 && i<Y1+N1+Y2+N2){
												vec[i] = -1; sumNeg ++;
											}
											else if(i>=Y1+N1+Y2+N2 && i<Y1+N1+Y2+N2+Y3){
												vec[i] = +1; sumPos ++;
											}
											else if(i>=Y1+N1+Y2+N2+Y3 && i<Y1+N1+Y2+N2+Y3+N3){
												vec[i] = -1; sumNeg ++;
											}
											else if(i>=Y1+N1+Y2+N2+Y3+N3 && i<Y1+N1+Y2+N2+Y3+N3+Y4){
												vec[i] = +1; sumPos ++;
											}
											else{
												vec[i] = -1; sumNeg ++;
											}
										}
										
//										System.out.println("poss: "+sumPos+"  negs "+sumNeg+ " total "+(sumPos+sumNeg));
										
										for (int i = 0; i < featLen; i++){
											vec[i] = (vec[i]>0)? vec[i]/sumPos : vec[i]/sumNeg ;
										}
										
										feat.add(vec);
										cnt++;
										
									}
								}
							}
						}
					}
				}
			}
		}		
		
		System.out.println("total " + cnt + " feats");
		
	}
	
	public int 			nrFilters(){
		
		return feat.size();
	
	}
	
	public void 		showFilters(){
		
		ImageStack viz_stk 		=  	new ImageStack(featLen, 1);
		for (int i = 0; i < feat.size(); i++) {
			FloatProcessor fp = new FloatProcessor(featLen, 1, feat.get(i));
			viz_stk.addSlice(fp);
		}
		
		ImagePlus viz_img =  new ImagePlus("features", viz_stk);
		
		viz_img.show();
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		viz_img.getCanvas().zoomIn(0, 0);
		
		ImageStack plot_stack 	= 	new ImageStack(800, 400);
		for (int i = 0; i < feat.size(); i++) {
			Plot p = new Plot("feat"+i, "angle", "filter");
			p.setSize(800, 400);
			p.setLimits(0, 2*Math.PI, -1, +1);
			p.setColor(Color.RED);
			p.addPoints(angProfile, feat.get(i), Plot.LINE);
			plot_stack.addSlice("feat"+i, p.getProcessor());
		}
		
		new ImagePlus("features", plot_stack).show();
		
		
		
	}

	public double[]		calculateProfileFeatures(double[] profile){
		
		double[] score = new double[feat.size()];
		
		for (int i = 0; i < feat.size(); i++) {
			
			score[i] = 0;
			
			for (int j = 0; j < score.length; j++) {
				
				score[i]  += feat.get(i)[j] * profile[j];
				
			}
			
		}
		
		return score;
		
	}
}

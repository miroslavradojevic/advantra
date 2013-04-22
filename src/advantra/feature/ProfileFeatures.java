package advantra.feature;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

import java.awt.Color;
import java.util.Vector;

public class ProfileFeatures {
	
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
	
	private double				angStep;
	private int					angStepLen;
	
	private float[] 			angProfile; // will be used to guide the construction of steps
	
	public ProfileFeatures(int featLen, double angStep){
		
		this.feat 		= new Vector<float[]>();
		this.featLen 	= featLen;
		
		this.angStep 	= angStep;
		this.angStepLen = (int)Math.ceil((angStep/(2*Math.PI))*featLen);
		
		this.angProfile = new float[featLen];
		for (int i = 0; i < featLen; i++) {
			this.angProfile[i] = (float) (i*2*Math.PI/featLen);
		}
		
	}
	
	public void createFeatures(){
		
		int 	sLim = 1;
		float 	sumNeg, sumPos;
		int 	cnt 	= 0;
		
		
		// 1x
		for (int s = 1; s <= sLim; s++) { 
			
			float[] vec = new float[featLen];
			sumNeg = sumPos = 0;
			
			for (int i = 0; i < featLen; i++){
				
				if(angProfile[i]<s*angStep){
					vec[i] = +1;
					sumPos ++;
				}
				else{
					vec[i] = -1;
					sumNeg ++;
				}
				
			}
			
			System.out.println("poss: "+sumPos+"  negs "+sumNeg+ " total "+(sumPos+sumNeg));
			
			for (int i = 0; i < featLen; i++){
				vec[i] = (vec[i]>0)? vec[i]/sumPos : vec[i]/sumNeg ;
			}
			
			feat.add(vec);
			cnt++;
			
		}
		
		//2x
		for (int s1 = 1; s1 <= sLim; s1++) { 
			for (double w1 = angStep; w1 < 2*Math.PI; w1+=angStep) {
				for (int s2 = 1; s2 <= sLim; s2++) {
					if(s1*angStep+w1+s2*angStep<2*Math.PI){
						
						float[] vec = new float[featLen];
						sumNeg = sumPos = 0;
						
						for (int i = 0; i < vec.length; i++) {
							if(angProfile[i]<s1*angStep){
								vec[i] = +1; sumPos ++;
							}
							else if(angProfile[i]>=s1*angStep && angProfile[i]<s1*angStep+w1){
								vec[i] = -1; sumNeg ++;
							}
							else if(angProfile[i]>=s1*angStep+w1 && angProfile[i]<s1*angStep+w1+s2*angStep){
								vec[i] = +1; sumPos ++;
							}
							else{
								vec[i] = -1; sumNeg ++;
							}
						}
						
						System.out.println("poss: "+sumPos+"  negs "+sumNeg+ " total "+(sumPos+sumNeg));
						
						for (int i = 0; i < featLen; i++){
							vec[i] = (vec[i]>0)? vec[i]/sumPos : vec[i]/sumNeg ;
						}
						
						feat.add(vec);
						cnt++;
						
					}
				}
			}
		}
		
		System.out.println("" + cnt + " feats");
		
//			float[] vec = new float[featLen];
//			sumNeg = sumPos = 0;
//			for (int i = 0; i < featLen; i++){
//				
//				if(angProfile[i]<sLim*angStep){
//					vec[i] = 1;
//					sumPos ++;
//				}
//				else if(angProfile[i]>=sLim*angStep && angProfile[i]<2*sLim*angStep) {
//					vec[i] = -1;
//					sumNeg ++;
//				}
//				else if(angProfile[i]>=2*sLim*angStep && angProfile[i]<sLim*angStep){
//					
//				}
//				else {
//					
//				}
//				
//			}
//			
//			for (int i = 0; i < angProfile.length; i++){
//				vec[i] = (angProfile[i]<s_limit*angStep)? vec[i]/sumPos : vec[i]/sumNeg ;
//			}
//			
//			feat.add(vec);
//			cnt++;
	}
	
	public int 			nrFeats(){
		
		return feat.size();
	
	}
	
	public void 		showFeatures(){
		
		ImageStack viz_stk 		=  	new ImageStack(featLen, featLen);
		
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
	
}

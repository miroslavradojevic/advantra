package advantra.feature;

import ij.ImagePlus;

import java.util.Vector;

public class CircularFilterSet {

	public int[] 				angScaleDegrees;
	public Vector<CircularFilterConfiguration> filts;
	
	public CircularFilterSet(int[] angular_scale_degrees){
		
		angScaleDegrees = new int[angular_scale_degrees.length];
		for (int i = 0; i < angular_scale_degrees.length; i++) {
			angScaleDegrees[i] = angular_scale_degrees[i];
		}
		
		filts = new Vector<CircularFilterConfiguration>();
		
		for (int takeScale = 0; takeScale < angScaleDegrees.length; takeScale++) {
			
			int currentScale = angScaleDegrees[takeScale];
			
			// 1x
			//filts.add(new CircularFilterConfiguration(currentScale, new int[]{360}));
			
			// 2x
			for (int alfa = currentScale; alfa < 360; alfa+=currentScale) {
				for (int beta = currentScale; beta < 360; beta+=currentScale) {
					if(alfa + beta == 360){
						
						// check if it exists already
						boolean exists = false;
						
						for (int i = 0; i < filts.size(); i++) {
							if(filts.get(i).nrPeaks==2){
								if(
										currentScale==filts.get(i).angResDeg &&
										alfa==filts.get(i).angBtwPeakDeg[0] &&
										beta==filts.get(i).angBtwPeakDeg[1]
												){
									exists = true;
								}
								if(
										currentScale==filts.get(i).angResDeg &&
										alfa==filts.get(i).angBtwPeakDeg[1] &&
										beta==filts.get(i).angBtwPeakDeg[0]
												){
									exists = true;
								}
							}
						}

						if(!exists){
							int[] tmp = new int[]{alfa, beta};
							filts.add(new CircularFilterConfiguration(currentScale, tmp));
						}
						
					}
				}
			}
			
			// 3x
			for (int alfa = currentScale; alfa < 360; alfa+=currentScale) {
				for (int beta = currentScale; beta < 360; beta+=currentScale) {
					for (int gamma = currentScale; gamma < 360; gamma+=currentScale) {
						if(alfa + beta + gamma == 360){
						
							// check if it exists already
							boolean exists = false;
						
							for (int i = 0; i < filts.size(); i++) {
								if(filts.get(i).nrPeaks==3){
									if(
											currentScale==filts.get(i).angResDeg &&
											alfa==filts.get(i).angBtwPeakDeg[0] &&
											beta==filts.get(i).angBtwPeakDeg[1] &&
											gamma==filts.get(i).angBtwPeakDeg[2]
													){
										exists = true;
									}
									if(
											currentScale==filts.get(i).angResDeg &&
											alfa==filts.get(i).angBtwPeakDeg[1] &&
											beta==filts.get(i).angBtwPeakDeg[2] &&
											gamma==filts.get(i).angBtwPeakDeg[0]
											){
										exists = true;
									}
									if(
											currentScale==filts.get(i).angResDeg &&
											alfa==filts.get(i).angBtwPeakDeg[2] &&
											beta==filts.get(i).angBtwPeakDeg[0] &&
											gamma==filts.get(i).angBtwPeakDeg[1]
											){
										exists = true;
									}
								}
							}

							if(!exists){
								int[] tmp = new int[]{alfa, beta, gamma};
								filts.add(new CircularFilterConfiguration(currentScale, tmp));
							}
						}
					}
				}
			}
			
			// 4x
			for (int alfa = currentScale; alfa < 360; alfa+=currentScale) {
				for (int beta = currentScale; beta < 360; beta+=currentScale) {
					for (int gamma = currentScale; gamma < 360; gamma+=currentScale) {
						for (int delta = currentScale; delta < 360; delta+=currentScale) {
							if(alfa + beta + gamma + delta == 360){
							
								// check if it exists already
								boolean exists = false;
							
								for (int i = 0; i < filts.size(); i++) {
									if(filts.get(i).nrPeaks==4){
										if(
												currentScale==filts.get(i).angResDeg &&
												alfa==filts.get(i).angBtwPeakDeg[0] &&
												beta==filts.get(i).angBtwPeakDeg[1] &&
												gamma==filts.get(i).angBtwPeakDeg[2] &&
												delta==filts.get(i).angBtwPeakDeg[3]
														){
											exists = true;
										}
										if(
												currentScale==filts.get(i).angResDeg &&
												alfa==filts.get(i).angBtwPeakDeg[1] &&
												beta==filts.get(i).angBtwPeakDeg[2] &&
												gamma==filts.get(i).angBtwPeakDeg[3] &&
												delta==filts.get(i).angBtwPeakDeg[0]
												){
											exists = true;
										}
										if(
												currentScale==filts.get(i).angResDeg &&
												alfa==filts.get(i).angBtwPeakDeg[2] &&
												beta==filts.get(i).angBtwPeakDeg[3] &&
												gamma==filts.get(i).angBtwPeakDeg[0] &&
												delta==filts.get(i).angBtwPeakDeg[1]
												){
											exists = true;
										}
										if(
												currentScale==filts.get(i).angResDeg &&
												alfa==filts.get(i).angBtwPeakDeg[3] &&
												beta==filts.get(i).angBtwPeakDeg[0] &&
												gamma==filts.get(i).angBtwPeakDeg[1] &&
												delta==filts.get(i).angBtwPeakDeg[2]
												){
											exists = true;
										}
									}
								}
							
	
								if(!exists){
									int[] tmp = new int[]{alfa, beta, gamma, delta};
									filts.add(new CircularFilterConfiguration(currentScale, tmp));
								}
							}
						}
					}
				}
			}
			
		} // scale loop
		
		System.out.println("total configs created: "+filts.size());
		
	}
	
	public double[] calculateScore(double[] angularProfile){
		double[] score = new double[filts.size()];
		
		for (int i = 0; i < filts.size(); i++) {
			
			score[i] = filts.get(i).calculateScore(angularProfile);
			
		}
		
		return score;
	}
	
	public void showConfigs(){
		
		for (int i = 0; i < filts.size(); i++) {
			
			CircularFilterConfiguration currentFilter = filts.get(i);
			
			String name;
			name = "R="+currentFilter.angResDeg;
			if(currentFilter.nrPeaks>=1){
				name+=",alfa="+currentFilter.angBtwPeakDeg[0];
			}
			if(currentFilter.nrPeaks>=2){
				name+=",beta="+currentFilter.angBtwPeakDeg[1];
			}
			if(currentFilter.nrPeaks>=3){
				name+=",gamma="+currentFilter.angBtwPeakDeg[2];
			}
			if(currentFilter.nrPeaks>=4){
				name+=",delta="+currentFilter.angBtwPeakDeg[3];
			}
			
			new ImagePlus(name, currentFilter.plot()).show();
			
		}
		
	}
	
	public void showConfigs(int idx){
		
			CircularFilterConfiguration currentFilter = filts.get(idx);
			
			String name;
			name = "R="+currentFilter.angResDeg;
			if(currentFilter.nrPeaks>=1){
				name+=",alfa="+currentFilter.angBtwPeakDeg[0];
			}
			if(currentFilter.nrPeaks>=2){
				name+=",beta="+currentFilter.angBtwPeakDeg[1];
			}
			if(currentFilter.nrPeaks>=3){
				name+=",gamma="+currentFilter.angBtwPeakDeg[2];
			}
			if(currentFilter.nrPeaks>=4){
				name+=",delta="+currentFilter.angBtwPeakDeg[3];
			}
			
			new ImagePlus(name, currentFilter.plot()).show();
		
	}
	
}
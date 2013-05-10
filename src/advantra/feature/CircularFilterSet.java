package advantra.feature;

import ij.IJ;
import ij.ImagePlus;

import java.util.Vector;

public class CircularFilterSet {

	public int[] 				angScaleDegrees;
	public int					minAngScaleDegrees;
	public Vector<CircularFilterConfiguration> filts;

	public CircularFilterSet(
            int[] angular_scale_degrees
    )
	{

		filts = new Vector<CircularFilterConfiguration>();

		angScaleDegrees = new int[angular_scale_degrees.length];
		minAngScaleDegrees = Integer.MAX_VALUE;
		for (int i = 0; i < angular_scale_degrees.length; i++) {
			angScaleDegrees[i] = angular_scale_degrees[i];
			if (angScaleDegrees[i]<minAngScaleDegrees){
				minAngScaleDegrees = angScaleDegrees[i];
			}
		}

		int nrScales = angular_scale_degrees.length;  // depends on the number of scales

        // make permutations of length 1 using given angles
		int nrPer1 = (int) Math.pow(nrScales, 1); // 1x
		int[][] per1 = new int[nrPer1][1];
		int c = 0;
		for (int k = 0; k < nrScales; k++){
			per1[c][0] = angular_scale_degrees[k];
			c++;
		}

		for (int i = 0; i < nrPer1; i++){
			if (360>per1[i][0]/2+per1[i][0]/2){
				// add it to the list of filters
				filts.add(new CircularFilterConfiguration(per1[i], new int[]{360}));
			}
		}


		int nrPer2 = (int) Math.pow(nrScales, 2); // 2x
		int[][] per2 = new int[nrPer2][2];
		c = 0;
		for (int k = 0; k < nrScales; k++){
			for (int l = 0; l < nrScales; l++){
				per2[c][0] = angular_scale_degrees[k];
				per2[c][1] = angular_scale_degrees[l];
				c++;
			}
		}

		for (int i = 0; i < nrPer2; i++){
			for (int d1 = minAngScaleDegrees; d1 < 360; d1+=minAngScaleDegrees/2){
				for (int d2 = minAngScaleDegrees; d2 < 360; d2+=minAngScaleDegrees/2){
					boolean isConfiguration = false;
					isConfiguration =
							(d1+d2==360) &&
							(d1>per2[i][0]/2+per2[i][1]/2) &&
							(d2>per2[i][1]/2+per2[i][0]/2);

					if(isConfiguration){

						boolean covered = false;
						// check if it exists so far in other rotations
						for (int k = 0; k < filts.size(); k++){
							if(filts.get(k).angResDeg.length==2){// if it's with 2 peaks
								if(
										(
											per2[i][0]==filts.get(k).angResDeg[0] &&
											d1        ==filts.get(k).angBtwPeakDeg[0] &&
											per2[i][1]==filts.get(k).angResDeg[1] &&
											d2        ==filts.get(k).angBtwPeakDeg[1]
										)
										||
												(
													per2[i][0]==filts.get(k).angResDeg[1] &&
													d1        ==filts.get(k).angBtwPeakDeg[1] &&
													per2[i][1]==filts.get(k).angResDeg[0] &&
													d2        ==filts.get(k).angBtwPeakDeg[0]
												)
										){
									covered = true;
								}
							}
						}

						if(!covered){
							filts.add(new CircularFilterConfiguration(per2[i], new int[]{d1, d2}));
						}

					}
				}
			}
		}


		int nrPer3 = (int) Math.pow(nrScales, 3); // 3x
		int[][] per3 = new int[nrPer3][3];
		c = 0;
		for (int k = 0; k < nrScales; k++){
			for (int l = 0; l < nrScales; l++){
				for (int m = 0; m < nrScales; m++){
					per3[c][0] = angular_scale_degrees[k];
					per3[c][1] = angular_scale_degrees[l];
					per3[c][2] = angular_scale_degrees[m];
					c++;
				}
			}
		}

		for (int i = 0; i < nrPer3; i++){
			for (int d1 = minAngScaleDegrees; d1 < 360; d1+=minAngScaleDegrees/2){
				for (int d2 = minAngScaleDegrees; d2 < 360; d2+=minAngScaleDegrees/2){
					for (int d3 = minAngScaleDegrees; d3 < 360; d3+=minAngScaleDegrees/2){

						boolean isConfiguration = false;
						isConfiguration =
								(d1+d2+d3==360) &&
								(d1>per3[i][0]/2+per3[i][1]/2) &&
								(d2>per3[i][1]/2+per3[i][2]/2) &&
								(d3>per3[i][2]/2+per3[i][0]/2);

						if(isConfiguration){

							boolean covered = false;
							// check if it exists so far in other rotations
							for (int k = 0; k < filts.size(); k++){
								if(filts.get(k).angResDeg.length==3){// if it's with 3 peaks
									if(
											(
											per3[i][0]==filts.get(k).angResDeg[0] &&
											d1        ==filts.get(k).angBtwPeakDeg[0] &&
											per3[i][1]==filts.get(k).angResDeg[1] &&
											d2        ==filts.get(k).angBtwPeakDeg[1]  &&
											per3[i][2]==filts.get(k).angResDeg[2] &&
											d3        ==filts.get(k).angBtwPeakDeg[2]
											)
											||
											(
											per3[i][0]==filts.get(k).angResDeg[1] &&
											d1        ==filts.get(k).angBtwPeakDeg[1] &&
											per3[i][1]==filts.get(k).angResDeg[2] &&
											d2        ==filts.get(k).angBtwPeakDeg[2]  &&
											per3[i][2]==filts.get(k).angResDeg[0] &&
											d3        ==filts.get(k).angBtwPeakDeg[0]
											)
											||
											(
											per3[i][0]==filts.get(k).angResDeg[2] &&
											d1        ==filts.get(k).angBtwPeakDeg[2] &&
											per3[i][1]==filts.get(k).angResDeg[0] &&
											d2        ==filts.get(k).angBtwPeakDeg[0]  &&
											per3[i][2]==filts.get(k).angResDeg[1] &&
											d3        ==filts.get(k).angBtwPeakDeg[1]
											)
									)
									{
										covered = true;
									}
								}
							}

							if(!covered){
								filts.add(new CircularFilterConfiguration(per3[i], new int[]{d1, d2, d3}));
							}
						}
					}
				}
			}
		}

		int nrPer4 = (int) Math.pow(nrScales, 4); // 4x
		int[][] per4 = new int[nrPer4][4];
		c = 0;
		for (int k = 0; k < nrScales; k++){
			for (int l = 0; l < nrScales; l++){
				for (int m = 0; m < nrScales; m++){
					for (int n = 0; n < nrScales; n++){
						per4[c][0] = angular_scale_degrees[k];
						per4[c][1] = angular_scale_degrees[l];
						per4[c][2] = angular_scale_degrees[m];
						per4[c][3] = angular_scale_degrees[n];
						c++;
					}
				}
			}
		}

		for (int i = 0; i < nrPer4; i++){
			for (int d1 = minAngScaleDegrees; d1 < 360; d1+=minAngScaleDegrees/2){
				for (int d2 = minAngScaleDegrees; d2 < 360; d2+=minAngScaleDegrees/2){
					for (int d3 = minAngScaleDegrees; d3 < 360; d3+=minAngScaleDegrees/2){
						for (int d4 = minAngScaleDegrees; d4 < 360; d4+=minAngScaleDegrees/2){

							boolean isConfiguration = false;
							isConfiguration =
									(d1+d2+d3+d4==360) &&
											(d1>per4[i][0]/2+per4[i][1]/2) &&
											(d2>per4[i][1]/2+per4[i][2]/2) &&
											(d3>per4[i][2]/2+per4[i][3]/2) &&
											(d4>per4[i][3]/2+per4[i][0]/2)
							;

							if(isConfiguration){

								boolean covered = false;
								// check if it exists so far in other rotations
								for (int k = 0; k < filts.size(); k++){
									if(filts.get(k).angResDeg.length==4){// if it's with 4 peaks
										if(
												(
												per4[i][0]==filts.get(k).angResDeg[0] &&
												d1        ==filts.get(k).angBtwPeakDeg[0] &&
												per4[i][1]==filts.get(k).angResDeg[1] &&
												d2        ==filts.get(k).angBtwPeakDeg[1]  &&
												per4[i][2]==filts.get(k).angResDeg[2] &&
												d3        ==filts.get(k).angBtwPeakDeg[2]  &&
												per4[i][3]==filts.get(k).angResDeg[3] &&
												d4        ==filts.get(k).angBtwPeakDeg[3]
												)
												||
												(
												per4[i][0]==filts.get(k).angResDeg[1] &&
												d1        ==filts.get(k).angBtwPeakDeg[1] &&
												per4[i][1]==filts.get(k).angResDeg[2] &&
												d2        ==filts.get(k).angBtwPeakDeg[2]  &&
												per4[i][2]==filts.get(k).angResDeg[3] &&
												d3        ==filts.get(k).angBtwPeakDeg[3]  &&
												per4[i][3]==filts.get(k).angResDeg[0] &&
												d4        ==filts.get(k).angBtwPeakDeg[0]
												)
												||
												(
												per4[i][0]==filts.get(k).angResDeg[2] &&
												d1        ==filts.get(k).angBtwPeakDeg[2] &&
												per4[i][1]==filts.get(k).angResDeg[3] &&
												d2        ==filts.get(k).angBtwPeakDeg[3]  &&
												per4[i][2]==filts.get(k).angResDeg[0] &&
												d3        ==filts.get(k).angBtwPeakDeg[0]  &&
												per4[i][3]==filts.get(k).angResDeg[1] &&
												d4        ==filts.get(k).angBtwPeakDeg[1]
												)
												||
												(
												per4[i][0]==filts.get(k).angResDeg[3] &&
												d1        ==filts.get(k).angBtwPeakDeg[3] &&
												per4[i][1]==filts.get(k).angResDeg[0] &&
												d2        ==filts.get(k).angBtwPeakDeg[0]  &&
												per4[i][2]==filts.get(k).angResDeg[1] &&
												d3        ==filts.get(k).angBtwPeakDeg[1]  &&
												per4[i][3]==filts.get(k).angResDeg[2] &&
												d4        ==filts.get(k).angBtwPeakDeg[2]
												)
										)
										{
											covered = true;
										}
									}
								}

								if(!covered){
									filts.add(new CircularFilterConfiguration(per4[i], new int[]{d1, d2, d3, d4}));
								}
							}


						}
					}
				}
			}
		}



//
//		for (int takeScale = 0; takeScale < angScaleDegrees.length; takeScale++) {
//
//			int currentScale = angScaleDegrees[takeScale];
//
//			// 1x
//			filts.add(new CircularFilterConfiguration(currentScale, new int[]{360}));
//
//			// 2x
//
//			for (int alfa = currentScale; alfa < 360; alfa+=currentScale/2) {
//				for (int beta = currentScale; beta < 360; beta+=currentScale/2) {
//					if(alfa + beta == 360){
//
//						// check if it exists already
//						boolean exists = false;
//
//						for (int i = 0; i < filts.size(); i++) {
//							if(filts.get(i).nrPeaks==2){
//								if(
//										currentScale==filts.get(i).angResDeg &&
//										alfa==filts.get(i).angBtwPeakDeg[0] &&
//										beta==filts.get(i).angBtwPeakDeg[1]
//												){
//									exists = true;
//								}
//								if(
//										currentScale==filts.get(i).angResDeg &&
//										alfa==filts.get(i).angBtwPeakDeg[1] &&
//										beta==filts.get(i).angBtwPeakDeg[0]
//												){
//									exists = true;
//								}
//							}
//						}
//
//						if(!exists){
//							int[] tmp = new int[]{alfa, beta};
//							filts.add(new CircularFilterConfiguration(currentScale, tmp));
//						}
//
//					}
//				}
//			}
//
//			// 3x
//			for (int alfa = currentScale; alfa < 360; alfa+=currentScale/2) {
//				for (int beta = currentScale; beta < 360; beta+=currentScale/2) {
//					for (int gamma = currentScale; gamma < 360; gamma+=currentScale/2) {
//						if(alfa + beta + gamma == 360){
//
//							// check if it exists already
//							boolean exists = false;
//
//							for (int i = 0; i < filts.size(); i++) {
//								if(filts.get(i).nrPeaks==3){
//									if(
//											currentScale==filts.get(i).angResDeg &&
//											alfa==filts.get(i).angBtwPeakDeg[0] &&
//											beta==filts.get(i).angBtwPeakDeg[1] &&
//											gamma==filts.get(i).angBtwPeakDeg[2]
//													){
//										exists = true;
//									}
//									if(
//											currentScale==filts.get(i).angResDeg &&
//											alfa==filts.get(i).angBtwPeakDeg[1] &&
//											beta==filts.get(i).angBtwPeakDeg[2] &&
//											gamma==filts.get(i).angBtwPeakDeg[0]
//											){
//										exists = true;
//									}
//									if(
//											currentScale==filts.get(i).angResDeg &&
//											alfa==filts.get(i).angBtwPeakDeg[2] &&
//											beta==filts.get(i).angBtwPeakDeg[0] &&
//											gamma==filts.get(i).angBtwPeakDeg[1]
//											){
//										exists = true;
//									}
//								}
//							}
//
//							if(!exists){
//								int[] tmp = new int[]{alfa, beta, gamma};
//								filts.add(new CircularFilterConfiguration(currentScale, tmp));
//							}
//						}
//					}
//				}
//			}

			// 4x
//			for (int alfa = currentScale; alfa < 360; alfa+=currentScale/2) {
//				for (int beta = currentScale; beta < 360; beta+=currentScale/2) {
//					for (int gamma = currentScale; gamma < 360; gamma+=currentScale/2) {
//						for (int delta = currentScale; delta < 360; delta+=currentScale/2) {
//							if(alfa + beta + gamma + delta == 360){
//
//								// check if it exists already
//								boolean exists = false;
//
//								for (int i = 0; i < filts.size(); i++) {
//									if(filts.get(i).nrPeaks==4){
//										if(
//												currentScale==filts.get(i).angResDeg &&
//												alfa==filts.get(i).angBtwPeakDeg[0] &&
//												beta==filts.get(i).angBtwPeakDeg[1] &&
//												gamma==filts.get(i).angBtwPeakDeg[2] &&
//												delta==filts.get(i).angBtwPeakDeg[3]
//														){
//											exists = true;
//										}
//										if(
//												currentScale==filts.get(i).angResDeg &&
//												alfa==filts.get(i).angBtwPeakDeg[1] &&
//												beta==filts.get(i).angBtwPeakDeg[2] &&
//												gamma==filts.get(i).angBtwPeakDeg[3] &&
//												delta==filts.get(i).angBtwPeakDeg[0]
//												){
//											exists = true;
//										}
//										if(
//												currentScale==filts.get(i).angResDeg &&
//												alfa==filts.get(i).angBtwPeakDeg[2] &&
//												beta==filts.get(i).angBtwPeakDeg[3] &&
//												gamma==filts.get(i).angBtwPeakDeg[0] &&
//												delta==filts.get(i).angBtwPeakDeg[1]
//												){
//											exists = true;
//										}
//										if(
//												currentScale==filts.get(i).angResDeg &&
//												alfa==filts.get(i).angBtwPeakDeg[3] &&
//												beta==filts.get(i).angBtwPeakDeg[0] &&
//												gamma==filts.get(i).angBtwPeakDeg[1] &&
//												delta==filts.get(i).angBtwPeakDeg[2]
//												){
//											exists = true;
//										}
//									}
//								}
//
//
//								if(!exists){
//									int[] tmp = new int[]{alfa, beta, gamma, delta};
//									filts.add(new CircularFilterConfiguration(currentScale, tmp));
//								}
//							}
//						}
//					}
//				}
//			}

//		} // scale loop

		//System.out.println("total configs created: "+filts.size());

	}

	public double[] calculateScore(double[] angularProfile)
	{

        double[] score = new double[filts.size()];

		for (int i = 0; i < filts.size(); i++) {

//			score[i] = filts.get(i).calculateScore(angularProfile);

		}

		return score;
	}

	public double[] calculateScore(float[] angularProfile)
	{

        double[] score = new double[filts.size()];

		for (int i = 0; i < filts.size(); i++) {

//			score[i] = filts.get(i).calculateScore(angularProfile);

		}

		return score;
	}

    public double[] calculateScore(double[] angularProfile, int[] choose_filters)
	{

        double[] score = new double[choose_filters.length];

        for (int i = 0; i < choose_filters.length; i++){

//            score[i] = filts.get(choose_filters[i]).calculateScore(angularProfile);

        }

        return score;

    }

	public void showConfigs()
	{

		for (int i = 0; i < filts.size(); i++) {

			CircularFilterConfiguration currentFilter = filts.get(i);

			String name;
			name = "CONF_"+i;
			if(currentFilter.angResDeg.length>=1){
				name += ",Pks="+currentFilter.angResDeg[0];
			}
			if(currentFilter.angResDeg.length>=2){
				name += ","+currentFilter.angResDeg[1];
			}
			if(currentFilter.angResDeg.length>=3){
				name += ","+currentFilter.angResDeg[2];
			}
			if(currentFilter.angResDeg.length>=4){
				name += ","+currentFilter.angResDeg[3];
			}

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

//			new ImagePlus(name, currentFilter.plotFilter()).show();

		}

	}

	public void showConfigs(int idx)
	{

			CircularFilterConfiguration currentFilter = filts.get(idx);

			String name;
			name = "CONF_"+idx+",R="+currentFilter.angResDeg;
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

//			new ImagePlus(name, currentFilter.plotFilter()).show();

	}

}
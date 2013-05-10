package advantra.feature;

import ij.IJ;
import ij.ImagePlus;

import java.util.Vector;

public class FilterSet {

	public int[] 				            angScaleDegrees;
	public int					            minAngScaleDegrees;

    public double[]                         radScale;
    public double                           minRadScale;

	public Vector<CircularConfiguration>    circConfs;
    public Vector<RadialConfiguration>      radlConfs;

	public FilterSet(
            int[]       angular_scale_degrees,
            double[]    radial_scale
    )
	{

        int nrScales;

        /*
        fill in circular configurations
         */

        angScaleDegrees = new int[angular_scale_degrees.length];
        minAngScaleDegrees = Integer.MAX_VALUE;
        for (int i = 0; i < angular_scale_degrees.length; i++) {
            angScaleDegrees[i] = angular_scale_degrees[i];
            if (angScaleDegrees[i]<minAngScaleDegrees){
                minAngScaleDegrees = angScaleDegrees[i];
            }
        }

        nrScales = angular_scale_degrees.length;


        circConfs = new Vector<CircularConfiguration>();

        // 1
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
				circConfs.add(new CircularConfiguration(per1[i], new int[]{360}));
			}
		}

        // 2
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
						for (int k = 0; k < circConfs.size(); k++){
							if(circConfs.get(k).angResDeg.length==2){// if it's with 2 peaks
								if(
										(
											per2[i][0]==circConfs.get(k).angResDeg[0] &&
											d1        ==circConfs.get(k).angBtwPeakDeg[0] &&
											per2[i][1]==circConfs.get(k).angResDeg[1] &&
											d2        ==circConfs.get(k).angBtwPeakDeg[1]
										)
										||
												(
													per2[i][0]==circConfs.get(k).angResDeg[1] &&
													d1        ==circConfs.get(k).angBtwPeakDeg[1] &&
													per2[i][1]==circConfs.get(k).angResDeg[0] &&
													d2        ==circConfs.get(k).angBtwPeakDeg[0]
												)
										){
									covered = true;
								}
							}
						}

						if(!covered){
                            circConfs.add(new CircularConfiguration(per2[i], new int[]{d1, d2}));
						}

					}
				}
			}
		}

        // 3
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
							for (int k = 0; k < circConfs.size(); k++){
								if(circConfs.get(k).angResDeg.length==3){// if it's with 3 peaks
									if(
											(
											per3[i][0]==circConfs.get(k).angResDeg[0] &&
											d1        ==circConfs.get(k).angBtwPeakDeg[0] &&
											per3[i][1]==circConfs.get(k).angResDeg[1] &&
											d2        ==circConfs.get(k).angBtwPeakDeg[1]  &&
											per3[i][2]==circConfs.get(k).angResDeg[2] &&
											d3        ==circConfs.get(k).angBtwPeakDeg[2]
											)
											||
											(
											per3[i][0]==circConfs.get(k).angResDeg[1] &&
											d1        ==circConfs.get(k).angBtwPeakDeg[1] &&
											per3[i][1]==circConfs.get(k).angResDeg[2] &&
											d2        ==circConfs.get(k).angBtwPeakDeg[2]  &&
											per3[i][2]==circConfs.get(k).angResDeg[0] &&
											d3        ==circConfs.get(k).angBtwPeakDeg[0]
											)
											||
											(
											per3[i][0]==circConfs.get(k).angResDeg[2] &&
											d1        ==circConfs.get(k).angBtwPeakDeg[2] &&
											per3[i][1]==circConfs.get(k).angResDeg[0] &&
											d2        ==circConfs.get(k).angBtwPeakDeg[0]  &&
											per3[i][2]==circConfs.get(k).angResDeg[1] &&
											d3        ==circConfs.get(k).angBtwPeakDeg[1]
											)
									)
									{
										covered = true;
									}
								}
							}

							if(!covered){
                                circConfs.add(new CircularConfiguration(per3[i], new int[]{d1, d2, d3}));
							}
						}
					}
				}
			}
		}
        // 4
		int nrPer4 = (int) Math.pow(nrScales, 4); // 4x
		int[][] per4 = new int[nrPer4][4];
		c = 0;
		for (int k = 0; k < nrScales; k++)
            for (int l = 0; l < nrScales; l++)
                for (int m = 0; m < nrScales; m++)
                    for (int n = 0; n < nrScales; n++) {
                        per4[c][0] = angular_scale_degrees[k];
                        per4[c][1] = angular_scale_degrees[l];
                        per4[c][2] = angular_scale_degrees[m];
                        per4[c][3] = angular_scale_degrees[n];
                        c++;
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
								for (int k = 0; k < circConfs.size(); k++){
									if(circConfs.get(k).angResDeg.length==4){// if it's with 4 peaks
										if(
												(
												per4[i][0]==circConfs.get(k).angResDeg[0] &&
												d1        ==circConfs.get(k).angBtwPeakDeg[0] &&
												per4[i][1]==circConfs.get(k).angResDeg[1] &&
												d2        ==circConfs.get(k).angBtwPeakDeg[1]  &&
												per4[i][2]==circConfs.get(k).angResDeg[2] &&
												d3        ==circConfs.get(k).angBtwPeakDeg[2]  &&
												per4[i][3]==circConfs.get(k).angResDeg[3] &&
												d4        ==circConfs.get(k).angBtwPeakDeg[3]
												)
												||
												(
												per4[i][0]==circConfs.get(k).angResDeg[1] &&
												d1        ==circConfs.get(k).angBtwPeakDeg[1] &&
												per4[i][1]==circConfs.get(k).angResDeg[2] &&
												d2        ==circConfs.get(k).angBtwPeakDeg[2]  &&
												per4[i][2]==circConfs.get(k).angResDeg[3] &&
												d3        ==circConfs.get(k).angBtwPeakDeg[3]  &&
												per4[i][3]==circConfs.get(k).angResDeg[0] &&
												d4        ==circConfs.get(k).angBtwPeakDeg[0]
												)
												||
												(
												per4[i][0]==circConfs.get(k).angResDeg[2] &&
												d1        ==circConfs.get(k).angBtwPeakDeg[2] &&
												per4[i][1]==circConfs.get(k).angResDeg[3] &&
												d2        ==circConfs.get(k).angBtwPeakDeg[3]  &&
												per4[i][2]==circConfs.get(k).angResDeg[0] &&
												d3        ==circConfs.get(k).angBtwPeakDeg[0]  &&
												per4[i][3]==circConfs.get(k).angResDeg[1] &&
												d4        ==circConfs.get(k).angBtwPeakDeg[1]
												)
												||
												(
												per4[i][0]==circConfs.get(k).angResDeg[3] &&
												d1        ==circConfs.get(k).angBtwPeakDeg[3] &&
												per4[i][1]==circConfs.get(k).angResDeg[0] &&
												d2        ==circConfs.get(k).angBtwPeakDeg[0]  &&
												per4[i][2]==circConfs.get(k).angResDeg[1] &&
												d3        ==circConfs.get(k).angBtwPeakDeg[1]  &&
												per4[i][3]==circConfs.get(k).angResDeg[2] &&
												d4        ==circConfs.get(k).angBtwPeakDeg[2]
												)
										)
										{
											covered = true;
										}
									}
								}

								if(!covered){
                                    circConfs.add(new CircularConfiguration(per4[i], new int[]{d1, d2, d3, d4}));
								}
							}


						}
					}
				}
			}
		}

        /*
        fill in radial configurations
        */

        radScale = new double[radial_scale.length];
        minRadScale = Double.MAX_VALUE;
        for (int i = 0; i < radial_scale.length; i++) {
            radScale[i] = radial_scale[i];
            if (radScale[i]<minRadScale){
                minRadScale = radScale[i];
            }
        }

        nrScales = radial_scale.length;

        radlConfs = new Vector<RadialConfiguration>();

        // 1
        int nrRep1 = (int) Math.pow(nrScales, 1); // 1x
        double[][] rep1 = new double[nrRep1][1];
        c = 0;
        for (int k = 0; k < nrScales; k++){
            rep1[c][0] = radial_scale[k];
            c++;
        }

        for (int i = 0; i < nrRep1; i++){
            if (rep1[i][0]<1){
                  radlConfs.add(new RadialConfiguration(rep1[i], new double[]{0}));
//                circConfs.add(new CircularConfiguration(per1[i], new int[]{360}));
            }
        }

        // 2
        // 3
        // 4





	}

    public void print()
    {
        for (int i = 0; i < circConfs.size(); i++){
            circConfs.get(i).print();
        }

        for (int i = 0; i < radlConfs.size(); i++){
                                     radlConfs.get(i).print();
        }

    }

	public void showConfigs()
	{

		for (int i = 0; i < circConfs.size(); i++) {

			CircularConfiguration currentFilter = circConfs.get(i);

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

			CircularConfiguration currentFilter = circConfs.get(idx);

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

    public double[] calculateScore(double[] angularProfile)
    {

        double[] score = new double[circConfs.size()];

        for (int i = 0; i < circConfs.size(); i++) {

//			score[i] = filts.get(i).calculateScore(); // TODO change this

        }

        return score;
    }

}

   /*
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

    */

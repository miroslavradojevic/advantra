package advantra.feature;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

import java.util.Vector;

public class FilterSet {

	public int[] 				            angScaleDegrees;
	public int					            minAngScaleDegrees;

    public double[]                         radScale;
    public double                           minRadScale;

	public Vector<CircularConfiguration>    circConfs;
    public Vector<RadialConfiguration>      radlConfs;

	public float[] score;

	public FilterSet(
            					int[]       angular_scale_degrees,
								double[]    angular_scale_radiuses,
            					double[]    radial_scale
    )
	{

        int nrScales;

        // temporary
        double innerRingRadiusForTesting = 0.0;

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

		int c;

//        // 1
//		int nrPer1 = (int) Math.pow(nrScales, 1); // 1x
//		int[][] per1 = new int[nrPer1][1];
//		c = 0;
//		for (int k = 0; k < nrScales; k++){
//			per1[c][0] = angular_scale_degrees[k];
//			c++;
//		}
//
//		for (int i = 0; i < nrPer1; i++){
//			if (360>per1[i][0]/2+per1[i][0]/2){
//				// add it to the list of filters
//                // add it for all radial scales
//                for (int k = 0; k < angular_scale_radiuses.length; k++){
//                    circConfs.add(new CircularConfiguration(per1[i], new int[]{360}, new double[]{0.0, angular_scale_radiuses[k]} ));
//                    if (angular_scale_radiuses[k]<0.99) circConfs.add(new CircularConfiguration(per1[i], new int[]{360}, new double[]{angular_scale_radiuses[k], 1.0} ));
//                }
//
//			}
//		}

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
                            for (int k = 0; k < angular_scale_radiuses.length-1; k++){
                                //circConfs.add(new CircularConfiguration(per2[i], new int[]{d1, d2},  new double[]{0.0, angular_scale_radiuses[k]}));
                                if (angular_scale_radiuses[k]<0.99)
                                    circConfs.add(new CircularConfiguration(per2[i], new int[]{d1, d2},  new double[]{angular_scale_radiuses[k], angular_scale_radiuses[k+1]}, innerRingRadiusForTesting));
                                int last_one = angular_scale_radiuses.length-1;
                                circConfs.add(new CircularConfiguration(per2[i], new int[]{d1, d2},  new double[]{angular_scale_radiuses[last_one], 1.0}, innerRingRadiusForTesting));
                            }
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
			for (int d1 = minAngScaleDegrees; d1 < 360; d1+=minAngScaleDegrees){
				for (int d2 = minAngScaleDegrees; d2 < 360; d2+=minAngScaleDegrees){
					for (int d3 = minAngScaleDegrees; d3 < 360; d3+=minAngScaleDegrees){

						boolean isConfiguration =
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
								// add it for all radial scales
								for (int k = 0; k < angular_scale_radiuses.length-1; k++){
									//circConfs.add(new CircularConfiguration(per3[i], new int[]{d1, d2, d3},  new double[]{0.0, angular_scale_radiuses[k]}));
                                    if (angular_scale_radiuses[k]<0.99)
                                        circConfs.add(new CircularConfiguration(per3[i], new int[]{d1, d2, d3},  new double[]{angular_scale_radiuses[k], angular_scale_radiuses[k+1]}, innerRingRadiusForTesting));

                                    int last_one = angular_scale_radiuses.length-1;
                                            circConfs.add(new CircularConfiguration(per3[i], new int[]{d1, d2, d3},  new double[]{angular_scale_radiuses[last_one], 1.0}, innerRingRadiusForTesting));
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

/*        // 1
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
            }
        }*/

        // 2
        // 3
        // 4

		score = new float[circConfs.size()+radlConfs.size()];

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

	public ImageStack plot(int N)
	{

		ImageStack viz = new ImageStack(N, N);

		for (int i = 0; i < circConfs.size(); i++) {

			CircularConfiguration currentFilter = circConfs.get(i);

			String name;
			name = "CCONF_"+i;

			for (int j = 0; j < currentFilter.nrPeaks; j++){
				name += ","+currentFilter.angResDeg[j];
			}

			for (int k = 0; k < currentFilter.nrPeaks; k++){
				name += ","+currentFilter.angBtwPeakDeg[k];
			}

			viz.addSlice(name, currentFilter.plot(N));

		}

		for (int i = 0; i < radlConfs.size(); i++) {

			RadialConfiguration currentFilter = radlConfs.get(i);

			String name;
			name = "RCONF_"+i;

			for (int j = 0; j < currentFilter.nrRings; j++){
				name += ","+currentFilter.ringRes[j];
			}

			viz.addSlice(name, currentFilter.plot(N).getProcessor(1));

		}

		return viz;

	}

	public ImageProcessor plotOne(
								int idx,
                                int N
	)
	{

		// they are aligned circ+radl
		// idx covers both

		if (idx<circConfs.size()){
			// take from circular configuration features
			return circConfs.get(idx).plot(N);
		}
		else{
			// take from radial configuration features
			return radlConfs.get(idx-circConfs.size()).plot(N).getProcessor(1);
		}
//			new ImagePlus(name, currentFilter.plot()).show();
//			new ImagePlus(name, currentFilter.plotFilter()).show();

	}

    public void score(
								float[] val,
								float[] ang,
								float[] rad
	)
    {
        for (int i = 0; i < circConfs.size(); i++) {

			score[i] = circConfs.get(i).score(val, ang, rad);

        }

		for (int i = circConfs.size(); i < (circConfs.size()+radlConfs.size()); i++){

			score[i] = radlConfs.get(i-circConfs.size()).score(val, rad);

		}
    }

	public float[] score(
								float[] val,
								float[] ang,
								float[] rad,
								int[] filt_idx
	)
	{
		// calculate score only on selected features and give it as a raw output
		float[] out = new float[filt_idx.length];
		int cnt = 0;
		for (int i = 0; i < filt_idx.length; i++){


			if (filt_idx[i]<circConfs.size()){
				// circular
				out[cnt] = circConfs.get(filt_idx[i]).score(val, ang, rad);

			}
			else{
				// radial
				out[cnt] = radlConfs.get(filt_idx[i]-circConfs.size()).score(val, rad);
			}

			cnt++;


		}

		return out;

	}

//	public void initFilter(
//								  int length
//	)
//	{
//
//		for (int i = 0; i < circConfs.size(); i++){
//
//			circConfs.get(i).initFilter(length);
//
//		}
//
//		for (int i = 0; i < radlConfs.size(); i++){
//
//			radlConfs.get(i).initFilter(length);
//
//		}
//
//	}

}
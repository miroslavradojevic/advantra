{

		// create detection image
		int w = _critpoint_det.length;
		int h = _critpoint_det[0].length;

        // sto
		byte[] t = new byte[w*h];
		for (int x=0; x<w; x++) {
			for (int y=0; y<h; y++) {
				float curr_on = _critpoint_det[x][y];
                // threshold score on critpoint detection
                if (curr_on>=output_membership_th) t[y * w + x] = (byte) 255;
			}
		}

		// take detections (binary image), find connected regions
		ByteProcessor bp = new ByteProcessor(w, h, t);
        ImagePlus imp_thresholded_scores = new ImagePlus("", bp);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(imp_thresholded_scores, true);
		conn_reg.run("");
		//conn_reg.showLabels().show();

        if (save_midresults) {
			if (_choose_type == CritpointRegion.RegionType.END) {
                IJ.saveAs(imp_thresholded_scores,                                   "Tiff", midresults_dir+"ends_thr_"+det_radius+".tif");
                IJ.saveAs(new ImagePlus("", new FloatProcessor(_critpoint_det)),    "Tiff", midresults_dir+"ends_sco_"+det_radius+".tif");
            }
			if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) {
                IJ.saveAs(imp_thresholded_scores,                                   "Tiff", midresults_dir+"jun_thr_"+det_radius+".tif");
                IJ.saveAs(new ImagePlus("", new FloatProcessor(_critpoint_det)),    "Tiff", midresults_dir+"jun_sco_"+det_radius+".tif");
            }
		}

		if (_choose_type== CritpointRegion.RegionType.END) {// DEBUG: just thresholding and connected components will give a lot of noisy regions
			// export binary image
//			ByteProcessor spots = new ByteProcessor(_critpoint_det.length, _critpoint_det[0].length);
//			System.out.println("Endpoints... what happened? total " + regs.size() + "regions" + minSize + " <> " +maxSize);
//			int cnt = 0;
//			for (int i=0; i<regs.size(); i++) {
//				if (true) { // regs.get(i).size()>=minSize && regs.get(i).size()<=maxSize
//					cnt++;
//					float C=0;        	// score
//					for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region
//						int xcoord = regs.get(i).get(aa)[1];
//						int ycoord = regs.get(i).get(aa)[0];
//						spots.set(xcoord, ycoord, (byte)255);
//						C += _critpoint_det[xcoord][ycoord];
//					}
//					C /= regs.get(i).size();   // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)
//					//System.out.println("region elements:  " + regs.get(i).size() + " : " + C );
//				}
//			}
		}

		//System.out.println("loop regions.. " + minSize + " till " + maxSize);

        ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

		 // will be filled up for each SALIENT branch
		 // will be filled up for each SALIENT branch


					int icoord = Profiler2D.xy2i[xcoord][ycoord];

					if (icoord!= -1) {  // peaks_i[icoord] is not null
						// element of the region is in foreground, take its thetas (if they exist)
						for (int peak_idx = 0; peak_idx < PeakExtractor2D.peaks_i[icoord].length; peak_idx++) {

							// iccord is the location index which is eqivalent to branch index
							// check if it exists and if it exists check whether the branch is on
							int curr_peak_i = PeakExtractor2D.peaks_i[icoord][peak_idx];
							boolean curr_peak_on = FuzzyDetector2D.branch_score[icoord][peak_idx]>output_membership_th;

							if (curr_peak_i!=-1 && curr_peak_i!=-2 && curr_peak_on) { // indexes of the spatial locations corresponding to peaks

								int peak_x = PeakExtractor2D.i2xy[curr_peak_i][0]; // PeakExtractor2D stores spatial location of the follow-up points
								int peak_y = PeakExtractor2D.i2xy[curr_peak_i][1];

								float[] unit_vxy = new float[]{peak_x-Cx, peak_y-Cy};
								float norm_vxy = (float) Math.sqrt(Math.pow(unit_vxy[0],2)+Math.pow(unit_vxy[1],2));
								unit_vxy[0] = unit_vxy[0] / norm_vxy;
								unit_vxy[1] = unit_vxy[1] / norm_vxy;

								vxy.add(unit_vxy);

							}

						}

					}


	}

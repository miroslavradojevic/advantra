
		
		
		
		if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
			found = true;
			tp++;
			annots[b] = true;
		}

			annots = new boolean[reader.nodes.size()];

			

			// END
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_endpoints);
			annots = new boolean[reader.nodes.size()];

			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.END) {

					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<reader.nodes.size(); b++) {

						float bx = reader.nodes.get(b)[reader.XCOORD];
						float by = reader.nodes.get(b)[reader.YCOORD];
						float br = reader.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					if (!found) fp++;  // detected but was not in the list of annotated ones

				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += tp + ",\t" + fp + ",\t" + fn + ",\t";




			// CRS
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_crosssections);
			annots = new boolean[reader.nodes.size()];

			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.CROSS) {

					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<reader.nodes.size(); b++) {

						float bx = reader.nodes.get(b)[reader.XCOORD];
						float by = reader.nodes.get(b)[reader.YCOORD];
						float br = reader.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					if (!found) fp++;  // detected but was not in the list of annotated ones

				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += tp + ",\t" + fp + ",\t" + fn + ",\t";





			// JUN (BIF+CRS)
			tp = 0; fp = 0; fn = 0;
			ReadSWC readerBIF = new ReadSWC(image_dir+image_gndtth_bifurcations);
			ReadSWC readerCRS = new ReadSWC(image_dir+image_gndtth_crosssections);
			annots = new boolean[readerBIF.nodes.size() + readerCRS.nodes.size()];

			for (int a=0; a<detected_regions.size(); a++) {
				if (detected_regions.get(a).type== CritpointRegion.RegionType.BIF || detected_regions.get(a).type== CritpointRegion.RegionType.CROSS) {

					boolean found = false;

					float ax = detected_regions.get(a).centroid[0];
					float ay = detected_regions.get(a).centroid[1];
					float ar = detected_regions.get(a).radius;

					for (int b=0; b<readerBIF.nodes.size(); b++) {

						float bx = readerBIF.nodes.get(b)[reader.XCOORD];
						float by = readerBIF.nodes.get(b)[reader.YCOORD];
						float br = readerBIF.nodes.get(b)[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}
					}

					for (int b = readerBIF.nodes.size(); b<readerBIF.nodes.size()+readerCRS.nodes.size(); b++) {

						float bx = readerCRS.nodes.get(b-readerBIF.nodes.size())[reader.XCOORD];
						float by = readerCRS.nodes.get(b-readerBIF.nodes.size())[reader.YCOORD];
						float br = readerCRS.nodes.get(b-readerBIF.nodes.size())[reader.RADIUS];

						if (circlesOverlap(ax, ay, ar, bx, by, br)) {// two roi circles overlap
							found = true;
							tp++;
							annots[b] = true;
						}

					}

					if (!found) fp++;  // detected but was not in the list of annotated ones

				}
			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;
			eval_string += tp + ",\t" + fp + ",\t" + fn;



//		}   // incomplete annotation
//		else
//			eval_string += "\"" + image_name + "\",\t" +
//								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +    // bif
//								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +    // end
//								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +    // crs
//								   Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN;             // jun (bif+crs)

		logWriter.println(eval_string);
		IJ.log("" + eval_string);

		logWriter.close();

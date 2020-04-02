/*
 * search folder given as argument, look for DET.* folders
 * need to evaluate each detection inside DET.* with respect to the 
 * given annotations folder - there can be more than one annotation
 * search annotation folder and it's subfolders for matching names
 * evaluate the detection, append one more line to eval.csv 
 * for each annotation
 * example call:
 
	java -jar ~/ImageJ/ij.jar -ijpath ~/ImageJ/plugins/ -macro eval.ijm parent_folder,annotations_folder

	arguments are current folder and the folder with annotations(A corresponds to plant dataset)
 
 */
 
args = getArgument;
print(getArgument);
if (args=="") exit("needs arguments: detection_folder(where the ijm is), annotation_folder,margin");
 
t = split(args, ",");
if (t.length!=3) exit("need 3 arguments");

det_dir = t[0];
ann_dir = t[1];
marg=t[2];

if (!File.isDirectory(det_dir)) exit(det_dir + " is not a dir");
if (!File.isDirectory(ann_dir)) exit(ann_dir + " is not a dir");

setBatchMode(true);

/* paths to all swcs in annotation folder and its subfolders */
count_swc = 0;
countSwcFiles(ann_dir);

//exit("found " + count_swc + " annotations: ");

if (count_swc<=0) exit("no annotations");

gndtth_paths = newArray(count_swc);
gndtth_tags = newArray(count_swc);
count_swc = 0;
storeSwcFiles(ann_dir, gndtth_paths, gndtth_tags);

print("" + count_swc + " SWCs found in annotations");

for (i=0; i<gndtth_paths.length; i++) print(i + " : " + gndtth_tags[i] + " --> " + gndtth_paths[i]);

ls = getFileList(det_dir);

cnt = 0;

for (i=0; i<ls.length; i++) {
        if (startsWith(ls[i], "REC.")) {
        
			ls1 = getFileList(det_dir + ls[i]);
			// delete all existing eval.csv
			for (j=0; j<ls1.length; j++) {
				if (endsWith(ls1[j], ".csv")) {
						dummy=File.delete(det_dir + ls[i] + ls1[j]); // because it returns 1 if file was deleted
				}
			}
			// do evaluation on .swc detection files
			for (j=0; j<ls1.length; j++) {
				
				//print(det_dir + ls[i] + ls1[j]);
				//print(replace(File.getName(det_dir + ls[i] + ls1[j]), ".det", ""));
				
				if (endsWith(ls1[j], ".swc")) {

						detection_image_path = det_dir + ls[i] + ls1[j];
						//print(detection_image_path);
						detection_image_name = replace(File.getName(det_dir + ls[i] + ls1[j]), ".swc", "");
						
						// search the corresopnding ground truth list
						for (k=0; k<gndtth_paths.length; k++) {
							
							gndtth_name = replace(File.getName(gndtth_paths[k]), ".swc", "");
							
							if (detection_image_name==gndtth_name) {
								
cnt++;

								print(detection_image_name+" A.K.A. "+detection_image_path+" \n\tevaluate with " + gndtth_tags[k] + " ----> " + gndtth_paths[k]);
		

// print("select=" + detection_image_path +
// 	" gndtth_path=" + gndtth_paths[k] +
// 	" gndtth_tag="  + gndtth_tags[k] + 
// 	" dst=" + marg);

run("Evaluator2D_reconstruction", 
	"select=" + detection_image_path +
	" gndtth_path=" + gndtth_paths[k] +
	" gndtth_tag="  + gndtth_tags[k] + 
	" dst=" + marg
	);


							}
							
						} 
						
				}
			}
			
        }
}

print("total comparisons : " + cnt);



function countSwcFiles(dir) {
     list = getFileList(dir);
     for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/"))
           countSwcFiles(dir+list[i]);
        else
			if (endsWith(dir+list[i],".swc")) count_swc++;
	 }
}

function storeSwcFiles(dir, gndtth_paths, gndtth_tags) {
     list = getFileList(dir);
     for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/"))
           storeSwcFiles(dir+list[i], gndtth_paths, gndtth_tags);
        else
			if (endsWith(dir+list[i],".swc")) {
				//print(dir+list[i]);
				gndtth_paths[count_swc] = dir+list[i];
				gndtth_tags[count_swc] = File.getName(File.getParent(dir+list[i]));
				count_swc++;
			}
	 }
}


print("DONE!");

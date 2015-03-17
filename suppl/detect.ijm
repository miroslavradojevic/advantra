/*
 * search the folder given as argument, seek for .tif files
 * apply detection parameters for each file
 */

print("\n************\nstalking...\n************");

neuron_diameter = newArray(6, 8);
percentile      = newArray(70, 80, 90);
correlation_boundary = newArray(0.65, 0.75, 0.85);
std_angle_deg  	= newArray(60, 80);
std_gcsstd_pix 	= newArray(2, 4);
ni 				= newArray(50, 100);
search_mode 	= newArray("SOMA_AWAY", "HIGH_CORR"); // "SOMA_CLOSE",

ns = 50;
zdist = 2;

main_folder = getArgument;

if (main_folder=="") exit ("Need argument: folder with .tifs to detect (choose current folder)");
if (!File.isDirectory(main_folder)) exit("Argument is not a folder!");

print("directory:\t" + main_folder);

//setBatchMode(true);

t_start = getTime();

/* list all files and detect with parameter grid */
images = getFileList(main_folder);

count = 0;

for (i=0; i<images.length; i++) {
        
        if (endsWith(images[i], ".tif")) { // it is image
        
			/* loop all the parameters */
			for (i1=0; i1<neuron_diameter.length; i1++) {
			for (i2=0; i2<percentile.length; i2++) {
			for (i3=0; i3<correlation_boundary.length; i3++) {
			for (i4=0; i4<std_angle_deg.length; i4++) {
			for (i5=0; i5<std_gcsstd_pix.length; i5++) {
			for (i6=0; i6<ni.length; i6++) {
			for (i7=0; i7<search_mode.length; i7++) {
				
				/* run for each combination of parameters */
				arg = 
					"select="           	+main_folder+images[i]+
					" neuron_diameter=" 	+toString(neuron_diameter[i1])+
					" percentile="  		+toString(percentile[i2])+
					" correlation_boundary="+toString(correlation_boundary[i3])+
					" std_angle_deg="		+toString(std_angle_deg[i4])+
					" std_gcsstd_pix="		+toString(std_gcsstd_pix[i5])+
					" ni="					+toString(ni[i6])+
					" ns="      			+toString(ns)+
					" zdist="				+toString(zdist)+
					" search_mode="			+search_mode[i7]
					;
				//print(arg);
				count++;
				print("\n************\n" + count + "/288\n************\n");
				run( "Neuron Stalker", arg);
				// run("Close All");

			}
			}
			}
			}
			}
			}
			}

        }
           
}

t_end = getTime();
print("elapsed " + ((t_end-t_start)/60000) + " min.");
exit();
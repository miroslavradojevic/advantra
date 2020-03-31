snr      = 4;

// generate configurations by varying (p1, p2, p3) triangle configurations (p1+p2+p3=1)
pmax     = 0.6;    // Dmin was hardcoded to 3, corresponds to lowest pi
pmin = (1-pmax)/2;
N = 3; 				// sampling

p1 = newArray(pmax, pmin, pmin);
p2 = newArray(pmin, pmax, pmin);
p3 = newArray(pmin, pmin, pmax);
p0 = newArray(0.333, 0.333, 0.333);

v1 = newArray(0,0,0);
a = 0;
for (i=0; i<v1.length; i++) { 
	v1[i] = p1[i] - p3[i];
	a += pow(v1[i],2);
}
a = sqrt(a);

v2 = newArray(0,0,0);
for (i=0; i<v2.length; i++) {
	v2[i] = p2[i] - p1[i];
}

for (i=0; i<v2.length; i++) {
	v1[i] = v1[i] / a;
	v2[i] = v2[i] / a;
}

L = a / pow(2,N);

root_dir = getDirectory("home");
print(root_dir);

//root_dir = "/scratch/mradojevic/critpoint_tests/synth_ext/"; // where the experiment is on the cluster
//print(root_dir);

out_dir = root_dir+"SNR"+snr+File.separator;
File.makeDirectory(out_dir);
print("exporting to: " + out_dir);

time = getTime();

p = newArray(0, 0, 0);
cnt = 0;
for (i=0; i<=pow(2, N); i++) {
  for (j=0; j<=i; j++) {
    
    for (ii=0; ii<3; ii++) p[ii] = p3[ii] + i * L * v1[ii] + j * L * v2[ii];
    
    arg = 	" SNR="+toString(snr)+
			" p1="+toString(p[0])+
			" p2="+toString(p[1])+
			" p3="+toString(p[2])+
			" nr_imgs=1"+
			" out_dir="+out_dir;
    
    print(p[0]+", "+p[1]+", "+p[2]);
    print(arg);
    
    run("Generate Junctions", arg);
    
    cnt++;
    
  }
}

print("DONE, "+cnt+" created, pmax="+pmax+", N="+N);
time = (getTime() - time) / 60000.0;
print ("" + time + " min.");
print("DONE");

// not used
function inTriangle(px, py, p0x, p0y, p1x, p1y, p2x, p2y) {
      
      Area = 1/2*(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);
      
      s = 1/(2*Area)*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py);
      t = 1/(2*Area)*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py);
      
      if (s>0 && t>0 && (1-s-t)>0) {
		return true;
	  }
	  else {
		return false;
	  }

}

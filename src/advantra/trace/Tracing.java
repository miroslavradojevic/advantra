package advantra.trace;

public interface Tracing {

	static int 		N 				= 100;					// maximum number of estimations used for tracing one branch
	static int		CPU_NR 			= 7;
	static int 		N_orientations 	= CPU_NR*10; 
	static int 		MS_PTS			= CPU_NR*30;
	
//	static int 		STEPS_BACKWARDS 		= 10;
	static double	RADIUS_SCALE_BORDER		= 1.1;
	
	static double 	radius_init 	= 1.0;
	static double 	radius_step		= 1.0;
	static double 	radius_limit	= 5.0;

	public static double 	k				= 2.0; 			// k*					radius is the total radius of the hypothesis
	public static double	jump_ahead		= 1.2; 			// jump_ahead*			radius is the jump ahead when tracing
	
	static double	radius_std 		= 1.0; // voxel
	static double	direction_std	= 0.9; // rad
	
	// mean-shift detection 3d
	static int 		extract_sphere_resolution 		= 32;
	
	int 			MS_PTS_TH 						= 10; // depends on sensitivity
	int 			MS_MAX_ITER 					= 100;
	double 			MS_EPS 							= 0.0001; 
	double 			MS_NEIGHBOUR 					= 8;  //degs
	double  		MS_NEIGHBOUR_RAD				= (MS_NEIGHBOUR/180)*Math.PI;
	double 			MS_ANGLE_RANGE_DEG 				= 30;
	double 			MS_ANGLE_RANGE_RAD 				= (MS_ANGLE_RANGE_DEG/180)*Math.PI;
	
	// mean-shift in 2d 
	static int 		extract_plane_resolution 		= 32;
	
	
	
}

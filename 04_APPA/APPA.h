
// FlowTable Structure
struct	patch_object {

	int pch;//the ID in descending order, start from 0

	int	patchID;
	int  landID;
	int	 inx;
	double lat;
	double lon;
	double z;
	double acc_area;

	//remained acc_area
	double re_acc_area;

	double runoff_coefficient;
	
	//this has been read 
	int    neigh_num{};
	int    neigh_ID[8]{};
	double neigh_gamma[8]{};//neighbor
	int	   neigh_pch[8]{};


	//need to be init here
	int    above_num{};
	int    above_ID[8]{};
	double above_gamma[8]{};//neighbor
	int	   above_pch[8]{};

	//new added data 
	int basin_inx{}; // ID of basins units
	int sthread{}; //sub-basin thread
	
	// channel layer and thread
	// e.g. 1000100 means the 1000 layer and the 100 thread
	int cthread{}; 
	int channel_inx{};//ID of channel units
	int layer{};//1000 the layer
	int layer_thread{};//100 the corresponding thread

	//state of channel_inx
	int channel_state{};//1 transient, 2 static 0 not valued
	int lock_state{};//0 unlocked, 1 locked 

	// patch 1 and 2 flow to 3, 1 and 2 are two outlets and 3 is a convergent point
	int outlet_flag{};
	int converge_flag{};

	//steam length
	double longest_stream_length{};

	//accumulated runtime
	int channel_acc;//original value 
	int hill_acc;//original value 
	
	int re_hill_acc;//remained value 

	int re_channel_acc;//remained value 
	int re_channel_acc_backup;//the backup is used to store the re_channel_acc while trying different parameters: tchm
};

// status information
struct status_information {

	int patch_num{};
	int hillslope_num{};
	int channel_num{};

	double DF{};
	int	   subbasin_num{};//sub-basins
	

	double TOV_layers[1000] = {};// maximum threads: 1000
	double LAYER = {};


	//for each try
	int     channel_unit_num{};//channel unit num
	double	TOV{}; //partitioned channel time in this layer
	double  TOi[1000]{}; // occupied time for each thread
	double  TAi[1000]{}; // available time for each thread

	int    rest_num{};//rest channel num



	double ESRch_layers[1000]{}; 

	double ESR{};//speed up ratio
	double ESRhs{};
	double ESRch{};
	double TMSRch{};
	double TMSRhs{};
	double TMSR{};

};


//Functions

struct patch_object* read_flow_table(char *flowtable_in, struct status_information* INF);

void init_hillslope_parallel(struct patch_object *patch, struct status_information*INF);
void init_channel_parallel(struct patch_object *patch, struct status_information *INF);


void hillslope_parallel_all(
	struct patch_object* patch,
	int available_threads,
	double xita_s_range[3],
	int D_range[3],
	struct status_information* INF,
	bool verbose
);
void channel_parallel_all(
	struct patch_object* patch,
	int available_threads,
	double xita_g,
	double Tcg,
	struct status_information* INF,
	bool verbose
);

void hillslope_parallel(struct patch_object* patch, double xita_s, int D, int threads, struct status_information* INF);
void channel_parallel(struct patch_object* patch, double Tchm,            int threads, struct status_information* INF);
void update_channel_parallel_within_layer(struct patch_object* patch,  int layer, int init_flag, struct status_information* INF);

void print_threads_map(struct patch_object* patch, char* patch_in, char* threads_out, int patch_num, int threads, int print_flag, int tail_flag);

void print_information(struct status_information* INF, int available_threads);

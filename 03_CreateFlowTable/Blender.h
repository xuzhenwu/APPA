#include <fstream>

//BLENDER.H
#ifndef Blender_H_
#define Blender_H_
#define MAXS   200	// size of filename string
#define LEN_GRASSHEADER 12
#define LEN_ARCHEADER	12 

#define MIN_RISE        0.001
#define DtoR	        0.01745329 
#define PI              3.1415926
#define THREADS			2

// STRUCTURE DEFINITIONS
struct ID_struct {
	int	patch;
};

struct	adj_struct {

	int	patchID;
	int  landID;
	int	 inx;
	double perimeter;
	double slope;
	double gamma;
	double z;
	int i_value;//xu i_value for the code 1-8 representing relative position 
	int DFA;//
	int DFB;//
	int DF;//A FINAL DISTRIBUTED FLAG

	double Dis;
	double Ste;
	int subbasin_id{};

	double DF_ratio;//to give in this facet
	double distribute_ratio;//in this tri and the ratio of this patch and the next one is calculated as 1-distribute_ratio
	double max_slope;//in this tri
	struct adj_struct *next;
};

struct	flow_struct	{
	int	patchID;
	int	ID_order;
	int area;
	int land;
	double x;
	double y;
	double z;
	double flna;
	double total_gamma;
	double acc_area;
	double road_dist;
	double slope;
	double internal_slope;
	double path_length;
	int	   inflow_cnt;
	int    num_adjacent;
	int	   stream_inx;
	int    num_acc;//xu
	int    real_adjacent;//xu
	double stream_gamma;

	int    acc_patch_num[8];//xu
	double acc_max_slope[8];
	int    acc_max_num[8];
	double acc_gamma_ratio[8];
	
	int subbasin_id{};

	double distance_to_stream;
	double elevation_above_stream;
	int number;

	int DF_num;

	double Dts;
	double Sts;
	double Dtr;
	double Str;
	double TWI;


	struct ID_struct  stream_ID;
	struct adj_struct *adj_list;

	struct adj_struct *adj_list_acc[8];//four potential acc area in this small region

	struct adj_struct *adj_ptr;

};


struct adj_struct  *check_list(int patchID, int max, adj_struct *list)
{
	int i, fnd;
	fnd = 0;
	i = 0;

	while ((fnd == 0) && (i < max)) {
		if ((list->patchID == patchID))
			fnd = 1;
		else {
			if (list->next == NULL)
				i = max;
			else {
				i += 1;
				list = list->next;
			}
		}
	}
	return(list);
}

class CReadArcInfoImages{
	public:
		int maxr;
		int maxc;
		double  cell;
		double  xllcorner, yllcorner;
		char    fndem[MAXS], fntable[MAXS], fnpatch[MAXS], fnslope[MAXS];
		char	fnstream[MAXS], fnroads[MAXS], fnflna[MAXS], fnsubbasin[MAXS];
		double  *dem = { nullptr }, *slope = { nullptr }, *flna = { nullptr };
		int     *patch = { nullptr }, *stream = { nullptr }, *roads = { nullptr };
		int     *patch_pch = { nullptr }, *subbasin{nullptr};
		void    initialize();
		void	input_prompt(int , int , char *, char *, char *, char *, char *, char *, char *, char *, int );

	private:
		void    input_header(int *maxr, int *maxc, char *fndem, int arc_flag);
		void	input_ascii_int(int *, char *, int, int, int);
		void	input_ascii_double(double *, char *, int, int, int);
};


class CReadCommandline{
	public:
		int		   f_flag, fl_flag, fh_flag, vflag;
		int		   s_flag, r_flag, sc_flag, arc_flag;
		int		   st_flag,algorithm_flag;
		int		   prefix_flag;
		
		//parallel groups
		int		   rcboun[THREADS+1];//threads+1 boun since it starts from 1
		int		   pchboun[THREADS + 1];//threads+1 boun since it starts from 1
	

		double     width;
		double     cell;
		char       *input_prefix;
		void       Commandline(int, char *[]);
		
		bool		acc;
		bool		dts;
		bool		sts;
		bool		dtr;
		bool		str;
		bool		tga;//total gamma
		bool		twi;

		bool		pch;

	private:
		void initialize();
};

class CCreatPitsFile{
public:
	std::ofstream   out1;
	std::ofstream   out2;
	// filenames for each image and file
	char name1[MAXS];
	char name2[MAXS];
	void creatPitsFile(char*,char*,std::ofstream&,std::ofstream&,char*);
private:
	void initialize();
};


class CreatFlowTable{
public:
	double    *dem , *slope, *flna ;
	int       *patch , *stream , *roads ;
	flow_struct	*flow_table;
	int       num_patches;
    void      initialize();
	void      add_roads(flow_struct *flow_table, int num_patches, double cell);
	void      adjust_pit(flow_struct *flow_table, int curr, int edge_inx, double edge_elev, double cell);
	void      compute_drainage_density(flow_struct *flow_table, int num_patches, double cell);
	void      compute_upslope_area(flow_struct *flow_table, int num_patches, int rflag, double cell);
	void      find_max_flna(flow_struct *flow_table, int curr, int *str_inx);
	void      find_min_flna(flow_struct *flow_table, int curr, int *str_inx);
	int       find_stream(flow_struct *flow_table, int curr, int *str_inx);
	double    find_top(flow_struct *flow_table, int curr, double pit_elev, int *num_pit,int *upslope_list, int *edge_inx);
	int	      in_list(int inx, int *list, int num_pit);
	double    path_lengths(flow_struct *flow_table, int curr, int *str_inx);
	void      remove_pits(flow_struct *flow_table, int num_patches, int sc_flag, double cell, std::ofstream & f1);
	void      route_roads_to_patches(flow_struct *flow_table, int num_patches, double cell, int fl_flag);
	int		  sort_flow_table(int *patch_pch, double *dem, int *patch, int maxr, int maxc);
	void      zero_flow_table(flow_struct *flow_table, int maxr, int maxc);
	void      print_topo_indices(int num_patches, struct flow_struct * flow_table, int *patch, CReadCommandline *commandline, double cell, char *input_prefix, double width, int maxr, int maxc, int algorithm_flag);
	void      print_flow_table(int num_patches, struct flow_struct * flow_table, int sc_flag, double cell, char *input_prefix, double width, int algorithm_flag);
	void      print_drain_stats(int num_patches, struct flow_struct *flow_table, char *input_prefix);
};

void      build_flow_table(flow_struct *flow_table, double *dem, double *slope, int *patch, int *stream, int *roads, int *subbasin,double *flna, int maxr, int maxc, int f_flag, int sc_flag, double cell, double xllcorner, double yllcorner,
				int*patch_pch, int *rcboun, int thread_inx);
int	      check_neighbours(int er, int ec, int *patch, int *patch_pch, double *dem, int *stream,int *subbasin, flow_struct *flow_entry, int num_adj, int maxr, int maxc, int sc_flag, double cell);
void	  compute_gamma(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag,int*pch_boun,int thread_inx);
void      compute_gamma_D8(flow_struct *flow_table, int num_patches, int sc_flag, double cell,int algorithm_flag, int*pch_boun, int thread_inx);
void      compute_gamma_MD8(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx);
void      compute_gamma_MFD_MD(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx);
void      compute_gamma_D_inf(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx);
void      compute_gamma_RMD_inf(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx);
int       find_patch(int num_patches, flow_struct *flow_table, int patchID);

#endif
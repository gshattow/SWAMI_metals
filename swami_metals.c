#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_sort.h>
//#include <gsl/gsl_integration.h>
#include <time.h>
//#undef DEF_MIN_OSX_VERSION 
//#define DEF_MIN_OSX_VERSION "10.6" 
#define directory "./" //"~/Testing/DiffuseGas/SWAMI_II/"
#define grid_dir "diffuse_files/"
#define sage_directory "./millennium_" //"~/windsage/output/results/millennium_"
#define sage_name "model_z"
#define metal_dir "metals_files/"
#define SIM "MM"
#define MODEL "W33"

// check box_side
// check xmin, ymin, zmin
// replace ALL type, e.g diffuse
// replace ALL SIM

int nbins = 100;
double pcls3 = 270.;
double npcls;
double box_side = 62.5;
double particle_mass = 8.61e8;
double T0 = 1.0e4;
double mstars_min = 0.0e10;

int lines;
int first = 0;
int last = 7;
int fnr;


int tot_n_trees = 0;
int tot_n_gals = 0;
int gals = 0;


int snapnum;
int numsnaps = 64;


double beta_HI(double temp);
double beta_HeI(double temp);
double beta_HeII(double temp);

double alpha_HII(double temp);
double alpha_HeII(double temp);
double alpha_HeIII(double temp);

double xi_HeII(double temp);

double n_HI(double temp, double n_H);
double n_HII(double temp, double n_H);
double n_HeI(double temp, double n_He);
double n_HeII(double temp, double n_He);
double n_HeIII(double temp, double n_He);

static double diffuse_delta[100][100][100] = {0.};
static double diffuse_temp[100][100][100] = {0.} ;
static double diffuse_vel[100][100][100][3]  = {0.};
static double diffuse_vdisp[100][100][100]  = {0.};
static double diffuse_H[100][100][100]  = {0.};

static double diffuse_HI[100][100][100]  = {0.};
static double diffuse_HII[100][100][100]  = {0.};
static double diffuse_HeI[100][100][100]  = {0.};
static double diffuse_HeII[100][100][100]  = {0.};
static double diffuse_HeIII[100][100][100]  = {0.};

static double diffuse_grid[100][100][100] = {0.};
static double diffuse_grid_old[100][100][100] = {0.};
static double bound_grid[100][100][100] = {0.};
static double halo_grid[100][100][100] = {0.};
static double bound_baryons[100][100][100] = {0.};
static double infall_grid[100][100][100] = {0.};

static double reject_mvir[100][100][100] = {0.};
static double reject_baryons[100][100][100] = {0.};
static double reject_metals[100][100][100] = {0.};

static double stars[100][100][100] = {0.};
static double coldgas[100][100][100] = {0.};
static double hotgas[100][100][100] = {0.};
static int WHIM[100][100][100] = {0.};

static double zstars[100][100][100] = {0.};
static double zcoldgas[100][100][100] = {0.};
static double zhotgas[100][100][100] = {0.};

static double metals_grid[100][100][100] = {0.};
static double gas_ejected[100][100][100] = {0.};
static double metals_ejected[100][100][100] = {0.};

void initialize_grids();
void read_Gamma(char transition[], double gammaarr[]);
void countdata();
int read_halos();
void grid_halos();
void read_diffuse_grid();
void fraction_flow();
void fraction_absorbed();
void calculate_baryons();

void reset_grids();
void write_metals();
void write_metals_redshift();
void write_bound_redshift();
void write_gas();
void timestamp();


double xmin = 0.; //5.8;
double ymin = 0.; //2.5;
double zmin = 0.; //12.5;

double z[64] = {127.000, 79.998, 50.000, 30.000, 19.916, 18.244, 16.725, 15.343, 14.086, 12.941, 11.897, 
	10.944, 10.073, 9.278, 8.550, 7.883, 7.272, 6.712, 6.197, 5.724, 5.289, 4.888, 4.520, 4.179, 3.866, 3.576, 
		3.308, 3.060, 2.831, 2.619, 2.422, 2.239, 2.070, 1.913, 1.766, 1.630, 1.504, 1.386, 1.276, 1.173, 1.078, 
			0.989, 0.905, 0.828, 0.755, 0.687, 0.624, 0.564, 0.509, 0.457, 0.408, 0.362, 0.320, 0.280, 0.242, 
				0.208, 0.175, 0.144, 0.116, 0.089, 0.064, 0.041, 0.020, 0.000};


double Gamma_HI[64] = {0.};
double Gamma_HeI[64] = {0.};
double Gamma_HeII[64] = {0.};

double tot_mstars[64] = {0.};
double tot_mmstars[64] = {0.};
double tot_mejected[64] = {0.};
double tot_mmejected[64] = {0.};
double tot_mhotgas[64] = {0.};
double tot_mmhotgas[64] = {0.};
double tot_mcoldgas[64] = {0.};
double tot_mmcoldgas[64] = {0.};
double tot_Hdiffuse[64] = {0.};
double tot_WHIM[64] = {0.};
double tot_mWHIM[64] = {0.};
double bound_z[64] = {0.};
double halos_z[64] = {0.};
double bound_baryons_z[64] = {0.};



double previous_mass = 0.;
double halo_mass = 0.;
double new_mass = 0.;
double bound_mass = 0.;

struct data
{
  int   Type;
  long long   GalaxyIndex;
  int   HaloIndex;
  int   FOFHaloIndex;
  int   TreeIndex;
  
  int   SnapNum;
  int   CentralGal;
  float CentralMvir;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  int   dT;

  // properties of subhalo at the last time this galaxy was a central galaaxy 
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;   
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;

  // baryonic reservoirs 
  float ColdGas;
  float StellarMass;
  float BulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;

  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float MetalsBulgeMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICS;

  // to calculate magnitudes
  float SfrDisk;
  float SfrBulge;
  float SfrDiskZ;
  float SfrBulgeZ;
  
  // misc 
  float DiskScaleRadius;
  float Cooling;
  float Heating;
  float LastMajorMerger;
  float OutflowRate;

  //infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;
};




static struct data G[100000];




int main()
{
	timestamp();
	npcls = pow(pcls3, 3);
	initialize_grids();
	char transition[5];
	strcpy(transition, "HI");
	read_Gamma(transition, Gamma_HI);
	strcpy(transition, "HeI");
	read_Gamma(transition, Gamma_HeI);
	strcpy(transition, "HeII");
	read_Gamma(transition, Gamma_HeII);

	for (snapnum = 63; snapnum < numsnaps; snapnum++)
	{

 		read_diffuse_grid();
		timestamp();

		countdata();
		printf("Reading the data for redshift z = %1.2f...\n", z[snapnum]);
		lines = read_halos();
		timestamp();
		grid_halos();
		timestamp();
 		calculate_baryons();
 		timestamp();
		write_metals();
		timestamp();
		write_gas();
		previous_mass = halo_mass;
 		reset_grids();
 		printf("\n\n");
	}
	write_metals_redshift();
	write_bound_redshift();
	timestamp();
	return 0;
}

void initialize_grids()
{
	int ii, jj, kk;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++)
			{
				diffuse_grid_old[ii][jj][kk] = npcls/(nbins*nbins*nbins/1.)*particle_mass;
				metals_ejected[ii][jj][kk] = 0.;
				gas_ejected[ii][jj][kk] = 0.;
				metals_grid[ii][jj][kk] = 0.;
			}
		}	
	}
}

void read_Gamma(char transition[], double gammaarr[])
{
	int iz = 0;
	char gammafile[200];
	double gamma;
	sprintf(gammafile, "%s_%s.dat", "Gamma", transition);

	FILE * fh;
	fh = fopen(gammafile, "rt");
	if (fh == NULL)
	{
		printf("Can't open gamma file %s!\n", gammafile);
		perror("Error");
		system("pwd");
		system("ls -laF Gamma_*");
	}
	else
	{
		printf("Opening gamma file: %s\n", gammafile);
	}
	for ( ; ; )
	{
		if (feof(fh))
		{
			break;
		}
		else
		{
			fscanf(fh, "%lf\n", &gamma);
			if (gamma > 0.) {gammaarr[iz] = gamma*1.0e12;}
			iz ++;
		}
	}
	fclose(fh);
	


}

void countdata()
{
	tot_n_gals = 0;
	tot_n_trees = 0;

	for (fnr = first; fnr < last + 1; fnr ++)
	{
		char infl[200];
		sprintf(infl, "millennium_%s/%s%0.3f_%d", MODEL, sage_name, z[snapnum], fnr);
		
		/* Open the file of data */
		FILE *fp;
		fp = fopen(infl, "rb");
		if (fp == NULL)
		{
			printf("Can't open file %s!\n", infl);
		}

		int numtrees;
		int numfile;
		fread(&numtrees, sizeof(int), 1, fp);
		fread(&numfile, sizeof(int), 1, fp);
		int galspertree[numtrees];
		fread(&galspertree, sizeof(int), numtrees, fp);

		tot_n_trees += numtrees;
		tot_n_gals += numfile;

		if (numfile == 0) {printf("File %d has no galaxies!\n", fnr);}


		fclose(fp);
		
	}
	printf("Total of %d galaxies found in %d trees.\n", tot_n_gals, tot_n_trees);

}

int read_halos()
{

	int numgal = 0;

	
	int offset = 0;
	tot_n_gals = 0;
	tot_n_trees = 0;
	double zej = 0.;

	for (fnr = first; fnr < last + 1; fnr ++)
	{

		char infl[200];
		sprintf(infl, "millennium_%s/%s%0.3f_%d", MODEL, sage_name, z[snapnum], fnr);

//		infile(infl);
		printf(" %d\t", fnr);
			

		FILE *fp;
		fp = fopen(infl, "rb");
		if (fp == NULL)
		{
			printf("Can't open file!\n");
		}
	
		int numtrees;
		int numfile;
		fread(&numtrees, sizeof(int), 1, fp);
		fread(&numfile, sizeof(int), 1, fp);
		int galspertree[numtrees];
		fread(&galspertree, sizeof(int), numtrees, fp);
		int len = numfile;
		tot_n_trees += numtrees;
		tot_n_gals += numfile;
		
		if (numfile == 0) {printf("File %d has no galaxies!\n", fnr);}
		else {printf("File %d has %d galaxies!\n", fnr, len);}


		struct data *GG = (struct data *) malloc(len * sizeof( struct data ));	
		if (GG == NULL) { printf("Uh Oh...\n"); }

		lines = 0;
		fread(GG, sizeof(struct data), len, fp);

		fclose(fp);
		
		int gg = 0;
		int file_off = 0;

		for (gg = 1; gg < len; gg ++)
		{
			G[offset + gg] = GG[gg];		
			file_off ++;
			gals ++;
			if (GG[gg].Type > 2) {printf("Galaxy %d has type %d.\n", offset+gg, GG[gg].Type);}
			if (GG[gg].MetalsEjectedMass > 0.) 
			{
//				printf("Galaxy %d has Mejected %.2f.\n", 
//					offset+gg, GG[gg].MetalsEjectedMass*1.0e10);
				zej += GG[gg].MetalsEjectedMass*1.0e10;
			}

		}

		free(GG);
		offset += file_off;
		
		printf("%d\t", offset);


	}

	lines = offset;


	printf("\n%d galaxies in mass range out of %d galaxies in %d trees total.\n", offset, tot_n_gals, tot_n_trees);
	printf("Galaxy 1: (Stars, CGas, Hotgas) = (%.2f, %.2f, %.2f)\n", 
		G[1].StellarMass, G[1].ColdGas, G[1].HotGas);
	printf("Total ejected metals = %.2f\n", zej);

	return lines;
	
	
}

void grid_halos()
{
	int gg = 0;
	int nhalos = 0;
	double total_mass = 0.;
	double total_metals = 0.;
	double total_bound = 0.;
	double mstars_max = 0.;
	double min_mstars = mstars_min; ///pow(1. + z[snapnum], 3);

	for (gg = 0; gg < lines; gg ++)
	{
		double xpos, ypos, zpos;
		double mvir, infall, otherb;
		double mstars, mmstars;
		double mejected, mmejected;
		double mhotgas, mmhotgas;
		double mcoldgas, mmcoldgas;
		int xbin, ybin, zbin;

		xpos = G[gg].Pos[0];
		ypos = G[gg].Pos[1];
		zpos = G[gg].Pos[2];
		
		xbin = (int)floor((xpos-xmin)/box_side*nbins);
		ybin = (int)floor((ypos-ymin)/box_side*nbins);
		zbin = (int)floor((zpos-zmin)/box_side*nbins);

		mvir = G[gg].Mvir*1.0e10;
		mstars = G[gg].StellarMass*1.0e10;
		mcoldgas = G[gg].ColdGas*1.0e10;
		mhotgas = G[gg].HotGas*1.0e10;
		mejected = G[gg].EjectedMass*1.0e10;
		infall = G[gg].infallMvir*1.0e10;
		otherb = G[gg].BlackHoleMass*1.0e10 + G[gg].ICS*1.0e10;
		
		mmstars = G[gg].MetalsStellarMass*1.0e10;
		mmcoldgas = G[gg].MetalsColdGas*1.0e10;
		mmhotgas = G[gg].MetalsHotGas*1.0e10;
		mmejected = G[gg].MetalsEjectedMass*1.0e10;


		double bb = mstars + mhotgas + mcoldgas + otherb;

		double temp = 14.4*diffuse_vdisp[xbin][ybin][zbin];
		
		halo_grid[xbin][ybin][zbin] += mvir;
		infall_grid[xbin][ybin][zbin] += infall;
		metals_ejected[xbin][ybin][zbin] += mmejected;
		gas_ejected[xbin][ybin][zbin] += mejected;

		stars[xbin][ybin][zbin] += mstars;
		zstars[xbin][ybin][zbin] +=mmstars;
		coldgas[xbin][ybin][zbin] += mcoldgas;
		zcoldgas[xbin][ybin][zbin] +=mmcoldgas;
		hotgas[xbin][ybin][zbin] += mhotgas;
		zhotgas[xbin][ybin][zbin] +=mmhotgas;

		bound_baryons[xbin][ybin][zbin] += bb;

		tot_mstars[snapnum] += mstars;
		tot_mmstars[snapnum] += mmstars;
		tot_mejected[snapnum] += mejected;
		tot_mmejected[snapnum] += mmejected;
		if (temp > 1.0e5) { 
			tot_WHIM[snapnum] += mejected;
			tot_mWHIM[snapnum] += mmejected;
		}
		tot_mhotgas[snapnum] += mhotgas;
		tot_mmhotgas[snapnum] += mmhotgas;
		tot_mcoldgas[snapnum] += mcoldgas;
		tot_mmcoldgas[snapnum] += mmcoldgas;

		nhalos ++;
		total_mass += mvir;
		total_metals += mmejected;
		total_bound += bb;
		if (mstars_max < mstars) {mstars_max = mstars;}

	}

	halo_mass = total_mass;
	halos_z[snapnum] = total_mass;
	bound_baryons_z[snapnum] = total_bound;
 	printf("%d halos read for snapshot %d. Largest = %lf 1e10 Msun\n", 
 		nhalos, snapnum, mstars_max/1.0e10);
 	printf("Total bound baryon fraction= %.1f/%.1f M_sun. (%f)\n", 
 		total_bound*1.0e-10, halo_mass*1.0e-10, total_bound/total_mass);
 	total_bound = 0.;
 	printf("Total ejected metal mass: %.1f e10 M_sun.\n", total_metals*1.0e-10);
}

void read_diffuse_grid()
{
	int Pcls = 0;
	double total_mass = 0.;
	double newly_bound_mass = 0.;
	double bound_pcl_mass = 0.;
	int empty = 0;
	int max_Nnb = 0.;
	double ave_vdisp = 0.;
	int min_Nnb = 1.;
	char grid_file[200];
	double bound_v_halo = 0.;
	double total_bound = 0.;
	sprintf(grid_file, "%s/%s%s_%03d", SIM, grid_dir, "diffuse", snapnum);

	FILE *fg;
	fg = fopen(grid_file, "r");
	if (fg == NULL)
	{
		printf("Can't open grid file %s!\n", grid_file);
	}
	else
	{
//		printf("Opening file: %s\n", grid_file);
	}
	for ( ; ; )
	{
		if (feof(fg))
		{
			break;
		}
		else
		{
			int xbin, ybin, zbin;
			int Npcls, Nnb, Nb;
			double xvel, yvel, zvel;
			double xdisp, ydisp, zdisp;
			fscanf(fg,"%d\t %d\t %d\t %d\t %d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", 
				&xbin, &ybin, &zbin, &Nnb, &Nb, &xvel, &yvel, &zvel, &xdisp, &ydisp, &zdisp);

			diffuse_grid[xbin][ybin][zbin] = Nnb*particle_mass;
			bound_grid[xbin][ybin][zbin] = Nb*particle_mass;

			diffuse_vel[xbin][ybin][zbin][0] = xvel;
			diffuse_vel[xbin][ybin][zbin][1] = yvel;
			diffuse_vel[xbin][ybin][zbin][2] = zvel;
			diffuse_vdisp[xbin][ybin][zbin] = (xdisp + ydisp + zdisp);

			total_bound += Nb*particle_mass; 
			bound_v_halo += Nb*particle_mass; 
			bound_v_halo -= halo_grid[xbin][ybin][zbin];
			Pcls += (Nnb + Nb);
			total_mass += Nnb*particle_mass;
			ave_vdisp +=diffuse_vdisp[xbin][ybin][zbin];
			if (Nnb == 0) {empty ++;}
		}
	}
	fclose(fg);
	ave_vdisp /= (1.0e6);
	bound_v_halo /= (1.0e16);
	bound_z[snapnum] = total_bound;

	printf("%d particles total (%.0f)^3.\n", Pcls, cbrt((float)Pcls)); 
	printf("Average vel disp = %f (T = %f).\n", ave_vdisp, ave_vdisp*14.4);
	printf("Ave bound v halo = %f (1e10 Msun).\n", bound_v_halo);
	printf("%d Empty cells (diffuse).\n", empty);

}

void calculate_baryons()
{
	int ii, jj, kk;
	double cell_mass = npcls*particle_mass/(nbins*nbins*nbins);
	double ave_n_H = 0;
	double min_d = 100.;
	double max_d = 0.;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++)
			{
				double delta = (diffuse_grid[ii][jj][kk] + bound_grid[ii][jj][kk])/cell_mass - 1.;
//				double delta = (diffuse_grid[ii][jj][kk] + gas_ejected[ii][jj][kk]/0.17)/cell_mass - 1.;
				if (min_d > delta) {min_d = delta;}
				if (max_d < delta) {max_d = delta;}
				double temp = 14.4*diffuse_vdisp[ii][jj][kk];
				if (temp < 3.) {temp = 3.;}
				if (temp > 1.0e5) {WHIM[ii][jj][kk] = 1;}
				else {WHIM[ii][jj][kk] = 0;}
				double n_H = 1.67*1.0e-7*pow(1. + z[snapnum], 3.)*(1. + delta);
				double n_He = 1./12. * n_H;
//				double volume_fraction = diffuse_grid[ii][jj][kk]/(diffuse_grid[ii][jj][kk] + bound_grid[ii][jj][kk])*0.1;
//				if (volume_fraction > 0.) {n_H = n_H*volume_fraction;}
//				else {n_H = 0.;}
				ave_n_H += n_H;
				diffuse_delta[ii][jj][kk] = delta;
				diffuse_temp[ii][jj][kk] = temp;
				diffuse_H[ii][jj][kk] = n_H*2.47e16;
				if (diffuse_H[ii][jj][kk] < bound_baryons[ii][jj][kk]) { diffuse_H[ii][jj][kk] = bound_baryons[ii][jj][kk]; }
				diffuse_HI[ii][jj][kk] = n_HI(temp, n_H)*2.47e16;
				diffuse_HII[ii][jj][kk] = n_HII(temp, n_H)*2.47e16;
				diffuse_HeI[ii][jj][kk] = n_HeI(temp, n_He)*2.47e16*4.;
				diffuse_HeII[ii][jj][kk] = n_HeII(temp, n_He)*2.47e16*4.;
				diffuse_HeIII[ii][jj][kk] = n_HeIII(temp, n_He)*2.47e16*4.;
				tot_Hdiffuse[snapnum] += (diffuse_H[ii][jj][kk] + reject_baryons[ii][jj][kk]);

				if (temp > 1.0e5) 
				{
					tot_WHIM[snapnum] += diffuse_H[ii][jj][kk] + reject_baryons[ii][jj][kk];
//					WHIM[ii][jj][kk] += n_H*2.47e16;
				}

			}
		}
	}
	ave_n_H = ave_n_H/pow(nbins, 3);
	printf("delta range: %f - %f.\n", min_d, max_d);
	printf("average n_H = %f e-6 (%f).\n", ave_n_H*1.0e6, ave_n_H*2.47e16*1.0e-10); 
}

void reset_grids()
{
	int ii, jj, kk;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++)
			{
				diffuse_grid_old[ii][jj][kk] = diffuse_grid[ii][jj][kk];
				if (diffuse_grid_old[ii][jj][kk] < 0.) {diffuse_grid_old[ii][jj][kk] = 0.;}
				diffuse_grid[ii][jj][kk] = 0.;
				diffuse_H[ii][jj][kk] = 0.;
				diffuse_HI[ii][jj][kk] = 0.;
				bound_grid[ii][jj][kk] = 0.;
				halo_grid[ii][jj][kk] = 0.;
				infall_grid[ii][jj][kk] = 0.;
				gas_ejected[ii][jj][kk] = 0.;
				metals_ejected[ii][jj][kk] = 0.;
				bound_baryons[ii][jj][kk] = 0.;
				reject_mvir[ii][jj][kk] = 0.;
				reject_baryons[ii][jj][kk] = 0.;
				reject_metals[ii][jj][kk] = 0.;
				stars[ii][jj][kk] = 0.;
				coldgas[ii][jj][kk] = 0.;
				hotgas[ii][jj][kk] = 0.;
				zstars[ii][jj][kk] = 0.;
				zcoldgas[ii][jj][kk] = 0.;
				zhotgas[ii][jj][kk] = 0.;
				WHIM[ii][jj][kk] = 0.;

			}
		}	
	}
	
}

void write_gas()
{
	char gas_file[200];


	sprintf(gas_file, "millennium_%s/%s%s_%03d", MODEL, metal_dir, "gas", snapnum);
 	FILE *fm;
 	fm = fopen(gas_file, "w");
 	if (fm == NULL) { printf("Can't open outfile %s!\n", gas_file); }

	int ii, jj, kk;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++)
			{
 				fprintf(fm, "%d\t %d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", ii, jj, kk, 
 					diffuse_HI[ii][jj][kk], diffuse_HII[ii][jj][kk], diffuse_HeI[ii][jj][kk], diffuse_HeII[ii][jj][kk],
 					diffuse_HeIII[ii][jj][kk], metals_ejected[ii][jj][kk], diffuse_temp[ii][jj][kk], bound_baryons[ii][jj][kk]);
			}
		}	
	}
	printf("Printed to file %s.\n", gas_file);
}

void write_metals()
{
	char metal_file[200];

	double gross_bound = 0.;
	double gross_metals = 0.;
	double gross_diffuse = 0.;
	double gross_diffuse_H = 0.;

	sprintf(metal_file, "millennium_%s/%s%s_%03d", MODEL, metal_dir, "metals", snapnum);
 	FILE *fm;
 	fm = fopen(metal_file, "w");
 	if (fm == NULL) { printf("Can't open outfile %s!\n", metal_file); }

	int ii, jj, kk;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++)
			{
				double metals = metals_ejected[ii][jj][kk];
 				fprintf(fm, "%d\t %d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %d\t %f\t %f\t %f\n", 
 					ii, jj, kk, 
 					stars[ii][jj][kk], zstars[ii][jj][kk], 
 					coldgas[ii][jj][kk], zcoldgas[ii][jj][kk],
 					hotgas[ii][jj][kk], zhotgas[ii][jj][kk], 
 					diffuse_H[ii][jj][kk] + gas_ejected[ii][jj][kk] - bound_baryons[ii][jj][kk],  
 					metals_ejected[ii][jj][kk], 
 					WHIM[ii][jj][kk], diffuse_temp[ii][jj][kk], 
 					bound_baryons[ii][jj][kk], diffuse_HI[ii][jj][kk]);
 				gross_metals += metals;
 				gross_diffuse += diffuse_grid[ii][jj][kk];
 				gross_diffuse_H += diffuse_H[ii][jj][kk];
 				gross_diffuse_H += gas_ejected[ii][jj][kk];
 				gross_bound += bound_grid[ii][jj][kk];
			}
		}	
	}
	printf("Printed to file %s.\n", metal_file);
	printf("Total diffuse metallicity = %f/%f = %f (%f Zsun).\n", gross_metals*1.0e-10, gross_diffuse_H*1.0e-10, gross_metals/gross_diffuse_H, (gross_metals/gross_diffuse_H)/0.02);
	printf("Total metals in diffuse matter = %f 1e10.\n", gross_metals*1.0e-10);
	printf("Total diffuse matter = %f 1e10.\n", gross_diffuse*1.0e-10);
	printf("Total bound matter = %f 1e10.\n", gross_bound*1.0e-10);
}

void write_metals_redshift()
{
	char metal_file[200];
	sprintf(metal_file, "millennium_%s/%s%s_%03d", MODEL, metal_dir, "metals_all_z", 0);
 	FILE *fm;
 	fm = fopen(metal_file, "w");
 	if (fm == NULL) { printf("Can't open outfile %s!\n", metal_file); }
	int zz;
	for (zz = 0; zz < numsnaps; zz ++)
	{
			fprintf(fm, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", 
				z[zz], tot_mstars[zz], tot_mmstars[zz], tot_mejected[zz], tot_mmejected[zz], tot_mhotgas[zz], 
				tot_mmhotgas[zz], tot_mcoldgas[zz], tot_mmcoldgas[zz], tot_Hdiffuse[zz], tot_WHIM[zz], tot_mWHIM[zz]);
	}

	printf("Stellar metallicity = %f (%f Zsun).\n", tot_mmstars[63]/tot_mstars[63], (tot_mmstars[63]/tot_mstars[63])/0.02);
	printf("Ejected metallicity = %f (%f Zsun).\n", tot_mmejected[63]/tot_mejected[63], (tot_mmejected[63]/tot_mejected[63])/0.02);
	printf("Hot Gas metallicity = %f (%f Zsun).\n", tot_mmhotgas[63]/tot_mhotgas[63], (tot_mmhotgas[63]/tot_mhotgas[63])/0.02);
	printf("Cold Gas metallicity = %f (%f Zsun).\n", tot_mmcoldgas[63]/tot_mcoldgas[63], (tot_mmcoldgas[63]/tot_mcoldgas[63])/0.02);
}

void write_bound_redshift()
{
	char metal_file[200];
	sprintf(metal_file, "millennium_%s/%s%s_%03d", MODEL, metal_dir, "bound_all_z", 0);
 	FILE *fm;
 	fm = fopen(metal_file, "w");
 	if (fm == NULL) { printf("Can't open outfile %s!\n", metal_file); }
	int zz;
	for (zz = 0; zz < numsnaps; zz ++)
	{
			fprintf(fm, "%f\t %f\t %f\t %f\n", z[zz], bound_z[zz], halos_z[zz], bound_baryons_z[zz]);
	}

}


void timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    printf("%s",asctime( localtime(&ltime) ) );
}

double beta_HI(double temp)
{
	return 5.85e-11 * pow(temp, 1./2.) / (1. + pow(temp/1.0e5, 1./2.)) * exp(-1.578e5/temp);
}
double beta_HeI(double temp)
{
	return 2.38e-11 * pow(temp, 1./2.) / (1. + pow(temp/1.0e5, 1./2.)) * exp(-2.853e5/temp);
}
double beta_HeII(double temp)
{
	return 5.68e-12 * pow(temp, 1./2.) / (1. + pow(temp/1.0e5, 1./2.)) * exp(-6.315e5/temp);
}

double alpha_HII(double temp)
{
	return 3.96e-13 * pow(temp/1.0e4, -0.7)/(1. + pow(temp/1.0e6, 0.7));
}
double alpha_HeII(double temp)
{
	return 4.31e-10 * pow(temp/1.0e4, -0.6353);
}
double alpha_HeIII(double temp)
{
	return 2.12e-12 * pow(temp/1.0e4, -0.7) / (1. + 0.379 * pow(temp/1.0e6, 0.7));
}

double xi_HeII(double temp)
{
	return 6.0e-10 * pow(temp/1.0e5, -1.5) * exp(-4.7e5/temp) * (1.0 + 0.3 * exp(-0.94e5/temp));
}

double n_HI(double temp, double n_H)
{
	return alpha_HII(temp) / (alpha_HII(temp) + beta_HI(temp) + Gamma_HI[snapnum]) * n_H;
}
double n_HII(double temp, double n_H)
{
	return n_H - alpha_HII(temp) / (alpha_HII(temp) + beta_HI(temp) + Gamma_HI[snapnum]) * n_H;
}
double n_HeI(double temp, double n_He)
{
	double ax_HeII = alpha_HeII(temp) + xi_HeII(temp);
	double a_HeII = alpha_HeII(temp);
	double a_HeIII = alpha_HeIII(temp);
	double b_HeI = beta_HeI(temp);
	double b_HeII = beta_HeII(temp);
	double G_HeI = Gamma_HeI[snapnum];
	double G_HeII = Gamma_HeII[snapnum];
	return ax_HeII * (b_HeI + a_HeIII) / 
	(ax_HeII * (a_HeIII + b_HeI + G_HeI) + (b_HeI + G_HeI) * (a_HeIII + b_HeII + G_HeII)) * n_He;
}
double n_HeII(double temp, double n_He)
{
	double ax_HeII = alpha_HeII(temp) + xi_HeII(temp);
	double a_HeII = alpha_HeII(temp);
	double a_HeIII = alpha_HeIII(temp);
	double b_HeI = beta_HeI(temp);
	double b_HeII = beta_HeII(temp);
	double G_HeI = Gamma_HeI[snapnum];
	double G_HeII = Gamma_HeII[snapnum];
	return  ((b_HeI + G_HeI) * (b_HeI - a_HeIII))/	(ax_HeII * a_HeIII - (b_HeI + G_HeI) * (ax_HeII  + a_HeIII + b_HeII + G_HeII)) * n_He;
}
double n_HeIII(double temp, double n_He)
{
	double ax_HeII = alpha_HeII(temp) + xi_HeII(temp);
	double a_HeII = alpha_HeII(temp);
	double a_HeIII = alpha_HeIII(temp);
	double b_HeI = beta_HeI(temp);
	double b_HeII = beta_HeII(temp);
	double G_HeI = Gamma_HeI[snapnum];
	double G_HeII = Gamma_HeII[snapnum];
	return (-b_HeI * (ax_HeII + b_HeI + G_HeI) + (b_HeI + G_HeI) * (ax_HeII + b_HeII + G_HeII))/
		(ax_HeII * (a_HeIII + b_HeI + G_HeI) + (b_HeI + G_HeI) * (a_HeIII + b_HeII + G_HeII)) * n_He;
}





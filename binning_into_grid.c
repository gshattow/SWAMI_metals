#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#define directory "/home/gshattow/projects/SWAMI_II/"
#define pcl_dir "/projects/p004_swin/abibiano/genev/minimill/"
#define grid_dir "MM/diffuse_files/"
#define file_name "snap"
#define unbound_file_name "unbound"

// check box_side
// check tot_pcls=X^3
// replace ALL type, e.g diffuse
// replace ALL SIM

int nbins = 100;
int pcls = 0;
int ub_pcls = 0;
float box_side = 62.5;
int tot_pcls = 270;

int fnr;
int snapnum;
int numsnaps = 64;


int were_unbound[100000000] = {0};
int now_unbound[100000000] = {1};
static int unbound[100][100][100];
static float vel[100][100][100][3];
static float vel_disp[100][100][100][3];
static int now_bound[100][100][100];
static int halos[100][100][100];

void bin_pcls();
int read_data();
int calculate_vdisp();
int unbound_ids();
char outfl[200];
void outfile();
void write_data();
void timestamp();

float xmin = 0.; //5.8;
float ymin = 0.; //2.5;
float zmin = 0.; //12.5;


int main()
{
	timestamp();
	for (snapnum = 63; snapnum < numsnaps; snapnum++)
	{
		int pp = unbound_ids();
		timestamp();
		int ub_particles = read_unbound();
		timestamp();
		int particles = read_data();
		timestamp();
		calculate_vdisp();
		timestamp();
		write_data();
		timestamp();
	}
	return 0;
}

int unbound_ids()
{

	int pp;
	for (pp = 0; pp < tot_pcls*tot_pcls*tot_pcls+1; pp ++)
	{
		were_unbound[pp] = now_unbound[pp];
		now_unbound[pp] = 0;

	}

	return pp;

}

int read_unbound()
{
		
	/* Open the file of data */
	char input_fname[200];
//	sprintf(input_fname, "%s%s%s_%03d", directory, pcl_dir, "diffuse_pcls", snapnum);
	sprintf(input_fname, "%s%s_%02d", pcl_dir, unbound_file_name, snapnum);
	printf("loading particle list... %s\n", input_fname);

	FILE *fp;
	fp = fopen(input_fname, "r");
	if (fp == NULL)
	{
		printf("Can't open file %s!\n", input_fname);
	}

	ub_pcls = 0;
	int maxID = 0;
	int minID = 5000;
	float x_max = 0.;
	float y_max = 0.;
	float z_max = 0.;

	for ( ; ; )
	{
		if (feof(fp))
		{
			break;
		}
		else
		{
			int dumID;
			float xpos, ypos, zpos;
			fscanf(fp,"%d,%f,%f,%f\n", 
				&dumID, &xpos, &ypos, &zpos);

			now_unbound[dumID] = 1;
			if (dumID < minID) {minID = dumID;}
			if (dumID > maxID) {maxID = dumID;}
			if (xpos > x_max) {x_max = xpos;}
			if (ypos > y_max) {y_max = ypos;}
			if (zpos > z_max) {z_max = zpos;}
			ub_pcls ++;
		}
	}
	fclose(fp);
	
	printf("minID, maxID, tot pcls: %d, %d, %d.\n", minID, maxID, tot_pcls*tot_pcls*tot_pcls); 
	printf("%d unbound particles.\n", ub_pcls); 
	printf("(xmax, ymax, zmax) = (%.2f, %.2f, %.2f)\n", x_max, y_max, z_max);
	
// 	int ii;
// 	for (ii = 0; ii < 100; ii ++)
// 	{
// 		printf("%d\t", now_unbound[ii]);
// 	}

	return ub_pcls;

}


int read_data()
{
	int nbins3 = nbins*nbins*nbins; 
		
	/* Open the file of data */
	char input_fname[200];
//	sprintf(input_fname, "%s%s%s_%03d", directory, pcl_dir, "diffuse_pcls", snapnum);
	sprintf(input_fname, "%s%s_%02d", pcl_dir, file_name, snapnum);
	printf("loading particle list... %s\n", input_fname);

	FILE *fp;
	fp = fopen(input_fname, "r");
	if (fp == NULL)
	{
		printf("Can't open file %s!\n", input_fname);
	}

	float x_max = 0.;
	float y_max = 0.;
	float z_max = 0.;
	pcls = 0;
	int maxID = 0;
	int minID = 5000;

	for ( ; ; )
	{
		if (feof(fp))
		{
			break;
		}
		else
		{
			int dumID, dumf;
			int xbin, ybin, zbin;
			float xpos, ypos, zpos;
			float xvel, yvel, zvel;
//			fscanf(fp,"%d\t %f\t %f\t %f\n", 
			fscanf(fp,"%d,%f,%f,%f,%f,%f,%f\n", 
				&dumID, &xpos, &ypos, &zpos, &xvel, &yvel, &zvel);

			xbin = (int)floor((xpos-xmin)/box_side*nbins);
			ybin = (int)floor((ypos-ymin)/box_side*nbins);
			zbin = (int)floor((zpos-zmin)/box_side*nbins);
			int index = xbin*nbins*nbins + ybin*nbins + zbin;


			if (now_unbound[dumID] == 1)
			{
				if (dumID < minID) {minID = dumID;}
				if (dumID > maxID) {maxID = dumID;}
				unbound[xbin][ybin][zbin] ++;
				vel[xbin][ybin][zbin][0] += xvel;
				vel[xbin][ybin][zbin][1] += yvel;
				vel[xbin][ybin][zbin][2] += zvel;
				
				if (xpos > x_max) {x_max = xpos;}
				if (ypos > y_max) {y_max = ypos;}
				if (zpos > z_max) {z_max = zpos;}
	
				pcls ++;
// 				if (were_unbound[dumID] == 0)
// 				{
// 					now_bound[xbin][ybin][zbin] --;
// 				}
				
			}
			else if (now_unbound[dumID] == 0) 
			{
				halos[xbin][ybin][zbin] ++;
				
// 				if (were_unbound[dumID] == 1) 
// 				{
// 					now_bound[xbin][ybin][zbin] ++;
// 				}
			}
		}
	}
	fclose(fp);
	
	printf("minID, maxID, tot pcls: %d, %d, %d.\n", minID, maxID, tot_pcls*tot_pcls*tot_pcls); 
	printf("\n%d particles in the file.\n", pcls);
	printf("(xmax, ymax, zmax) = (%.2f, %.2f, %.2f)\n", x_max, y_max, z_max);
	
	int ii, jj, kk, ll;
	int index = 0;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++ )
			{
				for (ll = 0; ll < 3; ll ++ )
				{
					vel[ii][jj][kk][ll] /= (float)unbound[ii][jj][kk];
					if (unbound[ii][jj][kk] == 0) {vel[ii][jj][kk][ll] = 0.;}
				}
			}
		}
	}


	return pcls;

}

int calculate_vdisp()
{
	int nbins3 = nbins*nbins*nbins; 
	/* Open the file of data */
	char input_fname[200];
//	sprintf(input_fname, "%s%s%s_%03d", directory, pcl_dir, "diffuse_pcls", snapnum);
	sprintf(input_fname, "%s%s_%02d", pcl_dir, file_name, snapnum);
	printf("calculating vdisp from particle list... %s\n", input_fname);

	FILE *fp;
	fp = fopen(input_fname, "r");
	if (fp == NULL)
	{
		printf("Can't open file %s!\n", input_fname);
	}


	for ( ; ; )
	{
		if (feof(fp))
		{
			break;
		}
		else
		{
			int dumID, dumf;
			int xbin, ybin, zbin;
			float xpos, ypos, zpos;
			float xvel, yvel, zvel;
//			fscanf(fp,"%d\t %f\t %f\t %f\n", 
			fscanf(fp,"%d,%f,%f,%f,%f,%f,%f\n", 
				&dumID, &xpos, &ypos, &zpos, &xvel, &yvel, &zvel);

			xbin = (int)floor((xpos-xmin)/box_side*nbins);
			ybin = (int)floor((ypos-ymin)/box_side*nbins);
			zbin = (int)floor((zpos-zmin)/box_side*nbins);


			if (now_unbound[dumID] == 1)
			{
				vel_disp[xbin][ybin][zbin][0] += pow(xvel - vel[xbin][ybin][zbin][0], 2);
				vel_disp[xbin][ybin][zbin][1] += pow(yvel - vel[xbin][ybin][zbin][1], 2);
				vel_disp[xbin][ybin][zbin][2] += pow(zvel - vel[xbin][ybin][zbin][2], 2);					
			}
		}
	}
	fclose(fp);
	
	
	int ii, jj, kk, ll;
	int index = 0;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++ )
			{
				for (ll = 0; ll < 3; ll ++ )
				{
					vel_disp[ii][jj][kk][ll] /= (float)unbound[ii][jj][kk];
					if (unbound[ii][jj][kk] == 0) {vel_disp[ii][jj][kk][ll] = 0.;}
				}
			}
		}
	}


	return pcls;

}


void write_data()
{
	char output_fname[200];
	sprintf(output_fname, "%s%s%s_%03d", directory, grid_dir, "diffuse", snapnum);


	FILE *fq;
	fq = fopen(output_fname, "w");
	if (fq == NULL) { printf("Can't open outfile %s!\n", output_fname); }
	int i = 0;
	
	int ii, jj, kk, ll;
	for (ii = 0; ii < nbins; ii ++)
	{
		for (jj = 0; jj < nbins; jj ++)
		{
			for (kk = 0; kk < nbins; kk ++)
			{
				fprintf(fq, "%d\t %d\t %d\t %d\t %d\t", 
					ii, jj, kk, unbound[ii][jj][kk], halos[ii][jj][kk]);
				for (ll = 0; ll < 3; ll ++)
				{
					fprintf(fq, "%f\t", vel[ii][jj][kk][ll]);
				}
				for (ll = 0; ll < 3; ll ++)
				{
					fprintf(fq, "%f\t", vel_disp[ii][jj][kk][ll]);
				}
				fprintf(fq, "\n");
				unbound[ii][jj][kk] = 0;
				halos[ii][jj][kk] = 0;
				for (ll = 0; ll < 3; ll ++ )
				{
					vel[ii][jj][kk][ll] = 0.;
					vel_disp[ii][jj][kk][ll] = 0.;
				}
				i ++;
			}
		}
	}
	printf("%d pos cells printed to file %s.\n", i, outfl);
	
	fclose(fq);

}

void timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    printf("%s",asctime( localtime(&ltime) ) );
}

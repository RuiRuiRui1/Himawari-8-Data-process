
#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include <vector>
#include <string.h>
#include <string>
#include <iostream>
#include <dirent.h>
#include <algorithm>
#include <ctime>
#include <thread>
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "gdalwarper.h"
using namespace std;

#define Pi 3.1415926
#define ALLPIXELNUM win_width*win_height
#define _10KMALLPIXELNUM win_width*win_height/25
#define MISSVALUE -32768
//for debugging
#define COUT_ALLVALUE(valuename,valuenum) for(int i=0;i<valuenum;i++){cout<<valuename[i]<<" ";}cout<<endl;
//define Resolution(unit:degree)
#define _2KM 0.02
#define _10KM 0.1


//filetype enum( NETCDF-->.nc ENVI***-->.img)
typedef enum FileType
{
	NETCDF =1,
	ENVIAHI =2,
	ENVIBT =3,
	ENVIANGLE=4,
	ENVIAHI10KM=5,
	ENVISOA10KM=6,
	ENVISZA10KM=7

}FILETYPE;

typedef enum ifCldmsk
{
	CLOUDMSK = 1,
	NOTCLOUDMSK = 0

}CLDMSK;


//rename filelist
typedef vector<string> FILELIST;
typedef vector<string>::iterator FILELIST_iter;

//file functions
bool FindFiles(vector<string>* FileNames, FILETYPE fileFormat, const char* path);
bool SaveImageToFile(const char *outputPath, int nWidth, int nHeight, float **data, const char* pszFormat, double adfGeoTransform[6], const char *pszSRS_WKT, int Bandnum);
bool SaveImageToFile(const char *outputPath, int nWidth, int nHeight, int **data, const char* pszFormat, double adfGeoTransform[6], const char *pszSRS_WKT, int Bandnum);
void FileDelect(const char* path);
void ReadIMG(const char *filePath,float **data, int bandnum);
void ReadIMG(string filePath, float **data, int bandnum);

//valuecopy functions---->copy value to all pixels
void valuecopy(float *input,float * output);
void valuecopy(int *input, int * output);
void valuecopy(int num, int *output);
void valuecopy(float num, float *output);
template<typename T>
void valuecopy(T *input, T *output);
template<typename T>
void valuecopy(T num,T *output);

//AOD workflow functions
void H8_read_cut(const char* inputDir, const char* outputTemp, float roi_xstart, float roi_ystart, float roi_xend, float roi_yend, vector<string>::iterator it, CLDMSK cldmsk);
void H8_cldmsk();
void H8_resample(const char* inputDir, const char* outputDir, float resulotion, int bandnum, vector<string>::iterator it);
void AOD_retrival();

//work functions
void cldmsk_land(float *TOAb1, float *BTb1, float *BTb4, float *BTb8, float *BTb9, int *msk);
void Absorb_Correction(int bandnum);

//input path
char input_dir[128]= "F:\\AHI\\";


//output path
char temp_dir[128]= "F:\\AHIout\\temp\\";
char AHI_TOA[128]= "F:\\AHIout\\toa\\";
char AHI_BT[128] = "F:\\AHIout\\bt\\";
char ANGLE_dir[128] = "F:\\AHIout\\angle\\";




//cloud mask bands ->TOA band1 and BT band1,4,8,9
float *TOAb1 = NULL;
float *BTb1 = NULL;
float *BTb4 = NULL;
float *BTb8 = NULL;
float *BTb9 = NULL;
//cloud mask result
int   *msk = NULL;

//subdatasets
float *SOA = NULL;//Ì«Ñô·½Î»½Ç(solar azimuth angle)
float *SZA = NULL;//Ì«ÑôÌì¶¥½Ç(solar zenith angle)
float *albedo = NULL;//AHI 6bands
float *tbb = NULL;//BT 10bands
float *SAA = NULL;//ÎÀÐÇ·½Î»½Ç(satellite azimuth angle)
float *SAZ = NULL;//ÎÀÐÇÌì¶¥½Ç(satellite zenith angle)

//absorb correction parameters
double *RTrans_H2O = NULL;
double *RTrans_O3 = NULL;
double *RTrans_other = NULL;



//cut window parameters
int win_width = 0;
int win_height = 0;
//time of data (eg. 20170401_0000)
string timestr;


int main(int argc,char* argv[])
{
	//clock for time spend
	clock_t start, finish;
	
	//roi parameters x->lon ¡ãE and y->lat ¡ãN
	float roi_xstart = 130;
	float roi_ystart = 60;
	float roi_xend = 140;
	float roi_yend = 50;

	//find files like temp
	
	//nc file list
	FILELIST H8_list;
	if (!FindFiles(&H8_list, NETCDF, input_dir))
	{
		cout << "NO .nc file exist!\n";
		return 0;
	}
	for (FILELIST_iter it = H8_list.begin(); it != H8_list.end(); it++)
	{
		
		cout << "AHI file:" << it->c_str() << endl;
		//read and cut ---->IDL P1+P2
		start = clock();
		cout << "----------H8 data read and cut-----(if cloud mask?)-----\n";
		H8_read_cut(input_dir, ANGLE_dir, roi_xstart, roi_ystart, roi_xend, roi_yend,it,CLOUDMSK);
	
		finish = clock();
		cout << "H8 read and cut time spent:" << (finish - start) << " ms" << endl << endl;

		//cloud mask ---->IDL P3
		/*start = clock();
		cout<<"----------H8 cloud mask----------\n";
		H8_cldmsk();
		cout<<"----------H8 cloud mask ok!----------";
		finish = clock();
		cout << "H8 cloud mask time spent:" << (finish - start) << " ms" << endl << endl;*/

		//all done
		CPLFree(SOA);
		CPLFree(SZA);
		CPLFree(SAA);
		CPLFree(SAZ);
		CPLFree(TOAb1);
		CPLFree(msk);
	
	}
	//AHI img file list
	FILELIST AHI_list;
	if (!FindFiles(&AHI_list, ENVIAHI, AHI_TOA))
	{
		cout<<"NO AHI file exist!\n";
		return 0;
	}
	for (FILELIST_iter it = AHI_list.begin(); it != AHI_list.end(); it++)
	{

		cout << "IMG file:" << it->c_str() << endl;

		start = clock();
		cout<<"----------H8 TOA resampling----------\n";
		H8_resample(AHI_TOA, AHI_TOA,_10KM,6,it);
		cout<<"Resolution:2KM---->10KM\n";
		cout<<"----------H8 TOA resampling ok!----------";
		finish = clock();
		cout << "H8 TOA resampling time spent:" << (finish - start) << " ms" << endl << endl;
	}
	//BT img file list
	FILELIST BT_list;
	if (!FindFiles(&BT_list, ENVIBT, AHI_BT))
	{
		cout << "NO BT file exist!\n";
		return 0;
	}
	for (FILELIST_iter it = BT_list.begin(); it != BT_list.end(); it++)
	{

		cout << "IMG file:" << it->c_str() << endl;

		start = clock();
		cout << "----------H8 BT resampling----------\n";
		H8_resample(AHI_BT, AHI_BT, _10KM, 10, it);
		cout << "Resolution:2KM---->10KM\n";
		cout << "----------H8 BT resampling ok!----------";
		finish = clock();
		cout << "H8 BT resampling time spent:" << (finish - start) << " ms" << endl << endl;
	}
	//ANGLE img file list
	FILELIST Angle_list;
	if (!FindFiles(&Angle_list, ENVIANGLE, ANGLE_dir))
	{
		cout << "NO BT file exist!\n";
		return 0;
	}
	for (FILELIST_iter it = Angle_list.begin(); it != Angle_list.end(); it++)
	{

		cout << "IMG file:" << it->c_str() << endl;

		start = clock();
		cout << "----------H8 ANGLE resampling----------\n";
		H8_resample(ANGLE_dir, ANGLE_dir, _10KM, 1, it);
		cout << "Resolution:2KM---->10KM\n";
		cout << "----------H8 ANGLE resampling ok!----------";
		finish = clock();
		cout << "H8 ANGLE resampling time spent:" << (finish - start) << " ms" << endl << endl;
	}

	//delect temp files
	FileDelect(temp_dir);
	cout<<"----------temp delect ok!----------\n";
	cout<<"###########   finished   ###########\n";
	getchar();
    return 0;
}

//inputDir and outputTemp need absolute path (eg. "F:\\AHI\\")
//roi_xstart etc. need upperleft and lowerright lat-->y lon-->x(eg.80,60,130,10 are 10-60¡ãN,80-130¡ãE)
//it--->fliename  cldmsk->if cloud mask (CLOUDMSK or NOTCLOUDEMSK)
void H8_read_cut(const char* inputDir, const char* outputTemp, float roi_xstart, float roi_ystart, float roi_xend, float roi_yend, vector<string>::iterator it, CLDMSK cldmsk)

{
		
	
	    //buffer-->full filename include path ;timestr-->data time (eg.20070401_0000);
	    char buffer[128];
	    strcpy_s(buffer, inputDir);
	    strcat_s(buffer, it->c_str());
		const char *pszSrcFile = buffer;
		
	     timestr = pszSrcFile;
		timestr=timestr.substr(14, 13);

		GDALAllRegister();

		//*pDataSet-->all subdatesets and metadata
		GDALDataset *pDataSet = (GDALDataset*)GDALOpen(pszSrcFile, GA_ReadOnly);
		if (pDataSet == NULL)
		{
			cout<<"File cannot open!";
		}
	  
		//get globle metadata
		const char * strupleft_lat = pDataSet->GetMetadataItem("upper_left_latitude");
		const char * strupleft_lon = pDataSet->GetMetadataItem("upper_left_longitude");
		const char * strunit_of_latlon = pDataSet->GetMetadataItem("grid_interval");
		const char * strlatnum = pDataSet->GetMetadataItem("line_number");
		const char * strlonnum = pDataSet->GetMetadataItem("pixel_number");
		float upleft_lat = atof(strupleft_lat);
		float upleft_lon = atof(strupleft_lon);
		float unit_of_latlon = atof(strunit_of_latlon);
		int latnum = atof(strlatnum);
		int lonnum = atof(strlonnum);

		//compute cut window
		roi_xstart = max(roi_xstart, upleft_lon);
		roi_ystart = min(roi_ystart, upleft_lat);
		roi_xend = min(roi_xend, upleft_lon+lonnum*unit_of_latlon);
		roi_yend = max(roi_yend, upleft_lat-latnum*unit_of_latlon);
		win_width = abs((roi_xend - roi_xstart) / unit_of_latlon);
		win_height= abs((roi_yend - roi_ystart) / unit_of_latlon);
		int win_upleft_lat=abs((roi_ystart- upleft_lat) / unit_of_latlon);
		int win_upleft_lon=abs((roi_xstart - upleft_lon) / unit_of_latlon);
		
		//set projection
		const char *pszSRS_WKT = NULL;

		double adfGeoTransform[6];
		pDataSet->GetGeoTransform(adfGeoTransform);
		adfGeoTransform[0] = roi_xstart;
		adfGeoTransform[3] = roi_ystart;
		adfGeoTransform[1] = 0;
		adfGeoTransform[5] = 0;
	    adfGeoTransform[2] = _2KM;   //ENVI img does not have rotational terms
		adfGeoTransform[4] = _2KM;

		//** papszSUBDATASETS-->all subdatasets
		char ** papszSUBDATASETS = GDALGetMetadata((GDALDatasetH)pDataSet, "SUBDATASETS");
		if (papszSUBDATASETS == NULL)
		{
			cout<<"Error: no subdataset!";
			return;
		}
		else
		{
			int iCount = CSLCount(papszSUBDATASETS);
			if (iCount <= 0)
			{
				GDALClose((GDALDriverH)pDataSet);
				return;
			}

			//read function
			GDALDataset * tmpdt;
			GDALRasterBand *poBand;
			//i for subdataset[num] (i=2*num)
			int i;

			
			//SOA: subdataset[3]
			i = 6;
			string tmpstr = string(papszSUBDATASETS[i]);
			tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);			
			const char *tmpc_str = tmpstr.c_str();
			tmpdt = (GDALDataset *)GDALOpen(tmpc_str, GA_ReadOnly);			
			poBand = tmpdt->GetRasterBand(1);

		

			SOA = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
			poBand->RasterIO(GF_Read, win_upleft_lon, win_upleft_lat, win_width, win_height, SOA, win_width, win_height, GDT_Float32, 0, 0);
			//change to float
			for (int i = 0; i < ALLPIXELNUM; i++)
			{
				if (SOA[i] == MISSVALUE)
				{
					continue;
				}
				SOA[i] *= 0.01f;

			}
		

			//SOZ: subdataset[4]
			i = 8;
			tmpstr = string(papszSUBDATASETS[i]);
			tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);
			
			tmpc_str = tmpstr.c_str();
			tmpdt = (GDALDataset *)GDALOpen(tmpc_str, GA_ReadOnly);
			
			poBand = tmpdt->GetRasterBand(1);		
			SZA = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
			float * miu = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);//for calculate
			poBand->RasterIO(GF_Read, win_upleft_lon, win_upleft_lat, win_width, win_height, SZA, win_width, win_height, GDT_Float32, 0, 0);
			for (int i = 0; i < ALLPIXELNUM; i++)//change to float
			{
				if (SZA[i] == MISSVALUE)
				{
					continue;
				}
				SZA[i] *= 0.01f;
				miu[i] = cos(SZA[i] * Pi / 180);

			}
		

			//SAA: subdataset[1]
			i = 2;
			tmpstr = string(papszSUBDATASETS[i]);
			tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);
			
			tmpc_str = tmpstr.c_str();
			tmpdt = (GDALDataset *)GDALOpen(tmpc_str, GA_ReadOnly);
			
			poBand = tmpdt->GetRasterBand(1);		
			SAA = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
			poBand->RasterIO(GF_Read, win_upleft_lon, win_upleft_lat, win_width, win_height, SAA, win_width, win_height, GDT_Float32, 0, 0);
			for (int i = 0; i < ALLPIXELNUM; i++)//change to float
			{
				if (SAA[i] == MISSVALUE)
				{
					continue;
				}
				SAA[i] *= 0.01;				

			}
		

			//SAZ: subdataset[2]
			i = 4;
			tmpstr = string(papszSUBDATASETS[i]);
			tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);
			
			tmpc_str = tmpstr.c_str();
			tmpdt = (GDALDataset *)GDALOpen(tmpc_str, GA_ReadOnly);
			
			poBand = tmpdt->GetRasterBand(1);
			SAZ = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
			poBand->RasterIO(GF_Read, win_upleft_lon, win_upleft_lat, win_width, win_height, SAZ, win_width, win_height, GDT_Float32, 0, 0);
			for (int i = 0; i < ALLPIXELNUM; i++)//change to float
			{
				if (SAZ[i] == MISSVALUE)
				{
					continue;
				}
				SAZ[i] *= 0.01;
			}
		
			
		    //create 6 bands TOA img
			GDALDataset *TOA;
			GDALRasterBand *toaBand;

			char **papszOptions = NULL;
			char tmpout[128];
			strcpy_s(tmpout, AHI_TOA);
			strcat_s(tmpout, "HIMAWARI8_AHI_");
			strcat_s(tmpout, timestr.c_str());
			strcat_s(tmpout, "_2KM.img");
			GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("ENVI");
			TOA = poDriver->Create(tmpout, win_width, win_height, 6, GDT_Float32, papszOptions);
						
			//read albedo (subdataset[5]---subdateset[10])and calculate TOA
			for (int j=10; j < 21;j+=2)
			{
				tmpstr = string(papszSUBDATASETS[j]);
				tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);
				
				tmpc_str = tmpstr.c_str();
				tmpdt = (GDALDataset *)GDALOpen(tmpc_str, GA_ReadOnly);				
				poBand = tmpdt->GetRasterBand(1);		
				albedo = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
				poBand->RasterIO(GF_Read, win_upleft_lon, win_upleft_lat, win_width, win_height, albedo, win_width, win_height, GDT_Float32, 0, 0);
				//calculate TOA
				for (int i = 0; i < ALLPIXELNUM; i++)
				{
					if (albedo[i] == MISSVALUE)
					{
						continue;
					}
					albedo[i] *= 0.0001;
				    albedo[i] = albedo[i]/miu[i];
				
				}
				if (j==10)//export TOA band1 for cloud mask
				{
					TOAb1= (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
					valuecopy(albedo, TOAb1);
					H8_cldmsk();
					/*thread trd (H8_cldmsk);
					trd.join();*/
					
				}		
				//TOA cloud mask				
				if (cldmsk == CLOUDMSK)
				{
					for (int i = 0; i < ALLPIXELNUM; i++)
					{
						albedo[i] = albedo[i] * msk[i];
					
					}

				}
				//write band 1-6
				toaBand =TOA->GetRasterBand(j/2-4);
				toaBand->RasterIO(GF_Write, 0, 0, win_width, win_height, albedo, win_width, win_height, GDT_Float32, 0, 0);			
				CPLFree(albedo);				
			}
			TOA->SetGeoTransform(adfGeoTransform);
			TOA->SetProjection(pszSRS_WKT);
			GDALClose(TOA);
			cout<<"*******TOA save ok!*******\n";

			//create 10 bands BT img

			//get parameters:sf-->scale factor; ao-->add offset; ms-->missing value
			const char * strsf = pDataSet->GetMetadataItem("tbb_07_scale_factor");
			const char * strao = pDataSet->GetMetadataItem("tbb_07_add_offset");
			const char * strms = pDataSet->GetMetadataItem("tbb_07_missing_value");
			float sf = atof(strsf);
		    float ao = atof(strao);
			float ms = atof(strms);

			GDALDataset *BT;
			GDALRasterBand *btBand;

			strcpy_s(tmpout, AHI_BT);
			strcat_s(tmpout, "HIMAWARI8_BT_");
			strcat_s(tmpout, timestr.c_str());
			strcat_s(tmpout, "_2KM.img");
			
			BT = poDriver->Create(tmpout, win_width, win_height, 10, GDT_Float32, papszOptions);
			//read tbb(subdataset[12]---subdateset[21]) and calculate BT
			for (int j = 24; j < 43; j += 2)
			{
				tmpstr = string(papszSUBDATASETS[j]);
				tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);
				
				tmpc_str = tmpstr.c_str();
				tmpdt = (GDALDataset *)GDALOpen(tmpc_str, GA_ReadOnly);
			
				poBand = tmpdt->GetRasterBand(1);
				tbb = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
				poBand->RasterIO(GF_Read, win_upleft_lon, win_upleft_lat, win_width, win_height, tbb, win_width, win_height, GDT_Float32, 0, 0);
				//calculate BT
				for (int i = 0; i < ALLPIXELNUM; i++)
				{
					if (tbb[i] == MISSVALUE)
					{
						continue;
					}
					tbb[i] = tbb[i]*sf+ao;
				}
				//BT band1,4,8,9 for cloud mask ( no use£¡£¡£¡£©
			/*	switch (j)
				{
				case 24:
				{
					BTb1 = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
					valuecopy(tbb, BTb1);
					break;
				}
				case 30:
				{
					BTb4 = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
					valuecopy(tbb, BTb4);
					break;
				}
				case 38:
				{
					BTb8 = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
					valuecopy(tbb, BTb8);
					break;
				}
				case 40:
				{
					BTb9 = (float*)CPLMalloc(sizeof(float)*ALLPIXELNUM);
					valuecopy(tbb,BTb9);
					break;
				}
				
				}
				*/
				//BT cloud mask
				if (cldmsk == CLOUDMSK)
				{
					for (int i = 0; i < ALLPIXELNUM; i++)
					{
						tbb[i] = tbb[i] * msk[i];
					}

				}
				//write band 1-10
				btBand = BT->GetRasterBand(j / 2 - 11);
				btBand->RasterIO(GF_Write, 0, 0, win_width, win_height, tbb, win_width, win_height, GDT_Float32, 0, 0);
				CPLFree(tbb);
			}
			BT->SetGeoTransform(adfGeoTransform);
			BT->SetProjection(pszSRS_WKT);
			GDALClose(BT);
			cout<<"*******BT save ok!*******\n";
		

		    //save function SOA SZA SAA SAZ
		

			//save SOA
			char out_SOA[128];
			strcpy_s(out_SOA, outputTemp);					
			strcat_s(out_SOA, "HIMAWARI8_");
			strcat_s(out_SOA, timestr.c_str());
			strcat_s(out_SOA, "_SOA_2KM.img");
			SaveImageToFile(out_SOA, win_width, win_height, &SOA, "ENVI", adfGeoTransform, pszSRS_WKT,1);
			
			cout<<"*******SOA save ok!*******\n";

			//save SZA
			char out_SZA[128];
			strcpy_s(out_SZA, outputTemp);
			strcat_s(out_SZA, "HIMAWARI8_");
			strcat_s(out_SZA, timestr.c_str());
			strcat_s(out_SZA, "_SZA_2KM.img");
			SaveImageToFile(out_SZA, win_width, win_height, &SZA, "ENVI", adfGeoTransform, pszSRS_WKT,1);
			
			cout<<"*******SZA save ok!*******\n";

			//save SAA
			char out_SAA[128];
			strcpy_s(out_SAA, outputTemp);		
			strcat_s(out_SAA, "HIMAWARI8_");
			strcat_s(out_SAA, timestr.c_str());
			strcat_s(out_SAA, "_sensor_azmuth_2KM.img");
			SaveImageToFile(out_SAA, win_width, win_height, &SAA, "ENVI", adfGeoTransform, pszSRS_WKT,1);
			
			cout<<"*******sensor_azimuth save ok!*******\n";

			//save SAZ
			char out_SAZ[128];
			strcpy_s(out_SAZ, outputTemp);
			strcat_s(out_SAZ, "HIMAWARI8_");
			strcat_s(out_SAZ, timestr.c_str());
			strcat_s(out_SAZ, "_sensor_zenith_2KM.img");
			SaveImageToFile(out_SAZ, win_width, win_height, &SAZ, "ENVI", adfGeoTransform, pszSRS_WKT,1);
			
			cout<<"*******sensor_zenith save ok!*******\n";
			GDALClose(tmpdt);
			CPLFree(miu);
		}	
		if (cldmsk==CLOUDMSK)
		{
			cout << "----------H8 read and cut ok!----(already cloud mask)------";
		}
		else
		{
			cout << "----------H8 read and cut ok!----(not cloud mask)------";
		}
		}
//get cloud mask data(*msk)
void H8_cldmsk()
{
	//*msk has two value 0-->cloud 1-->clear
	
	msk = (int *)CPLMalloc(sizeof(int)*ALLPIXELNUM);

	cldmsk_land(TOAb1, BTb1, BTb4, BTb8, BTb9, msk);

	cout << "******************************************" << endl;
	//count cloud and clear pixels
	int count = 0;
	for (int i = 0; i<ALLPIXELNUM; i++)
	{	
		if (msk[i] == 1)
		{
			count++;
		}	
	}

	cout << count<<"---->clear pixels"<<endl<<ALLPIXELNUM-count<<"---->cloud pixels"<<endl;
	//save cloud mask img
	char outcld[128];
	strcpy_s(outcld, temp_dir);	
	strcat_s(outcld, "cldmsk.img");
	SaveImageToFile(outcld, win_width, win_height, &msk, "ENVI", NULL, NULL,1);
	//cloud mask TOA

	
	cout << "******************************************" << endl;
}
//not finish yet
void H8_absorb_corect()
{
	for (int i = 1; i < 7; i++)
	{
		Absorb_Correction(i);
		//use Trans_H2O etc. and save .img
		CPLFree(RTrans_H2O);
		CPLFree(RTrans_O3);
		CPLFree(RTrans_other);

	}
}
//resulotion (eg._2KM or _10KM) 
void H8_resample(const char* inputDir,const char* outputDir,float resulotion,int bandnum,vector<string>::iterator it)
{   
	
	double adfGeoTransform[6];

	//set resamling method 
	GDALRasterIOExtraArg psExtrArg;
	INIT_RASTERIO_EXTRA_ARG(psExtrArg);		
	psExtrArg.eResampleAlg = GRIORA_Average;


	GDALAllRegister();

	//buffer-->full filename include path ;
	char buffer[128];
	strcpy_s(buffer, inputDir);
	strcat_s(buffer, it->c_str());
	const char *pszSrcFile = buffer;
	//*pDataSet-->all subdatesets and metadata
	GDALDataset *pDataSet = (GDALDataset*)GDALOpen(pszSrcFile, GA_ReadOnly);
	if (pDataSet == NULL)
	{
		cout<<"File cannot open!";
		return;
	}

	//set new projection
	pDataSet->GetGeoTransform(adfGeoTransform);
	adfGeoTransform[2] = resulotion;
	adfGeoTransform[4] = resulotion;
	adfGeoTransform[1] = 0;
	adfGeoTransform[5] = 0;


	int wid= pDataSet->GetRasterXSize();
	int hgt= pDataSet->GetRasterYSize();
	int newwid = wid/resulotion*_2KM;
	int newhgt = hgt/resulotion*_2KM;

	

	GDALDataset *newdataset;
	GDALRasterBand *newBand;
	char **papszOptions = NULL;
	char tmpout[128];
	strcpy_s(tmpout, outputDir);
	strcat_s(tmpout, it->substr(0, it->rfind("2KM")).c_str());
	cout << it->substr(0, it->rfind("_2KM")).c_str();
	strcat_s(tmpout, "10KM.img");
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("ENVI");
	newdataset = poDriver->Create(tmpout, newwid, newhgt, bandnum, GDT_Float32, papszOptions);
	float *data = (float*)CPLMalloc(sizeof(float)*newwid*newhgt);
	GDALRasterBand *poBand;
	for (int i = 1; i <= bandnum; i++)
	{
	
		poBand = pDataSet->GetRasterBand(i);
		
		poBand->RasterIO(GF_Read, 0, 0, wid, hgt, data, newwid, newhgt, GDT_Float32, 0, 0, &psExtrArg);	
		newBand = newdataset->GetRasterBand(i);	
		newBand->RasterIO(GF_Write, 0, 0, newwid, newhgt, data, newwid, newhgt, GDT_Float32, 0, 0);
		for (int i = 0; i < newwid*newhgt; i++)
		{
			 data[i]=0;
		}
				
	}
	newdataset->SetGeoTransform(adfGeoTransform);
	CPLFree(data);
	GDALClose(pDataSet);
	GDALClose(newdataset);

}

//retrival
void AOD_retrival()
{
	GDALAllRegister();
	//look up table AOD parameters(const)
	double AODARR_470[7] = {0,0.3378,0.6657,1.2882,2.4398,3.6622,6.1082};
	double AODARR_550[7] = {0,0.25,0.5,1.0,2.0,3.0,5.0};
	double AODARR_660[7] = {0,0.1876,0.3788,0.7768,1.6332,2.4478,4.076};
	double AODARR_2250[7] = {0,0.0446,0.0803,0.1443,0.2913,0.4283,0.6963};
	//prior file path
	string month = "04";
	char DEM_file[128] = "F:\\Data-for-AOD-Retrieval\\DEM\\China_dem.img";
	char MCD12_file[128] = "F:\\Data-for-AOD-Retrieval\\MCD12\\MCD12C1.A2012001.051.resize.img";
	char RatioBRN_file[128];
	string tmpRatioBRN_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\BRN_Mon" + month + "_10km.img";
	strcpy_s(RatioBRN_file, tmpRatioBRN_file.c_str());
	char RatioRRN_file[128];
	string tmpRatioRRN_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\RRN_Mon" + month + "_10km.img";
	strcpy_s(RatioRRN_file, tmpRatioRRN_file.c_str());
	char iso660_file[128];
	string tmpiso660_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\MCD43C2." + month + "_band1_P1_mean_10km.img";
	strcpy_s(iso660_file, tmpiso660_file.c_str());
	char vol660_file[128];
	string tmpvol660_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\MCD43C2." + month + "_band1_P2_mean_10km.img";
	strcpy_s(vol660_file, tmpvol660_file.c_str());
	char geo660_file[128];
	string tmpgeo660_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\MCD43C2." + month + "_band1_P3_mean_10km.img";
	strcpy_s(geo660_file, tmpgeo660_file.c_str());
	char iso470_file[128];
	string tmpiso470_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\MCD43C2." + month + "_band3_P1_mean_10km.img";
	strcpy_s(iso470_file, tmpiso470_file.c_str());
	char vol470_file[128];
	string tmpvol470_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\MCD43C2." + month + "_band3_P2_mean_10km.img";
	strcpy_s(vol470_file, tmpvol470_file.c_str());
	char geo470_file[128];
	string tmpgeo470_file = "F:\\Data-for-AOD-Retrieval\\prior-data\\MCD43C2." + month + "_band3_P3_mean_10km.img";
	strcpy_s(geo470_file, tmpgeo470_file.c_str());


	//prior data read
	float *DEM = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM );
	float *MCD12 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *RatioBRN = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *RatioRRN = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *iso660 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *vol660 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *geo660 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *iso470 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *vol470 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *geo470 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	ReadIMG(DEM_file,&DEM,1);
	ReadIMG(MCD12_file, &MCD12,1);
	ReadIMG(RatioBRN_file, &RatioBRN,1);
	ReadIMG(RatioRRN_file, &RatioRRN,1);
	ReadIMG(iso660_file, &iso660,1);
	ReadIMG(vol660_file, &vol660,1);
	ReadIMG(geo660_file, &geo660,1);
	ReadIMG(iso470_file, &iso470,1);
	ReadIMG(vol470_file, &vol470,1);
	ReadIMG(geo470_file, &geo470,1);
	//toa and angle data path
	FILELIST TOALIST;
	if (!FindFiles(&TOALIST, ENVIAHI10KM, AHI_TOA))
	{
		cout << "NO TOA10km file exist!\n";		
	}
	FILELIST SOALIST;
	if (!FindFiles(&SOALIST, ENVISOA10KM, ANGLE_dir))
	{
		cout << "NO SOA10km file exist!\n";
	}
	FILELIST SZALIST;
	if (!FindFiles(&SZALIST, ENVISZA10KM, ANGLE_dir))
	{
		cout << "NO SZA10km file exist!\n";
	}
	char VAA_file[128];
	string tmpVAA_file = "F:\\AHIout\\angle\\HIMAWARI8_20170401_0000_sensor_azmuth_10KM.img";
	strcpy_s(VAA_file, tmpVAA_file.c_str());
	char VZA_file[128];
	string tmpVZA_file = "F:\\AHIout\\angle\\HIMAWARI8_20170401_0000_sensor_zenith_10KM.img";
	strcpy_s(VZA_file, tmpVZA_file.c_str());
	//toa and angle data read
	float *VAA = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *VZA = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAa = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAb = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAc = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAd = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAe = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAf = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAg = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SOAh = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAa = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAb = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAc = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAd = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAe = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAf = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAg = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *SZAh = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470a = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470b = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470c = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470d = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470e = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470f = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470g = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA470h = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660a = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660b = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660c = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660d = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660e = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660f = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660g = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA660h = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250a = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250b = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250c = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250d = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250e = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250f = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250g = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *TOA2250h = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *AOD470 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *AOD550 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *AOD660 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *ERROR470 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *ERROR550 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *ERROR660 = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *AOD470b = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *AOD550b = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);
	float *AOD660b = (float*)CPLMalloc(sizeof(float)*_10KMALLPIXELNUM);

	ReadIMG(VAA_file,&VAA,1);
	ReadIMG(VZA_file, &VZA, 1);
	for (int k = 30; k > 0; k--)
	{
		string reffile = TOALIST[k];
		string reffile2 = TOALIST[k+4];
		string reffile3 = TOALIST[k-4];
		string reffile4 = TOALIST[k+3];
		string reffile5 = TOALIST[k-3];
		string reffile6 = TOALIST[k+2];
		string reffile7 = TOALIST[k-2];
		string reffile8 = TOALIST[k+1];

		string SOAfilea = SOALIST[k];
		string SOAfileb = SOALIST[k + 4];
		string SOAfilec = SOALIST[k - 4];
		string SOAfiled = SOALIST[k + 3];
		string SOAfilee = SOALIST[k - 3];
		string SOAfilef = SOALIST[k + 2];
		string SOAfileg = SOALIST[k - 2];
		string SOAfileh = SOALIST[k + 1];

		string SZAfilea = SZALIST[k];
		string SZAfileb = SZALIST[k + 4];
		string SZAfilec = SZALIST[k - 4];
		string SZAfiled = SZALIST[k + 3];
		string SZAfilee = SZALIST[k - 3];
		string SZAfilef = SZALIST[k + 2];
		string SZAfileg = SZALIST[k - 2];
		string SZAfileh = SZALIST[k + 1];


		ReadIMG(SOAfilea, &SOAa,1);
		ReadIMG(SOAfileb, &SOAb, 1);
		ReadIMG(SOAfilec, &SOAc, 1);
		ReadIMG(SOAfiled, &SOAd, 1);
		ReadIMG(SOAfilee, &SOAe, 1);
		ReadIMG(SOAfilef, &SOAf, 1);
		ReadIMG(SOAfileg, &SOAg, 1);
		ReadIMG(SOAfileh, &SOAh, 1);

		ReadIMG(SZAfilea, &SZAa, 1);
		ReadIMG(SZAfileb, &SZAb, 1);
		ReadIMG(SZAfilec, &SZAc, 1);
		ReadIMG(SZAfiled, &SZAd, 1);
		ReadIMG(SZAfilee, &SZAe, 1);
		ReadIMG(SZAfilef, &SZAf, 1);
		ReadIMG(SZAfileg, &SZAg, 1);
		ReadIMG(SZAfileh, &SZAh, 1);

		ReadIMG(reffile, &TOA470a, 1);
		ReadIMG(reffile2, &TOA470b, 1);
		ReadIMG(reffile3, &TOA470c, 1);
		ReadIMG(reffile4, &TOA470d, 1);
		ReadIMG(reffile5, &TOA470e, 1);
		ReadIMG(reffile6, &TOA470f, 1);
		ReadIMG(reffile7, &TOA470g, 1);
		ReadIMG(reffile8, &TOA470h, 1);
		ReadIMG(reffile, &TOA660a, 3);
		ReadIMG(reffile2, &TOA660b, 3);
		ReadIMG(reffile3, &TOA660c, 3);
		ReadIMG(reffile4, &TOA660d, 3);
		ReadIMG(reffile5, &TOA660e, 3);
		ReadIMG(reffile6, &TOA660f, 3);
		ReadIMG(reffile7, &TOA660g, 3);
		ReadIMG(reffile8, &TOA660h, 3);
		ReadIMG(reffile, &TOA2250a, 6);
		ReadIMG(reffile2, &TOA2250b, 6);
		ReadIMG(reffile3, &TOA2250c, 6);
		ReadIMG(reffile4, &TOA2250d, 6);
		ReadIMG(reffile5, &TOA2250e, 6);
		ReadIMG(reffile6, &TOA2250f, 6);
		ReadIMG(reffile7, &TOA2250g, 6);
		ReadIMG(reffile8, &TOA2250h, 6);

		for (int i = 0; i < _10KMALLPIXELNUM; i++)
		{
			//toa skip
			if(TOA2250a[i]<0|| TOA2250a[i]>0.6||TOA470a[i]<0||MCD12[i]<=0||SZAa[i]>75||VZA[i]>75)
			{
				continue;
			}
			//prior skip
			if (iso470[i]<=0||vol470[i]<=0||geo470[i]<=0||iso660[i]<=0||vol660[i]<=0||geo660[i]<=0)
			{
				continue;
			}
			float TOAR470a = 0;
			float TOAR660a = 0;
			float TOAR2250a = 0;
			float TOAR470b = 0;
			float TOAR660b = 0;
			float TOAR2250b = 0;
			float SOA1 = 0;
			float SZA1 = 0;
			float SOA2 = 0;
			float SZA2 = 0;
			float VOA1 = 0;
			float VZA1 = 0;
			float RR = 0;
			float RB = 0;
			//choose clear data for retrival
			float TOAarr[7] = { TOA2250b[i],TOA2250c[i],TOA2250d[i],TOA2250e[i] ,TOA2250f[i] ,TOA2250g[i] ,TOA2250h[i] };
			int index = 0;
			for (int j = 0; j < 7; j++)
			{
				if(TOAarr[j]<=0)
				{
					continue;
					index = -1;
				}
				else
				{
					index = j;
					break;
				}
			}
		
			if (index == -1)
			{
				continue;
			}
			else if (index == 0)
			{
				TOAR470b = TOA470b[i];
				TOAR660b = TOA660b[i];
			    TOAR2250b = TOA2250b[i];
			    SOA2 = SOAb[i];
				SZA2 = SZAb[i];
			}
			else if (index == 1)
			{
				TOAR470b = TOA470c[i];
				TOAR660b = TOA660c[i];
				TOAR2250b = TOA2250c[i];
				SOA2 = SOAc[i];
				SZA2 = SZAc[i];
			}
			else if (index == 2)
			{
				TOAR470b = TOA470d[i];
				TOAR660b = TOA660d[i];
				TOAR2250b = TOA2250d[i];
				SOA2 = SOAd[i];
				SZA2 = SZAd[i];
			}
			else if (index == 3)
			{
				TOAR470b = TOA470e[i];
				TOAR660b = TOA660e[i];
				TOAR2250b = TOA2250e[i];
				SOA2 = SOAe[i];
				SZA2 = SZAe[i];
			}
			else if (index == 4)
			{
				TOAR470b = TOA470f[i];
				TOAR660b = TOA660f[i];
				TOAR2250b = TOA2250f[i];
				SOA2 = SOAf[i];
				SZA2 = SZAf[i];
			}
			else if (index == 5)
			{
				TOAR470b = TOA470g[i];
				TOAR660b = TOA660g[i];
				TOAR2250b = TOA2250g[i];
				SOA2 = SOAg[i];
				SZA2 = SZAg[i];
			}
			else if (index == 6)
			{
				TOAR470b = TOA470h[i];
				TOAR660b = TOA660h[i];
				TOAR2250b = TOA2250h[i];
				SOA2 = SOAh[i];
				SZA2 = SZAh[i];
			}
			SOA1 = SOAa[i];
			SZA1 = SZAa[i];
			VOA1 = VAA[i];
			VZA1 = VZA[i];
			RR = RatioRRN[i];
			RB = RatioBRN[i];
			TOAR470a = TOA470a[i];
			TOAR660a = TOA660a[i];
			TOAR2250a = TOA2250a[i];

			float RAAngle1 = abs(SOA1-VOA1-180);
			if (RAAngle1 > 360)
			{
				RAAngle1 -= 360;
			}
			else if (RAAngle1 > 180 && RAAngle1 < 360)
			{
				RAAngle1 = 360 - RAAngle1;
			}
			float RAAngle2 = abs(SOA2 - VOA1 - 180);
			if (RAAngle2 > 360)
			{
				RAAngle2 -= 360;
			}
			else if (RAAngle2 > 180 && RAAngle2 < 360)
			{
				RAAngle2 = 360 - RAAngle2;
			}
			float SCAT_Angle1 = -cos(SZA1*Pi / 180)*cos(VZA1*Pi / 180) + sin(SZA1*Pi / 180)*sin(VZA1*Pi / 180)*cos(RAAngle1*Pi/180);
			SCAT_Angle1 = acos(SCAT_Angle1)*180/Pi;
			float SCAT_Angle2 = -cos(SZA2*Pi / 180)*cos(VZA1*Pi / 180) + sin(SZA2*Pi / 180)*sin(VZA1*Pi / 180)*cos(RAAngle2*Pi/180);
			SCAT_Angle2 = acos(SCAT_Angle2)*180/Pi;

			if (TOAR2250a <= 0.6&&TOAR2250b <= 0.6)
			{
				Inter_Angle(SZA1,VZA1,RAAngle1,AODARR_470);
				Inter_Angle(SZA1, VZA1, RAAngle1, AODARR_550);
				Inter_Angle(SZA1, VZA1, RAAngle1, AODARR_660);
				Inter_Angle(SZA1, VZA1, RAAngle1, AODARR_2250);
				//second image
				Inter_Angle(SZA2, VZA1, RAAngle2, AODARR_470);
				Inter_Angle(SZA2, VZA1, RAAngle2, AODARR_550);
				Inter_Angle(SZA2, VZA1, RAAngle2, AODARR_660);
				Inter_Angle(SZA2, VZA1, RAAngle2, AODARR_2250);
				
			}
		}

	}


	CPLFree(DEM);
	CPLFree(MCD12);
	CPLFree(RatioBRN);
	CPLFree(RatioRRN);
	CPLFree(iso660);
	CPLFree(vol660);
	CPLFree(geo660);
	CPLFree(iso470);
	CPLFree(vol470);
	CPLFree(geo470);
	CPLFree(VAA);
	CPLFree(VZA);
	CPLFree(SOAa);
	CPLFree(SOAb);
	CPLFree(SOAc);
	CPLFree(SOAd);
	CPLFree(SOAe);
	CPLFree(SOAf);
	CPLFree(SOAg);
	CPLFree(SOAh);
	CPLFree(SZAa);
	CPLFree(SZAb);
	CPLFree(SZAc);
	CPLFree(SZAd);
	CPLFree(SZAe);
	CPLFree(SZAf);
	CPLFree(SZAg);
	CPLFree(SZAh);
	CPLFree(TOA470a);
	CPLFree(TOA470b);
	CPLFree(TOA470c);
	CPLFree(TOA470d);
	CPLFree(TOA470e);
	CPLFree(TOA470f);
	CPLFree(TOA470g);
	CPLFree(TOA470h);
	CPLFree(TOA660a);
	CPLFree(TOA660b);
	CPLFree(TOA660c);
	CPLFree(TOA660d);
	CPLFree(TOA660e);
	CPLFree(TOA660f);
	CPLFree(TOA660g);
	CPLFree(TOA660h);
	CPLFree(TOA2250a);
	CPLFree(TOA2250b);
	CPLFree(TOA2250c);
	CPLFree(TOA2250d);
	CPLFree(TOA2250e);
	CPLFree(TOA2250f);
	CPLFree(TOA2250g);
	CPLFree(TOA2250h);
	CPLFree(AOD470);
	CPLFree(AOD550);
	CPLFree(AOD660);
	CPLFree(ERROR470);
	CPLFree(ERROR550);
	CPLFree(ERROR660);
	CPLFree(AOD470b);
	CPLFree(AOD550b);
	CPLFree(AOD660b);


}
//retrival functions 
//get atm_ref,fdn,T,sbar,fd
void Inter_Angle(float SZA,float VZA,float RAAngle,double aodarr[7])
{
	int SZAARR[11] = { 0,12,24,36,48,54,60,66,72,78,84};
	int VZAARR[15] = { 0,6,12,18,24,30,36,42,48,54,60,66,72,78,84};
	int RAAARR[15] = { 0,6,12,18,24,30,36,42,48,54,60,66,72,78,84 };	
	int SZAindex=0;
	int VZAindex=0;
	int RAAindex=0;
	for (int i = 0; i < 11; i++)
	{
		if (SZAARR[i]<SZA&&SZAARR[i + 1]>SZA)
		{
			SZAindex = i;
		}
	}
	for (int i = 0; i < 15; i++)
	{
		if (VZAARR[i]<VZA&&VZAARR[i + 1]>VZA)
		{
			VZAindex = i;
			
		}
		if (RAAARR[i]<RAAngle&&RAAARR[i + 1]>RAAngle)
		{
			RAAindex = i;
			
		}
	
	}




}
void Calculate_SRef(float TOAR2250,float ARR2250,float Rred,float Rblue,float SCANT_ANGLE)
{
	int y644 = 0;
	int y466 = 0;
	float slope644 = 0;
	float slope466 = 0;
	if (Rred > 0 && Rblue > 0)
	{
		slope644 = Rred;
		slope466 = Rblue;
	}
	else
	{
		slope644 = 0.00027*SCANT_ANGLE + 0.05651;
		slope466 = -2.663055*0.00001*pow(SCANT_ANGLE,2) + 8.592420*0.001*SCANT_ANGLE - 0.3671062;
	}
	   

}

//2 overload for int and float data[]
//outputPath need absolute path
//pszFormat (eg. "ENVI") ;adfGeoTransform[6] refers 6 parameters of projection;pszSRS_WKT
bool SaveImageToFile(const char *outputPath, int nWidth, int nHeight, float **data, const char* pszFormat, double adfGeoTransform[6], const char *pszSRS_WKT, int Bandnum)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if (poDriver == NULL) {
		cout<<"The driver is not supported!";
		return false;
	}

	GDALDataset *poDstDS;
	char **papszOptions = NULL;
	poDstDS = poDriver->Create(outputPath, nWidth, nHeight, Bandnum, GDT_Float32, papszOptions);

	//set map projection information
	poDstDS->SetProjection(pszSRS_WKT);
	poDstDS->SetGeoTransform(adfGeoTransform);
	for (int i = 1; i <= Bandnum; i++)
	{	
	poDstDS->RasterIO(GF_Write, 0, 0, nWidth, nHeight, *data+(i-1)*ALLPIXELNUM, nWidth, nHeight, GDT_Float32, Bandnum, NULL, 0, 0, 0);
    }
	GDALClose(poDstDS);
	return true;
}
bool SaveImageToFile(const char *outputPath, int nWidth, int nHeight, int **data, const char* pszFormat, double adfGeoTransform[6], const char *pszSRS_WKT,int Bandnum)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if (poDriver == NULL) {
		cout<<"The driver is not supported!";
		return false;
	}

	GDALDataset *poDstDS;
	GDALRasterBand *poBand;
	char **papszOptions = NULL;
	poDstDS = poDriver->Create(outputPath, nWidth, nHeight, Bandnum, GDT_Int32, papszOptions);

	//set map projection information
	
	poBand = poDstDS->GetRasterBand(Bandnum);

	poBand->RasterIO(GF_Write, 0, 0, nWidth, nHeight, *data+(Bandnum- 1)*ALLPIXELNUM, nWidth, nHeight, GDT_Int32,0, 0);
	GDALClose(poDstDS);
	return true;
}
//2 overload read img and save data (char[] and string)
void ReadIMG(const char *filePath, float **data,int bandnum)
{
	GDALAllRegister();
	//buffer-->full filename include path ;
	char buffer[128];
	strcpy_s(buffer, filePath);
	const char *pszSrcFile = buffer;
	//*pDataSet-->all subdatesets and metadata
	GDALDataset *pDataSet = (GDALDataset*)GDALOpen(pszSrcFile, GA_ReadOnly);
	if (pDataSet == NULL)
	{
		cout << "File cannot open!";
		return;
	}	
	GDALRasterBand *poBand;
	int wid = pDataSet->GetRasterXSize();
	int hgt = pDataSet->GetRasterYSize();
	poBand = pDataSet->GetRasterBand(bandnum);
	poBand->RasterIO(GF_Read, 0, 0, wid, hgt, data, wid, hgt, GDT_Float32, 0, 0, NULL);
	GDALClose(pDataSet);
}
void ReadIMG(string filePath, float **data, int bandnum)
{
	GDALAllRegister();
	//buffer-->full filename include path ;
	char buffer[128];
	strcpy_s(buffer, filePath.c_str());
	const char *pszSrcFile = buffer;
	//*pDataSet-->all subdatesets and metadata
	GDALDataset *pDataSet = (GDALDataset*)GDALOpen(pszSrcFile, GA_ReadOnly);
	if (pDataSet == NULL)
	{
		cout << "File cannot open!";
		return;
	}
	GDALRasterBand *poBand;
	int wid = pDataSet->GetRasterXSize();
	int hgt = pDataSet->GetRasterYSize();
	poBand = pDataSet->GetRasterBand(bandnum);
	poBand->RasterIO(GF_Read, 0, 0, wid, hgt, data, wid, hgt, GDT_Float32, 0, 0, NULL);
	GDALClose(pDataSet);
}
//FileNames is filelist; fileFormat (eg. netCDF) ;path is folder path
bool FindFiles(vector<string>* FileNames, FILETYPE fileFormat, const char* path)
{
	DIR *dir_p;
	struct dirent *direntp;
	if ((dir_p = opendir(path)) == NULL)
	{
		cout<<"dir ERROR !!\n";
		return false;
	}
	else
	{
		while ((direntp = readdir(dir_p)) != NULL)
		{

			switch (fileFormat)
			{
			case NETCDF:
			{
				if (direntp->d_namlen == 44)//match stadard Himawari .nc file
				{

					FileNames->push_back(direntp->d_name);
				}
				break;
			}
			case ENVIAHI:
			{
				if (direntp->d_namlen == 35 && direntp->d_name[34] == 'g')//match ENVI AHI 2KM.img file not hdr
				{
					FileNames->push_back(direntp->d_name);
				}
				break;
			}
			case ENVIBT:
			{
				if (direntp->d_namlen == 34 && direntp->d_name[33] == 'g')//match  ENVI BT .img file
				{

					FileNames->push_back(direntp->d_name);
				}
				break;
			}
			case ENVIANGLE:
			{
				if (direntp->d_namlen == 35 && direntp->d_name[34] == 'g')//match  ENVI SOA\SOZ .img file
				{

					FileNames->push_back(direntp->d_name);
				}
				else if (direntp->d_namlen == 45 && direntp->d_name[44] == 'g')//match  ENVI sensor azimuth\zenith .img file
				{
					FileNames->push_back(direntp->d_name);
				}
				break;
			}
			case ENVIAHI10KM:
			{
				if (direntp->d_namlen == 36 && direntp->d_name[35] == 'g')//match ENVI AHI 10KM.img file not hdr
				{
					FileNames->push_back(direntp->d_name);
				}
				break;
			}
			case ENVISOA10KM:
			{

				if (direntp->d_namlen == 36 && direntp->d_name[35] == 'g'&&direntp->d_name[25] == 'O')//match  ENVI 10KM SOA .img file
				{

					FileNames->push_back(direntp->d_name);
				}
			}
			case ENVISZA10KM:
			{

				if (direntp->d_namlen == 36 && direntp->d_name[35] == 'g'&&direntp->d_name[25] == 'Z')//match  ENVI 10KM SZA .img file
				{

					FileNames->push_back(direntp->d_name);
				}
			}
			
		   }
		  }
	
		
		sort(FileNames->begin(), FileNames->end());
		closedir(dir_p);
	}
	return true;
}
//delect all files in folder path
void FileDelect(const char* path)
{
	DIR *dir_p;
	struct dirent *direntp;
	if ((dir_p = opendir(path)) == NULL)
	{
		cout<<"dir ERROR !!\n";
		return ;
	}
	else
	{
		char filename[100];
		while ((direntp = readdir(dir_p)) != NULL)
		{
		
			strcpy_s(filename,path);		
			strcat_s(filename, direntp->d_name);
			remove(filename);
			
		}
	}
	closedir(dir_p);
}

//rewrite from IDL by YLK
void cldmsk_land(float *TOAb1, float *BTb1, float *BTb4, float *BTb8, float *BTb9,int *msk)
{
	//only use TOAb1 and output msk
	GDALAllRegister();

	//DESCRIPTION:Generate cloud mask over land using spatial variability of 0.47 (>0.0025) as well as absolute value of 0.47 um >(0.4) 
	//1-->clear  0-->cloud
	const float sptstdvar = 0.0025;
    const float absval = 0.4;
	int * cldmsk=NULL;
	cldmsk = (int *)CPLMalloc(sizeof(int)*ALLPIXELNUM);
	//all set 1-->clear
	valuecopy(1,cldmsk);

	//3*3 grid calculate TOA cloud mask
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		//ignore boundary (first and last line(sample) )
		if(i<win_width||i>win_width*(win_height-1)||i%win_width==0||(i+1)%win_width==0)
		{
			cldmsk[i] = 0;
			continue;
		}		
		else
		{
			
			/*  a b c      ----->3*3 grid
		        d e f
                g h k    */

			float a, b, c, d, e, f, g, h, k=0;
			float avg = 0;
			float stdvar = 0;
			a = TOAb1[i - win_width - 1];
			b = TOAb1[i - win_width];
			c = TOAb1[i - win_width + 1];
			d = TOAb1[i - 1];
			e = TOAb1[i];
			f = TOAb1[i + 1];
			g = TOAb1[i + win_width - 1];
			h = TOAb1[i + win_width];
			k = TOAb1[i + win_width + 1];
			
			if(a<=0&&b<=0&&c<=0&&d<=0&&e<=0&&f<=0&&g<=0&&h<=0&&k<=0)
			{ 
				cldmsk[i] = 0;
				continue;
			}
			else if (e>= absval)//not land
			{
				cldmsk[i] = 0;
				continue;
			}
			avg = (a+b+c+d+e+f+g+h+k)/9;
			stdvar = sqrt((pow(a - avg, 2) + pow(b - avg, 2) + pow(c - avg, 2) + pow(d - avg, 2) + pow(f - avg, 2) + pow(g - avg, 2) + pow(h - avg, 2) + pow(k - avg, 2))/8);
			float newTOAb1 = stdvar*avg/3;
			if (newTOAb1 >= sptstdvar&&stdvar >= 0.0075)
			{
				cldmsk[i] = 0;
			}
		
			
		}	
	}
	//export to msk
	valuecopy(cldmsk,msk);
	
	CPLFree(BTb1);
	CPLFree(BTb4);
	CPLFree(BTb8);
	CPLFree(BTb9);
	CPLFree(cldmsk);

}
//absorb correction TOA band(bandnum<=6)
void Absorb_Correction(int bandnum)
{
	int bdn = bandnum - 1;
	//H20 and O3 coefficient 
	double K1_H2O[6] = { -9.85,-7.91,-5.6,-5.07,-6.8,-3.98 };
	double K2_H2O[6] = { 1.23, 1.00, 0.94, 0.877, 1.03, 0.886 };
	double K3_H2O[6] = { -0.116, -0.0129, -0.0178, -0.024, -0.00429, -0.0256 };
	double K1_O3[6] = { -1.14e-4, 5.18e-6, 1.16e-4, 2.8e-7, 1.19e-7, 6.29e-7 };
	double K2_O3[6] = { 8.69e-6, 9.50e-5, 7.32e-5, 2.36e-6, 5.17e-26, 7.03e-8 };
	double Opt_H2O[6] = { 8.0e-5, 5.0e-4, 5.11e-3, 8.61e-3,1.62e-3,2.53e-2 };
	double Opt_O3[6] = { 2.90e-3, 3.26e-2,2.52e-2, 8.10e-4,0.0,2.0e-5 };
	double Tao_Other[6] = {1.25e-3, 9.5e-4,3.91e-3, 2.0e-5,9.98e-3,1.63e-2 };

    //read globle average img
	double * Pwater = NULL;
	Pwater = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);
	//Pwater *= 0.1;
	double * Pozone = NULL;
	Pozone = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);
	//Pozone = Pozone*100000/2.1414;
	
	//malloc
	RTrans_H2O = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);
	RTrans_O3 = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);
	RTrans_other = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);

	//compute
	double * GH2O = NULL;
	GH2O = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);
	double * GO3 = NULL;
	GO3 = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);
	double * Gother = NULL;
	Gother = (double *)CPLMalloc(sizeof(double)*ALLPIXELNUM);

	

	float a1[3] = { 268.45,0.0311,0.4567 };
	float a2[3] = { 0.5,0.1,0.07 };
	float a3[3] = { 115.42,92.471,96.484 };
	float a4[3] = { -3.2922,-1.3814,-1.6970 };
	
	int value_w, value_o, nvalue_w, nvalue_o=0;
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		if(SZA[i]<=0||SAZ[i]<=0||TOAb1[i]<0)
		{
			continue;
		}
		else
		{
			GO3[i] = 1 / cos((SZA[i] * Pi / 180) + pow(a1[0] * SZA[i], a2[0])*pow(a3[0] - SZA[i], a4[0])) + 1 / cos((SAZ[i] * Pi / 180) + pow(a1[0] * SAZ[i], a2[0])*pow(a3[0] - SAZ[i], a4[0]));
			GH2O[i] = 1 / cos((SZA[i] * Pi / 180) + pow(a1[1] * SZA[i], a2[1])*pow(a3[1] - SZA[i], a4[1])) + 1 / cos((SAZ[i] * Pi / 180) + pow(a1[1] * SAZ[i], a2[1])*pow(a3[1] - SAZ[i], a4[1]));
			Gother[i] = 1 / cos((SZA[i] * Pi / 180) + pow(a1[2] * SZA[i], a2[2])*pow(a3[2] - SZA[i], a4[2])) + 1 / cos((SAZ[i] * Pi / 180) + pow(a1[2] * SAZ[i], a2[2])*pow(a3[2] - SAZ[i], a4[2]));
		  if(Pwater[i]>0&& GH2O[i]>0)
		  { 
			  value_w = 1;
			  nvalue_w = 0;
		  }
		  else
		  {
			  value_w = 0;
			  nvalue_w = 1;
		  }
		  if (Pozone[i]>0 && GO3[i]>0)
		  {
			  value_o = 1;
			  nvalue_o = 0;
		  }
		  else
		  {
			  value_o = 0;
			  nvalue_o = 1;
		  }
		  double logcon = log(Pwater[i] * GH2O[i]);
		  double longcon2 = logcon*logcon;
		  double exponent = Pozone[i] * GO3[i];
		 
		  RTrans_H2O[i] = exp(exp(K1_H2O[bdn]+K2_H2O[bdn]*logcon+K3_H2O[bdn]*longcon2))*value_w+exp(Opt_H2O[bdn]*GH2O[i])*nvalue_w;
		  RTrans_O3[i] = exp(K1_O3[bdn] + K2_O3[bdn] * exponent)*value_o + exp(Opt_O3[bdn]*GO3[i])*nvalue_o;
		  RTrans_other[i] = exp(Tao_Other[bdn]*Gother[i]);
		}
	}

	CPLFree(Pwater);
	CPLFree(Pozone);
	CPLFree(GH2O);
	CPLFree(GO3);
	CPLFree(Gother);

}
//four overload two usage  1:copy value from array to array (int and float) 2:array initialization -->all set num(int)
void valuecopy(float *input, float *output)
{
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		output[i] = input[i];
	}
}
void valuecopy(int *input, int *output)
{
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		output[i] = input[i];
	}
}
void valuecopy(int num, int *output)
{
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		output[i] = num;
	}
}
void valuecopy(float num, float *output)
{
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		output[i] = num;
	}
}
//for extention types (double,char,string etc.)
template<typename T>
void valuecopy(T *input,T *output)
{
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		output[i] = input[i];
	}
}
template<typename T>
void valuecopy(T num, T *output)
{
	for (int i = 0; i < ALLPIXELNUM; i++)
	{
		output[i] = num;
	}
}

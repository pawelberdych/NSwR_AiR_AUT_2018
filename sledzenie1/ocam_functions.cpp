/*------------------------------------------------------------------------------
   Example code that shows the use of the 'cam2world" and 'world2cam" functions
   Shows also how to undistort images into perspective or panoramic images
   Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich
   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
------------------------------------------------------------------------------*/

#include "ocam_functions.h"

//------------------------------------------------------------------------------
int get_ocam_model(struct ocam_model *myocam_model, char *filename)
{
 double *pol        = myocam_model->pol;
 double *invpol     = myocam_model->invpol;
 double *xc         = &(myocam_model->xc);
 double *yc         = &(myocam_model->yc);
 double *c          = &(myocam_model->c);
 double *d          = &(myocam_model->d);
 double *e          = &(myocam_model->e);
 int    *width      = &(myocam_model->width);
 int    *height     = &(myocam_model->height);
 int *length_pol    = &(myocam_model->length_pol);
 int *length_invpol = &(myocam_model->length_invpol);
 FILE *f;
 char buf[CMV_MAX_BUF];
 int i;

 //Open file
 if(!(f=fopen(filename,"r")))
 {
   printf("File %s cannot be opened\n", filename);
   return -1;
 }

 //Read polynomial coefficients
 fgets(buf,CMV_MAX_BUF,f);
 fscanf(f,"\n");
 fscanf(f,"%d", length_pol);
 for (i = 0; i < *length_pol; i++)
 {
     fscanf(f," %lf",&pol[i]);
 }

 //Read inverse polynomial coefficients
 fscanf(f,"\n");
 fgets(buf,CMV_MAX_BUF,f);
 fscanf(f,"\n");
 fscanf(f,"%d", length_invpol);
 for (i = 0; i < *length_invpol; i++)
 {
     fscanf(f," %lf",&invpol[i]);
 }

 //Read center coordinates
 fscanf(f,"\n");
 fgets(buf,CMV_MAX_BUF,f);
 fscanf(f,"\n");
 fscanf(f,"%lf %lf\n", xc, yc);

 //Read affine coefficients
 fgets(buf,CMV_MAX_BUF,f);
 fscanf(f,"\n");
 fscanf(f,"%lf %lf %lf\n", c,d,e);

 //Read image size
 fgets(buf,CMV_MAX_BUF,f);
 fscanf(f,"\n");
 fscanf(f,"%d %d", height, width);

 fclose(f);
 return 0;
}

//------------------------------------------------------------------------------
void cam2world(double point3D[3], double point2D[2], struct ocam_model *myocam_model)
{
 double *pol    = myocam_model->pol;
 double xc      = (myocam_model->xc);
 double yc      = (myocam_model->yc);
 double c       = (myocam_model->c);
 double d       = (myocam_model->d);
 double e       = (myocam_model->e);
 int length_pol = (myocam_model->length_pol);
 double invdet  = 1/(c-d*e); // 1/det(A), where A = [c,d;e,1] as in the Matlab file

 double xp = invdet*(    (point2D[0] - xc) - d*(point2D[1] - yc) );
 double yp = invdet*( -e*(point2D[0] - xc) + c*(point2D[1] - yc) );

 double r   = sqrt(  xp*xp + yp*yp ); //distance [pixels] of  the point from the image center
 double zp  = pol[0];
 double r_i = 1;
 int i;

 for (i = 1; i < length_pol; i++)
 {
   r_i *= r;
   zp  += r_i*pol[i];
 }

 //normalize to unit norm
 double invnorm = 1/sqrt( xp*xp + yp*yp + zp*zp );

 point3D[0] = invnorm*xp;
 point3D[1] = invnorm*yp;
 point3D[2] = invnorm*zp;
}

//------------------------------------------------------------------------------
void world2cam(double point2D[2], double point3D[3], struct ocam_model *myocam_model)
{
 double *invpol     = myocam_model->invpol;
 double xc          = (myocam_model->xc);
 double yc          = (myocam_model->yc);
 double c           = (myocam_model->c);
 double d           = (myocam_model->d);
 double e           = (myocam_model->e);
 int    width       = (myocam_model->width);
 int    height      = (myocam_model->height);
 int length_invpol  = (myocam_model->length_invpol);
 double norm        = sqrt(point3D[0]*point3D[0] + point3D[1]*point3D[1]);
 double theta       = atan(point3D[2]/norm);
 double t, t_i;
 double rho, x, y;
 double invnorm;
 int i;

  if (norm != 0)
  {
    invnorm = 1/norm;
    t  = theta;
    rho = invpol[0];
    t_i = 1;

    for (i = 1; i < length_invpol; i++)
    {
      t_i *= t;
      rho += t_i*invpol[i];
    }

    x = point3D[0]*invnorm*rho;
    y = point3D[1]*invnorm*rho;

    point2D[0] = x*c + y*d + xc;
    point2D[1] = x*e + y   + yc;
  }
  else
  {
    point2D[0] = xc;
    point2D[1] = yc;
  }
}

//------------------------------------------------------------------------------
void create_panoramic_undistortion_LUT ( Mat mapx, Mat mapy, float Rmin, float Rmax, float xc, float yc )
{
     int i, j;
     float theta;
     int width = mapx.cols;
     int height = mapx.rows;
    // float *data_mapx = mapx->data.fl;
    // float *data_mapy = mapy->data.fl;
     float rho;

     for (i=0; i<height; i++)
         for (j=0; j<width; j++)
         {
             theta = -((float)j)/width*2*CV_PI; // Note, if you would like to flip the image, just inverte the sign of theta
             rho   = Rmax - (Rmax-Rmin)/height*i;
             mapx.at<float>(i,j) = yc + rho*sin(theta); //in OpenCV "x" is the
             mapy.at<float>(i,j) = xc + rho*cos(theta);
         }
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void Ovu2point2d(double point2d[2], double pointVU[2], ocam_model *myocam_model)
{
	//pobranie wlasciwosci z kalibracji
	double xc = (myocam_model->xc);
	double yc = (myocam_model->yc);
	double c = (myocam_model->c);
	double d = (myocam_model->d);
	double e = (myocam_model->e);

	//Przeksztalcenie wspolrzednych ze srodkiem w ukladzie centralnym na wspolrzedne macierzowe
	point2d[0] = pointVU[0] * c + pointVU[1] * d + xc;
	point2d[1] = pointVU[0] * e + pointVU[1] + yc;
}

void Opoint2d2vu(double pointVU[2], double point2d[2], ocam_model *myocam_model)
{
	//pobranie wlasciwosci z kalibracji
	double xc = (myocam_model->xc);
	double yc = (myocam_model->yc);
	double c = (myocam_model->c);
	double d = (myocam_model->d);
	double e = (myocam_model->e);

	//Przeksztalcenie wspolrzednych macierzowych na wspolrzedne ze srodkiem w ukladzie centralnym
	double invdet = 1 / (c - d*e); // 1/det(A), where A = [c,d;e,1] as in the Matlab file
	pointVU[0] = invdet*((point2d[0] - xc) - d*(point2d[1] - yc));
	pointVU[1] = invdet*(-e*(point2d[0] - xc) + c*(point2d[1] - yc));

}

void O3d2vu(double pointVU[2], double point3d[3], ocam_model *myocam_model)
{
	double *invpol = myocam_model->invpol;

	int   width = (myocam_model->width);
	int   height = (myocam_model->height);
	int length_invpol = (myocam_model->length_invpol);
	double norm = sqrt(point3d[0] * point3d[0] + point3d[1] * point3d[1]);
	double theta = atan(point3d[2] / norm);
	double t, t_i;
	double rho, x, y;
	double invnorm;
	int i;

	if (norm != 0)
	{
		invnorm = 1 / norm;
		t = theta;
		rho = invpol[0];
		t_i = 1;

		for (i = 1; i < length_invpol; i++)
		{
			t_i *= t;
			rho += t_i*invpol[i];
		}

		pointVU[0] = point3d[0] * invnorm*rho;
		pointVU[1] = point3d[1] * invnorm*rho;
	}
	else
	{
		pointVU[0] = 0;
		pointVU[1] = 0;
	}
}

void Ovu2d3(double point3d[3], double pointVU[2], ocam_model *myocam_model)
{
	double *pol = myocam_model->pol;

	int length_pol = (myocam_model->length_pol);

	double xp = pointVU[0];
	double yp = pointVU[1];

	double r = sqrt(xp*xp + yp*yp); //distance [pixels] of  the point from the image center
	double zp = pol[0];
	double r_i = 1;
	int i;

	for (i = 1; i < length_pol; i++)
	{
		r_i *= r;
		zp += r_i*pol[i];
	}

	//normalize to unit norm
	double invnorm = 1 / sqrt(xp*xp + yp*yp + zp*zp);

	point3d[0] = invnorm*xp;
	point3d[1] = invnorm*yp;
	point3d[2] = invnorm*zp;
}

void Ocreate_panoramic(Mat mapx, Mat mapy, float Rmin, float Rmax, ocam_model *myocam_model)
{
	int licznik = 0; //licznik ---------------------/
					 //std::cout << "Rozpoczecie kreowania panoramy" <<std::endl; //liczni ---------------------/
	int i, j, k;
	float theta;
	int width = mapx.cols;	// szerokosc panoramy
	int height = mapx.rows; // wysokosc panoramy
	double wysokosc = 0;				// aktualna wysokosc dla zwyklych point3d
	double skalowanie_wysokosci = 0;	// skalowanie wysokosci
	double norm1 = 0;
	double norm2 = 0;
	double scal = 0;

	double * tab_z = new double[height];// tablica ze wspolrzednymi z dla kazdej wysokosci o dlugosci tablicy = height
	double x, y;							// wspolrzedne x,y dla point3d dla kazdej wysokosci
	double point2d_map[2];				// pkt ktory zostaje zmapowany na xmap i ymap
	double point3d_map[3];

	double point2d_max[2];				// punkt dla i=0 (max height)
	double point2d_min[2];				// punkt dla i=height_pan (min height)
	double pointVU_max[2];				// punkt dla i=0 (max height)
	double pointVU_min[2];				// punkt dla i=height_pan (min height)
	double point3d_rmax[3];				// punkt point3d dla max height
	double point3d_rmin[3];				// punkt point3d dla min height
										//std::cout << "Kreowanie panoramy 1" <<std::endl; //liczni ---------------------/

										//Obliczenia panoramy
										//petla (j = 0 ---> theta = 0)  ; (j = width_pan ---> theta = 2*pi)
	for (j = 0; j < width; j++)
	{

		// Dla zmiany obrotu panoramy (flip lewo-prawo) zmieniamy znak dla theta
		theta = -((float)j) / width * 2 * CV_PI; // kat theta dla danej kolumny na panoramie 
		pointVU_max[0] = Rmax*cos(theta);
		pointVU_max[1] = Rmax*sin(theta);
		pointVU_min[0] = Rmin*cos(theta);
		pointVU_min[1] = Rmin*sin(theta);
		/*if(licznik == 0)
		{
		std::cout << pointVU_min[0] << " x " << pointVU_min[1] << std::endl;
		std::cout << pointVU_max[0] << " x " << pointVU_max[1] << std::endl;
		}*/
		Ovu2d3(point3d_rmax, pointVU_max, myocam_model);
		Ovu2d3(point3d_rmin, pointVU_min, myocam_model);
		norm1 = sqrt(point3d_rmax[0] * point3d_rmax[0] + point3d_rmax[1] * point3d_rmax[1]); // dlugosc norm point3d_rmax(dlugosc w plaszczyznie x,y)
		norm2 = sqrt(point3d_rmin[0] * point3d_rmin[0] + point3d_rmin[1] * point3d_rmin[1]); // dlugosc norm point4d_rmin(dlugosc w plaszczyznie x,y)
		scal = norm1 / norm2;	// wspolczynnik skalujacy oba wektory do rownej dlugosci znormalizowanej (dlugosc w plaszczyznie x,y)
		point3d_rmin[0] = point3d_rmin[0] * scal;  // skalowanie wektora point3d_rmin
		point3d_rmin[1] = point3d_rmin[1] * scal;
		point3d_rmin[2] = point3d_rmin[2] * scal;

		//przeskalowanie punktow point3d do uzyskania kolumny o wysokosci danej panoramy %%    
		wysokosc = point3d_rmax[2] - point3d_rmin[2]; // roznica pomiedzy z dla obu punktow
		skalowanie_wysokosci = height / wysokosc; // skalowanie ktore rozszerzy z dla obu punktow
		point3d_rmax[0] = point3d_rmax[0] * skalowanie_wysokosci; // pkt point3d dla maksymalnej wysokosci
		point3d_rmax[1] = point3d_rmax[1] * skalowanie_wysokosci;
		point3d_rmax[2] = point3d_rmax[2] * skalowanie_wysokosci;
		point3d_rmin[0] = point3d_rmin[0] * skalowanie_wysokosci;	// pkt point3d dla minimalnej wysokosci
		point3d_rmin[1] = point3d_rmin[1] * skalowanie_wysokosci;
		point3d_rmin[2] = point3d_rmin[2] * skalowanie_wysokosci;
		/*if(licznik == 0)
		{
		std::cout << point3d_rmax[0] << " x " << point3d_rmax[1] << " x " << point3d_rmax[2] << std::endl;
		std::cout << point3d_rmin[0] << " x " << point3d_rmin[1] << " x " << point3d_rmin[2] << std::endl;
		}*/
		//Obliczenie punktow point3d w ilosci (n = height) %%
		//petla obliczajca wspolczynniki z dla kazdego point3d(x,y,z) dla
		//kazdego i - wysokosci na panoramie 
		for (k = 0; k < height; k++)
		{
			tab_z[k] = point3d_rmax[2] - (k);
		}
		point3d_map[0] = point3d_rmax[0]; // wspolczynniki x wszystkich point3d
		point3d_map[1] = point3d_rmax[1]; // wspolczynniki y wszystkich point3d

										  //Petla obliczajca point2d dla kazdego i - wysokosci na panoramie %%
										  // i = 0 ->(max height_pan) i = height_pan -> (min height_pan)
		for (i = 0; i < height; i++)
		{
			point3d_map[2] = tab_z[i];
			world2cam(point2d_map, point3d_map, myocam_model); // wywolanie funkcji skalujacej point3d na point2d
															   /* if(licznik == 0)
															   {
															   std::cout << point3d_map[0] << " x " << point3d_map[1] << " x " << point3d_map[2] << std::endl;
															   std::cout << point2d_map[0] << " x " << point2d_map[1] << std::endl;
															   }*/
			mapx.at<float>(i, j) = (float)point2d_map[1];	// zapisanie przeskalowania do mapx
			mapy.at<float>(i, j) = (float)point2d_map[0];	// zapisanie przeskalowania do mapy
			licznik++;
		}
		licznik++;


	}
	//   std::cout << "Zakonczone kreowanie panoramy" <<std::endl; //liczni ---------------------/
	delete[] tab_z;
}

void Ocreate_projection(Mat mapx, Mat mapy, double f, double fi, double theta, ocam_model *myocam_model)
{
	int licznik = 0; //licznik ---------------------/
					 //std::cout << "Rozpoczecie kreowania projection plane" <<std::endl; //liczni ---------------------/
	int i, j;
	int width = mapx.cols;	// szerokosc obrazu
	int height = mapx.rows; // wysokosc obrazu
	double point3D_Pc[3];
	double point3D_Pmap[3];
	double point2D_Pmap[2];
	double l = 0;
	double lx = 0;
	double Rkat_max = 0;
	double Rkat = 0;
	double Rtan = 0;
	///////////////////////////////////////////
	double Rtheta = (theta / 360.0) * 2 * CV_PI;
	double Rfi = (fi / 360.0) * 2 * CV_PI;
	point3D_Pc[0] = -cos(Rtheta)*cos(Rfi)*f;
	point3D_Pc[1] = -sin(Rtheta)*cos(Rfi)*f;
	point3D_Pc[2] = sin(Rfi)*f;
	//std::cout << "Punkt srodkowy plaszczyzny obrazu wirtualnej kamery: " << point3D_Pc[0] << " x " << point3D_Pc[1] << " x " << point3D_Pc[2] << std::endl;
	//////////////////////////////////////////
	//Jezeli theta = 0 i fi = 0
	/*
	if(fi == 0 && theta == 0)
	{
	std::cout << "Kat theta = 0 i fi = 0. Wykonywana jest uproszczona funkcja" std::<< endl
	for(i = 0 ; i<height ; i++)
	{
	for(j =0 ; j <width ; j++)
	{
	point2D_Pmap[0] =
	point2D_Pmap[1] =
	mapx.at<float>(i, j) = (float)point2D_Pmap[1];
	mapy.at<float>(i, j) = (float)point2D_Pmap[0];
	licznik++;
	}
	}
	}
	else
	{
	*/
	for (i = 0; i < height; i++)
	{
		l = cos(Rfi)*f - sin(Rfi)*(height / 2 - i);
		point3D_Pmap[2] = point3D_Pc[2] + cos(Rfi)*(height / 2 - i);

		for (j = 0; j < width; j++)
		{
			Rtan = (width / 2.0 - j) / l;
			Rkat = atan(Rtan);
			lx = l / cos(Rkat);
			point3D_Pmap[0] = -cos(Rtheta + Rkat)*lx;
			point3D_Pmap[1] = -sin(Rtheta + Rkat)*lx;
			/*if(licznik == 0)
			{
			std::cout << "Dla lewego gornego rogu i=0 oraz j=0 , punkt point3D_Pmap to:" << std::endl;
			std::cout << point3D_Pmap[0] << " x " << point3D_Pmap[1] << " x " << point3D_Pmap[2] << std::endl;
			}
			if(licznik == 100)
			{
			std::cout << "Dla punktu i=0 oraz j=100 , punkt point3D_Pmap to:" << std::endl;
			std::cout << point3D_Pmap[0] << " x " << point3D_Pmap[1] << " x " << point3D_Pmap[2] << std::endl;
			}*/
			world2cam(point2D_Pmap, point3D_Pmap, myocam_model);
			/*if((i == height/2 && j == width/2))
			{
			std::cout << "Dla Pc - i=height/2 oraz j=width/2 , punkt point2D_Pmap to:" << std::endl;
			std::cout << point2D_Pmap[0] << " x " << point2D_Pmap[1] << std::endl;
			}*/
			mapx.at<float>(i, j) = (float)point2D_Pmap[1];
			mapy.at<float>(i, j) = (float)point2D_Pmap[0];
			licznik++;
		} // for j

	} // for i
	  //} //else
	  //std::cout << "Zakonczone kreowanie projection plane" <<std::endl; //liczni ---------------------/

}

void Ocreate_panoramic_central(Mat mapx, Mat mapy, float Rprocent, float Rmax_new, ocam_model *myocam_model)
{
	int licznik = 0; //licznik ---------------------/
					 //std::cout << "Rozpoczecie kreowania panoramy " <<std::endl; //liczni ---------------------/
	int i, j, k;
	float theta;
	int width = mapx.cols;	// szerokosc panoramy
	int height = mapx.rows; // wysokosc panoramy
	double wysokosc = 0;				// aktualna wysokosc dla zwyklych point3d
	double skalowanie_wysokosci = 0;	// skalowanie wysokosci
	double norm1 = 0;
	double norm2 = 0;
	double scal = 0;
	double Rsrednie = 0;
	double Rmax = 0, Rmin = 0;

	double * tab_z = new double[height];// tablica ze wspolrzednymi z dla kazdej wysokosci o dlugosci tablicy = height
	double x, y;					// wspolrzedne x,y dla point3d dla kazdej wysokosci

	double point2d_map[2];				// pkt ktory zostaje zmapowany na xmap i ymap
	double point3d_map[3];

	double point2d_max[2];				// punkt dla i=0 (max height)
	double point2d_min[2];				// punkt dla i=height_pan (min height)
	double pointVU_max[2];				// punkt dla i=0 (max height)
	double pointVU_min[2];				// punkt dla i=height_pan (min height)
	double point3d_rmax[3];				// punkt point3d dla max height
	double point3d_rmin[3];				// punkt point3d dla min height
										//std::cout << "Kreowanie panoramy 1" <<std::endl; //liczni ---------------------/


										// Test i wyznaczenie Rmax i Rmin
	double pointVU_test1[2];
	double pointVU_test2[2];
	double point3d_test1[3];
	point3d_test1[0] = 400;
	point3d_test1[1] = 0;
	point3d_test1[2] = 0;

	double point3d_test2[3];
	point3d_test2[0] = 1800;
	point3d_test2[1] = 0;
	point3d_test2[2] = 0;
	O3d2vu(pointVU_test1, point3d_test1, myocam_model);
	O3d2vu(pointVU_test2, point3d_test2, myocam_model);

	//std::cout << "Punkty VU test1 i test2 maja wspolrzedne po kolei:" << std::endl;
	//std::cout << pointVU_test1[0] << " x " << pointVU_test1[1] << std::endl;
	//std::cout << pointVU_test2[0] << " x " << pointVU_test2[1] << std::endl;  

	Rsrednie = pointVU_test1[0];
	Rmax = (Rprocent / 100.0)*(Rmax_new - Rsrednie) + Rsrednie;
	pointVU_max[0] = Rmax;
	pointVU_max[1] = 0;
	Ovu2d3(point3d_rmax, pointVU_max, myocam_model);

	//std::cout << "Punkt 3d_rmax:" << std::endl;
	//std::cout << point3d_rmax[0] << " x " << point3d_rmax[1] << " x " << point3d_rmax[2] << std::endl;

	point3d_rmin[0] = point3d_rmax[0];
	point3d_rmin[1] = point3d_rmax[1];
	point3d_rmin[2] = -point3d_rmax[2];

	//std::cout << "Punkty 3d_rmin:" << std::endl;
	//std::cout << point3d_rmin[0] << " x " << point3d_rmin[1] << " x " << point3d_rmin[2] << std::endl;
	O3d2vu(pointVU_min, point3d_rmin, myocam_model);
	Rmin = pointVU_min[0];
	//Obliczenia panoramy
	//petla (j = 0 ---> theta = 0)  ; (j = width_pan ---> theta = 2*pi)
	for (j = 0; j < width; j++)
	{

		// Dla zmiany obrotu panoramy (flip lewo-prawo) zmieniamy znak dla theta
		theta = -((float)j) / width * 2 * CV_PI; // kat theta dla danej kolumny na panoramie 
		pointVU_max[0] = Rmax*cos(theta);
		pointVU_max[1] = Rmax*sin(theta);
		pointVU_min[0] = Rmin*cos(theta);
		pointVU_min[1] = Rmin*sin(theta);
		/*if(licznik == 0)
		{
		std::cout << pointVU_min[0] << " x " << pointVU_min[1] << std::endl;
		std::cout << pointVU_max[0] << " x " << pointVU_max[1] << std::endl;
		}*/
		Ovu2d3(point3d_rmax, pointVU_max, myocam_model);

		Ovu2d3(point3d_rmin, pointVU_min, myocam_model);
		norm1 = sqrt(point3d_rmax[0] * point3d_rmax[0] + point3d_rmax[1] * point3d_rmax[1]); // dlugosc norm point3d_rmax(dlugosc w plaszczyznie x,y)
		norm2 = sqrt(point3d_rmin[0] * point3d_rmin[0] + point3d_rmin[1] * point3d_rmin[1]); // dlugosc norm point4d_rmin(dlugosc w plaszczyznie x,y)

		scal = norm1 / norm2;	// wspolczynnik skalujacy oba wektory do rownej dlugosci znormalizowanej (dlugosc w plaszczyznie x,y)
		point3d_rmin[0] = point3d_rmin[0] * scal;  // skalowanie wektora point3d_rmin

		point3d_rmin[1] = point3d_rmin[1] * scal;
		point3d_rmin[2] = point3d_rmin[2] * scal;

		//przeskalowanie punktow point3d do uzyskania kolumny o wysokosci danej panoramy %%    
		wysokosc = point3d_rmax[2] - point3d_rmin[2]; // roznica pomiedzy z dla obu punktow
		skalowanie_wysokosci = height / wysokosc; // skalowanie ktore rozszerzy z dla obu punktow
		point3d_rmax[0] = point3d_rmax[0] * skalowanie_wysokosci; // pkt point3d dla maksymalnej wysokosci
		point3d_rmax[1] = point3d_rmax[1] * skalowanie_wysokosci;
		point3d_rmax[2] = point3d_rmax[2] * skalowanie_wysokosci;
		point3d_rmin[0] = point3d_rmin[0] * skalowanie_wysokosci;	// pkt point3d dla minimalnej wysokosci
		point3d_rmin[1] = point3d_rmin[1] * skalowanie_wysokosci;
		point3d_rmin[2] = point3d_rmin[2] * skalowanie_wysokosci;
		/*if(licznik == 0)
		{
		std::cout << point3d_rmax[0] << " x " << point3d_rmax[1] << " x " << point3d_rmax[2] << std::endl;
		std::cout << point3d_rmin[0] << " x " << point3d_rmin[1] << " x " << point3d_rmin[2] << std::endl;
		}*/
		//Obliczenie punktow point3d w ilosci (n = height) %%
		//petla obliczajca wspolczynniki z dla kazdego point3d(x,y,z) dla
		//kazdego i - wysokosci na panoramie 

		for (k = 0; k < height; k++)
		{
			tab_z[k] = point3d_rmax[2] - (k);
		}

		point3d_map[0] = point3d_rmax[0]; // wspolczynniki x wszystkich point3d
		point3d_map[1] = point3d_rmax[1]; // wspolczynniki y wszystkich point3d

										  //Petla obliczajca point2d dla kazdego i - wysokosci na panoramie %%
										  // i = 0 ->(max height_pan) i = height_pan -> (min height_pan)

		for (i = 0; i < height; i++)
		{
			point3d_map[2] = tab_z[i];
			world2cam(point2d_map, point3d_map, myocam_model); // wywolanie funkcji skalujacej point3d na point2d

															   /*if(licznik == 0)
															   {
															   std::cout << point3d_map[0] << " x " << point3d_map[1] << " x " << point3d_map[2] << std::endl;
															   std::cout << point2d_map[0] << " x " << point2d_map[1] << std::endl;

															   }*/
			mapx.at<float>(i, j) = (float)point2d_map[1];	// zapisanie przeskalowania do mapx
			mapy.at<float>(i, j) = (float)point2d_map[0];	// zapisanie przeskalowania do mapy
			licznik++;

		}
		licznik++;



	}
	//std::cout << "Zakonczone kreowanie panoramy" <<std::endl; //liczni ---------------------/
	delete[] tab_z;
}




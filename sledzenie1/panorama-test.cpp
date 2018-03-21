/*------------------------------------------------------------------------------
   Example code that shows the use of the 'cam2world" and 'world2cam" functions
   Shows also how to undistort images into perspective or panoramic images
   
   NOTE, IF YOU WANT TO SPEED UP THE REMAP FUNCTION I STRONGLY RECOMMEND TO INSTALL
   INTELL IPP LIBRARIES ( http://software.intel.com/en-us/intel-ipp/ )
   YOU JUST NEED TO INSTALL IT AND INCLUDE ipp.h IN YOUR PROGRAM
   
   Copyright (C) 2009 DAVIDE SCARAMUZZA, ETH Zurich  
   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
------------------------------------------------------------------------------*/

#include "panorama.h"

long GetFileSize(std::string filename)
{
	struct stat stat_buf;
	int rc = stat(filename.c_str(), &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
} 

void createPanorama()
{
	struct ocam_model o_cata; // our ocam_models for the fisheye and catadioptric cameras
	get_ocam_model(&o_cata, "./calib_results.txt");

	/*
	string folder = "C:/Users/Natalia/Desktop/STUDIA/In¿ynierka/Œledzenie obiektów/Praca/Praca/klatki";
	vector<String> filenames;

	glob(folder, filenames);
	cout << filenames.size() << endl;

	
	Mat myImage;
	IplImage *image2;

	for (int i = 0; i < filenames.size(); ++i)
	{
		//image = cvLoadImage(filenames[i]);
		cout << "klatka " << filenames[i] << endl;
		int c = cv::waitKey(1);
	}
	*/

	IplImage *image = cvLoadImage("C:/Users/Natalia/Desktop/STUDIA/In¿ynierka/Œledzenie obiektów/Praca/Praca/kalibracja.jpg");      // source image 2  
	CvSize size_pan_image = cvSize(1200, 400);        // size of the undistorted panoramic image
	IplImage *dst_pan = cvCreateImage(size_pan_image, 8, 3);    // undistorted panoramic image

	CvMat* mapx_pan = cvCreateMat(dst_pan->height, dst_pan->width, CV_32FC1);
	CvMat* mapy_pan = cvCreateMat(dst_pan->height, dst_pan->width, CV_32FC1);

	Mat mapx = cvarrToMat(mapx_pan);
	Mat mapy = cvarrToMat(mapy_pan);

	
	float Rmax = 350;  // the maximum radius of the region you would like to undistort into a panorama
	float Rmin = 120;   // the minimum radius of the region you would like to undistort into a panorama  
	create_panoramic_undistortion_LUT(mapx, mapy, Rmin, Rmax, o_cata.xc, o_cata.yc);

	//Ocreate_panoramic_central(mapx, mapy, 50, 90, &o_cata);
	
	//cvFlip(image, 0, 0);

	cvRemap(image, dst_pan, mapx_pan, mapy_pan, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll(0));

	//cvFlip(dst_pan, 0, 1);
	cvSaveImage("panorama.jpg", dst_pan);
	printf("\nImage %s saved\n", "panorama.jpg");
	

	cvReleaseImage(&image);
	cvReleaseImage(&dst_pan);
	cvReleaseMat(&mapx_pan);
	cvReleaseMat(&mapy_pan);

}

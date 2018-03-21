#include "OpticalFlow.h"

OpticalFlow::OpticalFlow()
{
	previous_frame_obj = new Frame();
	current_frame_obj = new Frame();
	
	previous_frame_obj->setDimensions();

	this->frame_number = 1;

	R = (previous_frame_obj->getFrameSizeY()) / 2;
	perimeter = 360;

	r_flow = Mat::zeros(R, perimeter, CV_64FC1);
	theta_flow = Mat::zeros(R, perimeter, CV_64FC1);

	//____________________________________________________________________________________________________________________________________________________
	//---Variables - "Tracking unknown moving targets on omnidirectional vision"----------------------
	/*
	magnitude_threshold = 1;

	alpha = 500;
	saved_u = Mat::zeros(R, perimeter, CV_64FC1);
	saved_v = Mat::zeros(R, perimeter, CV_64FC1);

	u = Mat::zeros(R, perimeter, CV_64FC1);
	v = Mat::zeros(R, perimeter, CV_64FC1);

	u_previous = Mat::zeros(R, perimeter, CV_64FC1);
	v_previous = Mat::zeros(R, perimeter, CV_64FC1);
	
	Ir = Mat::zeros(R, perimeter, CV_64FC1);
	Itheta = Mat::zeros(R, perimeter, CV_64FC1);
	It = Mat::zeros(R, perimeter, CV_64FC1);
	*/
	//-----------------
	mask = bsMask((previous_frame_obj->getFrameCenterX()), (previous_frame_obj->getFrameCenterY()), (previous_frame_obj->getFrameSizeX()), (previous_frame_obj->getFrameSizeY()) ); //utworzenie maski (macierz 0 i 1 powtórzonych 3x (3 kana³y)dla ka¿dego pixela)

	cart_2_polar_matrix = Mat(R, 360, CV_16UC2); //macierz [Rx360] 16-bit unsigned z 2 kana³ami, R kolumn, 360 wierszy i ka¿dy element w wierszu ma dwie wartoœci (2 kana³y)
	
	for (int theta = 0; theta < cart_2_polar_matrix.cols; theta++)
	{
		for (int r = 0; r<cart_2_polar_matrix.rows; r++)
		{
		//	Vec2s point_pixel = Vec2s(x_coordinate(r,theta) + (previous_frame_obj->getFrameCenterX()), y_coordinate(r,theta) + (previous_frame_obj->getFrameCenterY())); //wspó³rzêdne punktu
			Vec2s point_pixel = Vec2s(r*cos(CV_PI*theta / 180) + (previous_frame_obj->getFrameCenterX()), r*sin(CV_PI*theta / 180) + (previous_frame_obj->getFrameCenterY())); //wspó³rzêdne punktu

			cart_2_polar_matrix.at<Vec2s>(r, theta) = point_pixel; //w polu [r,theta] wpisz wpó³rzêdne pixela
			//cout << "Theta: " << theta << ",  r: " << r << "  Wspolrzedne: "<<point_pixel<<endl;
		}
	}
	//cout << cart_2_polar_matrix << endl;
	

}

OpticalFlow::~OpticalFlow()
{
}

Mat OpticalFlow::bsMask(int x, int y, int width, int height)
{
	Mat mask = Mat(height, width, CV_8UC3); //macierz [height x width] 8-bitowa unsigned int z 3 kana³ami
	for (int i = 0; i<mask.rows; i++)
	{
		for (int j = 0; j<mask.cols; j++)
		{
			if (((i - y)*(i - y) + (j - x)*(j - x))>R*R) //sprawdzenie czy znajduje sie poza okrêgiem 
			{
				mask.at<Vec3b>(i, j) = Vec3b(0, 0, 0); //je¿eli tak - 0
			}
			else
			{
				mask.at<Vec3b>(i, j) = Vec3b(1, 1, 1); //je¿eli nie - 1
			}

		}
	}
	//cout << mask << endl;
	return mask;
}

Mat OpticalFlow::transformOmniImage(Mat src_image)											//transformacja obrazu omni
{
	Mat image_cropped = mask.mul(src_image);										//utworzenie macierzy z pomno¿enia mask * src_image o wymiarze [x_mask x y_src_image] 
	Mat image_transformed = cart2PolarTransformation(image_cropped);
	return image_transformed;
}

Mat OpticalFlow::cart2PolarTransformation(Mat src_image)
{
	Mat image_transformed = Mat::zeros(cart_2_polar_matrix.rows, cart_2_polar_matrix.cols, CV_8UC3); //macierz zer o wymiarze [Rx360]

	for (int r = 0; r<image_transformed.rows; r++)
	{
		for (int theta = 0; theta<image_transformed.cols; theta++)
		{
			Vec2s pixel_coordinates; //vector containing coordinates (x,y) on the source_image of pixel to insert into image_transformed
			pixel_coordinates = cart_2_polar_matrix.at<Vec2s>(r, theta); //przypisanie zmiennej wspó³rzednych 
			image_transformed.at<Vec3b>(r, theta) = src_image.at<Vec3b>(pixel_coordinates[1], pixel_coordinates[0]);
			//cout <<"r: "<<r<<"  theta: "<<theta<<"  Image: "<< image_transformed.at<Vec3b>(r,theta) << endl;
		}
	}
	//cout << image_transformed << endl;
	return image_transformed;

}

void OpticalFlow::getPreviousAndCurrentFrame()
{
	this->previous_frame_obj->changePath(frame_number - 1);
	this->current_frame_obj->changePath(frame_number);
	
	this->previous_frame_obj->setMatrix();
	this->current_frame_obj->setMatrix();

	this->previous_frame_obj->setDimensions();
	this->current_frame_obj->setDimensions();

	cout << "Frame number   " << frame_number << endl;
}

int OpticalFlow::r_coordinate(int x, int y)
{
	int r ;
	r = sqrt((pow(x, 2)) + (pow(y, 2)));
	return r;
}

int OpticalFlow::theta_coordinate(int x, int y)
{
	int theta;
	theta = atan2(y, x) * 180 / CV_PI;
	if (theta < 0)
	{
		theta = 360 + theta;
	}
	return theta;
}

/*
int OpticalFlow::x_coordinate(int r, int theta)
{
int x;
x = r*cos(CV_PI*theta / 180);
return x;
}

int OpticalFlow::y_coordinate(int r, int theta)
{
int y;
y = r*sin(CV_PI*theta / 180);
return y;
}

*/

//____________________________________________________________________________________________________________________________________________________
//---Functions - "Tracking unknown moving targets on omnidirectional vision"----------------------

/*
void OpticalFlow::preprocessing(Mat src_image) //wydobycie krawêdzi przedmiotów
{
	Mat grad;
	Mat abs_Ir, abs_Itheta;
	GaussianBlur(src_image, src_image, Size(3, 3), 1.0);
	Sobel(src_image, Itheta_int, CV_16S, 1, 0, 3, 1, 0, BORDER_DEFAULT);
	Sobel(src_image, Ir_int, CV_16S, 0, 1, 3, 1, 0, BORDER_DEFAULT);
	
	//cout << "----------SOBEL---------------" << endl;
	//cout << "Itheta_int = " << Itheta_int << endl;
	//cout << "Ir_int = " << Ir_int << endl;
	//cout << "It_int = " << It_int << endl;

	//convertScaleAbs(Itheta_int, abs_Itheta); //wartoœc bezwzglêdna
	//convertScaleAbs(Ir_int, abs_Ir);
	//addWeighted(abs_Itheta, 0.5, abs_Ir, 0.5, 0, grad); //blending two images

	//cout << "Abs_Ir = " << abs_Ir << endl;
	//cout << "Abs_Itheta = " << abs_Itheta << endl;


	//imshow("grad", grad);
	//waitKey(0);
	

}


Mat OpticalFlow::calculate()
{
	
	Mat  original, cart_original;

	previous_frame_obj->getMatrix().copyTo(cart_original);
	previous_frame_obj->getPolarMatrix().copyTo(original);


	double check_vector = (0, 0);
	Mat current_frame_gray, previous_frame_gray;
	
	cvtColor(current_frame_obj->getMatrix(), current_frame_gray, COLOR_BGR2GRAY); //konwersja barwy na szary
	cvtColor(previous_frame_obj->getMatrix(), previous_frame_gray, COLOR_BGR2GRAY);

	preprocessing(current_frame_gray); //wyznaczenie gradientów Ir i Itheta
	
	It_int = current_frame_gray - previous_frame_gray;//wyznaczenie gradientu It - ró¿nica klatek current i previous

	Itheta_int.convertTo(Itheta, CV_64FC1);
	Ir_int.convertTo(Ir, CV_64FC1);
	It_int.convertTo(It, CV_64FC1);
	
	int x_1 = cart_original.cols / 2;
	int y_1 = cart_original.rows / 2;

	
	//cout << "-------------PO---------------" << endl;
	//cout << "Itheta = " << Itheta << endl;
	//cout << "Ir = " << Ir << endl;
	//cout << "It = " << It << endl;
	//cout << "-------------------------" << endl;
	//cout << "Itheta_int = " << Itheta_int << endl;
	//cout << "Ir_int = " << Ir_int << endl;
	//cout << "It_int = " << It_int << endl;
	

	
	//utworzenie macierzy
	Mat kernel_u = (Mat_<double>(3, 3) << 1, 1, 1, 1, 0, 1, 1, 1, 1) / 8; //normalizacja ~~ œrednia (8 bo bierzemy 8 pikseli)
	Mat kernel_v = (Mat_<double>(3, 3) << 0, 1, 0, 1, 0, 1, 0, 1, 0) / 4;
	Mat u_filtered, v_filtered;
	

	if (frame_number == 1)
	{
		
		saved_u = Mat::zeros(R, perimeter, CV_64FC1) ;
		saved_v = Mat::zeros(R, perimeter, CV_64FC1) ;
	
	}

	//zastosowanie filtru z powy¿szymi macierzami
	filter2D(saved_u, u_filtered, CV_64F, kernel_u);
	filter2D(saved_v, v_filtered, CV_64F, kernel_v);
	
	//u = 1 / 8 * u;
	//v = 1 / 4 * v;//normalizacja
	
	//cout << "-------------------------" << endl;
	//cout << "u = " << u << endl;
	
	

	for (int r = 0; r < u.rows; r++)
	{
		for (int theta = 0; theta < u.cols; theta++)
		{
			double u_val = u.at<double>(r, theta);
			double v_val = v.at<double>(r, theta);

			double Ir_val = Ir.at<double>(r, theta);
			double It_val = It.at<double>(r, theta);
			double Itheta_val = Itheta.at<double>(r, theta);
			double u_avg = u_filtered.at<double>(r, theta);
			double v_avg = v_filtered.at<double>(r, theta);

			u_val = u_avg - Ir_val * (Ir_val * u_avg + Itheta_val*v_avg + It_val) / (alpha*alpha + Ir_val*Ir_val + It_val*It_val);
			v_val = v_avg - Itheta_val * (Ir_val * u_avg + Itheta_val*v_avg + It_val) / (alpha*alpha + Ir_val*Ir_val + It_val*It_val);

			u.at<double>(r, theta) = u_val * 100;
			v.at<double>(r, theta) = v_val * 100;


		}
	}


	for (int r = 0; r < original.rows; r += 3)
	{
		for (int theta = 0; theta < original.cols; theta += 3)
		{
			int x_begin = r*cos(CV_PI*theta / 180) + x_1;
			int y_begin = r*sin(CV_PI*theta / 180) + y_1;


			int x_end = cvRound((r + v.at<double>(r, theta))*cos(CV_PI*(theta + u.at<double>(r, theta)) / 180) + x_1);
			int y_end = cvRound((r + v.at<double>(r, theta))*sin(CV_PI*(theta + u.at<double>(r, theta)) / 180) + y_1);

			// draw line at flow direction
			line(cart_original, Point(x_begin, y_begin), Point(x_end, y_end), Scalar(255, 0, 0));

			// draw initial point
			//circle(cart_original, Point(x_begin, y_begin), 1, Scalar(0, 0, 0), -1);


		}

	}

	saved_u = u;
	saved_v = v;

	//cout << u << endl;
	
	//imshow("u ", u_previous);
	//waitKey(0);
	//u_previous = Mat::zeros(R, perimeter, CV_64FC1);
	//frame_number++;
	//current_frame_gray.copyTo(previous_frame_gray);

			
			//if (u.at<double>(r, theta) != 0 && v.at<double>(r, theta) != 0)
			//{
			//	int x = r*cos(CV_PI*theta / 180) + img_center_x;
			//	int y = r*sin(CV_PI*theta / 180) + img_center_y;
			//	cout <<"["<< x << " ," << y<<"]" << endl;
			//
			//	cout << "Ir = " << Ir_val << endl;
			//	cout << "It = " << It_val << endl;
			//	cout << "Itheta = " << Itheta_val << endl;
			//	cout << "u dla r = " << r << " theta = " << theta << "   u = " << u.at<double>(r, theta) << endl;
			//	cout << "v dla r = " << r << " theta = " << theta << "   v = " << v.at<double>(r, theta) << endl;
			//
			//	cout << "____________________________________________" << endl;
			//}
			
			
		}

	}
	

	//cout << "u \n" << u << endl;
	//cout << "\n u_fil \n" << u_filtered << endl;
	//frame_number++;
	return cart_original;
}
*/

//_____________________________________________________________________________________________________________________________________________________
//--- Functions - Farneback Optical Flow----------------------------------------------------------

Mat OpticalFlow::showFarnebackOpticalFlow()
{
	//---Farneback Optical flow----------------------------------------------------

	Mat flow, original, cart_original;

	Mat current_frame_gray, previous_frame_gray;

	previous_frame_obj->getMatrix().copyTo(cart_original);


	//---PolarMatrix----------------------------------------------------------------

	previous_frame_obj->getPolarMatrix().copyTo(original);

	cvtColor(current_frame_obj->getPolarMatrix(), current_frame_gray, COLOR_BGR2GRAY); //konwersja barwy na szary
	cvtColor(previous_frame_obj->getPolarMatrix(), previous_frame_gray, COLOR_BGR2GRAY);

	UMat  flowUmat, prevgray;

	// calculate optical flow 
	calcOpticalFlowFarneback(previous_frame_gray, current_frame_gray, flowUmat, 0.4, 1, 12, 2, 8, 1.2, 0);
	
	// copy Umat container to standard Mat
	flowUmat.copyTo(flow);

	int x_1 = cart_original.cols / 2;
	int y_1 = cart_original.rows / 2;

	for (int r = 0; r < original.rows; r +=3)
	{
		for (int theta = 0; theta < original.cols; theta +=3)
		{

			// get the flow from y, x position * 10 for better visibility
			Point2f flowatxy = flow.at<Point2f>(r, theta) * 5; //flow(r,theta)[0].x=theta,  flow(r,theta)[1].y=r

			if (flowatxy.x < 1.2 && flowatxy.y < 1.2)
			{
				flowatxy.x = 0;
				flowatxy.y = 0;
			}
			 
			int x_begin = r*cos(CV_PI*theta / 180) + x_1;
			int y_begin = r*sin(CV_PI*theta / 180) + y_1;

				//int x_end = cvRound((r + r_flow.at<double>(r, theta))*cos(CV_PI*(theta + theta_flow.at<double>(r, theta)) / 180) + x_1);
				//int y_end = cvRound((r + r_flow.at<double>(r, theta))*sin(CV_PI*(theta + theta_flow.at<double>(r, theta)) / 180) + y_1);

			int x_end = cvRound((r + flowatxy.y)*cos(CV_PI*(theta + flowatxy.x) / 180) + x_1 );
			int y_end = cvRound((r + flowatxy.y)*sin(CV_PI*(theta + flowatxy.x) / 180) + y_1);

				// draw line at flow direction
				line(cart_original, Point(x_begin, y_begin), Point(x_end, y_end), Scalar(255, 0, 0));

				// draw initial point
				circle(cart_original, Point(x_begin, y_begin), 1, Scalar(0, 0, 0), -1);
				

				//line(original, Point(theta, r), Point(theta+flowatxy.x, r+flowatxy.y), Scalar(255, 0, 0));
				//circle(original, Point(theta, r), 1, Scalar(0, 0, 0), -1);

			
		}

	}

	// fill previous image again
	current_frame_gray.copyTo(previous_frame_gray);

	//frame_number++;
	//return original;
	return cart_original;

	
}

void OpticalFlow::calculateFarnebackOpticalFlow()
{
	//---Farneback Optical flow----------------------------------------------------

	Mat flow, original, cart_original;

	Mat current_frame_gray, previous_frame_gray;

	previous_frame_obj->getMatrix().copyTo(cart_original);

	//---PolarMatrix----------------------------------------------------------------

	previous_frame_obj->getPolarMatrix().copyTo(original);

	cvtColor(current_frame_obj->getPolarMatrix(), current_frame_gray, COLOR_BGR2GRAY);		//konwersja barwy na szary
	cvtColor(previous_frame_obj->getPolarMatrix(), previous_frame_gray, COLOR_BGR2GRAY);

	UMat  flowUmat, prevgray;

	// calculate optical flow 
	calcOpticalFlowFarneback(previous_frame_gray, current_frame_gray, flowUmat, 0.4, 1, 12, 2, 8, 1.2, 0);

	// copy Umat container to standard Mat
	flowUmat.copyTo(flow);

	for (int r = 0; r < original.rows; r++)
	{
		for (int theta = 0; theta < original.cols; theta++)
		{
			const Point2f flowatxy = flow.at<Point2f>(r, theta); // get the flow from y, x position 
			r_flow.at<double>(r, theta) = flowatxy.y;
			theta_flow.at<double>(r, theta) = flowatxy.x;	
		}

	}

	current_frame_gray.copyTo(previous_frame_gray);

}

void OpticalFlow::searchForMovingPixels()
{
	int i = 0;
	all_moving_pixels.push_back(movingPixel());

	int x_center = previous_frame_obj->getFrameCenterX();
	int y_center = previous_frame_obj->getFrameCenterY();

	for (int r = 0; r < current_frame_obj->getPolarMatrix().rows; r++)
	{
		for (int theta = 0; theta < current_frame_obj->getPolarMatrix().cols; theta++)
		{
			if (abs(r_flow.at<double>(r, theta)) > 0.328 || abs(theta_flow.at<double>(r, theta)) > 0.328)
			{
				all_moving_pixels[i].x = r*cos(CV_PI*theta / 180) + x_center;
				all_moving_pixels[i].y = r*sin(CV_PI*theta / 180) + y_center;
				all_moving_pixels[i].r = r;
				all_moving_pixels[i].theta = theta;
				all_moving_pixels[i].r_motion = r_flow.at<double>(r, theta);
				all_moving_pixels[i].theta_motion = theta_flow.at<double>(r, theta);
				i++;
				all_moving_pixels.push_back(movingPixel());
			}
		}
	}
	cout << "Moving pixels: " << all_moving_pixels.size() << endl;


	/*//draw all moving pixels
	Mat temp;
	current_frame_obj->getMatrix().copyTo(temp);
	for (int k = 0; k < all_moving_pixels.size(); k++)
	{
		int r1, t1;
		r1 = all_moving_pixels[k].r;
		t1 = all_moving_pixels[k].theta;
		int x1, y1;
		x1 = r1*cos(CV_PI*t1 / 180) + (previous_frame_obj->getFrameCenterX());
		y1 = r1*sin(CV_PI*t1 / 180) + (previous_frame_obj->getFrameCenterY());
		circle(temp, Point(x1, y1), 1, Scalar(0, 255, 0), -1);
	}
	imshow("d", temp);
	waitKey(0);
	
	*/

	
	//all_moving_pixels.clear();

	//for (int j = 0; j < all_moving_pixels.size(); j++)
		//cout << "r="<< all_moving_pixels[j].r <<"  theta="<<all_moving_pixels[j].theta<<"  r_flow="<<all_moving_pixels[j].r_motion<<"  theta_flow="<<all_moving_pixels[j].theta_motion<< endl;
	
}

//_____________________________________________________________________________________________________________________________________________________
//---Unused Functions----------------------------------------------------------------------------------------
/*
Mat OpticalFlow::calculateMagnitude()
{
	Mat magnitude = u.mul(u) + v.mul(v);
	Mat magnitude_thresholded;
	sqrt(magnitude, magnitude);

	magnitude *= 255;
	magnitude.convertTo(magnitude, CV_8UC1);

	threshold(magnitude, magnitude_thresholded, magnitude_threshold, 255, THRESH_BINARY);

	return magnitude_thresholded;
}

void OpticalFlow::nextframe()
{
	frame_number++;
}

Mat calculateIt(Mat frame_previous, Mat frame_current)
{
	Mat newIt;
	newIt = frame_current - frame_previous;
	return newIt;
}

void performOpticalFlow(Mat frame_previous, Mat frame_current, Mat frame_next, OpticalFlow object_frame_previous, OpticalFlow object_frame_current)
{
	Mat omni_image_previous, omni_image_current, cos, current_previous, next_current, a, b, c;
	a = sobel1(frame_previous);
	b = sobel1(frame_current);
	c = sobel1(frame_next);
	current_previous = b - a;
	next_current = c - b;

	//current_next = frame_current - frame_next;
	//cout << "PREVIOUS \n" << previous_current << endl;
	//cout << "CURRENT \n" << current_next << endl;
	//convertScaleAbs(previous_current, previous_current);
	//imshow("previous_current", previous_current);
	//imshow("current_next", current_next);
	//waitKey(0);

	//omni_image_previous = object_frame_previous.transformOmniImage(frame_previous); //wype³nienie macierzy Rx360 wartosciami BGR z danej klatki
	//omni_image_current = object_frame_current.transformOmniImage(frame_current);
	//object_frame_current.calculate(omni_image_previous, omni_image_current);

	//cos = polar2CartTraformation(frame_current, omni_image_current);
	//cout << cos << endl;

	//cout << "PREVIOUS \n" << omni_image_previous << endl;
	//cout << "CURRENT \n" << omni_image_current << endl;

}

Mat polar2CartTraformation(Mat cart_source, Mat polar_transformed) //przekazuje macierz [x,y] wype³nion¹ wartoœciami BGR z macierzy wspó³rzêdnych biegunowych
{
	Mat transformed_cart_matrix = Mat::zeros(cart_source.rows, cart_source.cols, CV_8UC3);
	int x_1 = cart_source.cols / 2;
	int y_1 = cart_source.rows / 2;

	for (int r = 0; r < polar_transformed.rows; r++)
	{
		for (int theta = 0; theta < polar_transformed.cols; theta++)
		{
			int x = r*cos(CV_PI*theta / 180) + x_1;
			int y = r*sin(CV_PI*theta / 180) + y_1;
			transformed_cart_matrix.at<Vec3b>(x, y) = polar_transformed.at<Vec3b>(r, theta);
		}
	}

	return transformed_cart_matrix;
}

Mat polarToCartVector(Mat cart_source, Mat polar_transformed)
{
	Mat transformed_cart_matrix = Mat::zeros(cart_source.rows, cart_source.cols, CV_8UC2);
	int x_1 = cart_source.cols / 2;
	int y_1 = cart_source.rows / 2;


	return transformed_cart_matrix;
}

Mat sobel1(Mat src_image)
{
	cvtColor(src_image, src_image, COLOR_BGR2GRAY);
	Mat Itheta,Ir,abs_Ir, abs_Itheta, grad;
	GaussianBlur(src_image, src_image, Size(3, 3), 1.0);
	Sobel(src_image, Itheta, CV_16S, 1, 0, 3, 1, 0, BORDER_DEFAULT);
	Sobel(src_image, Ir, CV_16S, 0, 1, 3, 1, 0, BORDER_DEFAULT);

	convertScaleAbs(Itheta, abs_Itheta); //wartoœc bezwzglêdna
	convertScaleAbs(Ir, abs_Ir);
	addWeighted(abs_Itheta, 0.5, abs_Ir, 0.5, 0, grad); //blending two images
	return grad;
}
*/
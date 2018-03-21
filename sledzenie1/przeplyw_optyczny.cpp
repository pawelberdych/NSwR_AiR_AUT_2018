#include "OpticalFlow.h"
#include "ParticleFilter.h"
//#include "TrackedObject.h"

vector<int> xs;
vector<int> ys;

int point_counter = 0;
Rect area;

void mouseCallback(int event, int x, int y, int flags, void* userdata)
{
	
	if (event == EVENT_LBUTTONDOWN)
	{
		
		if (point_counter < 1)
		{
			xs.push_back(x);
			ys.push_back(y);

			point_counter++;
			
		}
		else
		{
			xs.push_back(x);
			ys.push_back(y);	

			area = Rect(Point(xs.at(0), ys.at(0)), Point(xs.at(1), ys.at(1))); //wyznaczenie prostokata
			point_counter++;
			destroyWindow("Click in two corners of object to track");
		}

	}
}

int main(int, char)
{
	

	
	OpticalFlow obj;	//utworzenie obiektu klasy OpticalFlow

	//createPanorama();	//utworzenie panoramy
	
	VideoWriter tracking("tracking.avi", CV_FOURCC('M', 'J', 'P', 'G'), 10, Size(1280, 720));
	
	//VideoWriter flow("flow.avi", CV_FOURCC('M', 'J', 'P', 'G'), 10, Size(1280, 720));

	clock_t start, stop;
	start = clock();

	while (obj.getFrameNumber() < obj.current_frame_obj->getIndex()) //15)
	{
		obj.current_frame_obj->getFramePaths();		//nadanie œcie¿ki obiektowi current_frame_obj klasy Frame (wykorzstanie metody klasy Film)
		obj.previous_frame_obj->getFramePaths();	//nadanie œcie¿ki obiektowi previous_frame_obj klasy Frame (wykorzstanie metody klasy Film)
		obj.getPreviousAndCurrentFrame();			//pobranie obrazów klatek previous i current

		//Mat fullImageHSV;
		//cvtColor(obj.previous_frame_obj->getMatrix(), fullImageHSV, CV_BGR2HSV);
		//imshow("HSV", fullImageHSV);
		//waitKey(0);
	
		
											
		if (obj.getFrameNumber() == 1)
		{
			namedWindow("Click in two corners of object to track");
			setMouseCallback("Click in two corners of object to track", mouseCallback, NULL);
			imshow("Click in two corners of object to track", obj.previous_frame_obj->getMatrix());		
			waitKey(0);
			///in that moment you should select top left and then bottom right corner of object to track
		}		
		
		
		
		obj.previous_frame_obj->setPolarMatrix(obj.transformOmniImage(obj.previous_frame_obj->getMatrix())); //utworzenie obrazu omni
		obj.current_frame_obj->setPolarMatrix(obj.transformOmniImage(obj.current_frame_obj->getMatrix()));

	
		//Mat cos = obj.calculate(); //---"Tracking unknown moving targets on omnidirectional vision"-----
	
		//Mat cos = obj.showFarnebackOpticalFlow(); //wygenerowanie i wyœwietlenie przep³ywu optycznego
		
		obj.calculateFarnebackOpticalFlow();	//wygenerowanie przep³ywu
		obj.searchForMovingPixels();	//wyodrêbnienie poruszaj¹cych siê pikseli
	
		ParticleFilter particle(obj.previous_frame_obj->getFrameSizeX(), obj.previous_frame_obj->getFrameSizeY(), area, obj.previous_frame_obj->getPolarMatrix(), obj.getAllMovingPixels());
		
		particle.createParticles(obj.current_frame_obj->getMatrix(), obj.previous_frame_obj->getFrameSizeY(), obj.getAllMovingPixels());
		particle.createParticles2(obj.current_frame_obj->getMatrix(), cvRound(obj.previous_frame_obj->getFrameSizeY() * 0.5) , obj.getAllMovingPixels());
		particle.createParticles3(obj.current_frame_obj->getMatrix(), cvRound(obj.previous_frame_obj->getFrameSizeY() * 0.25), obj.getAllMovingPixels());
		
		area = particle.getFoundObject(); //przypisanie jako obiekt wybranej cz¹steczki
		Mat object = particle.showTracking(obj.current_frame_obj->getMatrix()); //pokazanie klatki z zaznaczonym obiektem

		obj.all_moving_pixels.clear(); //wyczyszczenie wektora all_moving_pixels
		particle.clearTargetVector();	//wyczyszczenie wektora target
		
		obj.nextframe();
		
		//flow.write(cos); //zapisanie klatki z przep³ywem optycznym do filmu
		
		tracking.write(object); //zapisanie klatki z obiektem do filmu
		
		
		
	}
	stop = clock();
	cout << "Czas wykonania: " << ((stop - start) / (double)CLOCKS_PER_SEC) <<" s"<< endl;
	
	//---Wyœwietlenie filmu---
	VideoCapture display_tracking("tracking.avi");

	if (!display_tracking.isOpened())  // isOpened() returns true if capturing has been initialized.
	{
		cout << "Cannot open the video file. \n";
		return -1;
	}

	double fps = display_tracking.get(CV_CAP_PROP_FPS);
	int odstep_miedzy_klatkami = 1000 / fps;
	namedWindow("Tracking", CV_WINDOW_AUTOSIZE);
	while (1)
	{
		Mat frame;

		if (!display_tracking.read(frame)) 
		{
			cout << "\n Cannot read the video file. \n";
			break;
		}

		imshow("Tracking", frame);
		waitKey(odstep_miedzy_klatkami);
		if (waitKey(30) == 27) // Wait for 'esc' key press to exit
		{
			break;
		}
	}

	system("pause");
	return 0;
}


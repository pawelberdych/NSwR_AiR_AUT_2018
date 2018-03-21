#include "OpticalFlow.h"

Film::Film()
{
	//video_path = "C:/Users/Natalia/Desktop/STUDIA/In¿ynierka/Œledzenie obiektów/Praca/Praca/przep³yw optyczny/film.mp4";
	//video_path = "C:\\NSwR\\Filmy\\tracking_ja_133_cz120.avi";
	video_path = "xxx.avi";
	//this->sequence_path = "C:\\Users\\Natalia\\Desktop\\STUDIA\\In¿ynierka\\Œledzenie obiektów\\Praca\\Praca\\ja";
	this->sequence_path = "C:\\NSwR\\sledzenie1\\sledzenie1\\klatki3";
	
	//createFrames();
	getFramePaths();
}

Film::~Film()
{

}

/*
void Film::createFrames() //tworzy klatki z wczytanego filmu i zapisuje je
{
	
	// odczytanie pliku 
	CvCapture* vid = cvCreateFileCapture(video_path);

	// odczytanie pierwszej klatki - niezbedne do prawidlowego odczytania wlasciwosci pliku
	// przy uzyciu funkcji cvGetCaptureProperty
	cvQueryFrame(vid);

	// odczytujemy z wlasciwosci pliku liczbe klatek na sekunde
	double fps = cvGetCaptureProperty(vid, CV_CAP_PROP_FPS);
	
	int i = 000; //numerator klatek
	char plik[50];

	// wyliczamy czas potrzebny do odtwarzania pliku z prawidlowa prêdkoscia
	int odstep_miedzy_klatkami = 1000 / fps;


	while (true)
	{
		// pobranie kolejnej klatki
		IplImage* klatka = cvQueryFrame(vid);

		// jezeli nie ma wiêcej klatek to przerwij
		if (klatka == 0)
			break;

		if (i>32000)//i = 000;
		{
			stringstream filename;

			filename <<  "klatka_" << setw(4) << std::setfill('0') << ++i << ".png";//przygotuj nazwê nastêpnego pliku png
			sprintf_s( plik, filename.str().c_str());
			cout << plik << endl;
		}
		
		if (!cvSaveImage(plik, klatka))//gdy zapis siê nie powiedzie( w œrodku) to
			cout << "Nie uda³o siê zapisaæ klatki do pliku - " << plik << endl;
	}

}
*/


void Film::createFrames() //tworzy klatki z wczytanego filmu i zapisuje je
{
	// odczytanie pliku avi
	CvCapture* vid = cvCreateFileCapture(video_path);

	// odczytanie pierwszej klatki - niezbedne do prawidlowego odczytania wlasciwosci pliku
	// przy uzyciu funkcji cvGetCaptureProperty
	cvQueryFrame(vid);

	// odczytujemy z wlasciwosci pliku liczbe klatek na sekunde
	double fps = cvGetCaptureProperty(vid, CV_CAP_PROP_FPS);
	int i = 000; //numerator klatek
	char plik[50];

	// wyliczamy czas potrzebny do odtwarzania pliku z prawidlowa prêdkoscia
	int odstep_miedzy_klatkami = 1000 / fps;


	while (true)
	{
		// pobranie kolejnej klatki
		IplImage* klatka = cvQueryFrame(vid);
		// jezeli nie ma wiêcej klatek to przerwij
		cout << klatka << endl;
		if (klatka == 0)
			break;

		if (i>32000)i = 000;
		{
			stringstream filename;

			filename << "klatka_" << setw(4) << std::setfill('0') << ++i << ".png";//przygotuj nazwê nastêpnego pliku png
			sprintf_s(plik, filename.str().c_str());
		}

		if (!cvSaveImage(plik, klatka))//gdy zapis siê nie powiedzie( w œrodku) to
			cout << "Nie uda³o siê zapisaæ klatki do pliku - " << plik << endl;
	}

}
void Film::getFramePaths() //pobiera wczeœniej utworzone klatki
{
	// Prepare a list of all images
	cv::glob(sequence_path, image_path, false);

	//for (int i = 0; i < image_path.size(); i++)
		//cout << image_path[i] << endl;
}

int Film::getIndex() 
{ 
	return image_path.size(); 
}
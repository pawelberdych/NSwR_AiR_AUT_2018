#pragma once
#ifndef FRAME
#define FRAME

#include "libraries.h"
#include "Film.h"

class Frame : public Film
{
private:
	string frame_path;		//œcie¿ka do danej klatki
	Mat frame_matrix;		//macierz kartezjañska [x,y] obrazu klatki
	Mat polar_frame_matrix;	//macierz biegunowa [360,R] klatki
	int frame_size_x;		//szerokoœæ klatki
	int frame_size_y;		//wysokoœæ klatki
	int frame_center_x;		//centrum obrazu - x
	int frame_center_y;		//centrum obrazu - y
	

public:
	
	Frame();
	~Frame();

	void changePath(int i);

	void setMatrix();
	Mat getMatrix() { return frame_matrix; }			//zwraca macierz kartezjañsk¹

	void setPolarMatrix(Mat matrix);
	Mat getPolarMatrix() { return polar_frame_matrix; }	//zwraca macierz biegunow¹ 
	
	void setDimensions();
	int getFrameSizeX() { return frame_size_x; }		//zwraca szerokoœæ klatki
	int getFrameSizeY() { return frame_size_y; }		//zwraca wysokoœc klatki
	int getFrameCenterX() { return frame_center_x; }	//zwraca œrodek klatki w osi X
	int getFrameCenterY() { return frame_center_y; }	//zwraca œrodek klatki w osi Y
};

#endif

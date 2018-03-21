#pragma once
#ifndef FILM
#define FILM

#include "libraries.h"

class Film
{
protected:
	const char *video_path;
	string sequence_path; 
	vector<String> image_path;
	int index;

public:
	Film();
	~Film();
	void createFrames();			//tworzy klatki z wczytanego filmu i zapisuje je
	void getFramePaths();			//pobiera klatki
	int getIndex();					//zwraca liczbe klatek
};

#endif

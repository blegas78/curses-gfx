#include "curses-gfx-handler.h"

void *cursesHandlerRender(void *data) {
	CursesGfxHandler* This = (CursesGfxHandler*)data;
	
	return NULL;
}


CursesGfxHandler::CursesGfxHandler() {
	threads = new pthread_t[1];
	totalThreads = 1;
	busyThreads = false;
}

CursesGfxHandler::~CursesGfxHandler() {
	delete [] threads;
}


void CursesGfxHandler::setTotalThreads(int numthreads) {
	if (!busyThreads) {
		threads = new pthread_t[numthreads];
		totalThreads = numthreads;
	}
}

#pragma once

#include <vector>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include <queue>
#include "stdafx.h"

class ImgCacher
{
public:
	ImgCacher(int width, int height, int endNum, BSTR dir);
	~ImgCacher();

	void loadImages();
	std::shared_ptr<float> getImages(int i);

	void stop(){ mBStop = true; }

private:
	int mWidth;
	int mHeight;
	BSTR mDir;
	std::queue<int> mFilesToRead;
	std::map<int, std::shared_ptr<float>> mFiles;
	bool mBStop;
	int mEndNum;
};


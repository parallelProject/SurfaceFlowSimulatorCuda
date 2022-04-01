#include "ImgCacher.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "shapefil.h"

#include "ComBase.h"

ImgCacher::ImgCacher(int width, int height, int endNum, BSTR dir)
{
	mWidth = width;
	mHeight = height;
	mDir = dir;
	mBStop = false;
	mEndNum = endNum;
	for (int i = 0; i < 50; i++) {
		mFilesToRead.push(i);
	}
	loadImages();
}

ImgCacher::~ImgCacher()
{

}

void ImgCacher::loadImages()
{
	CComBase combase;

	std::thread t([&](){
		while (!mBStop) {
			if (mFilesToRead.size() == 0) {
				Sleep(50);
				continue;
			}
			int index = mFilesToRead.front();
			int endindex = mFilesToRead.back();
			if (endindex < mEndNum - 1) {
				++endindex;
				mFilesToRead.push(endindex);
			}
			mFilesToRead.pop();

			float* val = new float[mWidth * mHeight];
			combase.OpenDirfile(mDir, index, mWidth, mHeight, val/*.get()*/);
			std::shared_ptr<float> valptr(val);
			mFiles[index] = valptr;
		}
	});
	t.detach();
}

std::shared_ptr<float> ImgCacher::getImages(int i)
{
	bool bWait = false;
	while (mFiles.count(i) == 0) {
		//printf("----------------------------waiting for read!------------------------------\n");
		Sleep(50);
	}
	std::shared_ptr<float> img = mFiles[i];
	mFiles.erase(i);
	return img;
}
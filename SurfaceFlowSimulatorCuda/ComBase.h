// ComBase.h: interface for the CComBase class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_COMBASE_H__C606D7F6_9F66_4F4C_B994_A70CAE97BD92__INCLUDED_)
#define AFX_COMBASE_H__C606D7F6_9F66_4F4C_B994_A70CAE97BD92__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gdal.h"
#include "gdal_priv.h"
#include "shapefil.h"

#include <vector>
#include <algorithm>

#define DBL_MAX         1.7976931348623158e+308 /* max value */

struct point    
{
    double x,y;

    point(){
    };

    point(double X,double Y):x(X),y(Y){
    }

    bool operator == (const point &a)const{
        return x==a.x && y==a.y;
    }

    bool operator != (const point &a)const{
        return x!=a.x || y!=a.y;
    }

    bool operator < (const point &a)const{
        return x==a.x ? y<a.y : x<a.x;
    }

       bool operator > (const point &a)const{
        return x==a.x ? y>a.y : x>a.x;
    }   
};

struct point3D
{
    double x,y,z;
    double val;  //按照val进行排序

    point3D(){
    };

    point3D(double X,double Y,double Z):x(X),y(Y),z(Z){
    }

    point3D(double X,double Y,double Z,double value):x(X),y(Y),z(Z),val(value){
    }

    bool operator == (const point3D &a)const{
        return x==a.x && y==a.y && z==a.z;
    }

    bool operator != (const point3D &a)const{
        return x!=a.x || y!=a.y ||z!=a.z ;
    }

    bool operator < (const point3D &a)const{
        return val<a.val ;
    }

       bool operator > (const point3D &a)const{
        return val>a.val ;
    }   


};

typedef std::vector<int>      intArray;
typedef std::vector<float>    floatArray;
typedef std::vector<point>    polyline;   //定义一条折线。
typedef std::vector<point3D>  polyline3D; //定义一条折线。

class CComBase  
{
public:
	CComBase();
	virtual ~CComBase();

public:
	//获得影像的后缀名。
	CString GetImgFormat(LPCTSTR lpstrFileName);

	//通过影像的后缀名得到对应GDAL能处理的类型名
	CString GetGDALFormat(CString imgFormat);

	//释放GDALDataset数据集。
	void  RELEASE(GDALDataset* pData);



    //判断输入的strDEMName是否有效，此处的imgName支持*.tif & *.img格式。
    BOOL  bImgNameVerf(CString strImgName,CString& imgFormat,CString& strFormat);

    //创建一个新的影像。
    BOOL  CreateNewImg(BSTR ImgName, int imgWidth,int imgHeight, double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal,  float *demZ);

    BOOL  CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal, float *demZ);

	BOOL  CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal, int *demZ);

    //创建一个新的影像(带投影信息)。
    BOOL  CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal, CString projRef, float *demZ);

	
    //打开一影像数据，获得该影像的高、宽信息。
    BOOL  OpenImg(BSTR ImgName, int& imgWidth, int& imgHeight);

	//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
	BOOL  OpenImg(CString strImgName, int& imgWidth, int& imgHeight);

	//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
	BOOL  OpenImg(CString strImgName,int imgWidth, int imgHeight, float *pBuffer);

	//打开一影像数据，获得该影像的信息。
	BOOL  OpenImg(BSTR ImgName,int& imgWidth, int& imgHeight,double& dx, double& dy, double& Xmin, double& Ymax);

    //打开一影像数据，获得该影像的信息。
    BOOL  OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, 
                  double& Xmin, double& Ymax, float *pBuffer);

	//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
	BOOL  OpenImg(CString strImgName,int imgWidth, int imgHeight,double& dx, double& dy, 
		          double& Xmin, double& Ymax, float *pBuffer);

    //打开一影像数据，获得该影像的信息(包括投影信息)。
    BOOL  OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, 
                  double& Xmin, double& Ymax, CString& projRef, float *pBuffer);

	//打开BMP影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
	BOOL  OpenBMP(BSTR ImgName, int& imgWidth, int& imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo);

	//打开BMP影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
	BOOL  OpenBMP(BSTR ImgName, int imgWidth, int imgHeight, BYTE *pBuffer);


	//判断点(cenX,cenY,cenZ)是否在三角形(xt,yt,zt)的边界上。判断一个点到线两个端点距离和是否等于线的距离。
	BOOL IsInEdge(double cenX,double cenY,double* xt,double* yt);

	//判断点(cenX,cenY,cenZ)是否在三角形(xt,yt,zt)的内(包括边界)。判断一个点到三角形各个顶点的面积和是否等于三角形面积。
	BOOL IsInTriangle(double cenX,double cenY,double* xt,double* yt);

	//已知三顶点坐标求三角形的面积。
	double AreaTrig(double x1, double y1, double x2, double y2, double x3, double y3);

	//打开一个目录InDir,返回该目录下的文件总数目。
	int  OpenDirfile(BSTR InDir);

	//打开一个目录InDir，获得其中第fileNum幅影像的数据信息，放入pVal中。imgWidth, imgHeight为该影像的宽度和高度。
	BOOL  OpenDirfile(BSTR InDir, int fileNum, int imgWidth, int imgHeight, float* pVal);

	//保存BMP影像数据。返回True说明保存成功，否则保存失败。
	BOOL  SaveBMP(CString strImgName, int imgWidth, int imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo, BYTE *pBuffer);
};

#endif // !defined(AFX_COMBASE_H__C606D7F6_9F66_4F4C_B994_A70CAE97BD92__INCLUDED_)

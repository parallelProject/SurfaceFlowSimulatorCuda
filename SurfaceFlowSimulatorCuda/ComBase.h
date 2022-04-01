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
    double val;  //����val��������

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
typedef std::vector<point>    polyline;   //����һ�����ߡ�
typedef std::vector<point3D>  polyline3D; //����һ�����ߡ�

class CComBase  
{
public:
	CComBase();
	virtual ~CComBase();

public:
	//���Ӱ��ĺ�׺����
	CString GetImgFormat(LPCTSTR lpstrFileName);

	//ͨ��Ӱ��ĺ�׺���õ���ӦGDAL�ܴ����������
	CString GetGDALFormat(CString imgFormat);

	//�ͷ�GDALDataset���ݼ���
	void  RELEASE(GDALDataset* pData);



    //�ж������strDEMName�Ƿ���Ч���˴���imgName֧��*.tif & *.img��ʽ��
    BOOL  bImgNameVerf(CString strImgName,CString& imgFormat,CString& strFormat);

    //����һ���µ�Ӱ��
    BOOL  CreateNewImg(BSTR ImgName, int imgWidth,int imgHeight, double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal,  float *demZ);

    BOOL  CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal, float *demZ);

	BOOL  CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal, int *demZ);

    //����һ���µ�Ӱ��(��ͶӰ��Ϣ)��
    BOOL  CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                       double dx, double dy, double invalidVal, CString projRef, float *demZ);

	
    //��һӰ�����ݣ���ø�Ӱ��ĸߡ�����Ϣ��
    BOOL  OpenImg(BSTR ImgName, int& imgWidth, int& imgHeight);

	//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
	BOOL  OpenImg(CString strImgName, int& imgWidth, int& imgHeight);

	//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
	BOOL  OpenImg(CString strImgName,int imgWidth, int imgHeight, float *pBuffer);

	//��һӰ�����ݣ���ø�Ӱ�����Ϣ��
	BOOL  OpenImg(BSTR ImgName,int& imgWidth, int& imgHeight,double& dx, double& dy, double& Xmin, double& Ymax);

    //��һӰ�����ݣ���ø�Ӱ�����Ϣ��
    BOOL  OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, 
                  double& Xmin, double& Ymax, float *pBuffer);

	//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
	BOOL  OpenImg(CString strImgName,int imgWidth, int imgHeight,double& dx, double& dy, 
		          double& Xmin, double& Ymax, float *pBuffer);

    //��һӰ�����ݣ���ø�Ӱ�����Ϣ(����ͶӰ��Ϣ)��
    BOOL  OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, 
                  double& Xmin, double& Ymax, CString& projRef, float *pBuffer);

	//��BMPӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
	BOOL  OpenBMP(BSTR ImgName, int& imgWidth, int& imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo);

	//��BMPӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
	BOOL  OpenBMP(BSTR ImgName, int imgWidth, int imgHeight, BYTE *pBuffer);


	//�жϵ�(cenX,cenY,cenZ)�Ƿ���������(xt,yt,zt)�ı߽��ϡ��ж�һ���㵽�������˵������Ƿ�����ߵľ��롣
	BOOL IsInEdge(double cenX,double cenY,double* xt,double* yt);

	//�жϵ�(cenX,cenY,cenZ)�Ƿ���������(xt,yt,zt)����(�����߽�)���ж�һ���㵽�����θ��������������Ƿ���������������
	BOOL IsInTriangle(double cenX,double cenY,double* xt,double* yt);

	//��֪�����������������ε������
	double AreaTrig(double x1, double y1, double x2, double y2, double x3, double y3);

	//��һ��Ŀ¼InDir,���ظ�Ŀ¼�µ��ļ�����Ŀ��
	int  OpenDirfile(BSTR InDir);

	//��һ��Ŀ¼InDir��������е�fileNum��Ӱ���������Ϣ������pVal�С�imgWidth, imgHeightΪ��Ӱ��Ŀ�Ⱥ͸߶ȡ�
	BOOL  OpenDirfile(BSTR InDir, int fileNum, int imgWidth, int imgHeight, float* pVal);

	//����BMPӰ�����ݡ�����True˵������ɹ������򱣴�ʧ�ܡ�
	BOOL  SaveBMP(CString strImgName, int imgWidth, int imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo, BYTE *pBuffer);
};

#endif // !defined(AFX_COMBASE_H__C606D7F6_9F66_4F4C_B994_A70CAE97BD92__INCLUDED_)

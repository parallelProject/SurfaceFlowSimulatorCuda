#if !defined( _dataType_h_INC_ )

#define _dataType_h_INC_


#include "TINClass.h"



struct DEMInfo {
	double x0,y0,dirAng;
	short nx,ny;
	float dx,dy;
}	;	



//��
struct POINT3D {
	double x,y,z;
	triPOINT *tinPt;
};

//�㼯
//��
// 
enum objTYPE {
	typePOINTSET,
	typePOLYLINE
};



struct geoOBJ {
	BYTE type;
	BYTE delFlag;
	BYTE color;
	short pSum;
	POINT3D *pts;
};



struct rectREGISTER	{
	int maxSum;
	int objSum;
	geoOBJ **objs;
} ;



//�ļ���ʽ

struct fileHEAD {
	char tag[32];
	long objSum;
	long fileSize;
	double xOffset, yOffset;
}	;

//�����ε���������
struct triVERTEXES {
	triPOINT *vertex[3];
};


//�����ָ�봮
struct InsertPOINT {
	triPOINT *insertpt;
	InsertPOINT *father;
	double  distance;//�ۼӾ���
} ;


//�ǵ�����
struct cornerPOINT {
	triPOINT *triPoint;
	InsertPOINT *insPoint;
	double  distance;//�ۼӾ���
} ;


/////////////////////////////



#endif
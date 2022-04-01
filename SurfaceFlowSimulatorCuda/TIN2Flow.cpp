// TIN2Flow.cpp: implementation of the TIN2Flow class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include "TIN2Flow.h"
#include "shapefil.h"
#include <fstream>
#include <math.h>
#include <algorithm>
#include "tin.h"
#include "time.h"
#include <float.h>

#include <string>
#include <cstringt.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define PI 3.1415926

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TIN2Flow::TIN2Flow()
{
}

TIN2Flow::~TIN2Flow()
{
}



//���ݽڵ�a�ͽڵ�b��xy�Ƚϴ�С��
static inline bool cmp_byxy(node* a,node* b)
{
    if(a->m_x!=b->m_x)
        return a->m_x<b->m_x;
    if(a->m_y!=b->m_y)
        return a->m_y<b->m_y;
    return a<b;
}


//ʵ��MapNodes���ڲ�����������F�ڵ��ӳ��ڵ㡣���ӳ��ڵ�Ϊ�Լ����ʾ����ɾ�����������������Ҫɾ����
static node** _MapNodes(node** F)
{
    node** L;
    (**F).m_flow = *F;
    for(L=F; *L!=NULL; ++L)
    {
        if( (**L).m_x!=(**F).m_x )
            break;
        if( (**L).m_y!=(**F).m_y )
            break;
        (**L).m_flow = *F;
    }
    return L;
}


//���������TIN��DEM����Flow������m_nodeA�С�
BOOL TIN2Flow::CreateFlow(CString fin, BSTR InDEM)
{
    int i,j,k,m=1;

	std::cout<<"��ʼ��ȡTIN!"<<std::endl;
	if( !ReadTin(fin,m_tinA,m_nodeA) )  //��TIN�ļ�
		return false;
	std::cout<<"��TIN���!"<<std::endl;

    MapNodes(m_tinA,m_nodeA);   //�ж��ظ��㣬�����ظ���ĵ�ַȡͬһ��ַ��
	std::cout<<"�ж��ظ������!"<<std::endl;

	std::cout << "��ʼ����Edge!" << std::endl;
    CreateEdges(m_tinA,m_edgeA); //����Edge
	std::cout << "����Edge���!" << std::endl;

    SetDegrees(m_nodeA);  //����m_degrees
	std::cout << "����m_degrees���!" << std::endl;

    j = m_nodeA.size();   //j��¼ԭʼTIN�ж������

	//��DEM���������ĵ���뵽m_nodeA
	if (!AddStartNode(m_tinA,InDEM, m_nodeA))
		return false;

    //��TIN�ڵĵ㿪ʼ��m_flow������m_nodeA
    for(i=j; i<(int)m_nodeA.size(); ++i)
    {
        node *nd; 
		node *nt;
        nd = m_nodeA[i];
        if( nd->m_isbound )
            continue;
        if( nd->m_flow!=NULL )
            continue;
        nt = nd->find_flow(); //�ҵ���ǰ���m_flow,����m_nodeA������ɾ����

        if( nt!=NULL )    //nt!=NULL���������m_nodeA
        {
            m_nodeA.push_back(nt);
            continue;
        }
        //nt==NULL���������m_nodeA
        nt = nd->m_flow;
        if(nt!=NULL && nt->m_flow==NULL)
            m_nodeA.push_back(nt);
    }
	std::cout<<"CreateFlow()finished!"<<std::endl;
	return true;
}

//��������shp�ļ��ж�ȡ������ ÿ�ĸ����ʾһ��������
BOOL TIN2Flow::ReadTin(CString fin,std::vector<tin>& tinA,std::vector<node*>& nodeA)
{
	//����TIN��ShapeFile��ʽ�ļ��е�ÿһ�������Σ��ߣ�
	SHPHandle shpIn=SHPOpen(fin,"rb");
	if (shpIn == NULL)
	{
		std::cout<<"Can not open the shpIn file!"<<std::endl;
		return false;
	}

    DBFHandle dbfIn = DBFOpen(fin,"rb");

	SHPInfo shpHeader;
	SHPGetInfo(shpIn,&shpHeader.nRecords,&shpHeader.nShapeType,shpHeader.adBoundsMin,shpHeader.adBoundsMax);

    //��ö����TIN�ļ���BBOX
	double minX,minY,maxX,maxY;
	minX = shpHeader.adBoundsMin[0]; minY = shpHeader.adBoundsMin[1];
	maxX = shpHeader.adBoundsMax[0]; maxY = shpHeader.adBoundsMax[1];

    tinA.reserve(shpHeader.nRecords);

	//����ÿһ�������Σ�������m_trigs�����С�
	for(int i  =0 ; i <shpHeader.nRecords; i++)
	{
		SHPObject *pshpObject = SHPReadObject(shpIn, i );
		if(pshpObject == NULL ) continue;

		if(pshpObject->nParts != 1)
        {
            std::cout<<"The line is not a simple line!"<<std::endl;
            SHPDestroyObject(pshpObject);
            continue;
        }

        tin a; node* n; double x,y,z;
		for(int k = 3; --k>=0; )
		{
			x = pshpObject->padfX[k];
			y = pshpObject->padfY[k];
			z = pshpObject->padfZ[k];
            n = new node(x,y,z);
            if( x<=minX || x>=maxX || y<=minY || y>=maxY )
                n->m_isbound = true;
            nodeA.push_back(n);     //n����m_nodeA��
            a.m_nodes[k] = n;       //����m_tinA��
		}
        a.m_slope = DBFReadDoubleAttribute(dbfIn,i,1);
		a.m_aspect = DBFReadDoubleAttribute(dbfIn,i,2);
        tinA.push_back(a);

        SHPDestroyObject(pshpObject);
	}

	SHPClose(shpIn);
    DBFClose(dbfIn);
	return true;
}


//��������ͬ�Ľڵ��Ϊͬһ��node����ָ��
void TIN2Flow::MapNodes(std::vector<tin>& tinA,std::vector<node*>& nodeA)
{
    int i; node** it;

    //����m_flow������m_flowֻ��Ϊ�˼�¼���nodeA��Ҫӳ��Ϊ˭
    std::sort(nodeA.begin(),nodeA.end(),cmp_byxy);
    nodeA.push_back(NULL);
    nodeA.pop_back();
	std::vector<node*>::iterator p=nodeA.begin();
	for(it=p._Ptr; it!=nodeA.end()._Ptr;)
    {
        it = _MapNodes(it);
    }
    //����m_flow���

    //��tinA�ĸ����ڵ��滻
    for(i=tinA.size(); --i>=0;)
    {
        it = tinA[i].m_nodes;
        it[0] = it[0]->m_flow;
        it[1] = it[1]->m_flow;
        it[2] = it[2]->m_flow;
    }

    //ɾ���Ѿ�ȷ�����ظ��ĵ㣬��m_flow�������Լ��ĵ�
    for(i=nodeA.size(); --i>=0;)
    {
        node* n = nodeA[i];
        if( n->m_flow!=n )
        {
            delete n;
            nodeA[i] = NULL;
        }
        else
        {
            n->m_flow = NULL;
        }
    }


		nodeA.erase(std::remove(nodeA.begin(),nodeA.end(),
			(node*)NULL), nodeA.end());
}


//�����ߣ�����m_edgeA
void TIN2Flow::CreateEdges(std::vector<tin>& tinA,std::vector<edge>& edgeA)
{
    edgeA.resize(tinA.size()*3);
    edge* E = &edgeA[0];
    for(int i=tinA.size(); --i>=0;)
    {
        tin& a = tinA[i];
        E = E->create(&a,a.m_nodes[0],a.m_nodes[1]);
        E = E->create(&a,a.m_nodes[1],a.m_nodes[2]);
        E = E->create(&a,a.m_nodes[2],a.m_nodes[0]);
    }

	std::vector<edge>::iterator p;
	for(p=edgeA.begin();p!=edgeA.end();++p)
	{
		if(p._Ptr==E){
			edgeA.erase(p, edgeA.end());
			break;
		}
	}
}


//����ƽ̹���m_degrees
void TIN2Flow::SetDegrees(std::vector<node*>& nodeA)
{
    std::vector<node*> A,N,B;
    int i,j,k;
    //�ҳ����е�ƽ̹��(�õ���Χ�ĵ�Ҫô��Zֵ������Ҫô������ͬ����������һ��Zֵ������ͬ)
    for(i=nodeA.size(); --i>=0;)
    {
        node* a = nodeA[i];

        a->neighbors(N);
        if( a->m_isbound )
            continue;
        for(j=N.size(); --j>=0;)   //�ж�Zֵ����С�ĵ�
        {
            if(N[j]->m_z<a->m_z)
                break;
        }
        if(j!=-1) continue;
        for(j=N.size(); --j>=0;)   //�ж�Zֵ��ͬ�ĵ�
        {
            if(N[j]->m_z==a->m_z)
                break;
        }
        if(j==-1) continue;
        A.push_back(a);            //A�����¼����ƽ̹��
    }

    for(i=A.size(); --i>=0;)       //��ʼֵ��m_degree=-1 ��ʾΪƽ̹�㣬m_degree==0��ʾ��ƽ̹��
    {
        A[i]->m_degree = -1;
    }

    for(k=1; !A.empty(); ++k)     //����A
    {   //kΪ��ƽ̹�����ȡ�
        N.clear(); B.clear();  
		//N��ʾδ�ҵ���ƽ̹��Ľڵ㣻B��ʾ���ҵ���ƽ̹��Ľڵ�
        //�ҵ������ڵ㣬�ҵ�����ͷŵ�B�δ�ҵ�����ͷŵ�N��
        for(i=A.size(); --i>=0; )
        {
            node* a = A[i];
            if( a->find_degree() )
                B.push_back(a);
            else
                N.push_back(a);
        }

        if( B.empty() )  //û���ҵ��µ�����㣬�˳�ѭ��
            break;

        for(i=B.size(); --i>=0;)  //�����ҵ������������
            B[i]->m_degree = k;

        A.swap(N);  //��������N��ֱ��NΪ��
    }
}


//��DEM���������ĵ���뵽m_nodeA������ÿһ�����ĵ��ҵ���Ӧ��TIN
BOOL TIN2Flow::AddStartNode(std::vector<tin>& tinA, BSTR InDEM, std::vector<node*>& nodeA)
{
	CComBase comBase;
	int imgWidth=0, imgHeight=0;
	double dx=0, dy=0, Xmin=0, Ymax=0;
    if (!comBase.OpenImg(InDEM, imgWidth, imgHeight, dx, dy, Xmin, Ymax))
        return false;

	for(int i=0; i<tinA.size(); ++i) 
    {
		double xt[3]= {tinA[i].m_nodes[0]->m_x,tinA[i].m_nodes[1]->m_x,tinA[i].m_nodes[2]->m_x};
		double yt[3]= {tinA[i].m_nodes[0]->m_y,tinA[i].m_nodes[1]->m_y,tinA[i].m_nodes[2]->m_y};
		double zt[3]= {tinA[i].m_nodes[0]->m_z,tinA[i].m_nodes[1]->m_z,tinA[i].m_nodes[2]->m_z};
		
		int col[3],row[3];   //��¼�����������ǵ��Ӧ�ĸ������к�
		int colmin=imgWidth, colmax=0, rowmin=imgHeight, rowmax=0;
		//ȷ�������ε������ǵ��Ӧ�ĸ����������߳�ֵ��
		for(int j=0;j<3;j++)
		{
 			col[j] = (int)((xt[j]-Xmin)/dx+0.5);
			row[j] = (int)((Ymax- yt[j])/dy+0.5);

			if (colmin>col[j]) colmin=col[j];
			if (colmax<col[j]) colmax=col[j];
			if (rowmin>row[j]) rowmin=row[j];
			if (rowmax<row[j]) rowmax=row[j];
		}

		//������ƽ��Z=aX+bY+c����������
		double a =0, b=0, c=0, t=0;
		t = (xt[0]-xt[1])*(yt[0]-yt[2])-(xt[0]-xt[2])*(yt[0]-yt[1]);
		if(t==0) //һ��ƽ�淽�̱��ʽ:Ax+By+Cz+D=0�� A=0��B=0��C=0���
			continue;

		//�����ε����㲻��һ�����ϵ����
		a = ((yt[0]-yt[2])*(zt[0]-zt[1])-(yt[0]-yt[1])*(zt[0]-zt[2]))/t;
		b = -((xt[0]-xt[2])*(zt[0]-zt[1])-(xt[0]-xt[1])*(zt[0]-zt[2]))/t;
		c = zt[0]-a*xt[0]-b*yt[0];
			
		for(int m=rowmin;m<=rowmax;m++)
		for(int n=colmin;n<=colmax;n++)
		{
			double cenX= Xmin+dx*(n+0.25);
			double cenY= Ymax-dy*(m+0.35);
			double cenZ= a*cenX+b*cenY+c ;

			if ((cenX==xt[0] && cenY==yt[0]) || (cenX==xt[1] && cenY==yt[1]) ||(cenX==xt[2] && cenY==yt[2]) ) //�������εĶ���
			{
				if(cenZ>zt[0] || cenZ>zt[1] || cenZ>zt[2])
				{
					nodeA.push_back(tinA[i].InPnt(cenX,cenY,cenZ));
				}
			}
			else if (comBase.IsInEdge(cenX,cenY,xt,yt)) //�������εı���
			{
				if(cenZ>zt[0] || cenZ>zt[1] || cenZ>zt[2])
				{
					nodeA.push_back(tinA[i].InPnt(cenX,cenY,cenZ));
				}
			}
			else if (comBase.IsInTriangle(cenX,cenY,xt,yt))	//�жϵ��Ƿ�����������		
			{
				nodeA.push_back(tinA[i].InPnt(cenX, cenY, cenZ));
			}

		}
    }
	return true;
}

//�����starPnts��ʼ�����е�FlowTrack�ģ���inteltimeʱ������Track������polylines������polylines�����괮��inteltime���ʱ������������ֵ
//strFirstDayDepthһ��¼���������һ���ˮ���ļ���������ʽ���õ���ˮ���뾶water flow depth
//inManingIndDataһ������ʽ���õ�����ϵ������N=0.04
//(outPntXData, outPntYData,outPntZData)��ָ����������ڵ�����ꡣ��������ڵ㽫���뵽������ÿһ����ˮ�ߵ����һ��
BOOL TIN2Flow::CalFlowPnts(double intelTimeData, std::vector<polyline3D>& polylines, CString strFirstDayDepth, double inManingIndData, float *pTWI, float *pL, float *pC, float *pS, double outPntXData, double outPntYData, double outPntZData)
{
	double slope = 0, velocity = 0;

	//��ȡˮ������
	CComBase comBase;
	int imgWidth = 0, imgHeight = 0;
	if (!comBase.OpenImg(strFirstDayDepth,imgWidth,imgHeight))
	{
		return FALSE;
	}

	//��ȡˮ��������Ԫֵ
	float* pBufIn = new float[imgWidth * imgHeight];
	double dx = 0, dy = 0, Xmin = 0, Ymax = 0;
	if (!comBase.OpenImg(strFirstDayDepth,imgWidth,imgHeight,dx,dy,Xmin,Ymax,pBufIn))
	{
		return FALSE;
	}

	for (long i = 0; i < m_nodeA.size(); ++i)
	{
		node* fn = m_nodeA[i];
		if (fn->m_tin == NULL)
			continue;

		long datetime = 0;   //��¼ʵ��ʱ�䣬����Ϊ��λ
		polyline3D pline;  //һ��FlowTrack

		velocity = 0;
		//�����һ����
		pline.push_back(point3D(fn->m_x, fn->m_y, fn->m_z));

		double lefttime = 0;
		for (; fn;) //ÿ�α�����һ���ڵ㣬ֱ��m_flow==NULL��������
		{
 			double pntX0 = fn->m_x;
			double pntY0 = fn->m_y;
			double pntZ0 = fn->m_z;
			fn = fn->m_flow;

			if (fn != NULL)
			{
				double pntX1 = fn->m_x;
				double pntY1 = fn->m_y;
				double pntZ1 = fn->m_z;

				//����ָ���߶ε�ˮ����Ϣ
				int nX0 = (pntX0 - Xmin) / dx;
				int nY0 = (Ymax - pntY0) / dy;
				double dp0 = (double)pBufIn[nY0*imgWidth + nX0];
				if (dp0 < 0) {
					dp0 = 0.232998;
				}
				int nX1 = (pntX1 - Xmin) / dx;
				int nY1 = (Ymax - pntY1) / dy;
				double dp1 = (double)pBufIn[nY1*imgWidth + nX1];
				if (dp1 < 0) {
					dp1 = 0.232998;
				}
				double dep = (dp0 + dp1) / 2;

				AddOneToPline(pline, pntX0, pntY0, pntZ0, pntX1, pntY1, pntZ1, velocity, lefttime, intelTimeData, dep, inManingIndData, Xmin, Ymax, imgWidth, dx, dy, pTWI, pL, pC, pS);
			}
		}


		//��һ����ˮ�ߵ�������������ڵ�,������ˮ�ߵ����һ����������ڵ�С��һ����ֵ��Χ�ż����������ڵ�
		point3D lastPnt = pline[pline.size() - 1];
		if (((lastPnt.x - outPntXData)*(lastPnt.x - outPntXData) + (lastPnt.y - outPntYData)*(lastPnt.y - outPntYData))<3200000)
			pline.push_back(point3D(outPntXData, outPntYData, outPntZData));

		polylines.push_back(pline);
	}

	delete[]pBufIn;  pBufIn = NULL;

	return true;
}

//��polylines������һ���߶�(������)֮������е�
//depth-������ʽ���õ���ˮ���뾶water flow depth
//inManingInd-������ʽ���õ�������ϵ��
void TIN2Flow::AddOneToPline(polyline3D& pline, double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double& lefttime, double intelTime, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS)
{
	// ����ö�·��ʱ�䡣��1��Ϊ��λ��ȷ�������ĸ���
	double t = CalOneLineTime(X0, Y0, Z0, X1, Y1, Z1, velocity, depth, inManingInd, Xmin, Ymax, imgWidth, dx, dy, pTWI, pL, pC, pS);

	////test
	//if (std::isinf(t)) {
	//	return;
	//}

	if (t <= 0)
		return;
	else if (t + lefttime < intelTime)
	{
		lefttime = t + lefttime;
		return;
	}

	if (lefttime > 0.00000)   //������һ���߶β���һ��ʱ�������߶θպ�һ���ʱ�䡣
	{
		double smallt = intelTime - lefttime;

		double X = X0 + (X1 - X0) * smallt / t;
		double Y = Y0 + (Y1 - Y0) * smallt / t;
		double Z = Z0 + (Z1 - Z0) * smallt / t;

		pline.push_back(point3D(X, Y, Z));
		X0 = X; Y0 = Y; Z0 = Z;
		t = t - smallt;
	}

	if (t < intelTime)
	{
		lefttime = t;
		return;
	}

	//�����м�Ľڵ�
	for (int i = intelTime; i<(int)t; i += intelTime)
	{
		double X = X0 + (X1 - X0) * i / t;
		double Y = Y0 + (Y1 - Y0) * i / t;
		double Z = Z0 + (Z1 - Z0) * i / t;
		pline.push_back(point3D(X, Y, Z));
		//Corrected by W.Q
		lefttime = t - (int)i;
	}

	////�����м�Ľڵ�
	////Modify Carol 20200304
	//for (int i = 1; (i * intelTime) < t; i++) {
	//	double X = X0 + (X1 - X0) * (i * intelTime / t);
	//	double Y = Y0 + (Y1 - Y0) * (i * intelTime / t);
	//	double Z = Z0 + (Z1 - Z0) * (i * intelTime / t);

	//	pline.push_back(point3D(X, Y, Z));
	//	lefttime = t - i * intelTime;
	//}
}

//����һ���߶�(������)֮���ͨ��ʱ��
//depth-������ʽ���õ���ˮ���뾶water flow depth
//inManingInd-������ʽ���õ�������ϵ��
double TIN2Flow::CalOneLineTime(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS)
{
	if (Z0 < 0) {
		Z0 = 151.590988;
	}
	if (Z1 < 0) {
		Z1 = 151.590988;
	}
	double dis = sqrt((X1 - X0)*(X1 - X0) + (Y1 - Y0)*(Y1 - Y0));

	//������ʽ
 	double M = 1 / inManingInd;					//inManingInd =0.04
	double R = depth;							//water flow depth
	double Mslope = fabs(Z1 - Z0) / dis;		//������ʽ���õ���slope
	Mslope = fabs(Mslope) > 1e-5 ? Mslope : 0;

	////����ָ���߶ε�slope
	//int nX0s = (X0 - Xmin) / dx;
	//int nY0s= (Ymax - Y0) / dy;
	//float Ms0 = *(pS + nY0s*imgWidth + nX0s);
	//int nX1s = (X1 - Xmin) / dx;
	//int nY1s = (Ymax - Y1) / dy;
	//float Ms1 = *(pS + nY1s*imgWidth + nX1s);
	//float Mslope = (Ms0 + Ms1) / 2;


	//����ָ���߶ε�TWI
	int nX0 = (X0 - Xmin) / dx;
	int nY0 = (Ymax - Y0) / dy;
	float Mtwi0 = *(pTWI + nY0*imgWidth + nX0);
	int nX1 = (X1 - Xmin) / dx;
	int nY1 = (Ymax -Y1) / dy;
	float Mtwi1 = *(pTWI + nY1*imgWidth + nX1);
	float Mtwi = (Mtwi0 + Mtwi1) / 2;

	//����ָ���߶ε�L
	float Ml0 = *(pL + nY0*imgWidth + nX0);
	float Ml1 = *(pL + nY1*imgWidth + nX1);
	float Ml = (Ml0 + Ml1) / 2;

	//����ָ���߶ε�C
	float Mc0 = *(pC + nY0*imgWidth + nX0);
	float Mc1 = *(pC + nY1*imgWidth + nX1);
	float Mc = (Mc0 + Mc1) / 2;

	////TWI,L,C��slope��Ȩ��ϵ��(ͬ����Ҫ)
	//double a = 1;
	//double b = 1;
	//double c = 1;
	//double d = 1;

	//TWI,L,C��slope��Ȩ��ϵ��(�������ֵΪ4.0776��CRΪ0.0287)
	double a = 0.0563;
	double b = 0.1310;
	double c = 0.2388;
	double d = 0.5738;

	////TWI,L,C��slope��Ȩ��ϵ��(�������ֵΪ4.2365��CRΪ0.0876)
	//double a = 0.0464;
	//double b = 0.1334;
	//double c = 0.2667;
	//double d = 0.5536;

	////TWI,L,C��slope��Ȩ��ϵ��(�������ֵΪ4.2365��CRΪ0.0876)
	//double a = 0.0464;
	//double b = 0.0324;
	//double c = 0.1222;
	//double d = 0.7990;

	//��������
	if (Mslope > 0) {			//���¶Ȳ���ƽ��ʱ�����¼����ٶ�
		if (a * Mtwi + b * Ml - c * fabs(Mc) + d * Mslope < 0) {
			velocity = M * pow(R, 2.0 / 3.0) * sqrt(Mslope);
		}
		else {
			velocity = M * pow(R, 2.0 / 3.0) * sqrt(a * Mtwi + b * Ml - c * fabs(Mc) + d * Mslope);
		}
	} else {							//���¶�Ϊƽ��ʱ���ٶȱ���֮ǰ���ٶ�
		velocity = velocity;
	}
	//std::cout << velocity << std::endl;
	double t = dis / velocity;		//����ö�·��ʱ�䣬��1��Ϊ��λ��ȷ�������ĸ���
	return t;
}

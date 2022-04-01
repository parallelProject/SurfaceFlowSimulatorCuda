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



//根据节点a和节点b的xy比较大小。
static inline bool cmp_byxy(node* a,node* b)
{
    if(a->m_x!=b->m_x)
        return a->m_x<b->m_x;
    if(a->m_y!=b->m_y)
        return a->m_y<b->m_y;
    return a<b;
}


//实现MapNodes的内部函数。设置F节点的映射节点。如果映射节点为自己则表示不必删除，其他的情况都需要删除。
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


//根据输入的TIN和DEM创建Flow，放入m_nodeA中。
BOOL TIN2Flow::CreateFlow(CString fin, BSTR InDEM)
{
    int i,j,k,m=1;

	std::cout<<"开始读取TIN!"<<std::endl;
	if( !ReadTin(fin,m_tinA,m_nodeA) )  //读TIN文件
		return false;
	std::cout<<"打开TIN完毕!"<<std::endl;

    MapNodes(m_tinA,m_nodeA);   //判断重复点，并将重复点的地址取同一地址。
	std::cout<<"判断重复点完毕!"<<std::endl;

	std::cout << "开始构造Edge!" << std::endl;
    CreateEdges(m_tinA,m_edgeA); //构造Edge
	std::cout << "构造Edge完毕!" << std::endl;

    SetDegrees(m_nodeA);  //设置m_degrees
	std::cout << "设置m_degrees完毕!" << std::endl;

    j = m_nodeA.size();   //j记录原始TIN中顶点个数

	//将DEM的所有中心点加入到m_nodeA
	if (!AddStartNode(m_tinA,InDEM, m_nodeA))
		return false;

    //从TIN内的点开始找m_flow，放入m_nodeA
    for(i=j; i<(int)m_nodeA.size(); ++i)
    {
        node *nd; 
		node *nt;
        nd = m_nodeA[i];
        if( nd->m_isbound )
            continue;
        if( nd->m_flow!=NULL )
            continue;
        nt = nd->find_flow(); //找到当前点的m_flow,放入m_nodeA，便于删除。

        if( nt!=NULL )    //nt!=NULL的情况放入m_nodeA
        {
            m_nodeA.push_back(nt);
            continue;
        }
        //nt==NULL的情况放入m_nodeA
        nt = nd->m_flow;
        if(nt!=NULL && nt->m_flow==NULL)
            m_nodeA.push_back(nt);
    }
	std::cout<<"CreateFlow()finished!"<<std::endl;
	return true;
}

//描述：从shp文件中读取三角形 每四个点表示一个三角形
BOOL TIN2Flow::ReadTin(CString fin,std::vector<tin>& tinA,std::vector<node*>& nodeA)
{
	//遍历TIN的ShapeFile格式文件中的每一个三角形（线）
	SHPHandle shpIn=SHPOpen(fin,"rb");
	if (shpIn == NULL)
	{
		std::cout<<"Can not open the shpIn file!"<<std::endl;
		return false;
	}

    DBFHandle dbfIn = DBFOpen(fin,"rb");

	SHPInfo shpHeader;
	SHPGetInfo(shpIn,&shpHeader.nRecords,&shpHeader.nShapeType,shpHeader.adBoundsMin,shpHeader.adBoundsMax);

    //获得读入的TIN文件的BBOX
	double minX,minY,maxX,maxY;
	minX = shpHeader.adBoundsMin[0]; minY = shpHeader.adBoundsMin[1];
	maxX = shpHeader.adBoundsMax[0]; maxY = shpHeader.adBoundsMax[1];

    tinA.reserve(shpHeader.nRecords);

	//遍历每一个三角形，并存入m_trigs数组中。
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
            nodeA.push_back(n);     //n放入m_nodeA中
            a.m_nodes[k] = n;       //放入m_tinA中
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


//将坐标相同的节点归为同一个node对象指针
void TIN2Flow::MapNodes(std::vector<tin>& tinA,std::vector<node*>& nodeA)
{
    int i; node** it;

    //设置m_flow，这里m_flow只是为了记录这个nodeA点要映射为谁
    std::sort(nodeA.begin(),nodeA.end(),cmp_byxy);
    nodeA.push_back(NULL);
    nodeA.pop_back();
	std::vector<node*>::iterator p=nodeA.begin();
	for(it=p._Ptr; it!=nodeA.end()._Ptr;)
    {
        it = _MapNodes(it);
    }
    //设置m_flow完毕

    //将tinA的各个节点替换
    for(i=tinA.size(); --i>=0;)
    {
        it = tinA[i].m_nodes;
        it[0] = it[0]->m_flow;
        it[1] = it[1]->m_flow;
        it[2] = it[2]->m_flow;
    }

    //删除已经确定有重复的点，即m_flow不等于自己的点
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


//创建边，放入m_edgeA
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


//设置平坦点的m_degrees
void TIN2Flow::SetDegrees(std::vector<node*>& nodeA)
{
    std::vector<node*> A,N,B;
    int i,j,k;
    //找出所有的平坦点(该点周围的点要么比Z值比它大，要么和它相同，但至少有一个Z值和它相同)
    for(i=nodeA.size(); --i>=0;)
    {
        node* a = nodeA[i];

        a->neighbors(N);
        if( a->m_isbound )
            continue;
        for(j=N.size(); --j>=0;)   //判断Z值比它小的点
        {
            if(N[j]->m_z<a->m_z)
                break;
        }
        if(j!=-1) continue;
        for(j=N.size(); --j>=0;)   //判断Z值相同的点
        {
            if(N[j]->m_z==a->m_z)
                break;
        }
        if(j==-1) continue;
        A.push_back(a);            //A里面记录的是平坦点
    }

    for(i=A.size(); --i>=0;)       //初始值：m_degree=-1 表示为平坦点，m_degree==0表示非平坦点
    {
        A[i]->m_degree = -1;
    }

    for(k=1; !A.empty(); ++k)     //遍历A
    {   //k为找平坦点的深度。
        N.clear(); B.clear();  
		//N表示未找到非平坦点的节点；B表示已找到非平坦点的节点
        //找到相连节点，找到流向就放到B里，未找到流向就放到N中
        for(i=A.size(); --i>=0; )
        {
            node* a = A[i];
            if( a->find_degree() )
                B.push_back(a);
            else
                N.push_back(a);
        }

        if( B.empty() )  //没有找到新的流向点，退出循环
            break;

        for(i=B.size(); --i>=0;)  //设置找到的流向点的深度
            B[i]->m_degree = k;

        A.swap(N);  //继续搜索N，直到N为空
    }
}


//将DEM的所有中心点加入到m_nodeA，并给每一个中心点找到对应的TIN
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
		
		int col[3],row[3];   //记录三角形三个角点对应的格网行列号
		int colmin=imgWidth, colmax=0, rowmin=imgHeight, rowmax=0;
		//确定三角形的三个角点对应的格网，并赋高程值。
		for(int j=0;j<3;j++)
		{
 			col[j] = (int)((xt[j]-Xmin)/dx+0.5);
			row[j] = (int)((Ymax- yt[j])/dy+0.5);

			if (colmin>col[j]) colmin=col[j];
			if (colmax<col[j]) colmax=col[j];
			if (rowmin>row[j]) rowmin=row[j];
			if (rowmax<row[j]) rowmax=row[j];
		}

		//三角形平面Z=aX+bY+c的三个参数
		double a =0, b=0, c=0, t=0;
		t = (xt[0]-xt[1])*(yt[0]-yt[2])-(xt[0]-xt[2])*(yt[0]-yt[1]);
		if(t==0) //一般平面方程表达式:Ax+By+Cz+D=0中 A=0或B=0或C=0情况
			continue;

		//三角形的三点不在一条线上的情况
		a = ((yt[0]-yt[2])*(zt[0]-zt[1])-(yt[0]-yt[1])*(zt[0]-zt[2]))/t;
		b = -((xt[0]-xt[2])*(zt[0]-zt[1])-(xt[0]-xt[1])*(zt[0]-zt[2]))/t;
		c = zt[0]-a*xt[0]-b*yt[0];
			
		for(int m=rowmin;m<=rowmax;m++)
		for(int n=colmin;n<=colmax;n++)
		{
			double cenX= Xmin+dx*(n+0.25);
			double cenY= Ymax-dy*(m+0.35);
			double cenZ= a*cenX+b*cenY+c ;

			if ((cenX==xt[0] && cenY==yt[0]) || (cenX==xt[1] && cenY==yt[1]) ||(cenX==xt[2] && cenY==yt[2]) ) //在三角形的顶点
			{
				if(cenZ>zt[0] || cenZ>zt[1] || cenZ>zt[2])
				{
					nodeA.push_back(tinA[i].InPnt(cenX,cenY,cenZ));
				}
			}
			else if (comBase.IsInEdge(cenX,cenY,xt,yt)) //在三角形的边上
			{
				if(cenZ>zt[0] || cenZ>zt[1] || cenZ>zt[2])
				{
					nodeA.push_back(tinA[i].InPnt(cenX,cenY,cenZ));
				}
			}
			else if (comBase.IsInTriangle(cenX,cenY,xt,yt))	//判断点是否在三角形内		
			{
				nodeA.push_back(tinA[i].InPnt(cenX, cenY, cenZ));
			}

		}
    }
	return true;
}

//计算从starPnts开始的所有的FlowTrack的，在inteltime时间间隔的Track流域线polylines。其中polylines的坐标串是inteltime间隔时间计算出的坐标值
//strFirstDayDepth一记录整个区域第一天的水深文件，曼宁公式中用到的水力半径water flow depth
//inManingIndData一曼宁公式中用到曼宁系数，如N=0.04
//(outPntXData, outPntYData,outPntZData)—指定的流域出口点的坐标。该流域出口点将加入到构建的每一条流水线的最后一点
BOOL TIN2Flow::CalFlowPnts(double intelTimeData, std::vector<polyline3D>& polylines, CString strFirstDayDepth, double inManingIndData, float *pTWI, float *pL, float *pC, float *pS, double outPntXData, double outPntYData, double outPntZData)
{
	double slope = 0, velocity = 0;

	//读取水深数据
	CComBase comBase;
	int imgWidth = 0, imgHeight = 0;
	if (!comBase.OpenImg(strFirstDayDepth,imgWidth,imgHeight))
	{
		return FALSE;
	}

	//获取水深数据像元值
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

		long datetime = 0;   //记录实际时间，以秒为单位
		polyline3D pline;  //一条FlowTrack

		velocity = 0;
		//加入第一个点
		pline.push_back(point3D(fn->m_x, fn->m_y, fn->m_z));

		double lefttime = 0;
		for (; fn;) //每次遍历下一个节点，直到m_flow==NULL，结束。
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

				//计算指定线段的水深信息
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


		//在一条流水线的最后加入流域出口点,仅当流水线的最后一点离流域出口点小于一定阈值范围才加入该流域出口点
		point3D lastPnt = pline[pline.size() - 1];
		if (((lastPnt.x - outPntXData)*(lastPnt.x - outPntXData) + (lastPnt.y - outPntYData)*(lastPnt.y - outPntYData))<3200000)
			pline.push_back(point3D(outPntXData, outPntYData, outPntZData));

		polylines.push_back(pline);
	}

	delete[]pBufIn;  pBufIn = NULL;

	return true;
}

//向polylines中增加一条线段(两个点)之间的所有点
//depth-曼宁公式中用到的水力半径water flow depth
//inManingInd-曼宁公式中用到的曼宁系数
void TIN2Flow::AddOneToPline(polyline3D& pline, double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double& lefttime, double intelTime, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS)
{
	// 计算该段路程时间。以1秒为单位，确定插入点的个数
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

	if (lefttime > 0.00000)   //计算上一段线段不足一秒时，到本线段刚好一秒的时间。
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

	//增加中间的节点
	for (int i = intelTime; i<(int)t; i += intelTime)
	{
		double X = X0 + (X1 - X0) * i / t;
		double Y = Y0 + (Y1 - Y0) * i / t;
		double Z = Z0 + (Z1 - Z0) * i / t;
		pline.push_back(point3D(X, Y, Z));
		//Corrected by W.Q
		lefttime = t - (int)i;
	}

	////增加中间的节点
	////Modify Carol 20200304
	//for (int i = 1; (i * intelTime) < t; i++) {
	//	double X = X0 + (X1 - X0) * (i * intelTime / t);
	//	double Y = Y0 + (Y1 - Y0) * (i * intelTime / t);
	//	double Z = Z0 + (Z1 - Z0) * (i * intelTime / t);

	//	pline.push_back(point3D(X, Y, Z));
	//	lefttime = t - i * intelTime;
	//}
}

//计算一条线段(两个点)之间的通过时间
//depth-曼宁公式中用到的水力半径water flow depth
//inManingInd-曼宁公式中用到的曼宁系数
double TIN2Flow::CalOneLineTime(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS)
{
	if (Z0 < 0) {
		Z0 = 151.590988;
	}
	if (Z1 < 0) {
		Z1 = 151.590988;
	}
	double dis = sqrt((X1 - X0)*(X1 - X0) + (Y1 - Y0)*(Y1 - Y0));

	//曼宁公式
 	double M = 1 / inManingInd;					//inManingInd =0.04
	double R = depth;							//water flow depth
	double Mslope = fabs(Z1 - Z0) / dis;		//曼宁公式中用到的slope
	Mslope = fabs(Mslope) > 1e-5 ? Mslope : 0;

	////计算指定线段的slope
	//int nX0s = (X0 - Xmin) / dx;
	//int nY0s= (Ymax - Y0) / dy;
	//float Ms0 = *(pS + nY0s*imgWidth + nX0s);
	//int nX1s = (X1 - Xmin) / dx;
	//int nY1s = (Ymax - Y1) / dy;
	//float Ms1 = *(pS + nY1s*imgWidth + nX1s);
	//float Mslope = (Ms0 + Ms1) / 2;


	//计算指定线段的TWI
	int nX0 = (X0 - Xmin) / dx;
	int nY0 = (Ymax - Y0) / dy;
	float Mtwi0 = *(pTWI + nY0*imgWidth + nX0);
	int nX1 = (X1 - Xmin) / dx;
	int nY1 = (Ymax -Y1) / dy;
	float Mtwi1 = *(pTWI + nY1*imgWidth + nX1);
	float Mtwi = (Mtwi0 + Mtwi1) / 2;

	//计算指定线段的L
	float Ml0 = *(pL + nY0*imgWidth + nX0);
	float Ml1 = *(pL + nY1*imgWidth + nX1);
	float Ml = (Ml0 + Ml1) / 2;

	//计算指定线段的C
	float Mc0 = *(pC + nY0*imgWidth + nX0);
	float Mc1 = *(pC + nY1*imgWidth + nX1);
	float Mc = (Mc0 + Mc1) / 2;

	////TWI,L,C和slope的权重系数(同等重要)
	//double a = 1;
	//double b = 1;
	//double c = 1;
	//double d = 1;

	//TWI,L,C和slope的权重系数(最大特征值为4.0776，CR为0.0287)
	double a = 0.0563;
	double b = 0.1310;
	double c = 0.2388;
	double d = 0.5738;

	////TWI,L,C和slope的权重系数(最大特征值为4.2365，CR为0.0876)
	//double a = 0.0464;
	//double b = 0.1334;
	//double c = 0.2667;
	//double d = 0.5536;

	////TWI,L,C和slope的权重系数(最大特征值为4.2365，CR为0.0876)
	//double a = 0.0464;
	//double b = 0.0324;
	//double c = 0.1222;
	//double d = 0.7990;

	//计算流速
	if (Mslope > 0) {			//当坡度不是平地时，重新计算速度
		if (a * Mtwi + b * Ml - c * fabs(Mc) + d * Mslope < 0) {
			velocity = M * pow(R, 2.0 / 3.0) * sqrt(Mslope);
		}
		else {
			velocity = M * pow(R, 2.0 / 3.0) * sqrt(a * Mtwi + b * Ml - c * fabs(Mc) + d * Mslope);
		}
	} else {							//当坡度为平地时，速度保持之前的速度
		velocity = velocity;
	}
	//std::cout << velocity << std::endl;
	double t = dis / velocity;		//计算该段路程时间，以1秒为单位，确定插入点的个数
	return t;
}

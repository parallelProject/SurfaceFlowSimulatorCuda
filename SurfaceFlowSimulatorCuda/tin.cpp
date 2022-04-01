// tin.cpp: implementation of the tin class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include "tin.h"
#include <float.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//某个点要么流向F，要么流向T。在F和T中寻找相对较好的流向
static inline node* _best(node* F,node* T)
{
    if( F->m_z==T->m_z )
        return F->m_degree<T->m_degree?F:T;
    return F->m_z<T->m_z?F:T;
}

//从source点要么流向target点，要么流向old点，找出其中流向更好的点。
static node* _best2(node* source,node* target,node* old)
{
    if( target->m_z>=source->m_z )
        return old;
    if( old==NULL )
        return target;
    double x,y,z,l,r;
    x = source->m_x;
    y = source->m_y;
    z = source->m_z;

    l = (target->m_x-x)*(target->m_x-x)+(target->m_y-y)*(target->m_y-y);
    l = -(target->m_z-z)/sqrt(l);  //计算梯度

    r = (old->m_x-x)*(old->m_x-x)+(old->m_y-y)*(old->m_y-y);
    r = -(old->m_z-z)/sqrt(r);     //计算梯度

    return l>r?target:old;
}

static inline node* _find(node** A,node* a,node* b)
{
    if(A[0]!=a && A[0]!=b) return A[0];
    if(A[1]!=a && A[1]!=b) return A[1];
    if(A[2]!=a && A[2]!=b) return A[2];
    return NULL;
}

//返回和该节点相连的所有节点A
void node::neighbors(std::vector<node*>& A)
{
    A.clear();
    for(int i=m_edges.size(); --i>=0;)
    {
        edge& e = *m_edges[i];
        if( e.m_from==this )
            A.push_back(e.m_to);
        if( e.m_to==this )
            A.push_back(e.m_from);
    }
}

//在m_edges中找到一个edge对象，from和to必须是n和this对象，找不到返回NULL。
edge* node::find_edge(node* n)
{
    for(int i=m_edges.size(); --i>=0;)
    {
        edge* e = m_edges[i];
        if( e->m_from==n )
            return e;
        if( e->m_to==n )
            return e;
    }
    return NULL;
}

//判断neighbors相连点中是否有非平坦点或已设置degree点。
BOOL node::find_degree()
{
    node* n;
    for(int i=m_edges.size(); --i>=0;)
    {
        edge* e = m_edges[i];
        n = e->m_from;
        if( n->m_z<=m_z && (int)(n->m_degree)!=-1 )
            return true;
        n = e->m_to;
        if( n->m_z<=m_z && (int)(n->m_degree)!=-1 )
            return true;
    }
    return false;
}

//找出m_to。如果m_to是已经存在的节点则返回NULL，
//否则返回m_to。找到的结果m_to放入m_flow中。
node* node::find_flow()
{
    int i; node* n;			

    //三角形内部情况
    if( m_tin!=NULL )
    {
        return m_tin->find_flow(this);
    }

    //在三角形边上的情况
    if( m_edge!=NULL )
    {
		//m_edge是整个区域的边界
        if( m_edge->m_right==NULL )  
            return NULL;

		//m_edge不是整个区域的边界
        tin* t = m_edge->find_flow();
        if( t!=NULL )
        {
            n = t->find_flow(this);
            if( m_flow!=NULL )
                return n;
        }
        m_flow = _best(m_edge->m_from,m_edge->m_to);
        return NULL;
    }

    ///////////在三角形角点的情况//////////////////

    //首先试图在相邻的三角形中流出
    for(i=m_edges.size(); --i>=0;)
    {
        edge* e = m_edges[i];
        n = e->m_left->find_flow(this);
        if( m_flow!=NULL )
            return n;
        if( e->m_right!=NULL )
        {
            n = e->m_right->find_flow(this);
            if( m_flow!=NULL )
                return n;
        }
    }

    //其次试图在相邻的节点中流出
    for(i=m_edges.size(); --i>=0;)
    {
        edge* e = m_edges[i];
        m_flow = _best2(this,e->m_from,m_flow);
        m_flow = _best2(this,e->m_to,m_flow);
    }
    if( m_flow!=NULL )
        return NULL;
    
    //最后根据m_degree寻找其流向
    for(i=m_edges.size(); --i>=0;)
    {
        edge* e = m_edges[i];
        n = e->m_from;
        if( n->m_z<=m_z && n->m_degree<m_degree )
        {
            m_flow = n;
            return NULL;
        }
        n = e->m_to;
        if( n->m_z<=m_z && n->m_degree<m_degree )
        {
            m_flow = n;
            return NULL;
        }
    }

    //寻找m_to失败 则this一定是汇聚点
    return NULL;
}

//根据三角形t，和三角形的节点a、b初始化edge
edge* edge::create(tin* t,node* a,node* b)
{
    std::vector<edge*>& A = a->m_edges;
    for(int i=A.size(); --i>=0;)
    {
        if(A[i]->m_from==b || A[i]->m_to==b)
        {
            A[i]->m_right = t;
            return this;
        }
    }
    m_from = a; m_to = b;
    m_left = t;
    a->m_edges.push_back(this);
    b->m_edges.push_back(this);
    return this+1;
}

//根据x,y和aspect找到和此边交点的参数。返回值=0表示from点，=1表示to点，在(0,1)表示内部点，否则无交点。
double edge::intersect(double x,double y,double a)
{
    const double pi = 3.1415926535;
    double dx,dy,x0,y0,ax,ay,t,tx,ty;
    a = pi*a/180;
    ax = sin(a); ay = cos(a);

    x0 = m_from->m_x; y0 = m_from->m_y;
    dx = m_to->m_x-x0; dy = m_to->m_y-y0;

	if(dy*ax-dx*ay == 0.0)	return -1;	//处理除数为0的情况

    t = ((x0-x)*ay-(y0-y)*ax)/(dy*ax-dx*ay);
    
    if( t<0 || t>1 ) return -1;

    tx = x0+dx*t; ty = y0+dy*t;
    dx = tx-x; dy = ty-y;
    if( dx*dx+dy*dy<FLT_EPSILON )
        return -1;

    if( dx*ax+dy*ay<0 )
        return -1;

    if( t<0.001 ) return 0;
    if( t>0.999 ) return 1;
    return t;
}

//返回为NULL，表示在该边流动，否则流向返回的三角形。
tin* edge::find_flow()
{
    if( m_right==NULL )
        return NULL;

    double x,y,z; node *tmp,*ttmp;
    x = (m_from->m_x+m_to->m_x)/2;
    y = (m_from->m_y+m_to->m_y)/2;
    z = (m_from->m_z+m_to->m_z)/2;
    tmp = new node(x,y,z);

    ttmp = m_left->find_flow(tmp);
    if( tmp->m_flow!=NULL )
    {
        delete tmp;
        delete ttmp;
        return m_left;
    }

    ttmp = m_right->find_flow(tmp);
    if( tmp->m_flow!=NULL )
    {
        delete tmp;
        delete ttmp;
        return m_right;
    }

    delete tmp;
    return NULL;
}

//描述：在三角形内部寻找点n的流向点
//返回：如果成功则n的m_flow将被赋予其流向，否则n->m_flow仍然为NULL
node* tin::find_flow(node* n)
{
    double t; 
	edge* e;

    if( m_slope==0 )      //坡度为0
        return NULL;

    for(int i=0; i<3; ++i)
    {
        e = m_nodes[i]->find_edge(m_nodes[(i+1)%3]);
        t = e->intersect(n->m_x,n->m_y,m_aspect);
        if( t<0 ) continue;
        if( t==0 )
        {
            n->m_flow = e->m_from;
            return NULL;
        }
        if( t==1 )
        {
            n->m_flow = e->m_to;
            return NULL;
        }
        n->m_flow = e->createnode(t);
        return n->m_flow;
    }
    return NULL;
}

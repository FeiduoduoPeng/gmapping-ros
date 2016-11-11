#ifndef SMMAP_H
#define SMMAP_H
#include <gmapping/grid/map.h>
#include <gmapping/grid/harray2d.h>
#include <gmapping/utils/point.h>
#define SIGHT_INC 1

namespace GMapping {

/*表示地图中的一个cell*/
struct PointAccumulator
{
	typedef point<float> FloatPoint;
	/* before 
	PointAccumulator(int i=-1): acc(0,0), n(0), visits(0){assert(i==-1);}
	*/
	/*after begin*/
	PointAccumulator(): acc(0,0), n(0), visits(0){}
	PointAccumulator(int i): acc(0,0), n(0), visits(0){assert(i==-1);}
	/*after end*/
    inline void update(bool value, const Point& p=Point(0,0));
	
	inline Point mean() const {return 1./n*Point(acc.x, acc.y);}
	
	/*返回被占用的概率*/
	inline operator double() const { return visits?(double)n*SIGHT_INC/(double)visits:-1; }
	
	inline void add(const PointAccumulator& p) {acc=acc+p.acc; n+=p.n; visits+=p.visits; }
	
	static const PointAccumulator& Unknown();
	
	static PointAccumulator* unknown_ptr;
	
	/*cell的位置*/
	FloatPoint acc;
	
	/*n表示被hit中的次数  visits表示访问的次数*/
	int n, visits;
	inline double entropy() const;
};

/*
@desc 更新某个Cell的状态
@value 表示是否被击中
@p     表示机击中的点
*/
void PointAccumulator::update(bool value, const Point& p){
	if (value) 
	{
		acc.x+= static_cast<float>(p.x);
		acc.y+= static_cast<float>(p.y); 
		n++; 
		visits+=SIGHT_INC;
	} else
		visits++;
}

/*
#@desc 该点的熵，表示该点被击中的概率对应的熵。
*/
double PointAccumulator::entropy() const{
	if (!visits)
		return -log(.5);
	if (n==visits || n==0)
		return 0;
	double x=(double)n*SIGHT_INC/(double)visits;
	return -( x*log(x)+ (1-x)*log(1-x) );
}


typedef Map<PointAccumulator,HierarchicalArray2D<PointAccumulator> > ScanMatcherMap;

};

#endif 

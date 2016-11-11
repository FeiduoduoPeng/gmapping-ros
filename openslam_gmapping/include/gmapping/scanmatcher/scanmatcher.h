#ifndef SCANMATCHER_H
#define SCANMATCHER_H

#include "icp.h"
#include "smmap.h"
#include <gmapping/utils/macro_params.h>
#include <gmapping/utils/stat.h>
#include <iostream>
#include <gmapping/utils/gvalues.h>
#define LASER_MAXBEAMS 2048

namespace GMapping {

class ScanMatcher{
	public:
		typedef Covariance3 CovarianceMatrix;   //协方差
		
		ScanMatcher();
		~ScanMatcher();
		
		/*
		icp优化求解位姿
		@param	pnew		求解的新的位姿
		@param  map			匹配的地图
		@param  P			老的位置
		@param	readings    激光数据  
		*/
		double icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		/*
		@desc 			根据旧的位姿和激光数据，通过跟地图匹配得到一个新的位姿
		@param pnew 	得到的最优的位姿
		@param map  	地图
		@param p    	老的位置
		@pram readings  激光数据
		*/
		double optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double optimize(OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double   registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		
		/*设置激光参数 激光束、角度值、激光位置*/
		void setLaserParameters
			(unsigned int beams, double* angles, const OrientedPoint& lpose);
		
		void setMatchingParameters
			(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations, double likelihoodSigma=1, unsigned int likelihoodSkip=0 );
		
		void invalidateActiveArea();
		
		void computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);

		inline double icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		inline double score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		inline unsigned int likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double likelihood(double& lmax, OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		
		double likelihood(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, Gaussian3& odometry, const double* readings, double gain=180.);
		
		/*获取激光束的角度数据*/
		inline const double* laserAngles() const { return m_laserAngles; }
		
		/*获取激光束的数量*/
		inline unsigned int laserBeams() const { return m_laserBeams; }
		
		static const double nullLikelihood;
	protected:
		//state of the matcher
		bool m_activeAreaComputed;
		
		/**laser parameters*/
		unsigned int m_laserBeams;														//激光束的数量
		double       m_laserAngles[LASER_MAXBEAMS];										//各个激光束的角度
		//OrientedPoint m_laserPose;
		PARAM_SET_GET(OrientedPoint, laserPose, protected, public, public)				//激光的位置
		PARAM_SET_GET(double, laserMaxRange, protected, public, public)					//激光的最大测距范围
		/**scan_matcher parameters*/
		PARAM_SET_GET(double, usableRange, protected, public, public)					//使用的激光的最大范围
		PARAM_SET_GET(double, gaussianSigma, protected, public, public)
		PARAM_SET_GET(double, likelihoodSigma, protected, public, public)
		PARAM_SET_GET(int,    kernelSize, protected, public, public)
		PARAM_SET_GET(double, optAngularDelta, protected, public, public)				//优化时的角度增量
		PARAM_SET_GET(double, optLinearDelta, protected, public, public)				//优化时的长度增量
		PARAM_SET_GET(unsigned int, optRecursiveIterations, protected, public, public)	//优化时的迭代次数
		PARAM_SET_GET(unsigned int, likelihoodSkip, protected, public, public)
		PARAM_SET_GET(double, llsamplerange, protected, public, public)
		PARAM_SET_GET(double, llsamplestep, protected, public, public)
		PARAM_SET_GET(double, lasamplerange, protected, public, public)
		PARAM_SET_GET(double, lasamplestep, protected, public, public)
		PARAM_SET_GET(bool, generateMap, protected, public, public)
		PARAM_SET_GET(double, enlargeStep, protected, public, public)
		PARAM_SET_GET(double, fullnessThreshold, protected, public, public)				//被认为是占用的阈值
		PARAM_SET_GET(double, angularOdometryReliability, protected, public, public)	//里程计的角度可靠性
		PARAM_SET_GET(double, linearOdometryReliability, protected, public, public)		//里程计的长度可靠性	
		PARAM_SET_GET(double, freeCellRatio, protected, public, public)					//free和occupany的阈值
		PARAM_SET_GET(unsigned int, initialBeamsSkip, protected, public, public)		//去掉初始的几个激光束的数量	

		// allocate this large array only once
		IntPoint* m_linePoints;
};

/*
@desc	icpOptimize中的一步，目前还很不完善，不能实现正常的功能，还需要做修改。
@param	pret		新的位姿
@param  map			地图
@param  p			旧的位姿
@param  readings	激光数据
*/
inline double ScanMatcher::icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	
	/*把激光雷达的坐标转换到世界坐标系下*/
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	
	unsigned int skip=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	std::list<PointPair> pairs;
	
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange||*r==0.0) continue;
		if (skip) continue;
		
		/*被激光击中的点的地图坐标系*/
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		
		
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		
		/*在对应的kernerSize中搜索*/
		bool found=false;
		Point bestMu(0.,0.);
		Point bestCell(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++)
		{
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;

				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold)
				{
					Point mu=phit-cell.mean();
					if (!found)
					{
						bestMu=mu;
						bestCell=cell.mean();
						found=true;
					}else
					{	
						if((mu*mu)<(bestMu*bestMu))
						{
							bestMu=mu;
							bestCell=cell.mean();
						}
					}
				}
		}
		if (found)
		{
			pairs.push_back(std::make_pair(phit, bestCell));
		}
	}
	
	OrientedPoint result(0,0,0);
	std::cerr << "result(" << pairs.size() << ")=" << result.x << " " << result.y << " " << result.theta << std::endl;
	pret.x=p.x+result.x;
	pret.y=p.y+result.y;
	pret.theta=p.theta+result.theta;
	pret.theta=atan2(sin(pret.theta), cos(pret.theta));
	return score(map, p, readings);
}

/*
@desc 		根据地图、机器人位置、激光雷达数据，计算出一个得分：原理为likelihood_field_range_finder_model
@map  		对应的地图
@p    		机器人对应的初始位置
@readings	激光雷达数据
*/
inline double ScanMatcher::score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	double s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	OrientedPoint lp=p;
	/*
	把激光雷达的坐标转换到世界坐标系
	先旋转到机器人坐标系，然后再转换到世界坐标系
	*/
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	
	/*map.getDelta表示地图分辨率*/
	unsigned int skip=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	
	
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (skip||*r>m_usableRange||*r==0.0) continue;
		
		/*被激光雷达击中的点 在地图坐标系中的坐标*/
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		
		/*该束激光的最后一个空闲坐标，即紧贴hitCell的freeCell*/
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
		
		/*phit 和 freeCell的距离*/
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		
		/*在kernelSize中搜索*/
		bool found=false;
		Point bestMu(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++)
		{
			/*根据已开始的关系 搜索phit的时候，也计算出来pfree*/
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			
			/*得到各自对应的Cell*/
			const PointAccumulator& cell=map.cell(pr);
			const PointAccumulator& fcell=map.cell(pf);
			/*
			(double)cell返回的是该cell被占用的概率
			*/
			if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold)
			{
				/*计算出被击中的点与对应的cell的距离*/
				Point mu=phit-cell.mean();
				if (!found)
				{
					bestMu=mu;
					found=true;
				}else
					bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
			}
		}
	    /*socre的计算公式exp(-d^2 / sigma))*/
		if (found)
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
	}
	return s;
}

/*
@desc 				根据地图、机器人位置、激光雷达数据，同时计算出一个得分和似然：原理为likelihood_field_range_finder_model
					计算出来的似然即为粒子的权重
@param s			得分
@param l			似然
@param map  		对应的地图
@param p    		机器人对应的初始位置
@param readings		激光雷达数据
*/
inline unsigned int ScanMatcher::likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	using namespace std;
	l=0;
	s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	
	/*
	把激光雷达的坐标转换到世界坐标系
	先旋转到机器人坐标系，然后再转换到世界坐标系
	*/
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	
	
	double noHit=nullLikelihood/(m_likelihoodSigma);
	unsigned int skip=0;
	unsigned int c=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		/*每隔m_likelihoodSkip个激光束 就跳过一个激光束*/
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange) continue;
		if (skip) continue;
		
		/*被激光击中的点*/
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		
		Point pfree=lp;
		pfree.x+=(*r-freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-freeDelta)*sin(lp.theta+*angle);
		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		
		/*在对应的kernerSize中搜索*/
		bool found=false;
		Point bestMu(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				
				/*如果cell(pr)被占用 而cell(pf)没有被占用 则说明找到了一个合法的点*/
				if (((double)cell )>m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold)
				{
					Point mu=phit-cell.mean();
					if (!found)
					{
						bestMu=mu;
						found=true;
					}
					else
					{
						bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
					}
				}
		}
	    /*计算得分*/
		if (found)
		{
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
			c++;
		}
		
		/*计算似然*/
		if (!skip)
		{
			double f=(-1./m_likelihoodSigma)*(bestMu*bestMu);
			l+=(found)?f:noHit;
		}
	}
	return c;
}

};

#endif

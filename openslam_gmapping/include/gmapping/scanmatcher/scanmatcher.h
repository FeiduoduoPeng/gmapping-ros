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
		typedef Covariance3 CovarianceMatrix;   //Э����
		
		ScanMatcher();
		~ScanMatcher();
		
		/*
		icp�Ż����λ��
		@param	pnew		�����µ�λ��
		@param  map			ƥ��ĵ�ͼ
		@param  P			�ϵ�λ��
		@param	readings    ��������  
		*/
		double icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		/*
		@desc 			���ݾɵ�λ�˺ͼ������ݣ�ͨ������ͼƥ��õ�һ���µ�λ��
		@param pnew 	�õ������ŵ�λ��
		@param map  	��ͼ
		@param p    	�ϵ�λ��
		@pram readings  ��������
		*/
		double optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double optimize(OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double   registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		
		/*���ü������ ���������Ƕ�ֵ������λ��*/
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
		
		/*��ȡ�������ĽǶ�����*/
		inline const double* laserAngles() const { return m_laserAngles; }
		
		/*��ȡ������������*/
		inline unsigned int laserBeams() const { return m_laserBeams; }
		
		static const double nullLikelihood;
	protected:
		//state of the matcher
		bool m_activeAreaComputed;
		
		/**laser parameters*/
		unsigned int m_laserBeams;														//������������
		double       m_laserAngles[LASER_MAXBEAMS];										//�����������ĽǶ�
		//OrientedPoint m_laserPose;
		PARAM_SET_GET(OrientedPoint, laserPose, protected, public, public)				//�����λ��
		PARAM_SET_GET(double, laserMaxRange, protected, public, public)					//���������෶Χ
		/**scan_matcher parameters*/
		PARAM_SET_GET(double, usableRange, protected, public, public)					//ʹ�õļ�������Χ
		PARAM_SET_GET(double, gaussianSigma, protected, public, public)
		PARAM_SET_GET(double, likelihoodSigma, protected, public, public)
		PARAM_SET_GET(int,    kernelSize, protected, public, public)
		PARAM_SET_GET(double, optAngularDelta, protected, public, public)				//�Ż�ʱ�ĽǶ�����
		PARAM_SET_GET(double, optLinearDelta, protected, public, public)				//�Ż�ʱ�ĳ�������
		PARAM_SET_GET(unsigned int, optRecursiveIterations, protected, public, public)	//�Ż�ʱ�ĵ�������
		PARAM_SET_GET(unsigned int, likelihoodSkip, protected, public, public)
		PARAM_SET_GET(double, llsamplerange, protected, public, public)
		PARAM_SET_GET(double, llsamplestep, protected, public, public)
		PARAM_SET_GET(double, lasamplerange, protected, public, public)
		PARAM_SET_GET(double, lasamplestep, protected, public, public)
		PARAM_SET_GET(bool, generateMap, protected, public, public)
		PARAM_SET_GET(double, enlargeStep, protected, public, public)
		PARAM_SET_GET(double, fullnessThreshold, protected, public, public)				//����Ϊ��ռ�õ���ֵ
		PARAM_SET_GET(double, angularOdometryReliability, protected, public, public)	//��̼ƵĽǶȿɿ���
		PARAM_SET_GET(double, linearOdometryReliability, protected, public, public)		//��̼Ƶĳ��ȿɿ���	
		PARAM_SET_GET(double, freeCellRatio, protected, public, public)					//free��occupany����ֵ
		PARAM_SET_GET(unsigned int, initialBeamsSkip, protected, public, public)		//ȥ����ʼ�ļ���������������	

		// allocate this large array only once
		IntPoint* m_linePoints;
};

/*
@desc	icpOptimize�е�һ����Ŀǰ���ܲ����ƣ�����ʵ�������Ĺ��ܣ�����Ҫ���޸ġ�
@param	pret		�µ�λ��
@param  map			��ͼ
@param  p			�ɵ�λ��
@param  readings	��������
*/
inline double ScanMatcher::icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	
	/*�Ѽ����״������ת������������ϵ��*/
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
		
		/*��������еĵ�ĵ�ͼ����ϵ*/
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		
		
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		
		/*�ڶ�Ӧ��kernerSize������*/
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
@desc 		���ݵ�ͼ��������λ�á������״����ݣ������һ���÷֣�ԭ��Ϊlikelihood_field_range_finder_model
@map  		��Ӧ�ĵ�ͼ
@p    		�����˶�Ӧ�ĳ�ʼλ��
@readings	�����״�����
*/
inline double ScanMatcher::score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	double s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	OrientedPoint lp=p;
	/*
	�Ѽ����״������ת������������ϵ
	����ת������������ϵ��Ȼ����ת������������ϵ
	*/
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	
	/*map.getDelta��ʾ��ͼ�ֱ���*/
	unsigned int skip=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	
	
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (skip||*r>m_usableRange||*r==0.0) continue;
		
		/*�������״���еĵ� �ڵ�ͼ����ϵ�е�����*/
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		
		/*������������һ���������꣬������hitCell��freeCell*/
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
		
		/*phit �� freeCell�ľ���*/
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		
		/*��kernelSize������*/
		bool found=false;
		Point bestMu(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++)
		{
			/*�����ѿ�ʼ�Ĺ�ϵ ����phit��ʱ��Ҳ�������pfree*/
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			
			/*�õ����Զ�Ӧ��Cell*/
			const PointAccumulator& cell=map.cell(pr);
			const PointAccumulator& fcell=map.cell(pf);
			/*
			(double)cell���ص��Ǹ�cell��ռ�õĸ���
			*/
			if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold)
			{
				/*����������еĵ����Ӧ��cell�ľ���*/
				Point mu=phit-cell.mean();
				if (!found)
				{
					bestMu=mu;
					found=true;
				}else
					bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
			}
		}
	    /*socre�ļ��㹫ʽexp(-d^2 / sigma))*/
		if (found)
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
	}
	return s;
}

/*
@desc 				���ݵ�ͼ��������λ�á������״����ݣ�ͬʱ�����һ���÷ֺ���Ȼ��ԭ��Ϊlikelihood_field_range_finder_model
					�����������Ȼ��Ϊ���ӵ�Ȩ��
@param s			�÷�
@param l			��Ȼ
@param map  		��Ӧ�ĵ�ͼ
@param p    		�����˶�Ӧ�ĳ�ʼλ��
@param readings		�����״�����
*/
inline unsigned int ScanMatcher::likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	using namespace std;
	l=0;
	s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	
	/*
	�Ѽ����״������ת������������ϵ
	����ת������������ϵ��Ȼ����ת������������ϵ
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
		/*ÿ��m_likelihoodSkip�������� ������һ��������*/
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange) continue;
		if (skip) continue;
		
		/*��������еĵ�*/
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		
		Point pfree=lp;
		pfree.x+=(*r-freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-freeDelta)*sin(lp.theta+*angle);
		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		
		/*�ڶ�Ӧ��kernerSize������*/
		bool found=false;
		Point bestMu(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				
				/*���cell(pr)��ռ�� ��cell(pf)û�б�ռ�� ��˵���ҵ���һ���Ϸ��ĵ�*/
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
	    /*����÷�*/
		if (found)
		{
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
			c++;
		}
		
		/*������Ȼ*/
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

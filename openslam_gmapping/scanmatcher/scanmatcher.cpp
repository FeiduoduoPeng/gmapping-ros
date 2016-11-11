#include <cstring>
#include <limits>
#include <list>
#include <iostream>

#include <gmapping/scanmatcher/scanmatcher.h>
#include "gridlinetraversal.h"
//#define GENERATE_MAPS

namespace GMapping {

using namespace std;

const double ScanMatcher::nullLikelihood=-.5;

ScanMatcher::ScanMatcher(): m_laserPose(0,0,0){
	//m_laserAngles=0;
	m_laserBeams=0;
	m_optRecursiveIterations=3;
	m_activeAreaComputed=false;

	// This  are the dafault settings for a grid map of 5 cm
	m_llsamplerange=0.01;
	m_llsamplestep=0.01;
	m_lasamplerange=0.005;
	m_lasamplestep=0.005;
	m_enlargeStep=10.;
	m_fullnessThreshold=0.1;
	m_angularOdometryReliability=0.;
	m_linearOdometryReliability=0.;
	m_freeCellRatio=sqrt(2.);
	m_initialBeamsSkip=0;
	
/*	
	// This  are the dafault settings for a grid map of 10 cm
	m_llsamplerange=0.1;
	m_llsamplestep=0.1;
	m_lasamplerange=0.02;
	m_lasamplestep=0.01;
*/	
	// This  are the dafault settings for a grid map of 20/25 cm
/*
	m_llsamplerange=0.2;
	m_llsamplestep=0.1;
	m_lasamplerange=0.02;
	m_lasamplestep=0.01;
	m_generateMap=false;
*/

   m_linePoints = new IntPoint[20000];
}

ScanMatcher::~ScanMatcher(){
	delete [] m_linePoints;
}

void ScanMatcher::invalidateActiveArea(){
	m_activeAreaComputed=false;
}

/*
@desc ������Ч����ֻ���ڴ����һ֡���ݵ�ʱ�򣬲Ż�����������
@param map			��ͼ
@param p			������λ��
@param readings		�����״�����
*/
void ScanMatcher::computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	
	/*�Ѿ�������ˣ���û��Ҫ������*/
	if (m_activeAreaComputed)
		return;
	
	/*�Ѽ����״��λ��ת������ͼ����ϵ*/
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	
	/*��ͼ�ķ�Χ*/
	Point min(map.map2world(0,0));
	Point max(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
	       
	/*��չ��ͼ�ķ�Χ*/
	if (lp.x<min.x) min.x=lp.x;
	if (lp.y<min.y) min.y=lp.y;
	if (lp.x>max.x) max.x=lp.x;
	if (lp.y>max.y) max.y=lp.y;
	
	/*determine the size of the area*/
	/*���ݼ���������չ��ͼ�ķ�Χ*/
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		/*ȥ���������ֵ*/
		if (*r>m_laserMaxRange||*r==0.0||isnan(*r)) continue;
		double d=*r>m_usableRange?m_usableRange:*r;
		
		/*��������е�λ��*/
		Point phit=lp;
		phit.x+=d*cos(lp.theta+*angle);
		phit.y+=d*sin(lp.theta+*angle);
		
		/*���䷶Χ*/
		if (phit.x<min.x) min.x=phit.x;
		if (phit.y<min.y) min.y=phit.y;
		if (phit.x>max.x) max.x=phit.x;
		if (phit.y>max.y) max.y=phit.y;
	}
	//min=min-Point(map.getDelta(),map.getDelta());
	//max=max+Point(map.getDelta(),map.getDelta());
	
	/*�����ͼ��Ҫ��չ*/
	if ( !map.isInside(min)	|| !map.isInside(max))
	{
		/*�õ�Ŀǰ��ͼ�Ĵ�С*/
		Point lmin(map.map2world(0,0));
		Point lmax(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
		//cerr << "CURRENT MAP " << lmin.x << " " << lmin.y << " " << lmax.x << " " << lmax.y << endl;
		//cerr << "BOUNDARY OVERRIDE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
		
		/*�����Ҫ���䣬��Ѷ�Ӧ��ά����չm_enlargeStep�Ĵ�С*/
		min.x=( min.x >= lmin.x )? lmin.x: min.x-m_enlargeStep;
		max.x=( max.x <= lmax.x )? lmax.x: max.x+m_enlargeStep;
		min.y=( min.y >= lmin.y )? lmin.y: min.y-m_enlargeStep;
		max.y=( max.y <= lmax.y )? lmax.y: max.y+m_enlargeStep;
		map.resize(min.x, min.y, max.x, max.y);
		//cerr << "RESIZE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
	}
	
	/*��ͼ����Ч����(��ͼ����ϵ)*/
	HierarchicalArray2D<PointAccumulator>::PointSet activeArea;
	/*allocate the active area*/
	angle=m_laserAngles+m_initialBeamsSkip;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		/*�����Ҫ���ɵ�ͼ*/
		if (m_generateMap)
		{
			double d=*r;
			if (d>m_laserMaxRange||d==0.0||isnan(d))
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
			
			/*�������������յ�*/
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p0=map.world2map(lp);
			IntPoint p1=map.world2map(phit);
			
			//IntPoint linePoints[20000] ;
			/*bresenham�㷨������·��*/
			GridLineTraversalLine line;
			line.points=m_linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
			
			/*���µ�ͼ*/
			for (int i=0; i<line.num_points-1; i++)
			{
				assert(map.isInside(m_linePoints[i]));
				activeArea.insert(map.storage().patchIndexes(m_linePoints[i]));
				assert(m_linePoints[i].x>=0 && m_linePoints[i].y>=0);
			}
			if (d<m_usableRange)
			{
				IntPoint cp=map.storage().patchIndexes(p1);
				assert(cp.x>=0 && cp.y>=0);
				activeArea.insert(cp);
			}
		} 
		else 
		{
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			IntPoint cp=map.storage().patchIndexes(p1);
			assert(cp.x>=0 && cp.y>=0);
			activeArea.insert(cp);
		}
	}
	
	map.storage().setActiveArea(activeArea, true);
	m_activeAreaComputed=true;
}

/*
@desc		��֪λ�õĸ���դ���ͼ�㷨(ʹ�õ�ģ��ΪCounting Model)
            Counting Model��ĳһ���㱻���ǵĸ���Ϊ��ռ�õĴ����뱻���ʵĴ����ı�ֵ
@map		��Ӧ�ĵ�ͼ
@p			�����˵�λ��
@reading    �����״������
*/
double ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	if (!m_activeAreaComputed)
		computeActiveArea(map, p, readings);
		
	//this operation replicates the cells that will be changed in the registration operation
	map.storage().allocActiveArea();
	
	/*�Ѽ����״��λ��ת������ͼ����ϵ*/
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	
	
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	/*esum��ʾ������ص����ӻ��߼���*/
	double esum=0;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		if (m_generateMap)
		{
			/*ȥ���Ƿ��ļ�����*/
			double d=*r;
			if (d>m_laserMaxRange||d==0.0||isnan(d))
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
			
			/*���ü��������еĵ�ĵ�ͼ����*/
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);
			//IntPoint linePoints[20000] ;
			
			/*bresenham�����㷨������ ����λ�úͱ�������е�λ��֮��Ŀ���λ��*/
			GridLineTraversalLine line;
			line.points=m_linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
			
			/*���¿���λ��*/
			for (int i=0; i<line.num_points-1; i++)
			{
				PointAccumulator& cell=map.cell(line.points[i]);
				/*����ǰ���صĸ���*/
				double e=-cell.entropy();       
				cell.update(false, Point(0,0));
				/*���ϸ��º���أ������º���ؼ�ȥ����ǰ���أ������θ��¶����صı仯*/
				e+=cell.entropy();
				esum+=e;
			}
			
			/*���±����е�λ��*/
			if (d<m_usableRange)
			{
				double e=-map.cell(p1).entropy();
				map.cell(p1).update(true, phit);
				e+=map.cell(p1).entropy();
				esum+=e;
			}
		} 
		else 
		{
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;
			/*�����еĵ�ĵ�ͼ����*/
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			
			/*���¶�Ӧ��cell��ֵ*/
			map.cell(p1).update(true,phit);
		}
	}
	return esum;
}

/*
@desc	���ݵ�ͼ���������ݡ�λ����icp��������һ�����ŵ�λ����(ĿǰicpStep��û����ȫд�� �������������ʱ�ǲ����õ�)
@param	pnew		�µ�����λ��
@param  map			��ͼ
@param	init		��ʼλ��
@param  readings	��������
*/

double ScanMatcher::icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	double currentScore;
	double sc=score(map, init, readings);;
	OrientedPoint start=init;
	pnew=init;
	
	/*һֱ�Ż�������û��������Ϊֹ*/
	int iterations=0;
	do{
		currentScore=sc;
		sc=icpStep(pnew, map, start, readings);
		start=pnew;
		iterations++;
	} while (sc>currentScore);
	
	cerr << "i="<< iterations << endl;
	return currentScore;
}

/*
@desc	���ݵ�ͼ���������ݡ�λ�˵������һ�����ŵ��µ�λ�˳���
@param	pnew		�µ�����λ��
@param  map			��ͼ
@param	init		��ʼλ��
@param  readings	��������
*/
double ScanMatcher::optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	double bestScore=-1;
	/*���㵱ǰλ�õĵ÷�*/
	OrientedPoint currentPose=init;
	double currentScore=score(map, currentPose, readings);
	
	/*����ʱ�Ĳ�������*/
	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
	
	/*��ȷ�����Ĵ���*/
	unsigned int refinement=0;
	
	/*�����ķ���*/
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
	int c_iterations=0;
	do{
		/*�����һ��(currentScore)���������һ��(bestScore)����п�������̫���ˣ�Ҫ������������ ������Ը�LM�е���*/
		if (bestScore>=currentScore)
		{
			refinement++;
			adelta*=.5;
			ldelta*=.5;
		}
		bestScore=currentScore;
		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

		/*��8����������һ��  �õ���8������������õ�һ��λ�˺Ͷ�Ӧ�ĵ÷�*/
		Move move=Front;
		do 
		{
			localPose=currentPose;
			switch(move)
			{
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
			}
			
			double odo_gain=1;
			if (m_angularOdometryReliability>0.)
			{
				double dth=init.theta-localPose.theta; 	dth=atan2(sin(dth), cos(dth)); 	dth*=dth;
				odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.)
			{
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
				double drho=dx*dx+dy*dy;
				odo_gain*=exp(-m_linearOdometryReliability*drho);
			}
			/*����÷�*/
			double localScore=odo_gain*score(map, localPose, readings);
			
			/*����÷ָ��ã������*/
			if (localScore>currentScore)
			{
				currentScore=localScore;
				bestLocalPose=localPose;
			}
			c_iterations++;
		} while(move!=Done);
		
		/* �ѵ�ǰλ������ΪĿǰ���ŵ�λ��  ���8��ֵ�������˵Ļ�����ô���ֵ�������*/
		currentPose=bestLocalPose;
	}while (currentScore>bestScore || refinement<m_optRecursiveIterations);
	
	/*��������λ�ú͵÷�*/
	pnew=currentPose;
	return bestScore;
}

struct ScoredMove{
	OrientedPoint pose;
	double score;
	double likelihood;
};

typedef std::list<ScoredMove> ScoredMoveList;

/*
@desc	���ݵ�ͼ���������ݡ�λ�˵������һ������λ�������������������ĺ�����ͬ���ǣ���������������λ�˲���һ���ض���cell�����Ǹ���λ�˵ļ�Ȩ�͡�
		��������������һ����������λ�˵���Ȼ��Ȼ�������Ȼ���м�Ȩ��(Ȩ�ؼ�Ϊ��Ȼ)��ȷ�������˵�λ��
@param	pnew		�µ�����λ��
@param  map			��ͼ
@param	init		��ʼλ��
@param  readings	��������
*/
double ScanMatcher::optimize(OrientedPoint& _mean, ScanMatcher::CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	ScoredMoveList moveList;
	double bestScore=-1;
	
	/*���㵱ǰλ�õĵ÷�score����Ȼlikehood*/
	OrientedPoint currentPose=init;
	ScoredMove sm={currentPose,0,0};
	unsigned int matched=likelihoodAndScore(sm.score, sm.likelihood, map, currentPose, readings);
	double currentScore=sm.score;
	moveList.push_back(sm);
	
	/*�ǶȺ����Ե���*/
	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
	
	unsigned int refinement=0;
	int count=0;
	
	/*8����������Ż����ҵ�һ�����ŵ�λ��*/
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
	do{
		if (bestScore>=currentScore){
			refinement++;
			adelta*=.5;
			ldelta*=.5;
		}
		bestScore=currentScore;
		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

		Move move=Front;
		do {
			localPose=currentPose;
			switch(move){
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
			}
			double localScore, localLikelihood;
			
			double odo_gain=1;
			if (m_angularOdometryReliability>0.){
				double dth=init.theta-localPose.theta; 	dth=atan2(sin(dth), cos(dth)); 	dth*=dth;
				odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.){
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
				double drho=dx*dx+dy*dy;
				odo_gain*=exp(-m_linearOdometryReliability*drho);
			}
			localScore=odo_gain*score(map, localPose, readings);
			//update the score
			count++;
			
			/*���¼���score �� localLikelihood*/
			matched=likelihoodAndScore(localScore, localLikelihood, map, localPose, readings);
			if (localScore>currentScore){
				currentScore=localScore;
				bestLocalPose=localPose;
			}
			
			sm.score=localScore;
			sm.likelihood=localLikelihood;//+log(odo_gain);
			sm.pose=localPose;
			moveList.push_back(sm);
			//update the move list
		} while(move!=Done);
		currentPose=bestLocalPose;
	}while (currentScore>bestScore || refinement<m_optRecursiveIterations);
	
	//normalize the likelihood
	/*��һ����Ȼ*/
	double lmin=1e9;
	double lmax=-1e9;
	/*��ȡ�����С��Ȼ*/
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmin=it->likelihood<lmin?it->likelihood:lmin;
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	
	/*ִ�й�һ������*/
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		it->likelihood=exp(it->likelihood-lmax);
	}
	
	//compute the mean
	/*�����ֵ �� ����*/
	OrientedPoint mean(0,0,0);
	double lacc=0;
	
	/*�����ֵ*/
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		lacc+=it->likelihood;
	}
	mean=mean*(1./lacc);
	
	/*���㷽��*/
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lacc, cov.xy/=lacc, cov.xt/=lacc, cov.yy/=lacc, cov.yt/=lacc, cov.tt/=lacc;
	
	_mean=currentPose;
	_cov=cov;
	return bestScore;
}

/*
@desc			���ü���Ĳ���
@param beams	������������
@param angles   ��������Ӧ�ĽǶ�ֵ
@param lpose    �����λ��
*/
void ScanMatcher::setLaserParameters
	(unsigned int beams, double* angles, const OrientedPoint& lpose){
	assert(beams<LASER_MAXBEAMS);
	m_laserPose=lpose;
	m_laserBeams=beams;
	//m_laserAngles=new double[beams];
	memcpy(m_laserAngles, angles, sizeof(double)*m_laserBeams);	
}
	

/*
����ĳһ�������˵�λ�õ���Ȼ�������Ȼ���ڻ�����λ�õ�һ����Χ��(-mllsamplerange~mllsamplerange  -mllsamplerange~mllsamplerange -mlasamplerange m_lasamplerange)����Ȼ
����ͨ�������Χ�ڵ���Ȼ�ֲ�������������˵�λ�õ�����ֵ�ͷ���
*/
double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	ScoredMoveList moveList;
	
	/*������������ڵ���Ȼ �� �÷�*/
	for (double xx=-m_llsamplerange; xx<=m_llsamplerange; xx+=m_llsamplestep)
	for (double yy=-m_llsamplerange; yy<=m_llsamplerange; yy+=m_llsamplestep)
	for (double tt=-m_lasamplerange; tt<=m_lasamplerange; tt+=m_lasamplestep)
	{
		
		OrientedPoint rp=p;
		rp.x+=xx;
		rp.y+=yy;
		rp.theta+=tt;
		
		ScoredMove sm;
		sm.pose=rp;
		
		likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
		moveList.push_back(sm);
	}
	
	//normalize the likelihood
    /*��һ����Ȼ*/
	double lmax=-1e9;
	double lcum=0;
	
	/*��ȡ������Ȼ*/
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	
	/*����lcum �� ÿ����ĸ���p(x)    ��һ����Ȼ*/
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		//it->likelihood=exp(it->likelihood-lmax);
		lcum+=exp(it->likelihood-lmax);
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}
	//����������������ڵ���Ȼ���Ѿ���һ����
	
	/*�����ֵ  ������λ�õľ�ֵ  E(x) = \sum_{x*p(x)}*/
	OrientedPoint mean(0,0,0);
	double s=0,c=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		s+=it->likelihood*sin(it->pose.theta);
		c+=it->likelihood*cos(it->pose.theta);
	}

	mean=mean*(1./lcum);
	s/=lcum;
	c/=lcum;
	mean.theta=atan2(s,c);
	
	/*���㷽�� ������λ�õķ���*/
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++)
	{
		
		/*����ÿ�����mean�ľ����ͽǶȲ�*/
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		
		/*����Э�������*/
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lcum, cov.xy/=lcum, cov.xt/=lcum, cov.yy/=lcum, cov.yt/=lcum, cov.tt/=lcum;
	
	_mean=mean;
	_cov=cov;
	_lmax=lmax;
	return log(lcum)+lmax;
}


/*������ĺ�����������һ�����ܣ���ͬ�������������������̼Ƶ���Ȼ*/
double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p,
	Gaussian3& odometry, const double* readings, double gain){
	ScoredMoveList moveList;
	
	
	for (double xx=-m_llsamplerange; xx<=m_llsamplerange; xx+=m_llsamplestep)
	for (double yy=-m_llsamplerange; yy<=m_llsamplerange; yy+=m_llsamplestep)
	for (double tt=-m_lasamplerange; tt<=m_lasamplerange; tt+=m_lasamplestep)
	{
		
		OrientedPoint rp=p;
		rp.x+=xx;
		rp.y+=yy;
		rp.theta+=tt;
		
		ScoredMove sm;
		sm.pose=rp;
		
		
		likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
		sm.likelihood+=odometry.eval(rp)/gain;
		assert(!isnan(sm.likelihood));
		moveList.push_back(sm);
	}
	
   //normalize the likelihood	��һ����Ȼ
  double lmax=-std::numeric_limits<double>::max();
	double lcum=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++)
	{
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++)
	{
		//it->likelihood=exp(it->likelihood-lmax);
		lcum+=exp(it->likelihood-lmax);
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}
	//��Ȼ��һ�����
	
	//�����ֵ
	OrientedPoint mean(0,0,0);
	double s=0,c=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		s+=it->likelihood*sin(it->pose.theta);
		c+=it->likelihood*cos(it->pose.theta);
	}
	mean=mean*(1./lcum);
	s/=lcum;
	c/=lcum;
	mean.theta=atan2(s,c);
	
	//���㷽��
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lcum, cov.xy/=lcum, cov.xt/=lcum, cov.yy/=lcum, cov.yt/=lcum, cov.tt/=lcum;
	
	_mean=mean;
	_cov=cov;
	_lmax=lmax;
	double v=log(lcum)+lmax;
	assert(!isnan(v));
	return v;
}

/*����ƥ��Ĳ���*/
void ScanMatcher::setMatchingParameters
	(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations,  double likelihoodSigma, unsigned int likelihoodSkip){	
	m_usableRange=urange;						//��������ʹ�÷�Χ
	m_laserMaxRange=range;						//�����������Χ
	m_kernelSize=kernsize;						//kernsize��Ҫ���ڼ���scoreʱ������Ĵ�С
	m_optLinearDelta=lopt;						//�Ż�ʱ�����Բ���
	m_optAngularDelta=aopt;						//�Ż�ʱ�ĽǶȲ���
	m_optRecursiveIterations=iterations;		//�Ż�ʱ�ĵ�������
	m_gaussianSigma=sigma;						//����socreʱ�ķ���
	m_likelihoodSigma=likelihoodSigma;			//������Ȼʱ�ķ���	
	m_likelihoodSkip=likelihoodSkip;			//������Ȼʱ�������ļ�����
}

};


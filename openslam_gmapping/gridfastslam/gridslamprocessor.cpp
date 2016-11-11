#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
#include <gmapping/gridfastslam/gridslamprocessor.h>

//#define MAP_CONSISTENCY_CHECK
//#define GENERATE_TRAJECTORIES

namespace GMapping {

const double m_distanceThresholdCheck = 20;
 
using namespace std;

  GridSlamProcessor::GridSlamProcessor(): m_infoStream(cout){
    
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;
  }
  
  GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor& gsp) 
    :last_update_time_(0.0), m_particles(gsp.m_particles), m_infoStream(cout){
    
    period_ = 5.0;

    m_obsSigmaGain=gsp.m_obsSigmaGain;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_minimumScore=gsp.m_minimumScore;
    
    m_beams=gsp.m_beams;
    m_indexes=gsp.m_indexes;
    m_motionModel=gsp.m_motionModel;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_matcher=gsp.m_matcher;
    
    m_count=gsp.m_count;
    m_readingCount=gsp.m_readingCount;
    m_lastPartPose=gsp.m_lastPartPose;
    m_pose=gsp.m_pose;
    m_odoPose=gsp.m_odoPose;
    m_linearDistance=gsp.m_linearDistance;
    m_angularDistance=gsp.m_angularDistance;
    m_neff=gsp.m_neff;
	
    cerr << "FILTER COPY CONSTRUCTOR" << endl;
    cerr << "m_odoPose=" << m_odoPose.x << " " <<m_odoPose.y << " " << m_odoPose.theta << endl;
    cerr << "m_lastPartPose=" << m_lastPartPose.x << " " <<m_lastPartPose.y << " " << m_lastPartPose.theta << endl;
    cerr << "m_linearDistance=" << m_linearDistance << endl;
    cerr << "m_angularDistance=" << m_linearDistance << endl;
    
		
    m_xmin=gsp.m_xmin;
    m_ymin=gsp.m_ymin;
    m_xmax=gsp.m_xmax;
    m_ymax=gsp.m_ymax;
    m_delta=gsp.m_delta;
    
    m_regScore=gsp.m_regScore;
    m_critScore=gsp.m_critScore;
    m_maxMove=gsp.m_maxMove;
    
    m_linearThresholdDistance=gsp.m_linearThresholdDistance;
    m_angularThresholdDistance=gsp.m_angularThresholdDistance;
    m_obsSigmaGain=gsp.m_obsSigmaGain;
    
#ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ <<  ": trajectories copy.... ";
#endif
    TNodeVector v=gsp.getTrajectories();
    for (unsigned int i=0; i<v.size(); i++){
		m_particles[i].node=v[i];
    }
#ifdef MAP_CONSISTENCY_CHECK
    cerr <<  "end" << endl;
#endif


    cerr  << "Tree: normalizing, resetting and propagating weights within copy construction/cloneing ..." ;
    updateTreeWeights(false);
    cerr  << ".done!" <<endl;
  }
  
  GridSlamProcessor::GridSlamProcessor(std::ostream& infoS): m_infoStream(infoS){
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;
	
  }

  GridSlamProcessor* GridSlamProcessor::clone() const {
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ << ": performing preclone_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      if (a1.m_reference){
		PointerMap::iterator f=pmap.find(a1.m_reference);
		if (f==pmap.end())
		  pmap.insert(make_pair(a1.m_reference, 1));
		else
		  f->second++;
	      }
	    }
	  }
	}
	cerr << __PRETTY_FUNCTION__ <<  ": Number of allocated chunks" << pmap.size() << endl;
	for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
	  assert(it->first->shares==(unsigned int)it->second);

	cerr << __PRETTY_FUNCTION__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
	GridSlamProcessor* cloned=new GridSlamProcessor(*this);
	
# ifdef MAP_CONSISTENCY_CHECK
	cerr << __PRETTY_FUNCTION__ <<  ": trajectories end" << endl;
	cerr << __PRETTY_FUNCTION__ << ": performing afterclone_fit_test" << endl;
	ParticleVector::const_iterator jt=cloned->m_particles.begin();
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const ScanMatcherMap& m2(jt->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
	  const HierarchicalArray2D<PointAccumulator>& h2(m2.storage());
	  jt++;
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      const autoptr< Array2D<PointAccumulator> >& a2(h2.m_cells[x][y]);
	      assert(a1.m_reference==a2.m_reference);
	      assert((!a1.m_reference) || !(a1.m_reference->shares%2));
	    }
	  }
	}
	cerr << __PRETTY_FUNCTION__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
	return cloned;
}
  
  GridSlamProcessor::~GridSlamProcessor(){
    cerr << __PRETTY_FUNCTION__ << ": Start" << endl;
    cerr << __PRETTY_FUNCTION__ << ": Deleting tree" << endl;
    for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
#ifdef TREE_CONSISTENCY_CHECK		
      TNode* node=it->node;
      while(node)
	node=node->parent;
      cerr << "@" << endl;
#endif
      if (it->node)
	delete it->node;
      //cout << "l=" << it->weight<< endl;
    }
    
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ << ": performing predestruction_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
    for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      const ScanMatcherMap& m1(it->map);
      const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
      for (int x=0; x<h1.getXSize(); x++){
	for (int y=0; y<h1.getYSize(); y++){
	  const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	  if (a1.m_reference){
	    PointerMap::iterator f=pmap.find(a1.m_reference);
	    if (f==pmap.end())
	      pmap.insert(make_pair(a1.m_reference, 1));
	    else
	      f->second++;
	  }
	}
      }
    }
    cerr << __PRETTY_FUNCTION__ << ": Number of allocated chunks" << pmap.size() << endl;
    for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
      assert(it->first->shares>=(unsigned int)it->second);
    cerr << __PRETTY_FUNCTION__ << ": SUCCESS, the error is somewhere else" << endl;
# endif
  }



  /*
  @desc ����ɨ��ƥ��Ĳ��� ���������Ķ��� ���Կ�m_matcher.setMatchingParameters��������
  urange
  range
  sigma
  kersize
  lopt
  aopt
  iterations
  */
  void GridSlamProcessor::setMatchingParameters (double urange, double range, double sigma, int kernsize, double lopt, double aopt, 
						 int iterations, double likelihoodSigma, double likelihoodGain, unsigned int likelihoodSkip){
    m_obsSigmaGain=likelihoodGain;
    m_matcher.setMatchingParameters(urange, range, sigma, kernsize, lopt, aopt, iterations, likelihoodSigma, likelihoodSkip);
  }
  
/*
@desc ������̼��˶�ģ�͵Ĳ���
@srr  �����˶���ɵ��������ķ���
@srt  �����˶���ɵĽǶ����ķ���
@str  ��ת�˶���ɵ��������ķ���
@stt  ��ת�˶���ɵĽǶ����ķ���
*/  
void GridSlamProcessor::setMotionModelParameters
(double srr, double srt, double str, double stt){
  m_motionModel.srr=srr;
  m_motionModel.srt=srt;
  m_motionModel.str=str;
  m_motionModel.stt=stt;	
}
  
  /*
  @desc 				���û����˸����˲�������ֵ
  @linear  				�߹������Ծ���
  @angular 				�߹��ĽǶȾ���
  @resampleThreshold	������ʱ��
  */
  void GridSlamProcessor::setUpdateDistances(double linear, double angular, double resampleThreshold){
    m_linearThresholdDistance=linear; 
    m_angularThresholdDistance=angular;
    m_resampleThreshold=resampleThreshold;	
  }
  
  //HERE STARTS THE BEEF

  GridSlamProcessor::Particle::Particle(const ScanMatcherMap& m):
    map(m), pose(0,0,0), weight(0), weightSum(0), gweight(0), previousIndex(0){
    node=0;
  }
  
  
  void GridSlamProcessor::setSensorMap(const SensorMap& smap)
  {
    
    /*
      Construct the angle table for the sensor
      
      FIXME For now detect the readings of only the front laser, and assume its pose is in the center of the robot 
    */
    
    SensorMap::const_iterator laser_it=smap.find(std::string("FLASER"));
    if (laser_it==smap.end())
	{
      cerr << "Attempting to load the new carmen log format" << endl;
      laser_it=smap.find(std::string("ROBOTLASER1"));
      assert(laser_it!=smap.end());
    }
    const RangeSensor* rangeSensor=dynamic_cast<const RangeSensor*>((laser_it->second));
    assert(rangeSensor && rangeSensor->beams().size());
    
    m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
    double* angles=new double[rangeSensor->beams().size()];
    for (unsigned int i=0; i<m_beams; i++)
	{
      angles[i]=rangeSensor->beams()[i].pose.theta;
    }
    m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
    delete [] angles;
  }
  
  /*
  @desc GridFastSLAM��ʼ��
  */
  void GridSlamProcessor::init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, OrientedPoint initialPose)
  {
  
	//���õ�ͼ��С�ͷֱ���
	m_xmin=xmin;
    m_ymin=ymin;
    m_xmax=xmax;
    m_ymax=ymax;
    m_delta=delta;
    
	/*����ÿ�����ӵĳ�ʼֵ*/
    m_particles.clear();
    TNode* node=new TNode(initialPose, 0, 0, 0);
    ScanMatcherMap lmap(Point(xmin+xmax, ymin+ymax)*.5, xmax-xmin, ymax-ymin, delta);
    for (unsigned int i=0; i<size; i++)
	{
      m_particles.push_back(Particle(lmap));
      m_particles.back().pose=initialPose;
      m_particles.back().previousPose=initialPose;
      m_particles.back().setWeight(0);
      m_particles.back().previousIndex=0;
      
		// this is not needed
		//		m_particles.back().node=new TNode(initialPose, 0, node, 0);

		// we use the root directly
		m_particles.back().node= node;
    }
    m_neff=(double)size;
    m_count=0;
    m_readingCount=0;
    m_linearDistance=m_angularDistance=0;
  }

  /*
  @desc ������ʵλ��  ��Ҫ���ڷ��档����̼Ƶ����Ϊ0ʱʹ�á�
  */
  void GridSlamProcessor::processTruePos(const OdometryReading& o)
  {
    const OdometrySensor* os=dynamic_cast<const OdometrySensor*>(o.getSensor());
    if (os && os->isIdeal() && m_outputStream){
      m_outputStream << setiosflags(ios::fixed) << setprecision(3);
      m_outputStream <<  "SIMULATOR_POS " <<  o.getPose().x << " " << o.getPose().y << " " ;
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << o.getPose().theta << " " <<  o.getTime() << endl;
    }
  }
  
  
  /*
  @desc ����һ֡�������� ������gmapping�㷨����Ҫ������
  ���������������������ĺ�����������̼�Ԥ�⡢����������¡����Ӳ����ȵȲ��衣
  
  ��Ҫ�������£�
  �����˶�ģ�͸�����̼Ʒֲ�
  ���������һ�ι۲������proposal�ֲ���
  ����proposal�ֲ�+�����״�������ȷ���������ӵ�Ȩ��
  �����ӽ����ز���
  
  
  */
  bool GridSlamProcessor::processScan(const RangeReading & reading, int adaptParticles)
  {
    /**retireve the position from the reading, and compute the odometry*/
    /*�õ��µ���̼Ƶ�λ��*/
	OrientedPoint relPose=reading.getPose();
    
	/*m_count��ʾ������������õĴ��� ����ǵ�0�ε���,�����е�λ�˶���һ����*/
	if (!m_count)
	{
      m_lastPartPose=m_odoPose=relPose;
    }
    
    //write the state of the reading and update all the particles using the motion model
	/*����ÿһ�����ӣ�������̼��˶�ģ���в������õ����ӵĳ�������λ��  ��һ����Ӧ��   ��̼Ƶĸ��� */
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
      OrientedPoint& pose(it->pose);
      pose=m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose);
    }

    //invoke the callback
	/*�ص�����  ʵ����ʲô��û��*/
    onOdometryUpdate();
    
    // accumulate the robot translation and rotation
	/*����������̼Ƶ����� ������������˵�����λ�ƺͽǶ�λ�Ƶ��ۻ�ֵ*/
    OrientedPoint move=relPose-m_odoPose;
    move.theta=atan2(sin(move.theta), cos(move.theta));
    m_linearDistance+=sqrt(move*move);
    m_angularDistance+=fabs(move.theta);
    
    // if the robot jumps throw a warning
	/*�����̼Ƶ��˶� ����˶������򱨴�*/
    if (m_linearDistance>m_distanceThresholdCheck){
      cerr << "***********************************************************************" << endl;
      cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
      cerr << "m_distanceThresholdCheck=" << m_distanceThresholdCheck << endl;
      cerr << "Old Odometry Pose= " << m_odoPose.x << " " << m_odoPose.y 
	   << " " <<m_odoPose.theta << endl;
      cerr << "New Odometry Pose (reported from observation)= " << relPose.x << " " << relPose.y 
	   << " " <<relPose.theta << endl;
      cerr << "***********************************************************************" << endl;
      cerr << "** The Odometry has a big jump here. This is probably a bug in the   **" << endl;
      cerr << "** odometry/laser input. We continue now, but the result is probably **" << endl;
      cerr << "** crap or can lead to a core dump since the map doesn't fit.... C&G **" << endl;
      cerr << "***********************************************************************" << endl;
    }
    
	//���� ���µ�λ�ø�ֵ���ɵ�λ��
    m_odoPose=relPose;
    
    bool processed=false;

    // process a scan only if the robot has traveled a given distance or a certain amount of time has elapsed
	/*ֻ�е��������߹�һ���ľ���  ���� ��ת��һ���ĽǶ�  ���߹�һ��ָ����ʱ��Ŵ���������*/
    if (! m_count 
	|| m_linearDistance>=m_linearThresholdDistance 
	|| m_angularDistance>=m_angularThresholdDistance
    || (period_ >= 0.0 && (reading.getTime() - last_update_time_) > period_))
	{
      last_update_time_ = reading.getTime();      

      //this is for converting the reading in a scan-matcher feedable form
      assert(reading.size()==m_beams);
      double * plainReading = new double[m_beams];
      for(unsigned int i=0; i<m_beams; i++)
	  {
		plainReading[i]=reading[i];
      }
	  /*����һ֡����*/
      RangeReading* reading_copy = 
              new RangeReading(reading.size(),
                               &(reading[0]),
                               static_cast<const RangeSensor*>(reading.getSensor()),
                               reading.getTime());

	  /*������ǵ�һ֡����*/
      if (m_count>0)
	  {
		/*Ϊÿ�����ӽ���scanMatch���������ÿ�����ӵ�����λ�ˣ�ͬʱ���������λ�˵ĵ÷ֺ���Ȼ  ��Ӧ��gmapping�����е��������һ�β�������proposal���㷨*/
		scanMatch(plainReading);
	   
	    //���� ����proposal�ĸ�������ˣ��������Ǽ���Ȩ��
	    onScanmatchUpdate();
	
		/*����������ӵ�Ȩ�� ���ں�����״��������ȷ���������ӵ�Ȩ��*/
	    updateTreeWeights(false);
				
		/*�����ز���  ����neff�Ĵ�С�������ز���*/
		resample(plainReading, adaptParticles, reading_copy);
	  }
	  /*����ǵ�һ֡��������*/
	  else 
	  {
		//����ǵ�һ֡���ݣ������ֱ�Ӽ���activeArea����Ϊ���ʱ�򣬶Ի����˵�λ���Ƿǳ�ȷ���ģ�����(0,0,0)
		for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
		{	
			m_matcher.invalidateActiveArea();
			m_matcher.computeActiveArea(it->map, it->pose, plainReading);
			m_matcher.registerScan(it->map, it->pose, plainReading);
	  
			TNode* node=new	TNode(it->pose, 0., it->node,  0);
			node->reading = reading_copy;
			it->node=node;
		}
     }
      //		cerr  << "Tree: normalizing, resetting and propagating weights at the end..." ;
	 //�ٴμ���Ȩ��
     updateTreeWeights(false);
      //		cerr  << ".done!" <<endl;
      
	delete [] plainReading;
	m_lastPartPose=m_odoPose; //update the past pose for the next iteration
	m_linearDistance=0;
	m_angularDistance=0;
	m_count++;
	processed=true;
      
      //keep ready for the next step
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
		it->previousPose=it->pose;
    }
      
    }
    m_readingCount++;
    return processed;
 }
  
  
  std::ofstream& GridSlamProcessor::outputStream(){
    return m_outputStream;
  }
  
  std::ostream& GridSlamProcessor::infoStream(){
    return m_infoStream;
  }
  
  
  
  /*�õ�Ȩ���������ӵ��±�*/
  int GridSlamProcessor::getBestParticleIndex() const{
    unsigned int bi=0;
    double bw=-std::numeric_limits<double>::max();
    for (unsigned int i=0; i<m_particles.size(); i++)
      if (bw<m_particles[i].weightSum)
	  {
		bw=m_particles[i].weightSum;
		bi=i;
      }
    return (int) bi;
  }

  void GridSlamProcessor::onScanmatchUpdate(){}
  void GridSlamProcessor::onResampleUpdate(){}
  void GridSlamProcessor::onOdometryUpdate(){}

  
};// end namespace





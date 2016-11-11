#ifndef MOTIONMODEL_H
#define MOTIONMODEL_H

#include <gmapping/utils/point.h>
#include <gmapping/utils/stat.h>
#include <gmapping/utils/macro_params.h>

namespace  GMapping { 

/*��̼��˶�ģ��*/
struct MotionModel{
	
	/*���㵱ǰ���� ��̼���Ϣ ��ȥ�µ�λ��*/
	OrientedPoint drawFromMotion(const OrientedPoint& p, double linearMove, double angularMove) const;
	
	/*���㵱ǰ���� ��һ�ε���̼���Ϣ  ��һ�ε���̼���Ϣ  ��������µ�λ��*/
	OrientedPoint drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const;
	
	Covariance3 gaussianApproximation(const OrientedPoint& pnew, const OrientedPoint& pold) const;
	double srr, str, srt, stt;
	
	/*
	@srr  �����˶���ɵ��������ķ���
	@srt  �����˶���ɵĽǶ����ķ���
	@str  ��ת�˶���ɵ��������ķ���
	@stt  ��ת�˶���ɵĽǶ����ķ���
	*/
};

};

#endif

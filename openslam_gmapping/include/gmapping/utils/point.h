#ifndef _POINT_H_
#define _POINT_H_
#include <assert.h>
#include <math.h>
#include <iostream>
#include "gvalues.h"

#define DEBUG_STREAM cerr << __PRETTY_FUNCTION__ << ":" //FIXME

namespace GMapping {

/*�ṹ���*/
template <class T>
struct point{
	inline point():x(0),y(0) {}
	inline point(T _x, T _y):x(_x),y(_y){}
	T x, y;
};

/*��� + �����������*/
template <class T>
inline point<T> operator+(const point<T>& p1, const point<T>& p2){
	return point<T>(p1.x+p2.x, p1.y+p2.y);
}

/*��� - �����������*/
template <class T>
inline point<T> operator - (const point<T> & p1, const point<T> & p2){
	return point<T>(p1.x-p2.x, p1.y-p2.y);
}

/*��� * �����������  ����һ����*/
template <class T>
inline point<T> operator * (const point<T>& p, const T& v){
	return point<T>(p.x*v, p.y*v);
}

/*��� * �����������  ����һ����*/
template <class T>
inline point<T> operator * (const T& v, const point<T>& p){
	return point<T>(p.x*v, p.y*v);
}

/*��� * �����������  point���ڻ�*/
template <class T>
inline T operator * (const point<T>& p1, const point<T>& p2){
	return p1.x*p2.x+p1.y*p2.y;
}


/*������ĵ�  ������λ��֮�⣬���нǶȣ��������Ϊ�����˵�λ��*/
template <class T, class A>
struct orientedpoint: public point<T>{
  inline orientedpoint() : point<T>(0,0), theta(0) {};
	inline orientedpoint(const point<T>& p);
	inline orientedpoint(T x, T y, A _theta): point<T>(x,y), theta(_theta){}
        inline void normalize();
	inline orientedpoint<T,A> rotate(A alpha){
		T s=sin(alpha), c=cos(alpha);
		A a=alpha+theta;
		a=atan2(sin(a),cos(a));
		return orientedpoint(
			c*this->x-s*this->y,
			s*this->x+c*this->y, 
			a);
	}
	A theta;
};

/*�Ƕȹ�һ�� ���ѽǶȻ���-PI~PI*/
template <class T, class A>
void orientedpoint<T,A>::normalize() {
  if (theta >= -M_PI && theta < M_PI)
    return;
  
  int multiplier = (int)(theta / (2*M_PI));
  theta = theta - multiplier*2*M_PI;
  if (theta >= M_PI)
    theta -= 2*M_PI;
  if (theta < -M_PI)
    theta += 2*M_PI;
}


/*���캯��*/
template <class T, class A>
orientedpoint<T,A>::orientedpoint(const point<T>& p){
	this->x=p.x;
	this->y=p.y;
	this->theta=0.;
}


/*λ�˵� + ������������*/
template <class T, class A>
orientedpoint<T,A> operator+(const orientedpoint<T,A>& p1, const orientedpoint<T,A>& p2){
	return orientedpoint<T,A>(p1.x+p2.x, p1.y+p2.y, p1.theta+p2.theta);
}

/*λ�˵� - ������������*/
template <class T, class A>
orientedpoint<T,A> operator - (const orientedpoint<T,A> & p1, const orientedpoint<T,A> & p2){
	return orientedpoint<T,A>(p1.x-p2.x, p1.y-p2.y, p1.theta-p2.theta);
}

/*λ�˵� * ������������  ����һ����*/
template <class T, class A>
orientedpoint<T,A> operator * (const orientedpoint<T,A>& p, const T& v){
	return orientedpoint<T,A>(p.x*v, p.y*v, p.theta*v);
}

/*λ�˵� * ������������  ����һ����*/
template <class T, class A>
orientedpoint<T,A> operator * (const T& v, const orientedpoint<T,A>& p){
	return orientedpoint<T,A>(p.x*v, p.y*v, p.theta*v);
}

/*
@desc ����λ�˵Ĳ�ֵ�������P1����P2Ϊԭ�������ϵ��������ꡣ
*/
template <class T, class A>
orientedpoint<T,A> absoluteDifference(const orientedpoint<T,A>& p1,const orientedpoint<T,A>& p2){
	orientedpoint<T,A> delta=p1-p2;
	delta.theta=atan2(sin(delta.theta), cos(delta.theta));
	double s=sin(p2.theta), c=cos(p2.theta);
	return orientedpoint<T,A>(c*delta.x+s*delta.y, 
	                         -s*delta.x+c*delta.y, delta.theta);
}

/*
@desc ����λ�˵ĺ� p2��ʾ�������ú�����ʾ��P1��λ���ϣ�����һ��P2������
*/
template <class T, class A>
orientedpoint<T,A> absoluteSum(const orientedpoint<T,A>& p1,const orientedpoint<T,A>& p2){
	double s=sin(p1.theta), c=cos(p1.theta);
	return orientedpoint<T,A>(c*p2.x-s*p2.y,
				  s*p2.x+c*p2.y, p2.theta) + p1;
}

/*һ��λ�˺�һ��λ�õĺ�*/
template <class T, class A>
point<T> absoluteSum(const orientedpoint<T,A>& p1,const point<T>& p2){
	double s=sin(p1.theta), c=cos(p1.theta);
	return point<T>(c*p2.x-s*p2.y, s*p2.x+c*p2.y) + (point<T>) p1;
}


/*��ıȽϺ��� һ����������*/
template <class T>
struct pointcomparator{
	bool operator ()(const point<T>& a, const point<T>& b) const {
		return a.x<b.x || (a.x==b.x && a.y<b.y);
	}	
};

/*�㹹�ɵĽǶȵıȽϺ���  һ����������*/
template <class T>
struct pointradialcomparator{
	point<T> origin;
	bool operator ()(const point<T>& a, const point<T>& b) const {
		point<T> delta1=a-origin;
		point<T> delta2=b-origin;
		return (atan2(delta1.y,delta1.x)<atan2(delta2.y,delta2.x));
	}	
};

/*�������������X������Y����*/
template <class T>
inline point<T> max(const point<T>& p1, const point<T>& p2){
	point<T> p=p1;
	p.x=p.x>p2.x?p.x:p2.x;
	p.y=p.y>p2.y?p.y:p2.y;
	return p;
}

/*�����������С��X����С��Y����*/
template <class T>
inline point<T> min(const point<T>& p1, const point<T>& p2){
	point<T> p=p1;
	p.x=p.x<p2.x?p.x:p2.x;
	p.y=p.y<p2.y?p.y:p2.y;
	return p;
}

/*
P1P2��ݷ�������Բ�ֵ
T1��ʾP1��T2��ʾP2
T3��ʾP1��P2֮���һ����
*/
template <class T, class F>
inline point<T> interpolate(const point<T>& p1,  const F& t1, const point<T>& p2, const F& t2, const F& t3){
	F gain=(t3-t1)/(t2-t1);
	point<T> p=p1+(p2-p1)*gain;
	return p;
}

/*
����λ��֮��Ĳ�ֵ
T1.T2.T3�Ķ������ϸ�����
*/
template <class T, class A, class F>
inline orientedpoint<T,A> 
interpolate(const orientedpoint<T,A>& p1,  const F& t1, const orientedpoint<T,A>& p2, const F& t2, const F& t3){
	F gain=(t3-t1)/(t2-t1);
	orientedpoint<T,A> p;
	p.x=p1.x+(p2.x-p1.x)*gain;
	p.y=p1.y+(p2.y-p1.y)*gain;
	double  s=sin(p1.theta)+sin(p2.theta)*gain,
		c=cos(p1.theta)+cos(p2.theta)*gain;
	p.theta=atan2(s,c);
	return p;
}


/*�������ŷʽ����*/
template <class T>
inline double euclidianDist(const point<T>& p1, const point<T>& p2){
  return hypot(p1.x-p2.x, p1.y-p2.y);
}

/*����λ�˵�ŷʽ����*/
template <class T, class A>
inline double euclidianDist(const orientedpoint<T,A>& p1, const orientedpoint<T,A>& p2){
  return hypot(p1.x-p2.x, p1.y-p2.y);
}

/*λ�˺͵��ŷʽ����*/
template <class T, class A>
inline double euclidianDist(const orientedpoint<T,A>& p1, const point<T>& p2){
  return hypot(p1.x-p2.x, p1.y-p2.y);
}

/*���λ�˵�ŷʽ����*/
template <class T, class A>
inline double euclidianDist(const point<T>& p1, const orientedpoint<T,A>& p2 ){
  return hypot(p1.x-p2.x, p1.y-p2.y);
}


/*������ͬ���͵��ض���ָ��*/
typedef point<int> IntPoint;
typedef point<double> Point;
typedef orientedpoint<double, double> OrientedPoint;

}; //end namespace

#endif

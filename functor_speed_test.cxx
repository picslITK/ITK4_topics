#include <iostream>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <time.h>
#include "itkVector.h"

//#define PRINT 0
typedef itk::Vector<float,1000> voxeltype;
# ifndef PRINT 
unsigned long  ITERATIONS=10000000;
# else 
unsigned long  ITERATIONS=5;
# endif

template<class T>
class myfunctor
{
public:
  virtual const T& operator()(const T& x) const { return x; }  
  virtual const T& operator()(      T& x) const { return x; }
  myfunctor(){}
  virtual ~myfunctor(){}
};


template<class T>
class myfunctor_id : public myfunctor<T> 
{
public:
  const T& operator()(const T& x) const { 
   return x; 
  } 
};

template<class T>
class myfunctor_idb : public myfunctor<T> 
{
public:
  myfunctor_idb(const T& x) : x(x) {}
  const T& operator()(const T& x) const { 
   return x; 
  } 
private:
  const T& x;
};

template<class T>
class myfunctor_add1  : public myfunctor<T>
{
public:
  const T& operator()( T& x) const { 
   x=x+1;
   return x; 
  } 
private:
  //const T& x;
};

template<class T>
class myfunctor_add1b  : public myfunctor<T>
{
public:
  myfunctor_add1b(const T& x) : x(x) { this->x=x; }
 const T& operator()( T& x) const { 
   x=x+1;
   return x; 
 } 
private:
  const T& x;
};

template<class T>
class myVariableFunctorHolder
{
public:
  void SetFunctor( myfunctor<T>* g )
  {
    f=g;
  }
 void runf( )  { 
   voxeltype b; b.Fill(0);
   for (unsigned int i=0; i<ITERATIONS; i++) {
     (*this->f)(b);
#ifdef PRINT
    std::cout << b << std::endl;
#endif
   }
 } 
 void runf2( )  { 
   voxeltype b; b.Fill(0);
   for (unsigned int i=0; i<ITERATIONS; i++) 
     b=b;
 } 
  
private :
  myfunctor<T>* f;

};

template<class T, class functype>
class myTemplatedFunctorHolder
{
public:
 void runf( )  { 
   voxeltype b; b.Fill(0);
   functype f;
   for (unsigned int i=0; i<ITERATIONS; i++) {
     f(b);
#ifdef PRINT
    std::cout << b << std::endl;
#endif
   }
   }
private :
};

void noop()
{
  //unsigned int ct=0;
  for (unsigned int i=0; i<ITERATIONS; i++) {
    //  ct=ct+1;
  }
}

void idtest()
{
  myfunctor_id<voxeltype> a;
  voxeltype b; b.Fill(0);
  for (unsigned int i=0; i<ITERATIONS; i++) {
     a(b); 
#ifdef PRINT
    std::cout << b << std::endl;
#endif
  }
  
}


void idtestTemplatedFunctor()
{
  voxeltype b; b.Fill(0);
  myfunctor_id<voxeltype> a;
  typedef myfunctor_id<voxeltype> voxidtype;
  typedef myfunctor_add1<voxeltype> voxidtype2;
  myTemplatedFunctorHolder<voxeltype,voxidtype2> myVFH;
  myVFH.runf();
}

void idtestVariableFunctor()
{
  voxeltype b; b.Fill(0);
  myfunctor_id<voxeltype> a;
  myfunctor_add1<voxeltype> a2;
  myVariableFunctorHolder<voxeltype> myVFH;
  myVFH.SetFunctor( &a2 );
  myVFH.runf();
}
void idtestVariableFunctor2()
{
  voxeltype b; b.Fill(0);
  myfunctor_id<voxeltype> a;
  myVariableFunctorHolder<voxeltype> myVFH;
  myVFH.SetFunctor( &a );
  myVFH.runf2();
}

void modtest()
{
  voxeltype b; b.Fill(0);
  for (unsigned int i=0; i<ITERATIONS; i++) {
    myfunctor_add1<voxeltype> a; a(b); 
#ifdef PRINT
    std::cout << b << std::endl;
#endif
  }
}


int main()
{
  /** This function compares the speed of an identity functor 
   *   to an  add-1 functor wrt a add-1 counting loop.  The functors 
   *   are templated on voxel type.  
   */
  // tests organized from fastest to slowest 
  clock_t t1=clock();
  idtest();
  clock_t t2=clock();

  clock_t t3=clock();
  idtestTemplatedFunctor();
  //  idtestVariableFunctor2();
 clock_t t4=clock();

  clock_t t5=clock();
  idtestVariableFunctor();
  clock_t t6=clock();

  double nooptesttime=(t2-t1); 
  double idtesttime=(t4-t3); 
  double modtesttime=(t6-t5); 
  std::cout<< " tem-id-time " << idtesttime << " var-id-time " << modtesttime << std::endl;
  std::cout<< " temidtime/idtime " << idtesttime/nooptesttime << " varidtime/idtime " << modtesttime/nooptesttime << std::endl;
}

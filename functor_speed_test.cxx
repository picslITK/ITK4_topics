#include <iostream>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <time.h>
#include "itkVector.h"

#define ADDER 0
//#define PRINT 0
typedef itk::Vector<float,10> voxeltype;
# ifndef PRINT 
unsigned long  ITERATIONS=1000000000;
# else 
unsigned long  ITERATIONS=5;
# endif

template<class T>
class myfunctor
{
public:
  virtual const T& operator()(      T& x) const { return x; }
  myfunctor(){}
  virtual ~myfunctor(){}
};

template<class T>
class myfunctor_id : public myfunctor<T> 
{
public:
  const T& operator()( T& x) const { 
   return x; 
  } 
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
class myVariableFunctorHolder
{
public:
  void SetFunctor( myfunctor<T>* g )
  {
    this->f=g;
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



void optestTemplatedFunctor()
{
  typedef myfunctor_id<voxeltype> voxoptype_id;
  typedef myfunctor_add1<voxeltype> voxoptype_add1;
#ifdef ADDER
  myTemplatedFunctorHolder<voxeltype,voxoptype_add1> myVFH;
#else 
  myTemplatedFunctorHolder<voxeltype,voxoptype_id> myVFH;
#endif;
  myVFH.runf();
}

void optestVariableFunctor()
{
  myfunctor_id<voxeltype> id_op;
  myfunctor_add1<voxeltype> add1_op;
  myVariableFunctorHolder<voxeltype> myVFH;
#ifdef ADDER
  myVFH.SetFunctor( &add1_op );
#else 
  myVFH.SetFunctor( &id_op );
#endif;
  myVFH.runf();
}

void modtest()
{
  voxeltype b; b.Fill(0);
  myfunctor_id<voxeltype> id_op;
  myfunctor_add1<voxeltype> add1_op;
  for (unsigned int i=0; i<ITERATIONS; i++) {
#ifdef ADDER
    add1_op(b); 
#else 
    id_op(b); 
#endif;
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
  voxeltype vox;
 #ifdef ADDER
  std::cout << " operation-add , vec-length " << vox.Size() << " its " << ITERATIONS << std::endl;
#else 
  std::cout << " operation-id , vec-length " << vox.Size() << " its " << ITERATIONS << std::endl;
#endif;
 // tests organized from fastest to slowest 
  clock_t t1=clock();
  modtest();
  clock_t t2=clock();

  clock_t t3=clock();
  optestTemplatedFunctor();
  clock_t t4=clock();

  clock_t t5=clock();
  optestVariableFunctor();
  clock_t t6=clock();

  double nooptesttime=(t2-t1); 
  double idtesttime=(t4-t3); 
  double modtesttime=(t6-t5); 
  std::cout<< " templated-op-time " << idtesttime << " varifunc-op-time " << modtesttime << std::endl;
  std::cout<< " templatedtime/basetime " << idtesttime/nooptesttime << " varoptime/basetime " << modtesttime/nooptesttime << std::endl;
}

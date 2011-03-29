/*
 * itkCompositiveAffineTest.cxx
 *
 *  Created on: Apr 23, 2011
 *      Author: songgang
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "itkMatrixOffsetTransformBase.h"
#include "itkCompositeTransform.h"
#include "itkCenteredCompositeTransform.h"
#include <itkTimeProbe.h>
#include "vcl_cmath.h"
#include "itkScaleTransform.h"
#include "itkTranslationTransform.h"
#include "itkShear2DTransform.h"
// #include "itkRigid2DTransform.h"
#include "itkRotate2DTransform.h"

// first test: test the MatrixOffsetTransformBase from the old ITK
// second test: test the scaling transform
// third test: test my own affine (quaternion) transform
// third test: test the itkCompositeTransform works


template<typename,typename> struct ty { };

template<class T, class U> struct X {
    void f() {
        f(ty<T, U>());
    }

private:
    template<typename Tx, typename Ux>
    void f(ty<Tx, Ux>) { } // generic

    template<typename Tx>
    void f(ty<Tx, char>) { } // "specialized"
};



template<typename T, int Dim>
class MyTransform
{
public:
     MyTransform(){std::cout << "general!" << std::endl;};
     void rotate();
    T m;
};

template <typename T>
class MyTransform <T, 2>{
public:
    MyTransform() {
        std::cout << "2!" << std::endl;
    }
    void rotate();
    int m;
};

template <typename T, int Dim>
void
MyTransform<T, Dim>
::rotate() {
    std::cout << " rotate genereal" << std::endl;
};

template <typename T>
void
MyTransform<T, 2>
::rotate() {
    std::cout << " rotate 3" << std::endl;
};

int main1(){
    MyTransform<bool, 2> m;
    MyTransform<bool, 3> n;
    MyTransform<bool, 4> k;
    k.rotate();
    m.rotate();
    n.rotate();
    return 0;
}





typedef std::vector<double> RawVectorType;
typedef std::vector< RawVectorType > RawPointListType;


template<int Dim, class TRawVector, class TRawPointList>
void generate_random_point_list(const int n, TRawPointList &x, TRawPointList &y,
        TRawVector &A, TRawVector &t, TRawVector &c);

template<int Dim, class TRawVector, class TRawPointList>
void test_one_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter);


template<int Dim, class TRawVector, class TRawPointList>
void test_composite_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter);


template<int Dim, class TRawVector, class TRawPointList>
void test_centered_composite_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter);

template<int Dim, class TRawVector, class TRawPointList>
void test_centered_composite_RSKT_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter);

template<class TRawPointList>
double max_distance_between_point_list(const TRawPointList &a, const TRawPointList &b){
    if (a.size() != b.size() ) {
        std::cerr << "in distance_between_point_list: ERROR not equal size "
                << a.size() << " != " << b.size() << std::endl;
        return -1.0;

    }

    if (a.size()==0) {
        std::cerr << "in distance_between_point_list: WARNING empty point list ! "
                << std::endl;
        return 0;
    }

    double dmax=0.0;
    int nb_pts = a.size();
    int dim = a[1].size();
    for(int i=0; i<nb_pts; i++){
        double d1sum = 0;
        for(int d=0; d<dim; d++)
            d1sum += (a[i][d]-b[i][d]) * (a[i][d]-b[i][d]);
        dmax = (dmax > d1sum) ? (dmax) : (d1sum);
    }
    dmax = vcl_sqrt(dmax);
    return dmax;
}

int main(int argc, char **argv){

    const int Dim = 2;
    int n = 100;
    int nb_iter = 2000;

    srand(time(NULL));
    //    srand(1);

    std::cout << "Test composite affine" << std::endl;


    //    return EXIT_SUCCESS;

    RawPointListType raw_x, raw_y;
    RawVectorType A, t, c;

    generate_random_point_list<Dim, RawVectorType, RawPointListType>(n, raw_x, raw_y, A, t, c);

    itk::TimeProbe t1;
    t1.Start();
    //    test_one_transform<Dim, RawVectorType, RawPointListType>(raw_x, raw_y, A, t, c, nb_iter);
    t1.Stop();

//    itk::TimeProbe t2;
//    t2.Start();
//    test_composite_transform<Dim, RawVectorType, RawPointListType>(raw_x, raw_y, A, t, c, nb_iter);
//    t2.Stop();

    itk::TimeProbe t3;
    t3.Start();
    test_centered_composite_transform<Dim, RawVectorType, RawPointListType>(raw_x, raw_y, A, t, c, nb_iter);
    t3.Stop();



    itk::TimeProbe t4;
    t4.Start();
    test_centered_composite_RSKT_transform<Dim, RawVectorType, RawPointListType>(raw_x, raw_y, A, t, c, nb_iter);
    t4.Stop();


    std::cout << "test_one_transform timing:" << t1.GetTotal() << " seconds" << std::endl;
//    std::cout << "test_composite_transform timing:" << t2.GetTotal() << " seconds" << std::endl;
    std::cout << "test_centered_composite_transform timing:" << t3.GetTotal() << " seconds" << std::endl;
    std::cout << "test_centered_composite_RSKT_transform timing:" << t4.GetTotal() << " seconds" << std::endl;

    return EXIT_SUCCESS;
}


template<int Dim, class TRawVector, class TRawPointList>
void generate_random_point_list(const int n, TRawPointList &x, TRawPointList &y,
        TRawVector &A1, TRawVector &t1, TRawVector &c1){

    //    int n = 5;
    //    const int Dim = 2;

    // float A[Dim*Dim] = { 1.2, 0.2, 1.1, -0.1 };

    // float c[Dim] = { 0.4, 0.6 };
    //    float c[Dim] = { 0.4, 0.6 };
    //    float t[Dim] = { 0.3, 0.2 };

    x.clear();
    y.clear();
    A1.clear();
    t1.clear();
    c1.clear();



//    float A[Dim*Dim] = { 0.5, -3.6, 0.2, -2.9 };
//    float A[Dim*Dim] = { 1.2, -0.8, 0.8, 0.6};
//    float A[Dim*Dim] = { 0.6, -0.8, 0.8, 0.6};
//    float A[Dim*Dim] = { 0.8, -0.8, 0.8, 1.2};
//    float A[Dim*Dim] = { 1, 0, 0, 1};
    float A[Dim*Dim] = { 0.4330,    0.5294,     0.6500,    2.29581};
//        float A[Dim*Dim] = { 2, -1.1, 0, 1};
    float c[Dim] = { 0.5, 0.6 };
//    float c[Dim]  = { 0, 0 };
    float t[Dim] = { 0.4, -1.2 };
//    float t[Dim] = { 0, 0 };

    A1.resize(Dim*Dim);
    for(int i=0; i<Dim*Dim; i++) A1[i] = A[i];

    c1.resize(Dim*Dim);
    for(int i=0; i<Dim; i++) c1[i] = c[i];

    t1.resize(Dim*Dim);
    for(int i=0; i<Dim; i++) t1[i] = t[i];


    // std::vector< std::vector<float> > x, y;
    x.resize(n);
    y.resize(n);

    for(int i=0; i<n; i++) {
        x[i].resize(Dim);
        y[i].resize(Dim);
    }

    for(int i=0; i<n; i++) {
        for(int d=0; d<Dim; d++) {
            x[i][d] = ( rand() % 1000 ) / 1000.0;
        }

        for(int d=0; d<Dim; d++) {
            y[i][d] = 0.0;
            for(int f=0; f<Dim; f++) {
                y[i][d] += A[d*Dim +f ] * (x[i][f]-c[f]);
            }
            y[i][d] += t[d] + c[d];
            y[i][d] += (rand() % 1000 ) / 1000.0 * 0.5; //noise
        }

    }

    std::cout << "A = ";
    for(int i=0; i<Dim*Dim; i++) std::cout<<A1[i]<<",";
    std::cout << std::endl;

    std::cout << "c = ";
    for(int i=0; i<Dim; i++) std::cout <<c1[i]<<",";
    std::cout << std::endl;

    std::cout << "t = ";
    for(int i=0; i<Dim; i++) std::cout <<t1[i]<<",";
    std::cout << std::endl;

    std::cout <<"x = \t==>\t y = " << std::endl;
    for(int i=0; i<n; i++){
        for(int d=0; d<Dim; d++){
            std::cout << x[i][d] << " ";
        }
        std::cout << "\t==>\t";
        for(int d=0; d<Dim; d++){
            std::cout << y[i][d] << " ";
        }
        std::cout << std::endl;

    }

}


template<int Dim, class TRawVector, class TRawPointList>
void test_one_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter){

    typedef itk::MatrixOffsetTransformBase<double, Dim, Dim> AffineType;

    typename AffineType::Pointer affine = AffineType::New();

    // start from identity transform
    affine->SetIdentity();

    typedef typename AffineType::ParametersType ParaType;

    ParaType c1(Dim);
    for(int i=0; i<Dim; i++) c1[i]=c[i];

    affine->SetFixedParameters(c1);

    std::cout << "initial parameters" << std::endl;
    std::cout << "c1=" << c1 << std::endl;
    std::cout << "Matrix A: " << affine->GetMatrix() << std::endl;
    std::cout << "Offset t: " << affine->GetOffset() << std::endl << std::endl;

    int nb_pts = x.size();

    const int kParaDim = AffineType::ParametersDimension;
    ParaType current_para(kParaDim), delta_para(kParaDim);
    current_para = affine->GetParameters();

    int cnt = 0;
    while(cnt <= nb_iter){

        delta_para.Fill(0);

        for(int i=0; i<nb_pts; i++){
            typedef typename AffineType::InputPointType InputPointType;
            typedef typename AffineType::OutputPointType OutputPointType;

            InputPointType ptx;
            OutputPointType pty, ptz;

            for(int d=0;d<Dim;d++){
                ptx[d]=x[i][d];
                pty[d]=y[i][d];
            }
            ptz = affine->TransformPoint(ptx);
            //            std::cout << "ptx:" << ptx
            //                    << "\t" << "pty:" << pty
            //                    << "\t" << "ptz:" << ptz << std::endl;

            typedef typename AffineType::JacobianType JacobianType;
            JacobianType jac(Dim, kParaDim);

            affine->GetJacobianWithRespectToParameters(ptx, jac);

            JacobianType jac2(Dim, kParaDim);
            jac2 = affine->GetMatrix() * jac;

            for(int d=0; d<Dim; d++){
                for(int k=0; k<kParaDim; k++){
                    delta_para[k] += (ptz[d] - pty[d]) * jac[d][k];
                }
            }
        }

        delta_para *= 1.0 / nb_pts ;
        current_para -= delta_para;
        affine->SetParameters(current_para);
        cnt++;

    }

    std::cout << "current para:" << affine->GetParameters() << std::endl;

    for(int i=0; i<nb_pts; i++){
        typedef typename AffineType::InputPointType InputPointType;
        typedef typename AffineType::OutputPointType OutputPointType;

        InputPointType ptx;
        OutputPointType pty, ptz;

        for(int d=0;d<Dim;d++){
            ptx[d]=x[i][d];
            pty[d]=y[i][d];
        }
        ptz = affine->TransformPoint(ptx);
        std::cout << "ptx:" << ptx
                << "\t" << "pty:" << pty
                << "\t" << "ptz(=?=pty):" << ptz << std::endl;
    }



}



template<int Dim, class TRawVector, class TRawPointList>
void test_composite_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter){

    typedef itk::MatrixOffsetTransformBase<double, Dim, Dim> AffineType;
    typedef itk::CompositeTransform<double, Dim> CompositeType;

    typename AffineType::Pointer affine = AffineType::New();
    affine->SetIdentity();

    typedef typename AffineType::ParametersType ParaType;
    ParaType c1(Dim);
    for(int i=0; i<Dim; i++) c1[i]=c[i];

    affine->SetFixedParameters(c1);

    // set random to affine
    ParaType aff_para(Dim*Dim+Dim);
    aff_para[0] = 0.2;
    aff_para[1] = 0.8;
    aff_para[3] = 0.8;
    aff_para[4] = -1.1;
    aff_para[5] = 0.2;
    aff_para[6] = 0.3;
    //    for(int i=0; i<aff_para.Size(); i++) aff_para[i] = (rand() % 100) / 50.0;
    affine->SetParameters(aff_para);




    typename CompositeType::Pointer comp = CompositeType::New();
    comp->AddTransform(affine);



    typename AffineType::Pointer affine2 = AffineType::New();
    affine2->SetIdentity();
    ParaType c2(Dim);
    for(int i=0; i<Dim; i++) c2[i]=c[i];
    affine2->SetFixedParameters(c2);


    // set random to affine2
    ParaType aff2_para(Dim*Dim+Dim);
    aff2_para[0] = -0.7;
    aff2_para[1] = -0.2;
    aff2_para[3] = -0.3;
    aff2_para[4] = 0.9;
    aff2_para[5] = 0.1;
    aff2_para[6] = 0.2;
    //    for(int i=0; i<aff2_para.Size(); i++) aff2_para[i] = (rand() % 100) / 50.0;
    affine2->SetParameters(aff2_para);


    comp->AddTransform(affine2);

    std::cout << "affine2: " << affine2->GetMatrix() << std::endl;


    std::cout << "initial parameters" << std::endl;
    std::cout << comp->GetParameters() << std::endl;
    std::cout << comp->GetFixedParameters() << std::endl;

    int nb_pts = x.size();

    const int kParaDim = comp->GetNumberOfParameters();
    ParaType current_para(kParaDim), delta_para(kParaDim);
    current_para = comp->GetParameters();

    int cnt = 0;
    while(cnt <= nb_iter){

        delta_para.Fill(0);

        // std::cout << "iter: " << cnt << "---------------" << std::endl;

        for(int i=0; i<nb_pts; i++){
            typedef typename CompositeType::InputPointType InputPointType;
            typedef typename CompositeType::OutputPointType OutputPointType;

            InputPointType ptx;
            OutputPointType pty, ptz;

            for(int d=0;d<Dim;d++){
                ptx[d]=x[i][d];
                pty[d]=y[i][d];
            }
            ptz = comp->TransformPoint(ptx);
            //            std::cout << "ptx:" << ptx
            //                    << "\t" << "pty:" << pty
            //                    << "\t" << "ptz:" << ptz << std::endl;

            typedef typename AffineType::JacobianType JacobianType;
            JacobianType jac(Dim, kParaDim);

            comp->GetJacobianWithRespectToParameters(ptx, jac);

            for(int d=0; d<Dim; d++){
                for(int k=0; k<kParaDim; k++){
                    delta_para[k] += (ptz[d] - pty[d]) * jac[d][k];
                }
            }
        }

        delta_para *= 1.0 / nb_pts * 0.5;
        current_para -= delta_para;
        comp->SetParameters(current_para);
        cnt++;

    }

    std::cout << "current para:" << comp->GetParameters() << std::endl;

    for(int i=0; i<nb_pts; i++){
        typedef typename AffineType::InputPointType InputPointType;
        typedef typename AffineType::OutputPointType OutputPointType;

        InputPointType ptx;
        OutputPointType pty, ptz;

        for(int d=0;d<Dim;d++){
            ptx[d]=x[i][d];
            pty[d]=y[i][d];
        }
        ptz = comp->TransformPoint(ptx);
        std::cout << "ptx:" << ptx
                << "\t" << "pty:" << pty
                << "\t" << "ptz(=?=pty):" << ptz << std::endl;
    }

    return;

}







template<int Dim, class TRawVector, class TRawPointList>
void test_centered_composite_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter){

    typedef itk::MatrixOffsetTransformBase<double, Dim, Dim> AffineType;
    typedef itk::CenteredCompositeTransform<double, Dim> CompositeType;

    typename AffineType::Pointer affine = AffineType::New();
    affine->SetIdentity();

    // set the zero to the affine center
    typedef typename CompositeType::ParametersType ParaType;
    ParaType c1(Dim);
    for(int i=0; i<Dim; i++) c1[i]=0; // c[i];

    affine->SetFixedParameters(c1);

    // set random to affine
    ParaType aff_para(Dim*Dim+Dim);
//    aff_para[0] = 0.2;
//    aff_para[1] = 0.8;
//    aff_para[3] = 0.8;
//    aff_para[4] = -1.1;
//    aff_para[5] = 0.2;
//    aff_para[6] = 0.3;
    aff_para[0] = 1;
    aff_para[1] = 0.;
    aff_para[2] = 0.;
    aff_para[3] = 1;
    aff_para[4] = 0;
    aff_para[5] = 0;

    //    for(int i=0; i<aff_para.Size(); i++) aff_para[i] = (rand() % 100) / 50.0;
    affine->SetParameters(aff_para);

//    float scale[6] = {4,4,4,4,0.5,0.5};
    float scale[6] = {1,1,1,1,0.5,0.5};



    typename CompositeType::Pointer comp = CompositeType::New();
    comp->AddTransform(affine);


    // set fixed parameter to zeros, which will also set center
    typename AffineType::Pointer affine2 = AffineType::New();
    affine2->SetIdentity();
    ParaType c2(Dim);
    for(int i=0; i<Dim; i++) c2[i]=0; //c[i];
    affine2->SetFixedParameters(c2);


    // set random to affine2
    ParaType aff2_para(Dim*Dim+Dim);
//    aff2_para[0] = -0.7;
//    aff2_para[1] = -0.2;
//    aff2_para[2] = -0.3;
//    aff2_para[3] = 0.9;
//    aff2_para[4] = 0.1;
//    aff2_para[5] = 0.2;
    aff2_para[0] = 1;
    aff2_para[1] = 0;
    aff2_para[2] = 0.;
    aff2_para[3] = 1;
    aff2_para[4] = 0;
    aff2_para[5] = 0;

    //    for(int i=0; i<aff2_para.Size(); i++) aff2_para[i] = (rand() % 100) / 50.0;
    affine2->SetParameters(aff2_para);


//    comp->AddTransform(affine2);


    typename CompositeType::InputPointType c0;
    for(int i=0; i<Dim; i++) c0[i]=c[i];
    comp->SetCenter(c0);


    std::cout << "center: " << c0 << std::endl;

    std::cout << "affine2: " << affine2->GetMatrix() << std::endl;


    std::cout << "initial parameters" << std::endl;
    std::cout << comp->GetParameters() << std::endl;

    std::cout << "initial fixed parameters" << std::endl;
    std::cout << comp->GetFixedParameters() << std::endl;

    int nb_pts = x.size();

    const int kParaDim = comp->GetNumberOfParameters();

    ParaType current_para(kParaDim), delta_para(kParaDim);
    current_para = comp->GetParameters();

    int cnt = 0;
    TRawPointList y1;
    y1.resize(x.size());

    while(cnt <= nb_iter){

        delta_para.Fill(0);

        // std::cout << "iter: " << cnt << "---------------" << std::endl;

        for(int i=0; i<nb_pts; i++){
            typedef typename CompositeType::InputPointType InputPointType;
            typedef typename CompositeType::OutputPointType OutputPointType;

            InputPointType ptx;
            OutputPointType pty, ptz;

            for(int d=0;d<Dim;d++){
                ptx[d]=x[i][d];
                pty[d]=y[i][d];
            }
            ptz = comp->TransformPoint(ptx);
            //            std::cout << "ptx:" << ptx
            //                    << "\t" << "pty:" << pty
            //                    << "\t" << "ptz:" << ptz << std::endl;

            typedef typename CompositeType::JacobianType JacobianType;
            JacobianType jac(Dim, kParaDim);

            comp->GetJacobianWithRespectToParameters(ptx, jac);

            for(int d=0; d<Dim; d++){
                for(int k=0; k<kParaDim; k++){
                    delta_para[k] += (ptz[d] - pty[d]) * jac[d][k];
                }
            }
        }

        // delta_para *= 1.0 / nb_pts * 0.5;
        delta_para *= 1.0 / nb_pts;
        for(int i=0; i<kParaDim; i++){
               delta_para[i] *= scale[i];
           }


        current_para -= delta_para;
        comp->SetParameters(current_para);
        cnt++;


        for(int i=0; i<nb_pts; i++){
            typedef typename AffineType::InputPointType InputPointType;
            typedef typename AffineType::OutputPointType OutputPointType;

            InputPointType ptx;
            OutputPointType pty, ptz;

            for(int d=0;d<Dim;d++){
                ptx[d]=x[i][d];
                pty[d]=y[i][d];
            }
            ptz = comp->TransformPoint(ptx);
            y1[i].resize(Dim);
            for(int d=0; d<Dim; d++) y1[i][d] = ptz[d];
        }

        double dmax = max_distance_between_point_list(y, y1);

        if (cnt % 100 == 1)
            std::cout << "iter [" << cnt
            << "]: max ||T(ptx)-pty|| = " << dmax << std::endl;

        if (cnt == nb_iter) std::cout << "tired iterations: iter " << cnt << std::endl;

        if (dmax <= 1e-12) {
            std::cout << "reaching convergence (dmax <1e-6)! "
                    << "iter [" << cnt
                    << "]: max ||T(ptx)-pty|| = " << dmax << std::endl;
            break;
        }

    }


    std::cout << "current para:" << comp->GetParameters() << std::endl;
    std::cout << "current fixed para: " << comp->GetFixedParameters() << std::endl;


    for(int i=0; i<nb_pts; i++){
        typedef typename AffineType::InputPointType InputPointType;
        typedef typename AffineType::OutputPointType OutputPointType;

        InputPointType ptx;
        OutputPointType pty, ptz;

        for(int d=0;d<Dim;d++){
            ptx[d]=x[i][d];
            pty[d]=y[i][d];
        }
        ptz = comp->TransformPoint(ptx);
//        std::cout << "ptx:" << ptx
//                << "\t" << "pty:" << pty
//                << "\t" << "T(ptx)(=?=pty):" << ptz << std::endl;

        y1[i].resize(Dim);
        for(int d=0; d<Dim; d++) y1[i][d] = ptz[d];
    }

    std::cout << "max ||T(ptx)-pty|| = "
            << max_distance_between_point_list(y, y1) << std::endl;

    return;

}


template<int Dim, class TRawVector, class TRawPointList>
void test_centered_composite_RSKT_transform( const TRawPointList &x, const TRawPointList &y,
        const TRawVector &A, const TRawVector &t, const TRawVector &c, int nb_iter){

    // typedef itk::MatrixOffsetTransformBase<double, Dim, Dim> AffineType;
    typedef itk::CenteredCompositeTransform<double, Dim> CompositeType;

    typedef itk::ScaleTransform<double, Dim> ScaleTransformType;

    typename ScaleTransformType::Pointer scale_transform = ScaleTransformType::New();
    scale_transform->SetIdentity();
    typename ScaleTransformType::ParametersType s1(Dim);
    for(int i=0; i<Dim; i++) s1[i]=1; //(rand() % 100) / 50.0 ; //
    scale_transform->SetParameters(s1);


    typedef itk::Shear2DTransform<double, Dim> Shear2DTransformType;

    typename Shear2DTransformType::Pointer shear_transform = Shear2DTransformType::New();
    typename Shear2DTransformType::ParametersType k1(1);
    for(int i=0; i<1; i++) k1[i]=0; //(rand() % 100) / 50.0 - 1.0 ; //
    shear_transform->SetParameters(k1);


//    typedef itk::Rigid2DTransform<double> Rotation2DTransformType;
     typedef itk::Rotate2DTransform<double> Rotation2DTransformType;
    typename Rotation2DTransformType::Pointer rotation_transform = Rotation2DTransformType::New();

    typename Rotation2DTransformType::ParametersType r1(1);
    //angle: (in radius)
    // r1.Fill(-0.9273);
    r1.Fill(0);
//    r1[0] = 0.9273;
//    r1.Fill(53.13);
    // r1[0] = (rand() % 100) / 200.0 - 0.5;


    rotation_transform->SetParameters(r1);
    std::cout << "rotation matrix: " << std:: endl << rotation_transform->GetMatrix() << std::endl;


    typedef itk::TranslationTransform<double, Dim> TranslationTransformType;

    typename TranslationTransformType::Pointer translation_transform = TranslationTransformType::New();
    typename TranslationTransformType::ParametersType t1(Dim);
    for(int i=0; i<Dim; i++) t1[i]=0; // (rand() % 100) / 50.0 - 1.0 ; //
    translation_transform->SetParameters(t1);

    std::cout << "r.center=" << rotation_transform->GetCenter() << std::endl;
    std::cout << "r.fixed=" << rotation_transform->GetFixedParameters() << std::endl;
    std::cout << "r.trans=" << rotation_transform->GetTranslation() << std::endl;


    typename CompositeType::Pointer comp = CompositeType::New();

    // add in the reverse order of applying to the point


    typename Rotation2DTransformType::Pointer rotation_transform2 = Rotation2DTransformType::New();
    rotation_transform2->SetParameters(r1);

    // note: conjecture: rotation can not be between scale and shear
    // A = R*S*K or A = S*K*R, but not A = S*R*K

    comp->AddTransform(translation_transform);

    comp->AddTransform(rotation_transform);
    comp->AddTransform(scale_transform);
//    comp->AddTransform(rotation_transform2);
    comp->AddTransform(shear_transform);



//    comp->AddTransform(rotation_transform);


    float scale[7] = {1,1,1,1,1,0.5,0.5};

//    float scale[6] = {8,4,4,4,0.5,0.5};
//    float scale[7] = {2,4,2,4,4,0.5,0.5};
//    float scale[7] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5};


    // comp->SetNthTransformToOptimize(0, true);
    // comp->SetNthTransformToOptimize(1, true);
    // comp->SetNthTransformToOptimize(2, true);
    // comp->SetNthTransformToOptimize(3, true);



    typename CompositeType::InputPointType c0;
    for(int i=0; i<Dim; i++) c0[i]=c[i];
    comp->SetCenter(c0);


    std::cout << "center: " << c0 << std::endl;

    std::cout << "initial rotation : " << rotation_transform->GetParameters() << std::endl;
    std::cout << "initial shear : " << shear_transform->GetParameters() << std::endl;
    std::cout << "initial scale : " << scale_transform->GetParameters() << std::endl;
    std::cout << "initial translation : " << translation_transform->GetParameters() << std::endl;


    std::cout << "initial composite parameters" << std::endl;
    std::cout << comp->GetParameters() << std::endl;

    std::cout << "initial composite fixed parameters" << std::endl;
    std::cout << comp->GetFixedParameters() << std::endl;

    int nb_pts = x.size();

    typedef typename CompositeType::ParametersType ParaType;
    const int kParaDim = comp->GetNumberOfParameters();
    ParaType current_para(kParaDim), delta_para(kParaDim);
    current_para = comp->GetParameters();

    int cnt = 0;
    TRawPointList y1;
    y1.resize(x.size());

    while(cnt <= nb_iter){

        delta_para.Fill(0);

        // std::cout << "iter: " << cnt << "---------------" << std::endl;

        for(int i=0; i<nb_pts; i++){
            typedef typename CompositeType::InputPointType InputPointType;
            typedef typename CompositeType::OutputPointType OutputPointType;

            InputPointType ptx;
            OutputPointType pty, ptz;

            for(int d=0;d<Dim;d++){
                ptx[d]=x[i][d];
                pty[d]=y[i][d];
            }
            ptz = comp->TransformPoint(ptx);

            typedef typename CompositeType::JacobianType JacobianType;
            JacobianType jac(Dim, kParaDim);

            comp->GetJacobianWithRespectToParameters(ptx, jac);

            // manipulate the jacobian w.r.t the translation of the rigid transform
            // k s r t: 1 + 2 + 3 (angle + translation) + 2
//            jac[0][4] = 0;
//            jac[0][5] = 0;
//            jac[1][4] = 0;
//            jac[1][5] = 0;

            // r t k s ==> s k t r: 2 : 1 : 2 : 3
//            jac[0][6]=jac[0][7]=jac[1][6]=jac[1][7] = 0;

            // t s k r ==> r k s t : 3 : 1 : 2 : 2
            // jac[0][1]=jac[0][2]=jac[1][1]=jac[1][2] = 0;

//            const int OFFSET=3;
//            jac[0][1+OFFSET]=jac[0][2+OFFSET]=jac[1][1+OFFSET]=jac[1][2+OFFSET] = 0;

            for(int d=0; d<Dim; d++){
                for(int k=0; k<kParaDim; k++){


                    delta_para[k] += (ptz[d] - pty[d]) * jac[d][k];
                }
            }
        }

        delta_para *= 1.0 / nb_pts * 0.5; // * 4.0 ; // 0.5;

        for(int i=0; i<kParaDim; i++){
            delta_para[i] *= scale[i];
        }


        current_para -= delta_para;
        comp->SetParameters(current_para);
        cnt++;


        for(int i=0; i<nb_pts; i++){
            typedef typename CompositeType::InputPointType InputPointType;
            typedef typename CompositeType::OutputPointType OutputPointType;

            InputPointType ptx;
            OutputPointType pty, ptz;

            for(int d=0;d<Dim;d++){
                ptx[d]=x[i][d];
                pty[d]=y[i][d];
            }
            ptz = comp->TransformPoint(ptx);
            y1[i].resize(Dim);
            for(int d=0; d<Dim; d++) y1[i][d] = ptz[d];
        }

        double dmax = max_distance_between_point_list(y, y1);

        if (cnt % 100 == 1)
            std::cout << "iter [" << cnt
            << "]: max ||T(ptx)-pty|| = " << dmax << std::endl;

        if (cnt == nb_iter) std::cout << "tired iterations: iter " << cnt << std::endl;

        if (dmax <= 1e-12) {
            std::cout << "reaching convergence (dmax <1e-6)! "
                    << "iter [" << cnt
                    << "]: max ||T(ptx)-pty|| = " << dmax << std::endl;
            break;
        }

    }


    std::cout << "current para:" << comp->GetParameters() << std::endl;
    std::cout << "current fixed para: " << comp->GetFixedParameters() << std::endl;
    std::cout << "current combined K " << comp->GetParameters()[0] * comp->GetParameters()[1] << std::endl;

    std::cout << "final rotation : " << rotation_transform->GetParameters() << std::endl;
    std::cout << "final shear : " << shear_transform->GetParameters() << std::endl;
    std::cout << "final scale : " << scale_transform->GetParameters() << std::endl;
    std::cout << "final translation : " << translation_transform->GetParameters() << std::endl;

    std::cout << "r.center=" << rotation_transform->GetCenter() << std::endl;
    std::cout << "r.fixed=" << rotation_transform->GetFixedParameters() << std::endl;
    std::cout << "r.trans=" << rotation_transform->GetTranslation() << std::endl;

    for(int i=0; i<nb_pts; i++){
        typedef typename CompositeType::InputPointType InputPointType;
        typedef typename CompositeType::OutputPointType OutputPointType;

        InputPointType ptx;
        OutputPointType pty, ptz;

        for(int d=0;d<Dim;d++){
            ptx[d]=x[i][d];
            pty[d]=y[i][d];
        }
        ptz = comp->TransformPoint(ptx);
//        std::cout << "ptx:" << ptx
//                << "\t" << "pty:" << pty
//                << "\t" << "T(ptx)(=?=pty):" << ptz << std::endl;

        y1[i].resize(Dim);
        for(int d=0; d<Dim; d++) y1[i][d] = ptz[d];
    }

    std::cout << "max ||T(ptx)-pty|| = "
            << max_distance_between_point_list(y, y1) << std::endl;


    std::cout << "r.number_of_local_parameter" << rotation_transform->GetNumberOfLocalParameters() << std::endl;

    return;

}

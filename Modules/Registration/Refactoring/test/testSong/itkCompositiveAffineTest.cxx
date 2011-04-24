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
#include <itkTimeProbe.h>

// first test: test the MatrixOffsetTransformBase from the old ITK
// second test: test the scaling transform
// third test: test my own affine (quaternion) transform
// third test: test the itkCompositeTransform works


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

int main(int argc, char **argv){

    const int Dim = 2;
    int n = 10;
    int nb_iter = 100;

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

    itk::TimeProbe t2;
    t2.Start();
    test_composite_transform<Dim, RawVectorType, RawPointListType>(raw_x, raw_y, A, t, c, nb_iter);
    t2.Stop();

    std::cout << "test_one_transform timing:" << t1.GetTotal() << " seconds" << std::endl;
    std::cout << "test_composite_transform timing:" << t2.GetTotal() << " seconds" << std::endl;

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



    float A[Dim*Dim] = { -1, 0, 0, 1 };
    float c[Dim] = { 0.5, 0.5 };
    float t[Dim] = { 1, 1 };

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

            affine->GetLocalJacobian(ptx, jac);

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

            comp->GetLocalJacobian(ptx, jac);

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

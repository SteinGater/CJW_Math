/*********************************************************************************************/
//MMMR 实验室
//作者： CJW
//基于指数映射的刚体机器人数学库
//在该库中，默认通用名称和功能
/*1. 所有矩阵采用数组形式，即T*，如果外部定义的是二维数组指针可以通过强制转换为T*进行操作
矩阵检索方式为：M[i][j]=M[i*行长度+j]
2. 所有矩阵的大小都是非动态，作为参数注意传入的大小
3. 函数名称R=SO(3)表示旋转矩阵,P表示位置，G=SE(3)表示位姿变换矩阵,V表示速度，A表示加速度，
Exp表示指数映射，Log表示对数映射，Ad和ad表示李代数的伴随变换
T表示转置，Inv表示逆，Det表示行列式
4. 旋量采用（s,sxr)形式：速度为（w,v)加速度为（a,av)，力矩为（f,m)，注意：输出与输入保持一致！！！
*/
/*********************************************************************************************/
#ifndef CJW_MATH
#define CJW_MATH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#define USE_EIGEN

#ifdef USE_EIGEN

#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/SVD"

using namespace Eigen;

#endif


/***************************定义常量************************/
#define    M_PI          3.14159265358979323846   // pi
#define    MATH_ZERO     0.000000001
#define    EXP_ZERO      0.001
#define    INF           100000000
#define    M_ComputeMAX  100

/***************************定义计算错误************************/
#define    M_RIGHT       0
#define    M_ERROR       1
/***************************定义最大矩阵大小************************/
#define    M_N_MAX        30
/***************************定义点是否在平面多边形内******************/
#define IN_POLYGON_OUT     -1
#define IN_POLYGON_SIDE     0
#define IN_POLYGON_IN     1


#define IN_POLYGON_CROSS   0
#define IN_POLYGON_ANGLE   1

/**************************重力加速度*********************************/
#define WORLD_GRAVITY 9.81
/********************************************************
  basis function
  ******************************************************/
template <typename TT>
class CJW_Math 
{
public:
    /*************Constructor****************************************************************/
    CJW_Math(void){}
    ~CJW_Math(void){}


    /*************矩阵显示***********************************************************/
    void Show(int mym,int myn,TT *A);//矩阵显示

    /*************三维矢量运算***********************************************************/
    TT MyVector3Dot(TT Vect1[3],TT Vect2[3]);//矢量内积
    void MyVector3Cross(TT Vect1[3],TT Vect2[3],TT out[3]);//矢量外积
    void MyVector3ToDisMatrix(TT vector[3],TT Matrix[9]);//3维矢量反对称运算
    void MyDisMatrixToVector3(TT Matrix[9],TT vector[3]);//3维矢量反对称运算
    void MyVector6ToDisMatrix(TT vector[6],TT Matrix[16]);//6维矢量反对称运算
    void MyDisMatrixToVector6(TT Matrix[16],TT vector[6]);//6维矢量反对称运算
    TT MyVectorDot(TT* Vector1,TT* Vector2,int VectorL);//向量对应乘法
    TT MyVectorNorm1(TT* Vector,int VectorL);//矢量1-范数
    TT MyVectorNorm2(TT* Vector,int VectorL);//矢量2-范数
    /*************基本矩阵运算***********************************************************/
    void MyMatrixCopy(int mn, TT* A,TT *Acopy);//矩阵拷贝
    void MyMatrixAdd(int mn,TT *A,TT *B,TT *out);//矩阵相加A+B=cout
    void MyMatrixSub(int mn,TT *A,TT *B,TT *out);//矩阵减法A-B=out
    void MyMatrixMultiply(int mym,int myn,int myp,TT *A,TT *B,TT *out);//矩阵乘法A*B=out
    void MyMatrixMultiplyAT(int mym,int myn,int myp,TT *A,TT *B,TT *out);//矩阵转置乘法AT*B=out
    void MyMatrixMultiplyBT(int mym,int myn,int myp,TT *A,TT *B,TT *out);//矩阵转置乘法A*BT=out
    void MyMatrixRorate(int mym,int myn,TT *A,TT *out);//矩阵转置
    void MyMatrixQR_GS(int mym,TT *A,TT *Q,TT *R);//方阵高斯的QR分解
    TT MyMatrixdet(int mym,TT* A);//N阶方阵行列式
    int MyMatrixEigen_Jacobi(int myn,TT *A,TT *Q,TT* E);//N阶对称方阵的特征向量和特征矩阵，使用Jacobian搜索方法
    int MyMatrixNormal_Jacobi(int mym,int myn,TT *A,TT *V,TT* E);//N阶对称方阵的特征向量和特征矩阵，使用Jacobian搜索方法
    int MyMatrixSide_Jacobi(int mym,int myn,TT *A,TT *V,TT* UE);//单边Jacobian分解，其中行<=列
    int MyMatrixSVD_Jacobi(int mym,int myn,TT *A,TT *U,TT* V,TT* E);//3矩阵SVD分解，计算出两个单位正交矩阵，使用Jacobian方法
    int MyMatrix3SVD_Jacobi(TT *A,TT *U,TT* V,TT* E);//3矩阵SVD分解，计算出两个单位正交矩阵，使用Jacobian方法
    //使用一种快速地推算法（存在局限性）
    int  MyMatrixInv(int myn,TT *A,TT *Inv);//N阶方阵逆矩阵//每行组合都能实现逆矩阵

    /*************高斯数值法求解Ax=b的解（最大阶数为10）***********************************************************/
    int MySolveGS(int mym,TT* A, TT* b,TT* x);

    /*****************************三阶矩阵运算*****************************/
    TT   MyMatrix3det(TT A[9]);//三阶矩阵行列式
    int  MyMatrix3Inv(TT A[9],TT out[9]);//三阶矩阵逆矩阵

    /*****************************旋转矩阵运算****************************/
    void MyRCompositionR(TT R1[9], TT R2[9], TT Rsum[9]);//旋转矩阵相乘R1*R2=Rsum
    void MyRCompositionw(TT theR[9],TT thew[3],TT outw[3]);//旋转矩阵乘矢量theR*thew=outw
    void MyRCompositionwAddP(TT theR[9],TT thew[3],TT addP[3],TT outw[3]);//旋转矩阵乘矢量theR*thew+addP=outw
    void MyRInv(TT R[9],TT RT[9]);//旋转矩阵转置=逆矩阵
    void MyRInvCompositionR(TT R1[9], TT R2[9], TT Rsum[9]);//旋转矩阵转置乘法R1T*R2=Rsum
    void MyRInvCompositionw(TT theR[9],TT thew[3],TT outw[3]);//旋转逆解矩阵乘矢量theRInv*thew=outw
    void MyRInvCompositionwAddP(TT theR[9],TT thew[3],TT addP[3],TT outw[3]);//旋转矩阵乘矢量theRInv*thew+addP=outw

    /*****************************位姿矩阵运算****************************/
    void MyGCompositionG(TT G1[16], TT G2[16], TT GG[16]);//位姿矩阵相乘G1*G2=GG
    void MyGCompositionP(TT G[16], TT P[3], TT outP[3]);//位姿矩阵乘法位置G1*P=outP
    void MyGInv(TT G[16],TT InvG[16]);//位姿矩阵逆解
    void MyGInvCompositionG(TT G1[16], TT G2[16], TT GG[16]);//位姿矩阵逆解乘法G1Inv*G2=GG
    void MyGInvCompositionP(TT G[16], TT P[3], TT outP[3]);//位姿矩阵逆解乘法位置G1Inv*P=outP
    void MyRPToG(TT R[9],TT P[3],TT G[16]);//姿态+位置=位姿矩阵
    void MyGToRP(TT G[16],TT R[9],TT P[3]);//位姿矩阵=姿态+位置
    void MyGToRInvP(TT G[16],TT RT[9],TT P[3]);//位姿矩阵=姿态逆解+位置

    /*****************************姿态矩阵 和 其他表达的转换****************************/
    void MyEulerZYXToR(TT w[3], TT R[9]);//欧拉角ZYX转位姿矩阵
    void MyRToEulerZYX(TT R[9], TT w[3]);//姿态矩阵转欧拉角ZYX
    void MyQuaternionToR(TT q[4],TT R[9]);//q为 x y z w 四元数转旋转矩阵

    /*****************************基于指数坐标表达****************************/
    void MyExponent3ToR(TT w[3], TT R[9]);//指数坐标转旋转矩阵
    void MyRToExponent3(TT R[9], TT w[3]);//旋转矩阵转指数坐标
    void MyExponent3ToQuaternion(TT w[3], TT q[4]);//指数坐标转四元数xyzw
    void MyQuaternionToExponent3(TT q[4],TT w[3]);//四元数xyzw转指数坐标
    void MyExponent4ToG(TT w[6], TT R[16]);//指数坐标转位姿矩阵
    void MyGToExponent4(TT R[16], TT w[6]);//位姿矩阵转指数坐标

    /****************************************伴随运算***************************************/
    void MyEXPAdg(TT G[16],TT Adg[36]);//伴随运算
    void MyEXPAdgScrew(TT G[16],TT S[6],TT AdgV[6]);//伴随运算乘旋量
    void MyEXPAdgT(TT G[16],TT AdgT[36]);//伴随转置
    void MyEXPAdgTScrew(TT G[16],TT S[6],TT AdgTV[6]);//伴随转置乘旋量
    void MyEXPAdgInv(TT G[16],TT AdgInv[36]);//伴随逆解
    void MyEXPAdgInvScrew(TT G[16],TT S[6],TT AdgInvV[6]);//伴随逆解乘旋量
    void MyEXPAdgInvT(TT G[16],TT AdgInvT[36]);//伴随逆解转置
    void MyEXPAdgInvTScrew(TT G[16],TT S[6],TT AdgInvTV[6]);//伴随逆解转置乘旋量
    void MyEXPad(TT w[6],TT ad[36]);//旋量伴随运算
    void MyEXPadT(TT w[6],TT ad[36]);//旋量伴随转置

    /********************指数坐标的微分运算****************************************************/
    void MydEXP3(TT w[3],TT dexp[9]);//指数坐标一阶导数
    void MydEXP3Inv(TT w[3],TT dexp[9]);//指数坐标一阶导数逆矩阵
    void MyddEXP3(TT w[3],TT dw[3],TT ddexp[9]);//指数坐标二阶导数
    void MyddEXP3Inv(TT w[3],TT dw[3],TT ddexp[9]);//指数坐标二阶导数逆矩阵
    void MydEXP4(TT w[6],TT dexp[36]);//指数坐标一阶导数
    void MydEXP4Inv(TT w[6],TT dexp[36]);//指数坐标一阶导数逆矩阵
    void MyddEXP4(TT w[6],TT dw[6],TT ddexp[36]);//指数坐标二阶导数
    void MyddEXP4Inv(TT w[6],TT dw[6],TT ddexp[36]);//指数坐标二阶导数逆矩阵

    /********************一些指数坐标的转换*************************************************/
    //物体的姿态和物体坐标系速度 转指数坐标系
    void MyRVToExponet3(TT R[9],TT w[3],TT kesi[3],TT dkesi[3]);
    //物体的姿态和物体坐标系速度，加速度 转指数坐标系
    void MyRVAToExponet3(TT R[9],TT w[3],TT dw[3],TT kesi[3],TT dkesi[3],TT ddkesi[3]);
    //指数坐标系 转物体的姿态和物体坐标系速度
    void MyExponet3ToRV(TT kesi[3],TT dkesi[3],TT R[9],TT w[3]);
    //指数坐标系 转物体的姿态和物体坐标系速度，加速度
    void MyExponet3ToRVA(TT kesi[3],TT dkesi[3],TT ddkesi[3],TT R[9],TT w[3],TT dw[3]);
    //物体的位姿和物体坐标系速度 转指数坐标系
    void MyGVToExponet4(TT G[16],TT w[6],TT kesi[6],TT dkesi[6]);
    //物体的位姿和物体坐标系速度，加速度 转指数坐标系
    void MyGVAToExponet4(TT G[16],TT w[6],TT dw[6],TT kesi[6],TT dkesi[6],TT ddkesi[6]);
    //指数坐标系 转物体的位姿和物体坐标系速度
    void MyExponet4ToGV(TT kesi[6],TT dkesi[6],TT G[16],TT w[6]);
    //指数坐标系 转物体的位姿和物体坐标系速度，加速度
    void MyExponet4ToGVA(TT kesi[6],TT dkesi[6],TT ddkesi[6],TT G[16],TT w[6],TT dw[6]);
    //SE(3)空间映射成SE（2）
    void MySE3ToSE2(TT G0[16],TT Gout[16]);
    void MySE3ToSE2(TT G0[16],TT V0[6],TT Gout[16],TT Vout[6]);

    /********************单刚体指数坐标控制*************************************************/
    void RigidBody_SO3_Pcontroller(TT Re[9],TT Ve[3],TT Rr[9],TT Kp[3],TT Vout[3]);
    void RigidBody_SE3_Pcontroller(TT Ge[16],TT Ve[6],TT Gr[16],TT Kp[6],TT Vout[6]);
    void RigidBody_SO3_PDcontroller(TT Re[9],TT Ve[3],TT Ae[3],TT Rr[9],TT Vr[3],TT Kp[3],TT Kd[3],TT Aout[3]);
    void RigidBody_SE3_PDcontroller(TT Ge[16],TT Ve[6],TT Ae[6],TT Gr[16],TT Vr[6],TT Kp[6],TT Kd[6],TT Aout[6]);

    /********************广义惯性矩阵的坐标转换****************************************************/
    //trans based the body M to the absolute coordinate M
    void MyInertiaBasedGroup(TT M[36],TT G[16],TT outM[36]);
    void MyInertiaDBasedGroup(TT M[36],TT dM[36],TT G[16],TT VS[6],TT outM[36],TT outdM[36]);

    /*******************************基于惯性中心的动力学方程*****************************************/
    void DynamicBaseCOI(TT BodyM[36],TT BodydM[36],TT BodyV[6],TT BodyA[6],TT outFref[6]);

    /*************************************串联指数积公式*************************************/
    void MyExponent3ToR(TT ww[3], TT th, TT R[9]);//轴-角 转旋转矩阵
    void MyExponent4ToG(TT ww[6], TT th, TT G[16]);//旋量-角度 转位姿矩阵
    int MySeriesR(TT Rout[9], int N, ...);//旋转矩阵连乘Rout=R1*...*RN,变参数为关节角速度
    int MySeriesG(TT Gout[16], int N, ...);//位姿矩阵连乘Gout=G1*...*GN,变参数为关节旋量
    //轴-角 转 旋转矩阵连乘,变参数为初始状态绝对坐标系的关节角速度
    int MySeriesExp3ToR(TT Rout[9], TT* th, int N,  ...);
    //旋量-角度 转 位姿矩阵连乘,变参数为初始状态绝对坐标系的关节旋量
    int MySeriesExp4ToG(TT Gout[16], TT* th,int N,  ...);
    //旋量-角度 转 位姿矩阵连乘的雅克比-根部坐标系,变参数为初始状态绝对坐标系的关节旋量
    int MySeriesScrewToJacobianS( TT* Jout, TT* th,int N, ...);
    //旋量-角度 转 位姿矩阵连乘的雅克比-末端坐标系,变参数为初始状态绝对坐标系的关节旋量
    int MySeriesScrewToJacobianB( TT* Jout, TT* th,int N, ...);
    //旋量-角度 转 位姿矩阵连乘的质心-根部坐标系,变参数为初始状态绝对坐标系的质心位置和关节旋量：m_p[0],screw[0],m_p[1],screw[1] ...
    int MySeriesScrewMassToCoM(TT* sum_mass,TT CoMout[3], TT* mass, TT* th,int N, ...);
    //旋量-角度 转 位姿矩阵连乘的广义惯性矩阵-根部坐标系,变参数为初始状态绝对坐标系的广义惯性矩阵和关节旋量：M[0],screw[0],M[1],screw[1] ...
    int MySeriesScrewInertiaToCoI(TT CoIout[36], TT* th, int N, ...);

    /*************************************************一些几何运算****************************************/
    //三位空间点P0到（PA-PB）组成的直线距离
    TT Distance_PointToLine(TT P0[3],TT PA[3],TT PB[3]);
    //空间内两条直线距离 A1-A2和B1-B2的距离
    TT Distanc_LineToLine(TT A1[3],TT A2[3],TT B1[3],TT B2[3]);
    //判断平面内 点 是否在 多边形(逆时针，点表示） 内部 //返回-1表示不在，0表示在边界，1表示在内部
    //特别的 pplygon 存储顶点信息 X1 Y1 Z1 X2 Y2 Z2。。。。。
    int Inpolygon(int type,TT P0[2],TT* polygon,int polyN);
    //判断平面内 点 距离 多边形(逆时针，点表示） 最近距离（内部大于0，外部小于0）
    //特别的 pplygon 存储顶点信息 X1 Y1 Z1 X2 Y2 Z2。。。。
    TT MinDis_Inpolygon(int type,TT P0[3],TT* polygon,int polyN);
    //叉乘算法计算点是否在多边形内部及距离，只适用于凸多边形
    TT TriPoint_Plane_Cross(TT p0[2],TT p1[2],TT p2[2]);//平面三点P1-P0差乘P2-P0
    int Inpolygon_Cross(TT P0[2],TT* polygon,int pointN);//内部 //返回-1表示不在，0表示在边界，1表示在内部
    TT MinDis_Inpolygon_Cross(TT P0[3],TT* polygon,int pointN);//最近距离（内部大于0，外部小于0）
    //内角合方法计算点是否在多边形内部及距离，只适用于所有多边形
    TT TriPoint_Plane_Angle(TT p0[2],TT p1[2],TT p2[2]);//平面三点夹角P1-P0与P2-P0
    int Inpolygon_Angle(TT P0[2],TT* polygon,int pointN);//内部 //返回-1表示不在，0表示在边界，1表示在内部
    TT MinDis_Inpolygon_Angle(TT P0[3],TT* polygon,int pointN);//最近距离（内部大于0，外部小于0）



protected:

private:


};


#endif

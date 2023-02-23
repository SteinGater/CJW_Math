#include "CJW_Math.h"


/*************矩阵显示***********************************************************/
template <typename TT>
void CJW_Math<TT>::Show(int mym,int myn,TT *A)//矩阵显示
{
    //printf("\n the Matrix is\n");
    for(int i=0;i<mym;i++)//row
    {
        for(int j=0;j<myn;j++)//column
        {
            printf("%.3f\t",A[i*myn+j]);
        }
        printf("\n");
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixCopy(int mn, TT* A,TT *Acopy)//矩阵拷贝
{
    for(int i=0;i<mn;i++)
    {
        Acopy[i]=A[i];
    }
}
/*************三维矢量运算***********************************************************/
template <typename TT>
TT CJW_Math<TT>::MyVector3Dot(TT Vect1[3],TT Vect2[3])//矢量内积
{
    return (Vect1[0]*Vect2[0]+Vect1[1]*Vect2[1]+Vect1[2]*Vect2[2]);
}
template <typename TT>
void CJW_Math<TT>::MyVector3Cross(TT Vect1[3],TT Vect2[3],TT out[3])//矢量外积
{
    out[0]=Vect1[1]*Vect2[2]-Vect1[2]*Vect2[1];
    out[1]=Vect1[2]*Vect2[0]-Vect1[0]*Vect2[2];
    out[2]=Vect1[0]*Vect2[1]-Vect1[1]*Vect2[0];
}
template <typename TT>
void CJW_Math<TT>::MyVector3ToDisMatrix(TT vector[3],TT Matrix[9])//3维矢量反对称运算
{
    Matrix[0]=0;                 Matrix[1]=-vector[2];          Matrix[2]=vector[1];
    Matrix[3]=vector[2];    Matrix[4]=0;                        Matrix[5]=-vector[0];
    Matrix[6]=-vector[1];   Matrix[7]=vector[0];           Matrix[8]=0;
}
template <typename TT>
void CJW_Math<TT>::MyDisMatrixToVector3(TT Matrix[9],TT vector[3])//3维矢量反对称运算
{
    vector[0]=Matrix[7];vector[1]=Matrix[2];vector[2]=Matrix[3];
}
template <typename TT>
void CJW_Math<TT>::MyVector6ToDisMatrix(TT vector[6],TT Matrix[16])//6维矢量反对称运算
{
    Matrix[0]=0;                 Matrix[1]=-vector[2];          Matrix[2]=vector[1];        Matrix[3]=vector[3];
    Matrix[4]=vector[2];    Matrix[5]=0;                        Matrix[6]=-vector[0];       Matrix[7]=vector[4];
    Matrix[8]=-vector[1];   Matrix[9]=vector[0];           Matrix[10]=0;                   Matrix[11]=vector[5];
    Matrix[12]=0;               Matrix[13]=0;                      Matrix[14]=0;                   Matrix[15]=0;
}
template <typename TT>
void CJW_Math<TT>::MyDisMatrixToVector6(TT Matrix[16],TT vector[6])//6维矢量反对称运算
{
    vector[0]=Matrix[9];vector[1]=Matrix[2];vector[2]=Matrix[4];
    vector[3]=Matrix[3];vector[4]=Matrix[7];vector[5]=Matrix[11];
}
template <typename TT>
TT CJW_Math<TT>::MyVectorDot(TT* Vector1,TT* Vector2,int VectorL)//向量对应乘法
{
    TT result=0;
    for(int i=0;i<VectorL;i++){result=result+Vector1[i]*Vector2[i];}
    return result;
}

template <typename TT>
TT CJW_Math<TT>::MyVectorNorm1(TT* Vector,int VectorL)//矢量1-范数
{
    TT norm0=0;
    for(int i=0;i<VectorL;i++){norm0=norm0+fabs(Vector[i]);}
    return norm0;
}
template <typename TT>
TT CJW_Math<TT>::MyVectorNorm2(TT* Vector,int VectorL)//矢量2-范数
{
    TT norm0=0;
    for(int i=0;i<VectorL;i++){norm0=norm0+Vector[i]*Vector[i];}
    return sqrt(norm0);
}
/*************基本矩阵运算***********************************************************/
template <typename TT>
void CJW_Math<TT>::MyMatrixAdd(int mn,TT *A,TT *B,TT *out)//矩阵相加A+B=cout
{
    //A is m*p,B is p*n, C is out with m*n;
    for(int i=0;i<mn;i++)//row
    {
        out[i]=A[i]+B[i];
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixSub(int mn,TT *A,TT *B,TT *out)//矩阵减法A-B=out
{
    //A is m*p,B is p*n, C is out with m*n;
    for(int i=0;i<mn;i++)//row
    {
        out[i]=A[i]-B[i];
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixMultiply(int mym,int myn,int myp,TT *A,TT *B,TT *out)//矩阵乘法A*B=out
{
    //A is m*p,B is p*n, C is out with m*n;
    for(int i=0;i<mym;i++)//row
    {
        for(int j=0;j<myn;j++)//column
        {
            TT element=0;
            for(int k=0;k<myp;k++)
            {
                element=element+A[i*myp+k]*B[k*myn+j];
            }
            out[i*myn+j]=element;
        }
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixMultiplyAT(int mym,int myn,int myp,TT *A,TT *B,TT *out)//矩阵转置乘法AT*B=out
{
    //A is m*p,B is p*n, C is out with m*n;
    for(int i=0;i<mym;i++)//row
    {
        for(int j=0;j<myn;j++)//column
        {
            TT element=0;
            for(int k=0;k<myp;k++)
            {
                element=element+A[k*mym+i]*B[k*myn+j];
            }
            out[i*myn+j]=element;
        }
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixMultiplyBT(int mym,int myn,int myp,TT *A,TT *B,TT *out)//矩阵转置乘法AT*B=out
{
    //A is m*p,B is p*n, C is out with m*n;
    for(int i=0;i<mym;i++)//row
    {
        for(int j=0;j<myn;j++)//column
        {
            TT element=0;
            for(int k=0;k<myp;k++)
            {
                element=element+A[i*myp+k]*B[j*myp+k];
            }
            out[i*myn+j]=element;
        }
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixRorate(int mym,int myn,TT *A,TT *out)//矩阵转置
{
    //A is m*n,out is n*m;
    for(int i=0;i<mym;i++)//row
    {
        for(int j=0;j<myn;j++)//column
        {
            out[j*mym+i]=A[i*myn+j];
        }
    }
}
template <typename TT>
void CJW_Math<TT>::MyMatrixQR_GS(int mym,TT *A,TT *Q,TT *R)//方阵高斯的QR分解
{
    int mymmym=mym*mym;
    //set the R={0};
    for(int i=0;i<mymmym;i++){R[i]=0;}
    //initial
    TT flag=0;
    for(int i=0;i<mym;i++){flag=flag+A[i*mym]*A[i*mym];}
    if(flag<=0){for(int i=0;i<mym;i++){Q[i*mym]=0;}}
    else{flag=sqrt(flag);R[0]=flag;for(int i=0;i<mym;i++){Q[i*mym]=A[i*mym]/flag;}}
    //开始叠加
    int then=0;
    for(int j=1;j<mym;j++)
    {
        for(int i=0;i<mym;i++)
        {
            then=i*mym+j;
            Q[then]=A[then];
        }
        for(int i=0;i<j;i++)
        {
            then=i*mym+j;
            R[then]=0;
            for(int s=0;s<mym;s++)
            {
                R[then]=R[then]+Q[s*mym+i]*Q[s*mym+j];
            }
            for(int s=0;s<mym;s++)
            {
                Q[s*mym+j]=Q[s*mym+j]-R[then]*Q[s*mym+i];
            }
        }
        flag=0;
        for(int i=0;i<mym;i++){flag=flag+Q[i*mym+j]*Q[i*mym+j];}
        if(flag>0)
        {
            flag=sqrt(flag);R[j*mym+j]=flag;
            for(int i=0;i<mym;i++){Q[i*mym+j]=Q[i*mym+j]/flag;}
        }
    }
}
template <typename TT>
TT CJW_Math<TT>::MyMatrixdet(int mym,TT* A)//N阶方阵行列式
{
    TT A_R[M_N_MAX*M_N_MAX]={0};
    //record the data
    for(int i=0;i<mym;i++)
    {
        int mymi=mym*i;
        for(int j=0;j<mym;j++) {A_R[mymi+j]=A[mymi+j];}
    }
    int ChangeN=0;
    //GS get the QR function
    for(int i=0;i<(mym-1);i++)
    {
        int mymi=mym*i;
        for(int j=i+1;j<mym;j++)
        {
            int mymj=mym*j;
            if(A_R[mymj+i]!=0)//can be GS
            {
                if(A_R[mymi+i]==0)//change the postion
                {
                    ChangeN++;
                    TT thetemp;
                    for(int thec=i;thec<mym;thec++)
                    {
                        thetemp=A_R[mymj+thec];
                        A_R[mymj+thec]=A_R[mymi+thec];
                        A_R[mymi+thec]=thetemp;
                    }
                }
                else//the norm GS
                {
                    TT TheK=A_R[mymi+i]/A_R[mymj+i];
                    A_R[mymj+i]=0;
                    for(int k=i+1;k<mym;k++) {A_R[mymj+k]=TheK*A_R[mymj+k]-A_R[mymi+k];}
                }
            }
        }
    }
    TT det0=1;if(ChangeN%2==1) {det0=-1;}
    for(int i=0;i<(mym*mym);i=i+mym+1){det0=det0*A_R[i];}
    return det0;
}
template <typename TT>
int CJW_Math<TT>::MyMatrixInv(int myn,TT *A,TT *Inv)//N阶方阵逆矩阵//每行组合都能实现逆矩阵
{
    int mynmyn=myn*myn;
    TT theD[M_N_MAX*M_N_MAX];
    for(int i=0;i<mynmyn;i++){theD[i]=A[i];}
    for(int i=0;i<mynmyn;i=i+myn+1){theD[i]=theD[i]-1;}
    TT InvIf=1+theD[0];
    if(fabs(InvIf)<=MATH_ZERO){return M_ERROR;}
    TT InvCk[M_N_MAX*M_N_MAX]={0};
    for(int i=0;i<mynmyn;i=i+myn+1){InvCk[i]=1;}
    for(int i=0;i<myn;i++){InvCk[i]=InvCk[i]-theD[i]/InvIf;}
    TT TempdC[M_N_MAX];TT TempCkk[M_N_MAX];
    TT TempCkdC[M_N_MAX*M_N_MAX];
    for(int i=1;i<myn;i++)
    {
        MyMatrixMultiply(1,myn,myn,&(theD[myn*i]),InvCk,TempdC);
        InvIf=1+TempdC[i];if(fabs(InvIf)<=MATH_ZERO) {return i+1;}
        for(int j=0;j<myn;j++){TempCkk[j]=InvCk[j*myn+i];}
        MyMatrixMultiply(myn,myn,1,TempCkk,TempdC,TempCkdC);
        for(int k=0;k<mynmyn;k++)
        {
            InvCk[k]=InvCk[k]-TempCkdC[k]/InvIf;
        }
    }
    for(int i=0;i<mynmyn;i++){Inv[i]=InvCk[i];}
    return M_RIGHT;
}

template <typename TT>
int CJW_Math<TT>::MyMatrixEigen_Jacobi(int myn,TT *A,TT *Q,TT* E)//N阶对称方阵的特征向量和特征矩阵，使用Jacobian搜索方法
{
    int mynmyn=myn*myn;
    TT AA[M_N_MAX*M_N_MAX]={0};
    TT UU[M_N_MAX*M_N_MAX]={0};
    //初始化
    for(int i=0;i<mynmyn;i++){Q[i]=0;AA[i]=A[i];E[i]=AA[i];}
    for(int i=0;i<mynmyn;i=i+myn+1){Q[i]=1;UU[i]=1;}
    //记录最大点
    int maxi=0;int maxj=0;int maxp=0;
    TT maxmax=INF;
    TT theth=0;
    TT tempU[M_N_MAX*M_N_MAX]={0};
    TT tempQ[M_N_MAX*M_N_MAX]={0};
    TT ComputeN=-1;
    while(ComputeN<M_ComputeMAX)
    {
        ComputeN++;
        //寻找最大点
        maxmax=0;
        for(int i=0;i<myn-1;i++)
        {
            for(int j=i+1;j<myn;j++)
            {
                int then=myn*i+j;
                if(fabs(AA[then])>fabs(maxmax))
                {
                    maxi=i;maxj=j;maxp=then;
                    maxmax=AA[then];
                }
            }
        }
        if(fabs(maxmax)<=MATH_ZERO){break;}
        //更新U
        for(int i=0;i<mynmyn;i++){UU[i]=0;}
        for(int i=0;i<mynmyn;i=i+myn+1){UU[i]=1;}
        theth=atan2(2*AA[maxp],AA[maxi*myn+maxi]-AA[maxj*myn+maxj])/2;
        UU[maxi*myn+maxi]=cos(theth);UU[maxp]=-sin(theth);
        UU[maxj*myn+maxi]=-UU[maxp];UU[maxj*myn+maxj]=UU[maxi*myn+maxi];
        MyMatrixMultiplyAT(myn,myn,myn,UU,AA,tempU);
        MyMatrixMultiply(myn,myn,myn,tempU,UU,E);
        MyMatrixMultiply(myn,myn,myn,Q,UU,tempQ);
        for(int i=0;i<mynmyn;i++){Q[i]=tempQ[i];AA[i]=E[i];}
    }
    return ComputeN;
}
template <typename TT>
int CJW_Math<TT>::MyMatrixNormal_Jacobi(int mym,int myn,TT *A,TT *V,TT* E)//行正交方阵的特征向量和特征矩阵，使用Jacobian搜索方法
{
    //初始化
    for(int i=0;i<mym*myn;i++){E[i]=A[i];}
    for(int i=0;i<myn*myn;i++){V[i]=0;}
    for(int i=0;i<myn*myn;i=i+myn+1){V[i]=1;}
    //记录最大点
    int maxi=0;int maxj=0;
    TT maxmax=INF;TT theth=0;
    int ComputeN=-1;
    while(ComputeN<M_ComputeMAX)
    {
        ComputeN++;
        //寻找最大点
        maxmax=0;
        for(int i=0;i<mym;i++)
        {
            for(int j=0;j<myn;j++)
            {
                if(i!=j)
                {
                    int then=myn*i+j;
                    if(fabs(E[then])>fabs(maxmax))
                    {
                        maxi=i;maxj=j;
                        maxmax=E[then];
                    }
                }
            }
        }
        if(fabs(maxmax)<=MATH_ZERO){break;}
        //计算旋转矩阵
        if(maxi>maxj){theth=atan2(-E[maxi*myn+maxj],E[maxi*myn+maxi]);}
        else{theth=atan2(E[maxi*myn+maxj],E[maxi*myn+maxi]);}
        TT cos0=cos(theth);TT sin0=sin(theth);
        TT Vij[2];
        for(int i=0;i<mym;i++)
        {
            int i0=i*myn;
            Vij[0]=E[i0+maxi]*cos0+E[i0+maxj]*sin0;
            Vij[1]=-E[i0+maxi]*sin0+E[i0+maxj]*cos0;
            E[i0+maxi]=Vij[0];
            E[i0+maxj]=Vij[1];
        }
        for(int i=0;i<myn;i++)
        {
            int i0=i*myn;
            Vij[0]=V[i0+maxi]*cos0+V[i0+maxj]*sin0;
            Vij[1]=-V[i0+maxi]*sin0+V[i0+maxj]*cos0;
            V[i0+maxi]=Vij[0];V[i0+maxj]=Vij[1];
        }
    }
    return ComputeN;
}
template <typename TT>
int CJW_Math<TT>::MyMatrixSide_Jacobi(int mym,int myn,TT *A,TT *V,TT* UE)//单边Jacobian分解，其中行<=列
{
    int mymmyn=mym*myn;int mynmyn=myn*myn;
    //初始化
    for(int i=0;i<mymmyn;i++){UE[i]=A[i];}
    for(int i=0;i<mynmyn;i++){V[i]=0;}
    for(int i=0;i<mynmyn;i=i+myn+1){V[i]=1;}
    //开始迭代
    TT ComputeN=0;int Flag=1;
    TT C1[M_N_MAX]={0};TT C2[M_N_MAX]={0};
    while((ComputeN++)<M_ComputeMAX)
    {
        Flag=1;
        for(int row1=0;row1<(myn-1);row1++)
        {
            for(int row2=row1+1;row2<myn;row2++)
            {
                //提取当前C1和C2
                for(int i=0;i<mym;i++){C1[i]=UE[i*myn+row1];C2[i]=UE[i*myn+row2];}
                //printf("C1 and C2 is \n");Show(1,mym,C1);Show(1,mym,C2);
                //判断两列是否正交
                TT Norm12=MyVectorDot(C1,C2,mym);
                if(fabs(Norm12)>MATH_ZERO)
                {
                    Flag=0;
                    TT Norm11=MyVectorDot(C1,C1,mym);
                    TT Norm22=MyVectorDot(C2,C2,mym);
                    //排序特征//可以没有
                    //计算旋转矩阵
                    TT th=atan2(Norm12*2,Norm11-Norm22)/2;
                    TT cos0=cos(th);TT sin0=sin(th);
                    TT Vij[2];
                    for(int i=0;i<mym;i++)
                    {
                        int i0=i*myn;
                        UE[i0+row1]=C1[i]*cos0+C2[i]*sin0;
                        UE[i0+row2]=-C1[i]*sin0+C2[i]*cos0;
                    }
                    for(int i=0;i<myn;i++)
                    {
                        int i1=i*myn+row1;int i2=i*myn+row2;
                        Vij[0]=V[i1]*cos0+V[i2]*sin0;
                        Vij[1]=-V[i1]*sin0+V[i2]*cos0;
                        V[i1]=Vij[0];V[i2]=Vij[1];
                    }
                }
            }
        }
        if(Flag) {break;}//形成单位矩阵
    }
    return ComputeN;
}
template <typename TT>
int CJW_Math<TT>::MyMatrixSVD_Jacobi(int mym,int myn,TT *A,TT *U,TT* V,TT* E)//一般矩阵SVD分解，计算出两个单位正交矩阵，使用Jacobian方法
{
    /***********混合算法***********/
    int error=M_RIGHT;
    TT AA[M_N_MAX*M_N_MAX]={0};
    MyMatrixMultiplyBT(mym,mym,myn,A,A,AA);
    error=MyMatrixEigen_Jacobi(mym,AA,U,E);
    if(error>=M_ComputeMAX){return error;}
    //for(int i=0;i<fmin(mym,myn);i++){E[i*myn+i]=sqrt(E[i*myn+i]);}
    TT UTA[M_N_MAX*M_N_MAX]={0};
    MyMatrixMultiplyAT(mym,myn,mym,U,A,UTA);
    error=MyMatrixNormal_Jacobi(mym,myn,UTA,V,E);
    if(error>=M_ComputeMAX){return error;}
    for(int i=0;i<myn;i++)
    {
        TT theflag=E[i*myn+i];
        if(theflag<0)
        {
            E[i*myn+i]=fabs(theflag);
            for(int j=0;j<myn;j++)
            {
                V[i*myn+j]=-V[i*myn+j];
            }
        }
    }
//    printf("the A is \n");Show(mym,myn,A);
//    printf("the AA is \n");Show(mym,mym,AA);
//    printf("the E is \n");Show(mym,myn,E);
//    printf("the V is \n");Show(myn,myn,V);
//    printf("the U is \n");Show(mym,mym,U);
    return error;
}
template <typename TT>
int CJW_Math<TT>::MyMatrix3SVD_Jacobi(TT *A,TT *U,TT* V,TT* E)//一般矩阵SVD分解，计算出两个单位正交矩阵，使用Jacobian方法
{
    /******单边算法**************************/
    int error=M_RIGHT;
    for(int i=0;i<9;i++){E[i]=0;U[i]=0;}
    TT EUT[9]={0};
    //单边运算
    error=MyMatrixSide_Jacobi(3,3,A,V,EUT);
    if(error>=M_ComputeMAX){return error;}
    int detflag=0;int detN=0;
    for(int i=0;i<3;i++)
    {
        TT theNorm0=0;for(int j=0;j<3;j++){theNorm0=theNorm0+EUT[j*3+i]*EUT[j*3+i];}
        if(theNorm0>MATH_ZERO)
        {
            detN++;
            theNorm0=sqrt(theNorm0);
            E[i*3+i]=theNorm0;
            for(int j=0;j<3;j++){U[j*3+i]=EUT[j*3+i]/theNorm0;}
        }
        else
        {
            detflag=i;
        }
    }
    if(detN<=1)
    {
        error=M_ERROR;
    }
    else if(detN==2)
    {
        int the1=(detflag+1)%3;int the2=(detflag+2)%3;
        TT aa[3]={U[the1],U[the1+3],U[the1+6]};
        TT bb[3]={U[the2],U[the2+3],U[the2+6]};
        TT temp1=-(aa[0]*bb[1]-bb[0]*aa[1]);
        TT temp2=aa[0]*bb[2]-aa[2]*bb[0];
        TT thebeta=atan2(temp1,temp2);
        TT temp3=aa[1]*cos(thebeta)+aa[2]*sin(thebeta);
        TT theerfa=atan2(-aa[0],temp3);
        U[detflag]=cos(theerfa);
        U[detflag+3]=sin(theerfa)*cos(thebeta);
        U[detflag+6]=sin(theerfa)*sin(thebeta);
    }
//    printf("the A is \n");Show(3,3,A);
//    printf("the V is \n");Show(3,3,V);
//    printf("the EUT is \n");Show(3,3,EUT);
//    printf("the U is \n");Show(3,3,U);
    return error;
}

/*************高斯数值法求解Ax=b的解（最大阶数为10）***********************************************************/
template <typename TT>
int CJW_Math<TT>::MySolveGS(int mym,TT* A, TT* b,TT* x)
{
    TT A_R[M_N_MAX*M_N_MAX]/*[mym*mym]*/={0};TT b_R[M_N_MAX]/*[mym]*/={0};
    //record the data
    for(int i=0;i<mym;i++)
    {
        b_R[i]=b[i];int mymi=mym*i;
        for(int j=0;j<mym;j++) {A_R[mymi+j]=A[mymi+j];}
    }
    //GS get the QR function
    for(int i=0;i<(mym-1);i++)
    {
        int mymi=mym*i;
        for(int j=i+1;j<mym;j++)
        {
            int mymj=mym*j;
            if(A_R[mymj+i]!=0)//can be GS
            {
                if(A_R[mymi+i]==0)//change the postion
                {
                    TT thetemp=b_R[j];
                    b_R[j]=b_R[i];b_R[i]=thetemp;
                    for(int thec=i;thec<mym;thec++)
                    {
                        thetemp=A_R[mymj+thec];
                        A_R[mymj+thec]=A_R[mymi+thec];
                        A_R[mymi+thec]=thetemp;
                    }
                }
                else//the norm GS
                {
                    TT TheK=A_R[mymi+i]/A_R[mymj+i];
                    b_R[j]=TheK*b_R[j]-b_R[i];
                    A_R[mymj+i]=0;
                    for(int k=i+1;k<mym;k++) {A_R[mymj+k]=TheK*A_R[mymj+k]-A_R[mymi+k];}
                }
            }
        }
    }
    //get the x
    for(int i=mym-1;i>=0;i--)
    {
        int mymi=mym*i;
        if(fabs(A_R[mymi+i])<MATH_ZERO){return M_ERROR;}
        TT thesum=b_R[i];
        for(int j=mym-1;j>i;j--){thesum=thesum-A_R[mymi+j]*x[j];}
        x[i]=thesum/A_R[mymi+i];
    }
    return M_RIGHT;
}
/*****************************三阶矩阵运算*****************************/
template <typename TT>
TT CJW_Math<TT>::MyMatrix3det(TT A[9])//三阶矩阵行列式
{
    TT tempdet=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    return tempdet;
}
template <typename TT>
int CJW_Math<TT>::MyMatrix3Inv(TT A[9],TT out[9])//三阶矩阵逆矩阵
{
    TT tempdet=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    if(fabs(tempdet)<MATH_ZERO) {return M_ERROR;}
    out[0]=(A[4]*A[8]-A[7]*A[5])/tempdet;
    out[1]=(A[7]*A[2]-A[1]*A[8])/tempdet;
    out[2]=(A[1]*A[5]-A[2]*A[4])/tempdet;
    out[3]=(A[5]*A[6]-A[8]*A[3])/tempdet;
    out[4]=(A[8]*A[0]-A[2]*A[6])/tempdet;
    out[5]=(A[2]*A[3]-A[0]*A[5])/tempdet;
    out[6]=(A[3]*A[7]-A[6]*A[4])/tempdet;
    out[7]=(A[6]*A[1]-A[0]*A[7])/tempdet;
    out[8]=(A[0]*A[4]-A[1]*A[3])/tempdet;
    return M_RIGHT;
}
/*****************************旋转矩阵运算****************************/
template <typename TT>
void CJW_Math<TT>::MyRCompositionR(TT R1[9], TT R2[9], TT Rsum[9])//旋转矩阵相乘R1*R2=Rsum
{
    Rsum[0] = R1[0] * R2[0] + R1[1] * R2[3] + R1[2] * R2[6];
    Rsum[1] = R1[0] * R2[1] + R1[1] * R2[4] + R1[2] * R2[7];
    Rsum[2] = R1[0] * R2[2] + R1[1] * R2[5] + R1[2] * R2[8];
    Rsum[3] = R1[3] * R2[0] + R1[4] * R2[3] + R1[5] * R2[6];
    Rsum[4] = R1[3] * R2[1] + R1[4] * R2[4] + R1[5] * R2[7];
    Rsum[5] = R1[3] * R2[2] + R1[4] * R2[5] + R1[5] * R2[8];
    Rsum[6] = R1[6] * R2[0] + R1[7] * R2[3] + R1[8] * R2[6];
    Rsum[7] = R1[6] * R2[1] + R1[7] * R2[4] + R1[8] * R2[7];
    Rsum[8] = R1[6] * R2[2] + R1[7] * R2[5] + R1[8] * R2[8];
}
template <typename TT>
void CJW_Math<TT>::MyRCompositionw(TT theR[9],TT thew[3],TT outw[3])//旋转矩阵乘矢量theR*thew=outw
{
    outw[0] = theR[0] * thew[0] + theR[1] * thew[1] + theR[2] * thew[2];
    outw[1] = theR[3] * thew[0] + theR[4] * thew[1] + theR[5] * thew[2];
    outw[2] = theR[6] * thew[0] + theR[7] * thew[1] + theR[8] * thew[2];
}
template <typename TT>
void CJW_Math<TT>::MyRCompositionwAddP(TT theR[9],TT thew[3],TT addP[3],TT outw[3])//旋转矩阵乘矢量theR*thew+addP=outw
{
    outw[0] = theR[0] * thew[0] + theR[1] * thew[1] + theR[2] * thew[2]+addP[0];
    outw[1] = theR[3] * thew[0] + theR[4] * thew[1] + theR[5] * thew[2]+addP[1];
    outw[2] = theR[6] * thew[0] + theR[7] * thew[1] + theR[8] * thew[2]+addP[2];
}
template <typename TT>
void CJW_Math<TT>::MyRInv(TT R[9],TT RT[9])//旋转矩阵转置=逆矩阵
{
    RT[0]=R[0];RT[1]=R[3];RT[2]=R[6];
    RT[3]=R[1];RT[4]=R[4];RT[5]=R[7];
    RT[6]=R[2];RT[7]=R[5];RT[8]=R[8];
}
template <typename TT>
void CJW_Math<TT>::MyRInvCompositionR(TT R1[9], TT R2[9], TT Rsum[9])//旋转矩阵转置乘法R1T*R2=Rsum
{
    Rsum[0] = R1[0] * R2[0] + R1[3] * R2[3] + R1[6] * R2[6];
    Rsum[1] = R1[0] * R2[1] + R1[3] * R2[4] + R1[6] * R2[7];
    Rsum[2] = R1[0] * R2[2] + R1[3] * R2[5] + R1[6] * R2[8];
    Rsum[3] = R1[1] * R2[0] + R1[4] * R2[3] + R1[7] * R2[6];
    Rsum[4] = R1[1] * R2[1] + R1[4] * R2[4] + R1[7] * R2[7];
    Rsum[5] = R1[1] * R2[2] + R1[4] * R2[5] + R1[7] * R2[8];
    Rsum[6] = R1[2] * R2[0] + R1[5] * R2[3] + R1[8] * R2[6];
    Rsum[7] = R1[2] * R2[1] + R1[5] * R2[4] + R1[8] * R2[7];
    Rsum[8] = R1[2] * R2[2] + R1[5] * R2[5] + R1[8] * R2[8];
}
template <typename TT>
void CJW_Math<TT>::MyRInvCompositionw(TT theR[9],TT thew[3],TT outw[3])//旋转逆解矩阵乘矢量theRInv*thew=outw
{
    outw[0] = theR[0] * thew[0] + theR[3] * thew[1] + theR[6] * thew[2];
    outw[1] = theR[1] * thew[0] + theR[4] * thew[1] + theR[7] * thew[2];
    outw[2] = theR[2] * thew[0] + theR[5] * thew[1] + theR[8] * thew[2];
}
template <typename TT>
void CJW_Math<TT>::MyRInvCompositionwAddP(TT theR[9],TT thew[3],TT addP[3],TT outw[3])//旋转矩阵乘矢量theRInv*thew+addP=outw
{
    outw[0] = theR[0] * thew[0] + theR[3] * thew[1] + theR[6] * thew[2]+addP[0];
    outw[1] = theR[1] * thew[0] + theR[4] * thew[1] + theR[7] * thew[2]+addP[1];
    outw[2] = theR[2] * thew[0] + theR[5] * thew[1] + theR[8] * thew[2]+addP[2];
}

/*****************************位姿矩阵运算****************************/
template <typename TT>
void CJW_Math<TT>::MyGCompositionG(TT G1[16], TT G2[16], TT GG[16])//位姿矩阵相乘G1*G2=GG
{
    GG[0] = G1[0] * G2[0] + G1[1] * G2[4] + G1[2] * G2[8];
    GG[1] = G1[0] * G2[1] + G1[1] * G2[5] + G1[2] * G2[9];
    GG[2] = G1[0] * G2[2] + G1[1] * G2[6] + G1[2] * G2[10];
    GG[3] = G1[0] * G2[3] + G1[1] * G2[7] + G1[2] * G2[11]+G1[3];
    GG[4] = G1[4] * G2[0] + G1[5] * G2[4] + G1[6] * G2[8];
    GG[5] = G1[4] * G2[1] + G1[5] * G2[5] + G1[6] * G2[9];
    GG[6] = G1[4] * G2[2] + G1[5] * G2[6] + G1[6] * G2[10];
    GG[7] = G1[4] * G2[3] + G1[5] * G2[7] + G1[6] * G2[11]+G1[7];
    GG[8] = G1[8] * G2[0] + G1[9] * G2[4] + G1[10]* G2[8];
    GG[9] = G1[8] * G2[1] + G1[9] * G2[5] + G1[10]* G2[9];
    GG[10]= G1[8] * G2[2] + G1[9] * G2[6] + G1[10]* G2[10];
    GG[11]= G1[8] * G2[3] + G1[9] * G2[7] + G1[10]* G2[11]+G1[11];
    GG[12]=0;GG[13]=0;GG[14]=0;GG[15]=1;
}
template <typename TT>
void CJW_Math<TT>::MyGCompositionP(TT G[16], TT P[3], TT outP[3])//位姿矩阵乘法位置G1*P=outP
{
    outP[0] = G[0] * P[0] + G[1] * P[1] + G[2] * P[2] + G[3];
    outP[1] = G[4] * P[0] + G[5] * P[1] + G[6] * P[2] + G[7];
    outP[2] = G[8] * P[0] + G[9] * P[1] + G[10]* P[2] + G[11];
}
template <typename TT>
void CJW_Math<TT>::MyGInv(TT G[16],TT InvG[16])//位姿矩阵逆解
{
    InvG[0]=G[0];InvG[1]=G[4];InvG[2]=G[8];     InvG[3]=-(G[0]*G[3]+G[4]*G[7]+G[8]*G[11]);
    InvG[4]=G[1];InvG[5]=G[5];InvG[6]=G[9];     InvG[7]=-(G[1]*G[3]+G[5]*G[7]+G[9]*G[11]);
    InvG[8]=G[2];InvG[9]=G[6];InvG[10]=G[10];   InvG[11]=-(G[2]*G[3]+G[6]*G[7]+G[10]*G[11]);
    InvG[12]=0;InvG[13]=0;InvG[14]=0;InvG[15]=1;
}
template <typename TT>
void CJW_Math<TT>::MyGInvCompositionG(TT G1[16], TT G2[16], TT GG[16])//位姿矩阵逆解乘法G1Inv*G2=GG
{
    GG[0] = G1[0] * G2[0] + G1[4] * G2[4] + G1[8] * G2[8];
    GG[1] = G1[0] * G2[1] + G1[4] * G2[5] + G1[8] * G2[9];
    GG[2] = G1[0] * G2[2] + G1[4] * G2[6] + G1[8] * G2[10];
    GG[3] = G1[0] * G2[3] + G1[4] * G2[7] + G1[8] * G2[11]-(G1[0]*G1[3]+G1[4]*G1[7]+G1[8]*G1[11]);
    GG[4] = G1[1] * G2[0] + G1[5] * G2[4] + G1[9] * G2[8];
    GG[5] = G1[1] * G2[1] + G1[5] * G2[5] + G1[9] * G2[9];
    GG[6] = G1[1] * G2[2] + G1[5] * G2[6] + G1[9] * G2[10];
    GG[7] = G1[1] * G2[3] + G1[5] * G2[7] + G1[9] * G2[11]-(G1[1]*G1[3]+G1[5]*G1[7]+G1[9]*G1[11]);
    GG[8] = G1[2] * G2[0] + G1[6] * G2[4] + G1[10]* G2[8];
    GG[9] = G1[2] * G2[1] + G1[6] * G2[5] + G1[10]* G2[9];
    GG[10]= G1[2] * G2[2] + G1[6] * G2[6] + G1[10]* G2[10];
    GG[11]= G1[2] * G2[3] + G1[6] * G2[7] + G1[10]* G2[11]-(G1[2]*G1[3]+G1[6]*G1[7]+G1[10]*G1[11]);
    GG[12]=0;GG[13]=0;GG[14]=0;GG[15]=1;
}
template <typename TT>
void CJW_Math<TT>::MyGInvCompositionP(TT G[16], TT P[3], TT outP[3])//位姿矩阵逆解乘法位置G1Inv*P=outP
{
    outP[0] = G[0] * P[0] + G[4] * P[1] + G[8] * P[2] - (G[0]*G[3]+G[4]*G[7]+G[8]*G[11]);
    outP[1] = G[1] * P[0] + G[5] * P[1] + G[9] * P[2] - (G[1]*G[3]+G[5]*G[7]+G[9]*G[11]);
    outP[2] = G[2] * P[0] + G[6] * P[1] + G[10]* P[2] - (G[2]*G[3]+G[6]*G[7]+G[10]*G[11]);
}
template <typename TT>
void CJW_Math<TT>::MyRPToG(TT R[9],TT P[3],TT G[16])//姿态+位置=位姿矩阵
{
    G[0]=R[0];  G[1]=R[1];  G[2]=R[2];  G[3]=P[0];
    G[4]=R[3];  G[5]=R[4];  G[6]=R[5];  G[7]=P[1];
    G[8]=R[6];  G[9]=R[7];  G[10]=R[8]; G[11]=P[2];
    G[12]=0;    G[13]=0;    G[14]=0;    G[15]=1;
}
template <typename TT>
void CJW_Math<TT>::MyGToRP(TT G[16],TT R[9],TT P[3])//位姿矩阵=姿态+位置
{
    R[0]=G[0];  R[1]=G[1];  R[2]=G[2];
    R[3]=G[4];  R[4]=G[5];  R[5]=G[6];
    R[6]=G[8];  R[7]=G[9];  R[8]=G[10];

    P[0]=G[3];  P[1]=G[7];  P[2]=G[11];
}
template <typename TT>
void CJW_Math<TT>::MyGToRInvP(TT G[16],TT RT[9],TT P[3])//位姿矩阵=姿态逆解+位置
{
    RT[0]=G[0];  RT[1]=G[4];  RT[2]=G[8];
    RT[3]=G[1];  RT[4]=G[5];  RT[5]=G[9];
    RT[6]=G[2];  RT[7]=G[6];  RT[8]=G[10];

    P[0]=G[3];  P[1]=G[7];  P[2]=G[11];
}

/*****************************姿态矩阵 和 其他表达的转换****************************/
template <typename TT>
void CJW_Math<TT>::MyEulerZYXToR(TT w[3], TT R[9])//欧拉角ZYX转位姿矩阵
{
    //05š®0409¡À010303000109š®
    R[0] = cos(w[1])*cos(w[0]);
    R[1] = sin(w[2])*sin(w[1])*cos(w[0]) - cos(w[2])*sin(w[0]);
    R[2] = cos(w[2])*sin(w[1])*cos(w[0]) + sin(w[2])*sin(w[0]);
    R[3] = cos(w[1])*sin(w[0]);
    R[4] = sin(w[2])*sin(w[1])*sin(w[0]) + cos(w[2])*cos(w[0]);
    R[5] = cos(w[2])*sin(w[1])*sin(w[0]) - sin(w[2])*cos(w[0]);
    R[6] = -sin(w[1]);
    R[7] = sin(w[2])*cos(w[1]);
    R[8] = cos(w[2])*cos(w[1]);
}
template <typename TT>
void CJW_Math<TT>::MyRToEulerZYX(TT R[9], TT w[3])//姿态矩阵转欧拉角ZYX
{
    //05š®0409ZYX(yaw,pitch,roll)
    if(R[6]>1)
    {
        w[1]=M_PI/2;
    }
    else if(R[6]<-1)
    {
        w[1]=-M_PI/2;
    }
    else
    {
        w[1] = -asin(R[6]);
    }
    w[2] = atan2(R[7], R[8]);
    w[0] = atan2(R[3], R[0]);
}
template <typename TT>
void CJW_Math<TT>::MyQuaternionToR(TT q[4], TT R[9]) //q为 x y z w 四元数转旋转矩阵
{
    if (q[3] >= 1)
    {
        R[0]=1;R[1]=0;R[2]=0;
        R[3]=0;R[4]=1;R[5]=0;
        R[6]=0;R[7]=0;R[8]=1;
    }
    else
    {
        R[0] = 1 - 2 * q[1] * q[1] - 2 * q[2] * q[2];
        R[1] = 2 * q[0] * q[1] - 2 * q[2] * q[3];
        R[2] = 2 * q[0] * q[2] + 2 * q[1] * q[3];
        R[3] = 2 * q[0] * q[1] + 2 * q[2] * q[3];
        R[4] = 1 - 2 * q[0] * q[0] - 2 * q[2] * q[2];
        R[5] = 2 * q[1] * q[2] - 2 * q[0] * q[3];
        R[6] = 2 * q[0] * q[2] - 2 * q[1] * q[3];
        R[7] = 2 * q[1] * q[2] + 2 * q[0] * q[3];
        R[8] = 1 - 2 * q[0] * q[0] - 2 * q[1] * q[1];
    }
}

/*****************************基于指数坐标表达****************************/
template <typename TT>
void CJW_Math<TT>::MyExponent3ToR(TT w[3], TT R[9])//指数坐标转旋转矩阵
{
    //求解变换矩阵
    TT wo = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if ((wo <EXP_ZERO) && (wo >-EXP_ZERO))
    {
        R[0] = 1; R[1] = 0; R[2] = 0;
        R[3] = 0; R[4] = 1; R[5] = 0;
        R[6] = 0; R[7] = 0; R[8] = 1;
    }
    else
    {
        TT ww[3] = { w[0] / wo,w[1] / wo,w[2] / wo };
        R[0] = 1 + (-ww[1] * ww[1] - ww[2] * ww[2])*(1 - cos(wo));
        R[1] = -ww[2] * sin(wo) + ww[0] * ww[1] * (1 - cos(wo));
        R[2] = ww[1] * sin(wo) + ww[0] * ww[2] * (1 - cos(wo));
        R[3] = ww[2] * sin(wo) + ww[0] * ww[1] * (1 - cos(wo));
        R[4] = 1 + (-ww[0] * ww[0] - ww[2] * ww[2])*(1 - cos(wo));
        R[5] = -ww[0] * sin(wo) + ww[1] * ww[2] * (1 - cos(wo));
        R[6] = -ww[1] * sin(wo) + ww[0] * ww[2] * (1 - cos(wo));
        R[7] = ww[0] * sin(wo) + ww[1] * ww[2] * (1 - cos(wo));
        R[8] = 1 + (-ww[0] * ww[0] - ww[1] * ww[1])*(1 - cos(wo));
    }
}
template <typename TT>
void CJW_Math<TT>::MyRToExponent3(TT R[9], TT w[3])//旋转矩阵=指数坐标
{
	//求解变换矩阵
	TT trR = R[0] + R[4] + R[8];
	if (trR > 3)
	{
		trR = 3;
	}
	else if (trR < -1)
	{
		trR = -1;
	}
	TT tro = acos((trR - 1) / 2);
	TT tmp_sino = sin(tro);
	if (fabs(tmp_sino) < EXP_ZERO)
	{
		//轴线绝对值
		TT wabs[3] = { sqrt((R[0] - trR) / 2),sqrt((R[4] - trR) / 2),sqrt((R[8] - trR) / 2) };
		if (wabs[0] < EXP_ZERO)
		{
			w[0] = 0;
			if (wabs[1] < EXP_ZERO)
			{
				w[1] = 0;
				w[2] = wabs[2] * tro;
			}
			else
			{
				w[1] = wabs[1] * tro;
				w[2] = R[5] / w[1] / 2 * tro;
			}
		}
		else
		{
			w[0] = wabs[0] * tro;
			if (wabs[1] < EXP_ZERO)
			{
				w[1] = 0;
				w[2] = R[2] / wabs[0] / 2 * tro;
			}
			else
			{
				w[1] = R[1] / 2 / wabs[0] * tro;
				w[2] = R[5] / w[1] / 2 * tro;
			}
		}
	}
	else
	{
		w[0] = (R[7] - R[5])*tro / 2 / tmp_sino;
		w[1] = (R[2] - R[6])*tro / 2 / tmp_sino;
		w[2] = (R[3] - R[1])*tro / 2 / tmp_sino;
	}
}
/*template <typename TT>
void CJW_Math<TT>::MyRToExponent3(TT R[9], TT w[3])//旋转矩阵=指数坐标
{
	//求解变换矩阵
	TT trR = R[0]+R[4]+R[8];
	if(trR>3)
	{
		trR=3;
	}
	else if(trR<-1)
	{
		trR=-1;
	}
	TT tro = acos((trR - 1) / 2);
	TT tmp_sino = sin(tro);
	if (fabs(tmp_sino) < EXP_ZERO)
	{
		if(tro<M_PI/2)
		{
			w[0] = 0; w[1] = 0; w[2] = 0;
		}
		else
		{
			w[0] = sqrt((R[0]-trR)/2)*tro;
			w[1] = sqrt((R[4]-trR)/2)*tro;
			w[2] = sqrt((R[8]-trR)/2)*tro;
		}
	}
	else
	{
		w[0] = (R[7] - R[5])*tro / 2 / tmp_sino;
		w[1] = (R[2] - R[6])*tro / 2 / tmp_sino;
		w[2] = (R[3] - R[1])*tro / 2 / tmp_sino;
	}
}*/
template <typename TT>
void CJW_Math<TT>::MyExponent3ToQuaternion(TT w[3], TT q[4])//指数坐标转四元数xyzw
{
    TT a = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if (a < 0.001)
    {
        q[0] = 0; q[1] = 0; q[2] = 0; q[3] = 1;
    }
    else
    {
        q[0] = w[0] * sin(a / 2) / a; q[1] = w[1] * sin(a / 2) / a; q[2] = w[2] * sin(a / 2) / a;
        q[3] = cos(a / 2);
    }
}
template <typename TT>
void CJW_Math<TT>::MyQuaternionToExponent3(TT q[4],TT w[3])//四元数xyzw转指数坐标
{
    if (q[3]>=1)
    {
        w[0] = 0; w[1] = 0; w[2] = 0;
    }
    else
    {
        TT a = acos(q[3]);
        w[0] = q[0] * 2 * a / sin(a); w[1] = q[1] * 2 * a / sin(a); w[2] = q[2] * 2 * a / sin(a);
    }
}
template <typename TT>
void CJW_Math<TT>::MyExponent4ToG(TT w[6], TT R[16])//指数坐标转位姿矩阵
{
    //求解变换矩阵
    TT wo = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if ((wo <EXP_ZERO) && (wo >-EXP_ZERO))
    {
        R[0] = 1; R[1] = 0; R[2] = 0;      R[3] = w[3];
        R[4] = 0; R[5] = 1; R[6] = 0;      R[7] = w[4];
        R[8] = 0; R[9] = 0; R[10] = 1;    R[11] = w[5];
    }
    else
    {
        TT ww[3] = { w[0] / wo,w[1] / wo,w[2] / wo};
        TT vv[3]={ w[3] / wo ,w[4] / wo,w[5] / wo};
        R[0] = 1 + (-ww[1] * ww[1] - ww[2] * ww[2])*(1 - cos(wo));
        R[1] = -ww[2] * sin(wo) + ww[0] * ww[1] * (1 - cos(wo));
        R[2] = ww[1] * sin(wo) + ww[0] * ww[2] * (1 - cos(wo));
        R[4] = ww[2] * sin(wo) + ww[0] * ww[1] * (1 - cos(wo));
        R[5] = 1 + (-ww[0] * ww[0] - ww[2] * ww[2])*(1 - cos(wo));
        R[6] = -ww[0] * sin(wo) + ww[1] * ww[2] * (1 - cos(wo));
        R[8] = -ww[1] * sin(wo) + ww[0] * ww[2] * (1 - cos(wo));
        R[9] = ww[0] * sin(wo) + ww[1] * ww[2] * (1 - cos(wo));
        R[10] = 1 + (-ww[0] * ww[0] - ww[1] * ww[1])*(1 - cos(wo));
        TT dotwv=MyVector3Dot(ww,vv);
        TT crosswv[3];MyVector3Cross(ww,vv,crosswv);
        R[3]=(1-R[0])*crosswv[0]+(-R[1])*crosswv[1]+(-R[2])*crosswv[2]+dotwv*w[0];
        R[7]=(-R[4])*crosswv[0]+(1-R[5])*crosswv[1]+(-R[6])*crosswv[2]+dotwv*w[1];
        R[11]=(-R[8])*crosswv[0]+(-R[9])*crosswv[1]+(1-R[10])*crosswv[2]+dotwv*w[2];
    }
    R[12] = 0; R[13] = 0; R[14] = 0; R[15] = 1;
}
template <typename TT>
void CJW_Math<TT>::MyGToExponent4(TT R[16], TT w[6])//位姿矩阵转指数坐标
{
	//求解变换矩阵
	TT trR = R[0] + R[5] + R[10];
	if (trR > 3)
	{
		trR = 3;
	}
	else if (trR < -1)
	{
		trR = -1;
	}
	TT tro = acos((trR - 1) / 2);
	if ((tro < EXP_ZERO) && (tro > -EXP_ZERO))
	{
		w[0] = 0; w[1] = 0; w[2] = 0;
		w[3] = R[3]; w[4] = R[7]; w[5] = R[11];
	}
	else
	{
		TT tmp_sino = sin(tro);
		TT ww[3];
		if (fabs(tmp_sino) < EXP_ZERO)
		{
			// 轴线绝对值
			TT wabs[3] = { sqrt((R[0] - trR) / 2), sqrt((R[5] - trR) / 2), sqrt((R[10] - trR) / 2) };
			if (wabs[0] < EXP_ZERO)
			{
				ww[0] = 0;
				if (wabs[1] < EXP_ZERO)
				{
					ww[1] = 0;
					ww[2] = wabs[2];
				}
				else
				{
					ww[1] = wabs[1];
					ww[2] = R[6] / ww[1] / 2;
				}
			}
			else
			{
				ww[0] = wabs[0];
				if (wabs[1] < EXP_ZERO)
				{
					ww[1] = 0;
					ww[2] = R[2] / wabs[0] / 2;
				}
				else
				{
					ww[1] = R[1] / 2 / wabs[0];
					ww[2] = R[6] / ww[1] / 2;
				}
			}
		}
		else
		{
			ww[0] = (R[9] - R[6]) / 2 / tmp_sino;
			ww[1] = (R[2] - R[8]) / 2 / tmp_sino;
			ww[2] = (R[4] - R[1]) / 2 / tmp_sino;
		}
		TT tempM[9];
		tempM[0] = -R[1] * ww[2] + R[2] * ww[1] + tro * ww[0] * ww[0];
		tempM[3] = (1 - R[5])*ww[2] + R[6] * ww[1] + tro * ww[1] * ww[0];
		tempM[6] = -R[9] * ww[2] + (R[10] - 1)*ww[1] + tro * ww[2] * ww[0];
		tempM[1] = (R[0] - 1)*ww[2] - R[2] * ww[0] + tro * ww[0] * ww[1];
		tempM[4] = R[4] * ww[2] - R[6] * ww[0] + tro * ww[1] * ww[1];
		tempM[7] = R[8] * ww[2] + (1 - R[10])*ww[0] + tro * ww[2] * ww[1];
		tempM[2] = (1 - R[0])*ww[1] + R[1] * ww[0] + tro * ww[0] * ww[2];
		tempM[5] = -R[4] * ww[1] + (R[5] - 1)*ww[0] + tro * ww[1] * ww[2];
		tempM[8] = -R[8] * ww[1] + R[9] * ww[0] + tro * ww[2] * ww[2];
		TT tempMInv[9];
		MyMatrix3Inv(tempM, tempMInv);
		TT vv[3];
		vv[0] = tempMInv[0] * R[3] + tempMInv[1] * R[7] + tempMInv[2] * R[11];
		vv[1] = tempMInv[3] * R[3] + tempMInv[4] * R[7] + tempMInv[5] * R[11];
		vv[2] = tempMInv[6] * R[3] + tempMInv[7] * R[7] + tempMInv[8] * R[11];
		/////////
		w[0] = ww[0] * tro; w[1] = ww[1] * tro; w[2] = ww[2] * tro;
		w[3] = vv[0] * tro; w[4] = vv[1] * tro; w[5] = vv[2] * tro;
	}
}

/****************************************伴随运算***************************************/
template <typename TT>
void CJW_Math<TT>::MyEXPAdg(TT G[16],TT Adg[36])//伴随运算
{
    Adg[0]=G[0];  Adg[1]=G[1];  Adg[2]=G[2];        Adg[3]=0;Adg[4]=0;Adg[5]=0;
    Adg[6]=G[4];  Adg[7]=G[5];  Adg[8]=G[6];        Adg[9]=0;Adg[10]=0;Adg[11]=0;
    Adg[12]=G[8];Adg[13]=G[9];Adg[14]=G[10];      Adg[15]=0;Adg[16]=0;Adg[17]=0;

    Adg[18]=-G[11]*G[4]+G[7]*G[8];Adg[19]=-G[11]*G[5]+G[7]*G[9]; Adg[20]=-G[11]*G[6]+G[7]*G[10];      Adg[21]=Adg[0];Adg[22]=Adg[1];Adg[23]=Adg[2];
    Adg[24]=G[11]*G[0]-G[3]*G[8]; Adg[25]=G[11]*G[1]-G[3]*G[9];  Adg[26]=G[11]*G[2]-G[3]*G[10];      Adg[27]=Adg[6];Adg[28]=Adg[7];Adg[29]=Adg[8];
    Adg[30]=-G[7]*G[0]+G[3]*G[4]; Adg[31]=-G[7]*G[1]+G[3]*G[5];  Adg[32]=-G[7]*G[2]+G[3]*G[6];      Adg[33]=Adg[12];Adg[34]=Adg[13];Adg[35]=Adg[14];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgScrew(TT G[16],TT S[6],TT AdgV[6])//伴随运算乘旋量
{
    AdgV[0]=G[0]*S[0]+G[1]*S[1]+G[2]*S[2];
    AdgV[1]=G[4]*S[0]+G[5]*S[1]+G[6]*S[2];
    AdgV[2]=G[8]*S[0]+G[9]*S[1]+G[10]*S[2];

    AdgV[3]=(-G[11]*G[4]+G[7]*G[8])*S[0]+(-G[11]*G[5]+G[7]*G[9])*S[1]+(-G[11]*G[6]+G[7]*G[10])*S[2]+G[0]*S[3]+G[1]*S[4]+G[2]*S[5];
    AdgV[4]=(G[11]*G[0]-G[3]*G[8])*S[0] +(G[11]*G[1]-G[3]*G[9])*S[1] +(G[11]*G[2]-G[3]*G[10])*S[2] +G[4]*S[3]+G[5]*S[4]+G[6]*S[5];
    AdgV[5]=(-G[7]*G[0]+G[3]*G[4])*S[0] +(-G[7]*G[1]+G[3]*G[5])*S[1] +(-G[7]*G[2]+G[3]*G[6])*S[2]  +G[8]*S[3]+G[9]*S[4]+G[10]*S[5];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgT(TT G[16],TT AdgT[36])//伴随转置
{
    AdgT[0]=G[0];  AdgT[1]=G[4];  AdgT[2]=G[8];        AdgT[3]=-G[11]*G[4]+G[7]*G[8];  AdgT[4]=G[11]*G[0]-G[3]*G[8];  AdgT[5]=-G[7]*G[0]+G[3]*G[4];
    AdgT[6]=G[1];  AdgT[7]=G[5];  AdgT[8]=G[9];        AdgT[9]=-G[11]*G[5]+G[7]*G[9];  AdgT[10]=G[11]*G[1]-G[3]*G[9]; AdgT[11]=-G[7]*G[1]+G[3]*G[5];
    AdgT[12]=G[2]; AdgT[13]=G[6]; AdgT[14]=G[10];      AdgT[15]=-G[11]*G[6]+G[7]*G[10];AdgT[16]=G[11]*G[2]-G[3]*G[10];AdgT[17]=-G[7]*G[2]+G[3]*G[6];

    AdgT[18]=0; AdgT[19]=0; AdgT[20]=0;      AdgT[21]=AdgT[0];AdgT[22]=AdgT[1];AdgT[23]=AdgT[2];
    AdgT[24]=0; AdgT[25]=0; AdgT[26]=0;      AdgT[27]=AdgT[6];AdgT[28]=AdgT[7];AdgT[29]=AdgT[8];
    AdgT[30]=0;AdgT[31]=0;AdgT[32]=0;      AdgT[33]=AdgT[12];AdgT[34]=AdgT[13];AdgT[35]=AdgT[14];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgTScrew(TT G[16],TT S[6],TT AdgTV[6])//伴随转置乘旋量
{
    AdgTV[0]=G[0]*S[0]+G[4]*S[1]+G[8]*S[2] +(-G[11]*G[4]+G[7]*G[8])*S[3] +(G[11]*G[0]-G[3]*G[8])*S[4] +(-G[7]*G[0]+G[3]*G[4])*S[5];
    AdgTV[1]=G[1]*S[0]+G[5]*S[1]+G[9]*S[2] +(-G[11]*G[5]+G[7]*G[9])*S[3] +(G[11]*G[1]-G[3]*G[9])*S[4] +(-G[7]*G[1]+G[3]*G[5])*S[5];
    AdgTV[2]=G[2]*S[0]+G[6]*S[1]+G[10]*S[2]+(-G[11]*G[6]+G[7]*G[10])*S[3]+(G[11]*G[2]-G[3]*G[10])*S[4]+(-G[7]*G[2]+G[3]*G[6])*S[5];

    AdgTV[3]=G[0]*S[3]+G[4]*S[4]+G[8]*S[5];
    AdgTV[4]=G[1]*S[3]+G[5]*S[4]+G[9]*S[5];
    AdgTV[5]=G[2]*S[3]+G[6]*S[4]+G[10]*S[5];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgInv(TT G[16],TT AdgInv[36])//伴随逆解
{
    AdgInv[0]=G[0];  AdgInv[1]=G[4];  AdgInv[2]=G[8];        AdgInv[3]=0;AdgInv[4]=0;AdgInv[5]=0;
    AdgInv[6]=G[1];  AdgInv[7]=G[5];  AdgInv[8]=G[9];        AdgInv[9]=0;AdgInv[10]=0;AdgInv[11]=0;
    AdgInv[12]=G[2]; AdgInv[13]=G[6]; AdgInv[14]=G[10];      AdgInv[15]=0;AdgInv[16]=0;AdgInv[17]=0;

    AdgInv[18]=-G[11]*G[4]+G[7]*G[8]; AdgInv[19]=G[11]*G[0]-G[3]*G[8]; AdgInv[20]=-G[7]*G[0]+G[3]*G[4];      AdgInv[21]=AdgInv[0];AdgInv[22]=AdgInv[1];AdgInv[23]=AdgInv[2];
    AdgInv[24]=-G[11]*G[5]+G[7]*G[9]; AdgInv[25]=G[11]*G[1]-G[3]*G[9]; AdgInv[26]=-G[7]*G[1]+G[3]*G[5];      AdgInv[27]=AdgInv[6];AdgInv[28]=AdgInv[7];AdgInv[29]=AdgInv[8];
    AdgInv[30]=-G[11]*G[6]+G[7]*G[10];AdgInv[31]=G[11]*G[2]-G[3]*G[10];AdgInv[32]=-G[7]*G[2]+G[3]*G[6];      AdgInv[33]=AdgInv[12];AdgInv[34]=AdgInv[13];AdgInv[35]=AdgInv[14];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgInvScrew(TT G[16],TT S[6],TT AdgInvV[6])//伴随逆解乘旋量
{
    AdgInvV[0]=G[0]*S[0]+G[4]*S[1]+G[8]*S[2];
    AdgInvV[1]=G[1]*S[0]+G[5]*S[1]+G[9]*S[2];
    AdgInvV[2]=G[2]*S[0]+G[6]*S[1]+G[10]*S[2];

    AdgInvV[3]=(-G[11]*G[4]+G[7]*G[8])*S[0] +(G[11]*G[0]-G[3]*G[8])*S[1] +(-G[7]*G[0]+G[3]*G[4])*S[2] +G[0]*S[3]+G[4]*S[4]+G[8]*S[5];
    AdgInvV[4]=(-G[11]*G[5]+G[7]*G[9])*S[0] +(G[11]*G[1]-G[3]*G[9])*S[1] +(-G[7]*G[1]+G[3]*G[5])*S[2] +G[1]*S[3]+G[5]*S[4]+G[9]*S[5];
    AdgInvV[5]=(-G[11]*G[6]+G[7]*G[10])*S[0]+(G[11]*G[2]-G[3]*G[10])*S[1]+(-G[7]*G[2]+G[3]*G[6])*S[2] +G[2]*S[3]+G[6]*S[4]+G[10]*S[5];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgInvT(TT G[16],TT AdgInvT[36])//伴随逆解转置
{
    AdgInvT[0]=G[0];  AdgInvT[1]=G[1];  AdgInvT[2]=G[2];        AdgInvT[3]=-G[11]*G[4]+G[7]*G[8];AdgInvT[4]=-G[11]*G[5]+G[7]*G[9];AdgInvT[5]=-G[11]*G[6]+G[7]*G[10];
    AdgInvT[6]=G[4];  AdgInvT[7]=G[5];  AdgInvT[8]=G[6];        AdgInvT[9]=G[11]*G[0]-G[3]*G[8]; AdgInvT[10]=G[11]*G[1]-G[3]*G[9];AdgInvT[11]=G[11]*G[2]-G[3]*G[10];
    AdgInvT[12]=G[8]; AdgInvT[13]=G[9]; AdgInvT[14]=G[10];      AdgInvT[15]=-G[7]*G[0]+G[3]*G[4];AdgInvT[16]=-G[7]*G[1]+G[3]*G[5];AdgInvT[17]=-G[7]*G[2]+G[3]*G[6];

    AdgInvT[18]=0;AdgInvT[19]=0;AdgInvT[20]=0;      AdgInvT[21]=AdgInvT[0];AdgInvT[22]=AdgInvT[1];AdgInvT[23]=AdgInvT[2];
    AdgInvT[24]=0;AdgInvT[25]=0;AdgInvT[26]=0;      AdgInvT[27]=AdgInvT[6];AdgInvT[28]=AdgInvT[7];AdgInvT[29]=AdgInvT[8];
    AdgInvT[30]=0;AdgInvT[31]=0;AdgInvT[32]=0;      AdgInvT[33]=AdgInvT[12];AdgInvT[34]=AdgInvT[13];AdgInvT[35]=AdgInvT[14];
}
template <typename TT>
void CJW_Math<TT>::MyEXPAdgInvTScrew(TT G[16],TT S[6],TT AdgInvTV[6])//伴随逆解转置乘旋量
{
    AdgInvTV[0]=G[0]*S[0]+G[1]*S[1]+G[2]*S[2] +(-G[11]*G[4]+G[7]*G[8])*S[3]+(-G[11]*G[5]+G[7]*G[9])*S[4]+(-G[11]*G[6]+G[7]*G[10])*S[5];
    AdgInvTV[1]=G[4]*S[0]+G[5]*S[1]+G[6]*S[2] +(G[11]*G[0]-G[3]*G[8])*S[3] +(G[11]*G[1]-G[3]*G[9])*S[4] +(G[11]*G[2]-G[3]*G[10])*S[5];
    AdgInvTV[2]=G[8]*S[0]+G[9]*S[1]+G[10]*S[2]+(-G[7]*G[0]+G[3]*G[4])*S[3] +(-G[7]*G[1]+G[3]*G[5])*S[4] +(-G[7]*G[2]+G[3]*G[6])*S[5];

    AdgInvTV[3]=G[0]*S[3]+G[1]*S[4]+G[2]*S[5];
    AdgInvTV[4]=G[4]*S[3]+G[5]*S[4]+G[6]*S[5];
    AdgInvTV[5]=G[8]*S[3]+G[9]*S[4]+G[10]*S[5];
}
template <typename TT>
void CJW_Math<TT>::MyEXPad(TT w[6],TT ad[36])//旋量伴随运算
{
    ad[0]=0;ad[1]=-w[2];ad[2]=w[1];               ad[3]=0;ad[4]=0;ad[5]=0;
    ad[6]=w[2];ad[7]=0;ad[8]=-w[0];               ad[9]=0;ad[10]=0;ad[11]=0;
    ad[12]=-w[1];ad[13]=w[0];ad[14]=0;        ad[15]=0;ad[16]=0;ad[17]=0;

    ad[18]=0;ad[19]=-w[5];ad[20]=w[4];               ad[21]=0;ad[22]=-w[2];ad[23]=w[1];
    ad[24]=w[5];ad[25]=0;ad[26]=-w[3];               ad[27]=w[2];ad[28]=0;ad[29]=-w[0];
    ad[30]=-w[4];ad[31]=w[3];ad[32]=0;               ad[33]=-w[1];ad[34]=w[0];ad[35]=0;
}
template <typename TT>
void CJW_Math<TT>::MyEXPadT(TT w[6],TT ad[36])//旋量伴随转置
{
    ad[0]=0;    ad[6]=-w[2];ad[12]=w[1];                ad[18]=0;ad[24]=0;   ad[30]=0;
    ad[1]=w[2]; ad[7]=0;    ad[13]=-w[0];               ad[19]=0;ad[25]=0;  ad[31]=0;
    ad[2]=-w[1];ad[8]=w[0]; ad[14]=0;                   ad[20]=0;ad[26]=0;  ad[32]=0;

    ad[3]=0;    ad[9]=-w[5];ad[15]=w[4];                ad[21]=0;       ad[27]=-w[2];   ad[33]=w[1];
    ad[4]=w[5]; ad[10]=0;   ad[16]=-w[3];               ad[22]=w[2];    ad[28]=0;       ad[34]=-w[0];
    ad[5]=-w[4];ad[11]=w[3];ad[17]=0;                   ad[23]=-w[1];   ad[29]=w[0];    ad[35]=0;
}
/********************指数坐标的微分运算****************************************************/
template <typename TT>
void CJW_Math<TT>::MydEXP3(TT w[3],TT dexp[9])//指数坐标一阶导数
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        dexp[0]=1;dexp[1]=0;dexp[2]=0;
        dexp[3]=0;dexp[4]=1;dexp[5]=0;
        dexp[6]=0;dexp[7]=0;dexp[8]=1;
    }
    else
    {
        TT aerf=sin(th)/th;TT beta=2*(1-cos(th))/th/th;
        //dexp=I+beta*dis_w/2+(1-aerf)/th/th*dis_w*dis_w;
        TT A1=beta/2;TT A2=(1-aerf)/th/th;
        dexp[0]=1           +A2*(-w[1]*w[1]-w[2]*w[2]);
        dexp[1]=0+A1*(-w[2])+A2*(w[0]*w[1]);
        dexp[2]=0+A1*( w[1])+A2*(w[0]*w[2]);
        dexp[3]=0+A1*( w[2])+A2*(w[1]*w[0]);
        dexp[4]=1           +A2*(-w[0]*w[0]-w[2]*w[2]);
        dexp[5]=0+A1*(-w[0])+A2*(w[1]*w[2]);
        dexp[6]=0+A1*(-w[1])+A2*(w[0]*w[2]);
        dexp[7]=0+A1*( w[0])+A2*(w[1]*w[2]);
        dexp[8]=1           +A2*(-w[0]*w[0]-w[1]*w[1]);
    }

}

template <typename TT>
void CJW_Math<TT>::MydEXP3Inv(TT w[3],TT dexp[9])//指数坐标一阶导数逆矩阵
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        dexp[0]=1;dexp[1]=0;dexp[2]=0;
        dexp[3]=0;dexp[4]=1;dexp[5]=0;
        dexp[6]=0;dexp[7]=0;dexp[8]=1;
    }
    else
    {
        TT gama=th*cos(th/2)/sin(th/2)/2;
        //dexp=I-dis_w/2+(1-gama)/th/th*dis_w*dis_w;
        TT A1=-1.0/2;TT A2=(1-gama)/th/th;
        dexp[0]=1           +A2*(-w[1]*w[1]-w[2]*w[2]);
        dexp[1]=0+A1*(-w[2])+A2*(w[0]*w[1]);
        dexp[2]=0+A1*( w[1])+A2*(w[0]*w[2]);
        dexp[3]=0+A1*( w[2])+A2*(w[1]*w[0]);
        dexp[4]=1           +A2*(-w[0]*w[0]-w[2]*w[2]);
        dexp[5]=0+A1*(-w[0])+A2*(w[1]*w[2]);
        dexp[6]=0+A1*(-w[1])+A2*(w[0]*w[2]);
        dexp[7]=0+A1*( w[0])+A2*(w[1]*w[2]);
        dexp[8]=1           +A2*(-w[0]*w[0]-w[1]*w[1]);
    }
}
template <typename TT>
void CJW_Math<TT>::MyddEXP3(TT w[3],TT dw[3],TT ddexp[9])//指数坐标二阶导数
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        ///
        ddexp[0]=0; ddexp[1]=-dw[2]/2; ddexp[2]=dw[1]/2;
        ddexp[3]=dw[2]/2; ddexp[4]=0; ddexp[5]=-dw[0]/2;
        ddexp[6]=-dw[1]/2;ddexp[7]=dw[0]/2;ddexp[8]=0;
    }
    else
    {
        TT aerf=sin(th)/th;TT beta=2*(1-cos(th))/th/th;
        ///////
        TT Twdw=w[0]*dw[0]+w[1]*dw[1]+w[2]*dw[2];
        TT C0=beta/2;TT C1=(1.0-aerf)/th/th;
        TT C2=((aerf-beta)/th/th)*Twdw; TT C3=((beta/2-3.0/th/th*(1-aerf))/th/th)*Twdw;
        ddexp[0] =          +C1*(-w[2]*dw[2]-w[1]*dw[1]-dw[2]*w[2]-dw[1]*w[1])             +C3*(-w[1]*w[1]-w[2]*w[2]);
        ddexp[1] =C0*(-dw[2])+C1*(w[1]*dw[0]+dw[1]*w[0])                      +C2*(-w[2]) +C3*(w[0]*w[1]);
        ddexp[2] =C0*( dw[1])+C1*(w[2]*dw[0]+dw[2]*w[0])                      +C2*( w[1]) +C3*(w[0]*w[2]);
        ddexp[3] =C0*( dw[2])+C1*(w[0]*dw[1]+dw[0]*w[1])                      +C2*( w[2]) +C3*(w[1]*w[0]);
        ddexp[4] =          +C1*(-w[2]*dw[2]-w[0]*dw[0]-dw[2]*w[2]-dw[0]*w[0])             +C3*(-w[0]*w[0]-w[2]*w[2]);
        ddexp[5] =C0*(-dw[0])+C1*(w[2]*dw[1]+dw[2]*w[1])                      +C2*(-w[0]) +C3*(w[1]*w[2]);
        ddexp[6]=C0*(-dw[1])+C1*(w[0]*dw[2]+dw[0]*w[2])                      +C2*(-w[1]) +C3*(w[0]*w[2]);
        ddexp[7]=C0*( dw[0])+C1*(w[1]*dw[2]+dw[1]*w[2])                      +C2*( w[0]) +C3*(w[1]*w[2]);
        ddexp[8]=          +C1*(-w[1]*dw[1]-w[0]*dw[0]-dw[1]*w[1]-w[0]*dw[0])             +C3*(-w[0]*w[0]-w[1]*w[1]);
    }
}
template <typename TT>
void CJW_Math<TT>::MyddEXP3Inv(TT w[3],TT dw[3],TT ddexp[9])//指数坐标二阶导数逆矩阵
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        ///
        ddexp[0]=0; ddexp[1]=dw[2]/2; ddexp[2]=-dw[1]/2;
        ddexp[3]=-dw[2]/2; ddexp[4]=0; ddexp[5]=dw[0]/2;
        ddexp[6]=dw[1]/2;ddexp[7]=-dw[0]/2;ddexp[8]=0;
    }
    else
    {
        TT aerf=sin(th)/th;TT beta=2*(1-cos(th))/th/th;TT gama=th*cos(th/2)/sin(th/2)/2;
        ///////
        TT Twdw=w[0]*dw[0]+w[1]*dw[1]+w[2]*dw[2];
        TT C0=-1.0/2;TT C1=(1.0-gama)/th/th;//TT C2=0;
        TT C3=((1/beta+gama-2)/th/th/th/th)*Twdw;
        ddexp[0] =          +C1*(-w[2]*dw[2]-w[1]*dw[1]-dw[2]*w[2]-dw[1]*w[1])+C3*(-w[1]*w[1]-w[2]*w[2]);
        ddexp[1] =C0*(-dw[2])+C1*(w[1]*dw[0]+dw[1]*w[0])                      +C3*(w[0]*w[1]);
        ddexp[2] =C0*( dw[1])+C1*(w[2]*dw[0]+dw[2]*w[0])                      +C3*(w[0]*w[2]);
        ddexp[3] =C0*( dw[2])+C1*(w[0]*dw[1]+dw[0]*w[1])                      +C3*(w[1]*w[0]);
        ddexp[4] =          +C1*(-w[2]*dw[2]-w[0]*dw[0]-dw[2]*w[2]-dw[0]*w[0])+C3*(-w[0]*w[0]-w[2]*w[2]);
        ddexp[5] =C0*(-dw[0])+C1*(w[2]*dw[1]+dw[2]*w[1])                      +C3*(w[1]*w[2]);
        ddexp[6]=C0*(-dw[1])+C1*(w[0]*dw[2]+dw[0]*w[2])                      +C3*(w[0]*w[2]);
        ddexp[7]=C0*( dw[0])+C1*(w[1]*dw[2]+dw[1]*w[2])                      +C3*(w[1]*w[2]);
        ddexp[8]=          +C1*(-w[1]*dw[1]-w[0]*dw[0]-dw[1]*w[1]-w[0]*dw[0])+C3*(-w[0]*w[0]-w[1]*w[1]);
    }
}
template <typename TT>
void CJW_Math<TT>::MydEXP4(TT w[6],TT dexp[36])//指数坐标一阶导数
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        dexp[0]=1; dexp[1]=0; dexp[2]=0;                dexp[3]=0; dexp[4]=0; dexp[5]=0;
        dexp[6]=0; dexp[7]=1; dexp[8]=0;                dexp[9]=0; dexp[10]=0;dexp[11]=0;
        dexp[12]=0;dexp[13]=0;dexp[14]=1;               dexp[15]=0;dexp[16]=0;dexp[17]=0;

        dexp[18]=0;dexp[19]=-w[2]/2;dexp[20]=w[1]/2;    dexp[21]=1;dexp[22]=0;dexp[23]=0;
        dexp[24]=w[2]/2;dexp[25]=0;dexp[26]=-w[0]/2;    dexp[27]=0;dexp[28]=1;dexp[29]=0;
        dexp[30]=-w[1]/2;dexp[31]=w[0]/2;dexp[32]=0;    dexp[33]=0;dexp[34]=0;dexp[35]=1;

    }
    else
    {
        TT aerf=sin(th)/th;TT beta=2*(1-cos(th))/th/th;
        //dexp=I+beta*dis_w/2+(1-aerf)/th/th*dis_w*dis_w;
        TT A1=beta/2;TT A2=(1-aerf)/th/th;
        dexp[0] =1           +A2*(-w[1]*w[1]-w[2]*w[2]);
        dexp[1] =0+A1*(-w[2])+A2*(w[0]*w[1]);
        dexp[2] =0+A1*( w[1])+A2*(w[0]*w[2]);
        dexp[6] =0+A1*( w[2])+A2*(w[1]*w[0]);
        dexp[7] =1           +A2*(-w[0]*w[0]-w[2]*w[2]);
        dexp[8] =0+A1*(-w[0])+A2*(w[1]*w[2]);
        dexp[12]=0+A1*(-w[1])+A2*(w[0]*w[2]);
        dexp[13]=0+A1*( w[0])+A2*(w[1]*w[2]);
        dexp[14]=1           +A2*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        dexp[3]=0; dexp[4]=0; dexp[5]=0;
        dexp[9]=0; dexp[10]=0;dexp[11]=0;
        dexp[15]=0;dexp[16]=0;dexp[17]=0;
        ///////
        TT Twv=w[0]*w[3]+w[1]*w[4]+w[2]*w[5];
        TT C0=beta/2;TT C1=(1.0-aerf)/th/th;
        TT C2=((aerf-beta)/th/th)*Twv; TT C3=((beta/2-3.0/th/th*(1-aerf))/th/th)*Twv;
        dexp[18]=          +C1*(-w[2]*w[5]-w[1]*w[4]-w[5]*w[2]-w[4]*w[1])             +C3*(-w[1]*w[1]-w[2]*w[2]);
        dexp[19]=C0*(-w[5])+C1*(w[1]*w[3]+w[4]*w[0])                      +C2*(-w[2]) +C3*(w[0]*w[1]);
        dexp[20]=C0*( w[4])+C1*(w[2]*w[3]+w[5]*w[0])                      +C2*( w[1]) +C3*(w[0]*w[2]);
        dexp[24]=C0*( w[5])+C1*(w[0]*w[4]+w[3]*w[1])                      +C2*( w[2]) +C3*(w[1]*w[0]);
        dexp[25]=          +C1*(-w[2]*w[5]-w[0]*w[3]-w[5]*w[2]-w[3]*w[0])             +C3*(-w[0]*w[0]-w[2]*w[2]);
        dexp[26]=C0*(-w[3])+C1*(w[2]*w[4]+w[5]*w[1])                      +C2*(-w[0]) +C3*(w[1]*w[2]);
        dexp[30]=C0*(-w[4])+C1*(w[0]*w[5]+w[3]*w[2])                      +C2*(-w[1]) +C3*(w[0]*w[2]);
        dexp[31]=C0*( w[3])+C1*(w[1]*w[5]+w[4]*w[2])                      +C2*( w[0]) +C3*(w[1]*w[2]);
        dexp[32]=          +C1*(-w[1]*w[4]-w[0]*w[3]-w[4]*w[1]-w[0]*w[3])             +C3*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        dexp[21]=dexp[0]; dexp[22]=dexp[1]; dexp[23]=dexp[2];
        dexp[27]=dexp[6]; dexp[28]=dexp[7]; dexp[29]=dexp[8];
        dexp[33]=dexp[12];dexp[34]=dexp[13];dexp[35]=dexp[14];
    }

}

template <typename TT>
void CJW_Math<TT>::MydEXP4Inv(TT w[6],TT dexp[36])//指数坐标一阶导数逆矩阵
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        dexp[0]=1; dexp[1]=0; dexp[2]=0;                dexp[3]=0; dexp[4]=0; dexp[5]=0;
        dexp[6]=0; dexp[7]=1; dexp[8]=0;                dexp[9]=0; dexp[10]=0;dexp[11]=0;
        dexp[12]=0;dexp[13]=0;dexp[14]=1;               dexp[15]=0;dexp[16]=0;dexp[17]=0;

        dexp[18]=0;dexp[19]=w[2]/2;dexp[20]=-w[1]/2;    dexp[21]=1;dexp[22]=0;dexp[23]=0;
        dexp[24]=-w[2]/2;dexp[25]=0;dexp[26]=w[0]/2;    dexp[27]=0;dexp[28]=1;dexp[29]=0;
        dexp[30]=w[1]/2;dexp[31]=-w[0]/2;dexp[32]=0;    dexp[33]=0;dexp[34]=0;dexp[35]=1;
    }
    else
    {
        TT gama=th*cos(th/2)/sin(th/2)/2;TT beta=2*(1-cos(th))/th/th;
        //dexp=I-dis_w/2+(1-gama)/th/th*dis_w*dis_w;
        TT A1=-1.0/2;TT A2=(1-gama)/th/th;
        dexp[0] =1           +A2*(-w[1]*w[1]-w[2]*w[2]);
        dexp[1] =0+A1*(-w[2])+A2*(w[0]*w[1]);
        dexp[2] =0+A1*( w[1])+A2*(w[0]*w[2]);
        dexp[6] =0+A1*( w[2])+A2*(w[1]*w[0]);
        dexp[7] =1           +A2*(-w[0]*w[0]-w[2]*w[2]);
        dexp[8] =0+A1*(-w[0])+A2*(w[1]*w[2]);
        dexp[12]=0+A1*(-w[1])+A2*(w[0]*w[2]);
        dexp[13]=0+A1*( w[0])+A2*(w[1]*w[2]);
        dexp[14]=1           +A2*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        dexp[3]=0; dexp[4]=0; dexp[5]=0;
        dexp[9]=0; dexp[10]=0;dexp[11]=0;
        dexp[15]=0;dexp[16]=0;dexp[17]=0;
        ///////
        TT Twv=w[0]*w[3]+w[1]*w[4]+w[2]*w[5];
        TT C0=-1.0/2;TT C1=(1.0-gama)/th/th;//TT C2=0;
        TT C3=((1/beta+gama-2)/th/th/th/th)*Twv;
        dexp[18]=          +C1*(-w[2]*w[5]-w[1]*w[4]-w[5]*w[2]-w[4]*w[1]) +C3*(-w[1]*w[1]-w[2]*w[2]);
        dexp[19]=C0*(-w[5])+C1*(w[1]*w[3]+w[4]*w[0])                      +C3*(w[0]*w[1]);
        dexp[20]=C0*( w[4])+C1*(w[2]*w[3]+w[5]*w[0])                      +C3*(w[0]*w[2]);
        dexp[24]=C0*( w[5])+C1*(w[0]*w[4]+w[3]*w[1])                      +C3*(w[1]*w[0]);
        dexp[25]=          +C1*(-w[2]*w[5]-w[0]*w[3]-w[5]*w[2]-w[3]*w[0]) +C3*(-w[0]*w[0]-w[2]*w[2]);
        dexp[26]=C0*(-w[3])+C1*(w[2]*w[4]+w[5]*w[1])                      +C3*(w[1]*w[2]);
        dexp[30]=C0*(-w[4])+C1*(w[0]*w[5]+w[3]*w[2])                      +C3*(w[0]*w[2]);
        dexp[31]=C0*( w[3])+C1*(w[1]*w[5]+w[4]*w[2])                      +C3*(w[1]*w[2]);
        dexp[32]=          +C1*(-w[1]*w[4]-w[0]*w[3]-w[4]*w[1]-w[0]*w[3]) +C3*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        dexp[21]=dexp[0]; dexp[22]=dexp[1]; dexp[23]=dexp[2];
        dexp[27]=dexp[6]; dexp[28]=dexp[7]; dexp[29]=dexp[8];
        dexp[33]=dexp[12];dexp[34]=dexp[13];dexp[35]=dexp[14];

    }
}
template <typename TT>
void CJW_Math<TT>::MyddEXP4(TT w[6],TT dw[6],TT ddexp[36])//指数坐标二阶导数
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        ///
        ddexp[0]=0; ddexp[1]=-dw[2]/2; ddexp[2]=dw[1]/2;
        ddexp[6]=dw[2]/2; ddexp[7]=0; ddexp[8]=-dw[0]/2;
        ddexp[12]=-dw[1]/2;ddexp[13]=dw[0]/2;ddexp[14]=0;
        ///////
        ddexp[3]=0; ddexp[4]=0; ddexp[5]=0;
        ddexp[9]=0; ddexp[10]=0;ddexp[11]=0;
        ddexp[15]=0;ddexp[16]=0;ddexp[17]=0;
        ///////
        ddexp[18]=          (-dw[2]*w[5]-dw[1]*w[4]-w[5]*dw[2]-w[4]*dw[1])/6;
        ddexp[19]=(-w[5])/2 +(dw[1]*w[3]+w[4]*dw[0])/6;
        ddexp[20]=(w[4])/2  +(dw[2]*w[3]+w[5]*dw[0])/6;
        ddexp[24]=(w[5])/2  +(dw[0]*w[4]+w[3]*dw[1])/6;
        ddexp[25]=          (-dw[2]*w[5]-dw[0]*w[3]-w[5]*dw[2]-w[3]*dw[0])/6;
        ddexp[26]=(-w[3])/2 +(dw[2]*w[4]+w[5]*dw[1])/6;
        ddexp[30]=(-w[4])/2 +(dw[0]*w[5]+w[3]*dw[2])/6;
        ddexp[31]=(w[3])/2  +(dw[1]*w[5]+w[4]*dw[2])/6;
        ddexp[32]=          (-dw[1]*w[4]-dw[0]*w[3]-w[4]*dw[1]-w[3]*dw[0])/6;
        ///////
        ddexp[21]=ddexp[0]; ddexp[22]=ddexp[1]; ddexp[23]=ddexp[2];
        ddexp[27]=ddexp[6]; ddexp[28]=ddexp[7]; ddexp[29]=ddexp[8];
        ddexp[33]=ddexp[12];ddexp[34]=ddexp[13];ddexp[35]=ddexp[14];

    }
    else
    {
        TT aerf=sin(th)/th;TT beta=2*(1-cos(th))/th/th;
        ///////
        TT Twdw=w[0]*dw[0]+w[1]*dw[1]+w[2]*dw[2];
        TT C0=beta/2;TT C1=(1.0-aerf)/th/th;
        TT C2=((aerf-beta)/th/th)*Twdw; TT C3=((beta/2-3.0/th/th*(1-aerf))/th/th)*Twdw;
        ddexp[0] =          +C1*(-w[2]*dw[2]-w[1]*dw[1]-dw[2]*w[2]-dw[1]*w[1])             +C3*(-w[1]*w[1]-w[2]*w[2]);
        ddexp[1] =C0*(-dw[2])+C1*(w[1]*dw[0]+dw[1]*w[0])                      +C2*(-w[2]) +C3*(w[0]*w[1]);
        ddexp[2] =C0*( dw[1])+C1*(w[2]*dw[0]+dw[2]*w[0])                      +C2*( w[1]) +C3*(w[0]*w[2]);
        ddexp[6] =C0*( dw[2])+C1*(w[0]*dw[1]+dw[0]*w[1])                      +C2*( w[2]) +C3*(w[1]*w[0]);
        ddexp[7] =          +C1*(-w[2]*dw[2]-w[0]*dw[0]-dw[2]*w[2]-dw[0]*w[0])             +C3*(-w[0]*w[0]-w[2]*w[2]);
        ddexp[8] =C0*(-dw[0])+C1*(w[2]*dw[1]+dw[2]*w[1])                      +C2*(-w[0]) +C3*(w[1]*w[2]);
        ddexp[12]=C0*(-dw[1])+C1*(w[0]*dw[2]+dw[0]*w[2])                      +C2*(-w[1]) +C3*(w[0]*w[2]);
        ddexp[13]=C0*( dw[0])+C1*(w[1]*dw[2]+dw[1]*w[2])                      +C2*( w[0]) +C3*(w[1]*w[2]);
        ddexp[14]=          +C1*(-w[1]*dw[1]-w[0]*dw[0]-dw[1]*w[1]-w[0]*dw[0])             +C3*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        ddexp[3]=0; ddexp[4]=0; ddexp[5]=0;
        ddexp[9]=0; ddexp[10]=0;ddexp[11]=0;
        ddexp[15]=0;ddexp[16]=0;ddexp[17]=0;
        ///////
        ////brief Twv
        C2=(aerf-beta)/th/th;C3=(beta/2-3.0/th/th*(1-aerf))/th/th;
        TT Twv=w[0]*w[3]+w[1]*w[4]+w[2]*w[5];
        TT Tdwv=dw[0]*w[3]+dw[1]*w[4]+dw[2]*w[5];TT Twdv=w[0]*dw[3]+w[1]*dw[4]+w[2]*dw[5];
        TT Snov=Twv*Twdw/th/th;
        TT Cw=(C1-beta/2)*Snov+C2*(Tdwv+Twdv-4*Snov);
        TT Cv=C2*Twdw;
        TT Cdw=C2*Twv;
        TT Cdv=beta/2;
        TT Cww=C2*Snov+C3*(Tdwv+Twdv-5*Snov);
        TT Cwdw=C3*Twv;
        TT Cwv=C3*Twdw;
        TT Cdwv=C1;
        TT Cwdv=C1;
        ddexp[18]=
                +Cwv*(-w[2]*w[5]-w[1]*w[4]-w[5]*w[2]-w[4]*w[1])+Cdwv*(-dw[2]*w[5]-dw[1]*w[4]-w[5]*dw[2]-w[4]*dw[1])+Cwdv*(-w[2]*dw[5]-w[1]*dw[4]-dw[5]*w[2]-dw[4]*w[1])
                +Cwdw*(-w[2]*dw[2]-w[1]*dw[1]-dw[2]*w[2]-dw[1]*w[1])+Cww*(-w[1]*w[1]-w[2]*w[2]);
        ddexp[19]=Cw*(-w[2])+Cv*(-w[5])+Cdw*(-dw[2])+Cdv*(-dw[5])
                +Cwv*(w[1]*w[3]+w[4]*w[0])+Cdwv*(dw[1]*w[3]+w[4]*dw[0])+Cwdv*(w[1]*dw[3]+dw[4]*w[0])
                +Cwdw*(w[1]*dw[0]+dw[1]*w[0])+Cww*(w[0]*w[1]);
        ddexp[20]=Cw*( w[1])+Cv*( w[4])+Cdw*( dw[1])+Cdv*( dw[4])
                +Cwv*(w[2]*w[3]+w[5]*w[0])+Cdwv*(dw[2]*w[3]+w[5]*dw[0])+Cwdv*(w[2]*dw[3]+dw[5]*w[0])
                +Cwdw*(w[2]*dw[0]+dw[2]*w[0])+Cww*(w[0]*w[2]);
        ddexp[24]=Cw*( w[2])+Cv*( w[5])+Cdw*( dw[2])+Cdv*( dw[5])
                +Cwv*(w[0]*w[4]+w[3]*w[1])+Cdwv*(dw[0]*w[4]+w[3]*dw[1])+Cwdv*(w[0]*dw[4]+dw[3]*w[1])
                +Cwdw*(w[0]*dw[1]+dw[0]*w[1])+Cww*(w[1]*w[0]);
        ddexp[25]=
                +Cwv*(-w[2]*w[5]-w[0]*w[3]-w[5]*w[2]-w[3]*w[0])+Cdwv*(-dw[2]*w[5]-dw[0]*w[3]-w[5]*dw[2]-w[3]*dw[0])+Cwdv*(-w[2]*dw[5]-w[0]*dw[3]-dw[5]*w[2]-dw[3]*w[0])
                +Cwdw*(-w[2]*dw[2]-w[0]*dw[0]-dw[2]*w[2]-dw[0]*w[0])+Cww*(-w[0]*w[0]-w[2]*w[2]);
        ddexp[26]=Cw*(-w[0])+Cv*(-w[3])+Cdw*(-dw[0])+Cdv*(-dw[3])
                +Cwv*(w[2]*w[4]+w[5]*w[1])+Cdwv*(dw[2]*w[4]+w[5]*dw[1])+Cwdv*(w[2]*dw[4]+dw[5]*w[1])
                +Cwdw*(w[2]*dw[1]+dw[2]*w[1])+Cww*(w[1]*w[2]);
        ddexp[30]=Cw*(-w[1])+Cv*(-w[4])+Cdw*(-dw[1])+Cdv*(-dw[4])
                +Cwv*(w[0]*w[5]+w[3]*w[2])+Cdwv*(dw[0]*w[5]+w[3]*dw[2])+Cwdv*(w[0]*dw[5]+dw[3]*w[2])
                +Cwdw*(w[0]*dw[2]+dw[0]*w[2])+Cww*(w[0]*w[2]);
        ddexp[31]=Cw*( w[0])+Cv*( w[3])+Cdw*( dw[0])+Cdv*( dw[3])
                +Cwv*(w[1]*w[5]+w[4]*w[2])+Cdwv*(dw[1]*w[5]+w[4]*dw[2])+Cwdv*(w[1]*dw[5]+dw[4]*w[2])
                +Cwdw*(w[1]*dw[2]+dw[1]*w[2])+Cww*(w[1]*w[2]);
        ddexp[32]=
                +Cwv*(-w[1]*w[4]-w[0]*w[3]-w[4]*w[1]-w[3]*w[0])+Cdwv*(-dw[1]*w[4]-dw[0]*w[3]-w[4]*dw[1]-w[3]*dw[0])+Cwdv*(-w[1]*dw[4]-w[0]*dw[3]-dw[4]*w[1]-dw[3]*w[0])
                +Cwdw*(-w[1]*dw[1]-w[0]*dw[0]-dw[1]*w[1]-dw[0]*w[0])+Cww*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        ddexp[21]=ddexp[0]; ddexp[22]=ddexp[1]; ddexp[23]=ddexp[2];
        ddexp[27]=ddexp[6]; ddexp[28]=ddexp[7]; ddexp[29]=ddexp[8];
        ddexp[33]=ddexp[12];ddexp[34]=ddexp[13];ddexp[35]=ddexp[14];
    }
}
template <typename TT>
void CJW_Math<TT>::MyddEXP4Inv(TT w[6],TT dw[6],TT ddexp[36])//指数坐标二阶导数逆矩阵
{
    TT th=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if((th<EXP_ZERO)&&(th>-EXP_ZERO))
    {
        ///
        ddexp[0]=0; ddexp[1]=dw[2]/2; ddexp[2]=-dw[1]/2;
        ddexp[6]=-dw[2]/2; ddexp[7]=0; ddexp[8]=dw[0]/2;
        ddexp[12]=dw[1]/2;ddexp[13]=-dw[0]/2;ddexp[14]=0;
        ///////
        ddexp[3]=0; ddexp[4]=0; ddexp[5]=0;
        ddexp[9]=0; ddexp[10]=0;ddexp[11]=0;
        ddexp[15]=0;ddexp[16]=0;ddexp[17]=0;
        ///////
        ddexp[18]=           (-dw[2]*w[5]-dw[1]*w[4]-w[5]*dw[2]-w[4]*dw[1])/12;
        ddexp[19]=-(-w[5])/2 +(dw[1]*w[3]+w[4]*dw[0])/6;
        ddexp[20]=-(w[4])/2  +(dw[2]*w[3]+w[5]*dw[0])/6;
        ddexp[24]=-(w[5])/2  +(dw[0]*w[4]+w[3]*dw[1])/6;
        ddexp[25]=           (-dw[2]*w[5]-dw[0]*w[3]-w[5]*dw[2]-w[3]*dw[0])/12;
        ddexp[26]=-(-w[3])/2 +(dw[2]*w[4]+w[5]*dw[1])/6;
        ddexp[30]=-(-w[4])/2 +(dw[0]*w[5]+w[3]*dw[2])/6;
        ddexp[31]=-(w[3])/2  +(dw[1]*w[5]+w[4]*dw[2])/6;
        ddexp[32]=           (-dw[1]*w[4]-dw[0]*w[3]-w[4]*dw[1]-w[3]*dw[0])/12;
        ///////
        ddexp[21]=ddexp[0]; ddexp[22]=ddexp[1]; ddexp[23]=ddexp[2];
        ddexp[27]=ddexp[6]; ddexp[28]=ddexp[7]; ddexp[29]=ddexp[8];
        ddexp[33]=ddexp[12];ddexp[34]=ddexp[13];ddexp[35]=ddexp[14];

    }
    else
    {
        TT aerf=sin(th)/th;TT beta=2*(1-cos(th))/th/th;TT gama=th*cos(th/2)/sin(th/2)/2;
        ///////
        TT Twdw=w[0]*dw[0]+w[1]*dw[1]+w[2]*dw[2];
        TT C0=-1.0/2;TT C1=(1.0-gama)/th/th;//TT C2=0;
        TT C3=((1/beta+gama-2)/th/th/th/th)*Twdw;
        ddexp[0] =          +C1*(-w[2]*dw[2]-w[1]*dw[1]-dw[2]*w[2]-dw[1]*w[1])+C3*(-w[1]*w[1]-w[2]*w[2]);
        ddexp[1] =C0*(-dw[2])+C1*(w[1]*dw[0]+dw[1]*w[0])                      +C3*(w[0]*w[1]);
        ddexp[2] =C0*( dw[1])+C1*(w[2]*dw[0]+dw[2]*w[0])                      +C3*(w[0]*w[2]);
        ddexp[6] =C0*( dw[2])+C1*(w[0]*dw[1]+dw[0]*w[1])                      +C3*(w[1]*w[0]);
        ddexp[7] =          +C1*(-w[2]*dw[2]-w[0]*dw[0]-dw[2]*w[2]-dw[0]*w[0])+C3*(-w[0]*w[0]-w[2]*w[2]);
        ddexp[8] =C0*(-dw[0])+C1*(w[2]*dw[1]+dw[2]*w[1])                      +C3*(w[1]*w[2]);
        ddexp[12]=C0*(-dw[1])+C1*(w[0]*dw[2]+dw[0]*w[2])                      +C3*(w[0]*w[2]);
        ddexp[13]=C0*( dw[0])+C1*(w[1]*dw[2]+dw[1]*w[2])                      +C3*(w[1]*w[2]);
        ddexp[14]=          +C1*(-w[1]*dw[1]-w[0]*dw[0]-dw[1]*w[1]-w[0]*dw[0])+C3*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        ddexp[3]=0; ddexp[4]=0; ddexp[5]=0;
        ddexp[9]=0; ddexp[10]=0;ddexp[11]=0;
        ddexp[15]=0;ddexp[16]=0;ddexp[17]=0;
        ///////
        /// \brief Twv
        C3=((1/beta+gama-2)/th/th/th/th);
        TT Twv=w[0]*w[3]+w[1]*w[4]+w[2]*w[5];
        TT Tdwv=dw[0]*w[3]+dw[1]*w[4]+dw[2]*w[5];TT Twdv=w[0]*dw[3]+w[1]*dw[4]+w[2]*dw[5];
        TT Snov=Twv*Twdw/th/th;
        TT Cw=0;
        TT Cv=0;
        TT Cdw=0;
        TT Cdv=-1.0/2;
        TT Cww=2*(1-gama/beta)/th/th/th/th*Snov+C3*(Tdwv+Twdv-3*Snov);
        TT Cwdw=C3*Twv;
        TT Cwv=C3*Twdw;
        TT Cdwv=C1;
        TT Cwdv=C1;
        ddexp[18]=
                +Cwv*(-w[2]*w[5]-w[1]*w[4]-w[5]*w[2]-w[4]*w[1])+Cdwv*(-dw[2]*w[5]-dw[1]*w[4]-w[5]*dw[2]-w[4]*dw[1])+Cwdv*(-w[2]*dw[5]-w[1]*dw[4]-dw[5]*w[2]-dw[4]*w[1])
                +Cwdw*(-w[2]*dw[2]-w[1]*dw[1]-dw[2]*w[2]-dw[1]*w[1])+Cww*(-w[1]*w[1]-w[2]*w[2]);
        ddexp[19]=Cw*(-w[2])+Cv*(-w[5])+Cdw*(-dw[2])+Cdv*(-dw[5])
                +Cwv*(w[1]*w[3]+w[4]*w[0])+Cdwv*(dw[1]*w[3]+w[4]*dw[0])+Cwdv*(w[1]*dw[3]+dw[4]*w[0])
                +Cwdw*(w[1]*dw[0]+dw[1]*w[0])+Cww*(w[0]*w[1]);
        ddexp[20]=Cw*( w[1])+Cv*( w[4])+Cdw*( dw[1])+Cdv*( dw[4])
                +Cwv*(w[2]*w[3]+w[5]*w[0])+Cdwv*(dw[2]*w[3]+w[5]*dw[0])+Cwdv*(w[2]*dw[3]+dw[5]*w[0])
                +Cwdw*(w[2]*dw[0]+dw[2]*w[0])+Cww*(w[0]*w[2]);
        ddexp[24]=Cw*( w[2])+Cv*( w[5])+Cdw*( dw[2])+Cdv*( dw[5])
                +Cwv*(w[0]*w[4]+w[3]*w[1])+Cdwv*(dw[0]*w[4]+w[3]*dw[1])+Cwdv*(w[0]*dw[4]+dw[3]*w[1])
                +Cwdw*(w[0]*dw[1]+dw[0]*w[1])+Cww*(w[1]*w[0]);
        ddexp[25]=
                +Cwv*(-w[2]*w[5]-w[0]*w[3]-w[5]*w[2]-w[3]*w[0])+Cdwv*(-dw[2]*w[5]-dw[0]*w[3]-w[5]*dw[2]-w[3]*dw[0])+Cwdv*(-w[2]*dw[5]-w[0]*dw[3]-dw[5]*w[2]-dw[3]*w[0])
                +Cwdw*(-w[2]*dw[2]-w[0]*dw[0]-dw[2]*w[2]-dw[0]*w[0])+Cww*(-w[0]*w[0]-w[2]*w[2]);
        ddexp[26]=Cw*(-w[0])+Cv*(-w[3])+Cdw*(-dw[0])+Cdv*(-dw[3])
                +Cwv*(w[2]*w[4]+w[5]*w[1])+Cdwv*(dw[2]*w[4]+w[5]*dw[1])+Cwdv*(w[2]*dw[4]+dw[5]*w[1])
                +Cwdw*(w[2]*dw[1]+dw[2]*w[1])+Cww*(w[1]*w[2]);
        ddexp[30]=Cw*(-w[1])+Cv*(-w[4])+Cdw*(-dw[1])+Cdv*(-dw[4])
                +Cwv*(w[0]*w[5]+w[3]*w[2])+Cdwv*(dw[0]*w[5]+w[3]*dw[2])+Cwdv*(w[0]*dw[5]+dw[3]*w[2])
                +Cwdw*(w[0]*dw[2]+dw[0]*w[2])+Cww*(w[0]*w[2]);
        ddexp[31]=Cw*( w[0])+Cv*( w[3])+Cdw*( dw[0])+Cdv*( dw[3])
                +Cwv*(w[1]*w[5]+w[4]*w[2])+Cdwv*(dw[1]*w[5]+w[4]*dw[2])+Cwdv*(w[1]*dw[5]+dw[4]*w[2])
                +Cwdw*(w[1]*dw[2]+dw[1]*w[2])+Cww*(w[1]*w[2]);
        ddexp[32]=
                +Cwv*(-w[1]*w[4]-w[0]*w[3]-w[4]*w[1]-w[3]*w[0])+Cdwv*(-dw[1]*w[4]-dw[0]*w[3]-w[4]*dw[1]-w[3]*dw[0])+Cwdv*(-w[1]*dw[4]-w[0]*dw[3]-dw[4]*w[1]-dw[3]*w[0])
                +Cwdw*(-w[1]*dw[1]-w[0]*dw[0]-dw[1]*w[1]-dw[0]*w[0])+Cww*(-w[0]*w[0]-w[1]*w[1]);
        ///////
        ddexp[21]=ddexp[0]; ddexp[22]=ddexp[1]; ddexp[23]=ddexp[2];
        ddexp[27]=ddexp[6]; ddexp[28]=ddexp[7]; ddexp[29]=ddexp[8];
        ddexp[33]=ddexp[12];ddexp[34]=ddexp[13];ddexp[35]=ddexp[14];
    }
}


/********************一些指数坐标的转换*************************************************/
//物体的姿态和物体坐标系速度 转指数坐标系
template <typename TT>
void CJW_Math<TT>::MyRVToExponet3(TT R[9],TT w[3],TT kesi[3],TT dkesi[3])
{
    //get kesi
    MyRToExponent3(R,kesi);
    //get dkesi
    TT kesi_1[3]={-kesi[0],-kesi[1],-kesi[2]};
    TT dexpInvkesi_1[9];MydEXP3Inv(kesi_1,dexpInvkesi_1);
    MyRCompositionw(dexpInvkesi_1,w,dkesi);
}
//物体的姿态和物体坐标系速度，加速度 转指数坐标系
template <typename TT>
void CJW_Math<TT>::MyRVAToExponet3(TT R[9],TT w[3],TT dw[3],TT kesi[3],TT dkesi[3],TT ddkesi[3])
{
    //get kesi
    MyRToExponent3(R,kesi);
    //get dkesi
    TT kesi_1[3]={-kesi[0],-kesi[1],-kesi[2]};
    TT dexpInvkesi_1[9];MydEXP3Inv(kesi_1,dexpInvkesi_1);
    MyRCompositionw(dexpInvkesi_1,w,dkesi);
    //get ddkesi
    TT dkesi_1[3]={-dkesi[0],-dkesi[1],-dkesi[2]};
    TT ddexpkesi_1[9];MyddEXP3(kesi_1,dkesi_1,ddexpkesi_1);
    TT temp1[3];MyRCompositionw(ddexpkesi_1,dkesi,temp1);
    TT out1[3]={dw[0]-temp1[0],dw[1]-temp1[1],dw[2]-temp1[2]};
    MyRCompositionw(dexpInvkesi_1,out1,ddkesi);
}
//指数坐标系 转物体的姿态和物体坐标系速度
template <typename TT>
void CJW_Math<TT>::MyExponet3ToRV(TT kesi[3],TT dkesi[3],TT R[9],TT w[3])
{
    //get R
    MyExponent3ToR(kesi,R);
    //get w
    TT kesi_1[3]={-kesi[0],-kesi[1],-kesi[2]};
    TT dexp_kesi_1[9];MydEXP3(kesi_1,dexp_kesi_1);
    MyRCompositionw(dexp_kesi_1,dkesi,w);
}
//指数坐标系 转物体的姿态和物体坐标系速度，加速度
template <typename TT>
void CJW_Math<TT>::MyExponet3ToRVA(TT kesi[3],TT dkesi[3],TT ddkesi[3],TT R[9],TT w[3],TT dw[3])
{
    //get R
    MyExponent3ToR(kesi,R);
    //get w
    TT kesi_1[3]={-kesi[0],-kesi[1],-kesi[2]};
    TT dexp_kesi_1[9];MydEXP3(kesi_1,dexp_kesi_1);
    MyRCompositionw(dexp_kesi_1,dkesi,w);
    //get dw
    TT dkesi_1[3]={-dkesi[0],-dkesi[1],-dkesi[2]};
    TT ddexp_kesi_1[9];MyddEXP3(kesi_1,dkesi_1,ddexp_kesi_1);
    TT out1[3];MyRCompositionw(dexp_kesi_1,ddkesi,out1);
    TT out2[3];MyRCompositionw(ddexp_kesi_1,dkesi,out2);
    dw[0]=out1[0]+out2[0];dw[1]=out1[1]+out2[1];dw[2]=out1[2]+out2[2];
}
//物体的位姿和物体坐标系速度 转指数坐标系
template <typename TT>
void CJW_Math<TT>::MyGVToExponet4(TT G[16],TT w[6],TT kesi[6],TT dkesi[6])
{
    //get kesi
    MyGToExponent4(G,kesi);
    //get dkesi
    TT kesi_1[6]={-kesi[0],-kesi[1],-kesi[2],-kesi[3],-kesi[4],-kesi[5]};
    TT dexpInvkesi_1[36];MydEXP4Inv(kesi_1,dexpInvkesi_1);
    MyMatrixMultiply(6,1,6,dexpInvkesi_1,w,dkesi);
}
//物体的位姿和物体坐标系速度，加速度 转指数坐标系
template <typename TT>
void CJW_Math<TT>::MyGVAToExponet4(TT G[16],TT w[6],TT dw[6],TT kesi[6],TT dkesi[6],TT ddkesi[6])
{
    //get kesi
    MyGToExponent4(G,kesi);
    //get dkesi
    TT kesi_1[6]={-kesi[0],-kesi[1],-kesi[2],-kesi[3],-kesi[4],-kesi[5]};
    TT dexpInvkesi_1[36];MydEXP4Inv(kesi_1,dexpInvkesi_1);
    MyMatrixMultiply(6,1,6,dexpInvkesi_1,w,dkesi);
    //get ddkesi
    TT dkesi_1[6]={-dkesi[0],-dkesi[1],-dkesi[2],-dkesi[3],-dkesi[4],-dkesi[5]};
    TT ddexpkesi_1[36];MyddEXP4(kesi_1,dkesi_1,ddexpkesi_1);
    TT temp1[6];MyMatrixMultiply(6,1,6,ddexpkesi_1,dkesi,temp1);
    TT out1[6]={dw[0]-temp1[0],dw[1]-temp1[1],dw[2]-temp1[2],dw[3]-temp1[3],dw[4]-temp1[4],dw[5]-temp1[5]};
    MyMatrixMultiply(6,1,6,dexpInvkesi_1,out1,ddkesi);
}
//指数坐标系 转物体的位姿和物体坐标系速度
template <typename TT>
void CJW_Math<TT>::MyExponet4ToGV(TT kesi[6],TT dkesi[6],TT G[16],TT w[6])
{
    //get G
    MyExponent4ToG(kesi,G);
    //get w
    TT kesi_1[6]={-kesi[0],-kesi[1],-kesi[2],-kesi[3],-kesi[4],-kesi[5]};
    TT dexp_kesi_1[36];MydEXP4(kesi_1,dexp_kesi_1);
    MyMatrixMultiply(6,1,6,dexp_kesi_1,dkesi,w);
}
//指数坐标系 转物体的位姿和物体坐标系速度，加速度
template <typename TT>
void CJW_Math<TT>::MyExponet4ToGVA(TT kesi[6],TT dkesi[6],TT ddkesi[6],TT G[16],TT w[6],TT dw[6])
{
    //get G
    MyExponent4ToG(kesi,G);
    //get w
    TT kesi_1[6]={-kesi[0],-kesi[1],-kesi[2],-kesi[3],-kesi[4],-kesi[5]};
    TT dexp_kesi_1[36];MydEXP4(kesi_1,dexp_kesi_1);
    MyMatrixMultiply(6,1,6,dexp_kesi_1,dkesi,w);
    //get dw
    TT dkesi_1[6]={-dkesi[0],-dkesi[1],-dkesi[2],-dkesi[3],-dkesi[4],-dkesi[5]};
    TT ddexp_kesi_1[36];MyddEXP4(kesi_1,dkesi_1,ddexp_kesi_1);
    TT out1[6];MyMatrixMultiply(6,1,6,dexp_kesi_1,ddkesi,out1);
    TT out2[6];MyMatrixMultiply(6,1,6,ddexp_kesi_1,dkesi,out2);
    dw[0]=out1[0]+out2[0];dw[1]=out1[1]+out2[1];dw[2]=out1[2]+out2[2];
    dw[3]=out1[3]+out2[3];dw[4]=out1[4]+out2[4];dw[5]=out1[5]+out2[5];
}

//SE(3)空间映射成SE（2）
template <typename TT>
void CJW_Math<TT>::MySE3ToSE2(TT G0[16],TT Gout[16])
{
    TT R[9]={G0[0],G0[1],G0[2],
                G0[4],G0[5],G0[6],
                G0[8],G0[9],G0[10]};
    TT RZYX[3]={0};MyRToEulerZYX(R,RZYX);
    Gout[0]=cos(RZYX[0]);Gout[1]=-sin(RZYX[0]);Gout[2]=0;Gout[3]=G0[3];
    Gout[4]=sin(RZYX[0]);Gout[5]=cos(RZYX[0]);Gout[6]=0;Gout[7]=G0[7];
    Gout[8]=0;Gout[9]=0;Gout[10]=1;Gout[11]=0;
    Gout[12]=0;Gout[13]=0;Gout[14]=0;Gout[15]=1;
}
template <typename TT>
void CJW_Math<TT>::MySE3ToSE2(TT G0[16],TT V0[6],TT Gout[16],TT Vout[6])
{
    TT R[9]={G0[0],G0[1],G0[2],
                G0[4],G0[5],G0[6],
                G0[8],G0[9],G0[10]};
    TT RZYX[3]={0};MyRToEulerZYX(R,RZYX);
    Gout[0]=cos(RZYX[0]);Gout[1]=-sin(RZYX[0]);Gout[2]=0;Gout[3]=G0[3];
    Gout[4]=sin(RZYX[0]);Gout[5]=cos(RZYX[0]);Gout[6]=0;Gout[7]=G0[7];
    Gout[8]=0;Gout[9]=0;Gout[10]=1;Gout[11]=0;
    Gout[12]=0;Gout[13]=0;Gout[14]=0;Gout[15]=1;
    Vout[0]=0;
    Vout[1]=0;
    Vout[2]=R[6]*V0[0]+R[7]*V0[1]+R[8]*V0[2];
    double absvxy[2]={R[0]*V0[3]+R[1]*V0[4]+R[2]*V0[4],
                      R[3]*V0[3]+R[4]*V0[4]+R[5]*V0[5]};
    Vout[3]=Gout[0]*absvxy[0]+Gout[4]*absvxy[0];
    Vout[4]=Gout[1]*absvxy[0]+Gout[5]*absvxy[1];
    Vout[5]=0;
}

/********************单刚体指数坐标PD控制*************************************************/
template <typename TT>
void CJW_Math<TT>::RigidBody_SO3_Pcontroller(TT Re[9],TT Ve[3],TT Rr[9],TT Kp[3],TT Vout[3])
{
    //TT tempvector6
    TT temp[36]={0};
    //get Sd
    TT Rd[9]={0};MyRInvCompositionR(Re,Rr,Rd);
    //Show(4,4,Gd);
    TT Sd[3]={0};MyRToExponent3(Rd,Sd);
    TT dSd[3]={0};
    dSd[0]=-Kp[0]*Sd[0];dSd[1]=-Kp[1]*Sd[1];dSd[2]=-Kp[2]*Sd[2];
    MydEXP3(dSd,temp);
    TT Vd[6]={0};MyRCompositionw(temp,dSd,Vd);
    TT VLL[6]={0};
    VLL[0]=Ve[0]+Vd[0];VLL[1]=Ve[1]+Vd[1];VLL[2]=Ve[2]+Vd[2];
    MyRInvCompositionw(Rd,VLL,Vout);
}
template <typename TT>
void CJW_Math<TT>::RigidBody_SE3_Pcontroller(TT Ge[16],TT Ve[6],TT Gr[16],TT Kp[6],TT Vout[6])
{
    //TT tempvector6
    TT temp[36]={0};
    //get Sd
    TT Gd[16]={0};MyGInvCompositionG(Ge,Gr,Gd);
    //Show(4,4,Gd);
    TT Sd[6]={0};MyGToExponent4(Gd,Sd);
    TT dSd[6]={0};
    dSd[0]=-Kp[0]*Sd[0];dSd[1]=-Kp[1]*Sd[1];dSd[2]=-Kp[2]*Sd[2];
    dSd[3]=-Kp[3]*Sd[3];dSd[4]=-Kp[4]*Sd[4];dSd[5]=-Kp[5]*Sd[5];
    MydEXP4(dSd,temp);
    TT Vd[6]={0};MyMatrixMultiply(6,1,6,temp,dSd,Vd);
    TT VLL[6]={0};
    VLL[0]=Ve[0]+Vd[0];VLL[1]=Ve[1]+Vd[1];VLL[2]=Ve[2]+Vd[2];
    VLL[3]=Ve[3]+Vd[3];VLL[4]=Ve[4]+Vd[4];VLL[5]=Ve[5]+Vd[5];
    MyEXPAdgInvScrew(Gd,VLL,Vout);
}
    
template <typename TT>
void CJW_Math<TT>::RigidBody_SO3_PDcontroller(TT Re[9],TT Ve[3],TT Ae[3],TT Rr[9],TT Vr[3],TT Kp[3],TT Kd[3],TT Aout[3])
{
    //TT tempvector6
    TT temp[36]={0};
    //get Sd
    TT Rd[9]={0};MyRInvCompositionR(Rr,Re,Rd);
    TT Sd[3]={0};MyRToExponent3(Rd,Sd);
    //get dSd
    MyRInvCompositionw(Rd,Vr,temp);
    TT Vd[3]={0};
    Vd[0]=Ve[0]-temp[0];Vd[1]=Ve[1]-temp[1];Vd[2]=Ve[2]-temp[2];
    TT Sd_1[6]={-Sd[0],-Sd[1],-Sd[2]};
    MydEXP3Inv(Sd_1,temp);
    TT dSd[3];MyRCompositionw(temp,Vd,dSd);
    //PD controller
    TT ddSd[3]={0};
    ddSd[0]=Kp[0]*Sd[0]+Kd[0]*dSd[0];ddSd[1]=Kp[1]*Sd[1]+Kd[1]*dSd[1];ddSd[2]=Kp[2]*Sd[2]+Kd[2]*dSd[2];
    //mapping to the Ar
    TT AA1[3]={0};
    MyVector3ToDisMatrix(Vd,temp);
    MyRCompositionw(temp,Vr,AA1);
    TT AA2[3]={0};
    TT dSd_1[3]={-dSd[0],-dSd[1],-dSd[2]};
    MyddEXP3(Sd_1,dSd_1,temp);
    MyRCompositionw(temp,dSd,AA2);
    TT AA3[3]={0};MydEXP3(Sd_1,temp);
    MyRCompositionw(temp,ddSd,AA3);
    TT ALL[3]={0};
    ALL[0]=Ae[0]-AA2[0]+AA3[0];ALL[1]=Ae[1]-AA2[1]+AA3[1];ALL[2]=Ae[2]-AA2[2]+AA3[2];
    TT ALL2[3]={0};MyRCompositionw(Rd,ALL,ALL2);
    Aout[0]=ALL2[0]+AA1[0];Aout[1]=ALL2[1]+AA1[1];Aout[2]=ALL2[2]+AA1[2];
}
template <typename TT>
void CJW_Math<TT>::RigidBody_SE3_PDcontroller(TT Ge[16],TT Ve[6],TT Ae[6],TT Gr[16],TT Vr[6],TT Kp[6],TT Kd[6],TT Aout[6])
{
    //TT tempvector6
    TT temp[36]={0};
    //get Sd
    TT Gd[16]={0};MyGInvCompositionG(Gr,Ge,Gd);
    //Show(4,4,Gd);
    TT Sd[6]={0};MyGToExponent4(Gd,Sd);
    //Show(1,6,Sd);
    //get dSd
    MyEXPAdgInvScrew(Gd,Vr,temp);
    TT Vd[6]={0};
    Vd[0]=Ve[0]-temp[0];Vd[1]=Ve[1]-temp[1];Vd[2]=Ve[2]-temp[2];
    Vd[3]=Ve[3]-temp[3];Vd[4]=Ve[4]-temp[4];Vd[5]=Ve[5]-temp[5];
    TT Sd_1[6]={-Sd[0],-Sd[1],-Sd[2],-Sd[3],-Sd[4],-Sd[5]};
    MydEXP4Inv(Sd_1,temp);
    TT dSd[6]={0};MyMatrixMultiply(6,1,6,temp,Vd,dSd);
    //PD controller
    TT ddSd[6]={0};
    ddSd[0]=Kp[0]*Sd[0]+Kd[0]*dSd[0];ddSd[1]=Kp[1]*Sd[1]+Kd[1]*dSd[1];ddSd[2]=Kp[2]*Sd[2]+Kd[2]*dSd[2];
    ddSd[3]=Kp[3]*Sd[3]+Kd[3]*dSd[3];ddSd[4]=Kp[4]*Sd[4]+Kd[4]*dSd[4];ddSd[5]=Kp[5]*Sd[5]+Kd[5]*dSd[5];
    //mapping to the Ar
    TT AA1[6]={0};
    MyEXPad(Vd,temp);
    MyMatrixMultiply(6,1,6,temp,Vr,AA1);
    TT AA2[6]={0};
    TT dSd_1[6]={-dSd[0],-dSd[1],-dSd[2],-dSd[3],-dSd[4],-dSd[5]};
    //Show(1,6,dSd_1);
    MyddEXP4(Sd_1,dSd_1,temp);
    //Show(6,6,temp);
    MyMatrixMultiply(6,1,6,temp,dSd,AA2);
    TT AA3[6]={0};MydEXP4(Sd_1,temp);
    MyMatrixMultiply(6,1,6,temp,ddSd,AA3);
    TT ALL[6]={0};
    ALL[0]=Ae[0]-AA2[0]+AA3[0];ALL[1]=Ae[1]-AA2[1]+AA3[1];ALL[2]=Ae[2]-AA2[2]+AA3[2];
    ALL[3]=Ae[3]-AA2[3]+AA3[3];ALL[4]=Ae[4]-AA2[4]+AA3[4];ALL[5]=Ae[5]-AA2[5]+AA3[5];
    TT ALL2[6]={0};MyEXPAdgScrew(Gd,ALL,ALL2);
    Aout[0]=ALL2[0]+AA1[0];Aout[1]=ALL2[1]+AA1[1];Aout[2]=ALL2[2]+AA1[2];
    Aout[3]=ALL2[3]+AA1[3];Aout[4]=ALL2[4]+AA1[4];Aout[5]=ALL2[5]+AA1[5];
}


/********************广义惯性矩阵的坐标转换****************************************************/
//trans based the body M to the absolute coordinate M
template <typename TT>
void CJW_Math<TT>::MyInertiaBasedGroup(TT M[36],TT G[16],TT outM[36])
{
    TT AdgInv[36]={0};MyEXPAdgInv(G,AdgInv);
    TT AdgInvT[36]={0};MyEXPAdgInvT(G,AdgInvT);
    TT temp[36]={0};MyMatrixMultiply(6,6,6,M,AdgInv,temp);
    MyMatrixMultiply(6,6,6,AdgInvT,temp,outM);
}
template <typename TT>
void CJW_Math<TT>::MyInertiaDBasedGroup(TT M[36],TT dM[36],TT G[16],TT VS[6],TT outM[36],TT outdM[36])
{
    //the M
    TT AdgInv[36]={0};MyEXPAdgInv(G,AdgInv);
    TT AdgInvT[36]={0};MyEXPAdgInvT(G,AdgInvT);
    TT temp0[36]={0};MyMatrixMultiply(6,6,6,M,AdgInv,temp0);
    MyMatrixMultiply(6,6,6,AdgInvT,temp0,outM);
    //the dM
    TT temp1[36]={0};MyMatrixMultiply(6,6,6,dM,AdgInv,temp1);
    TT tempdM[36]={0};MyMatrixMultiply(6,6,6,AdgInvT,temp1,tempdM);
    TT adv[36]={0};MyEXPad(VS,adv);
    TT advT[36]={0};MyEXPadT(VS,advT);
    TT tempdL[36]={0};MyMatrixMultiply(6,6,6,advT,outM,tempdL);
    TT tempdR[36]={0};MyMatrixMultiply(6,6,6,outM,adv,tempdR);
    for(int i=0;i<36;i++){outdM[i]=-tempdL[i]+tempdM[i]-tempdR[i];}
}

/*******************************基于惯性中心的动力学方程*****************************************/
template <typename TT>
void CJW_Math<TT>::DynamicBaseCOI(TT BodyM[36],TT BodydM[36],TT BodyV[6],TT BodyA[6],TT outFref[6])
{
    TT tempA1[6]={0};MyMatrixMultiply(6,1,6,BodyM,BodyA,tempA1);
    TT tempA2[6]={0};MyMatrixMultiply(6,1,6,BodydM,BodyV,tempA2);
    TT temp3[36]={0};MyEXPadT(BodyV,temp3);
    TT temp4[36]={0};MyMatrixMultiply(6,6,6,temp3,BodyM,temp4);
    TT tempA3[6]={0};MyMatrixMultiply(6,1,6,temp4,BodyV,tempA3);
    for(int i=0;i<6;i++)
    {
        outFref[i]=tempA1[i]+tempA2[i]-tempA3[i];//with the COI compen
        //outFref[i]=tempA1[i];//wihtout COI
    }
}


/*************************************************一些几何运算****************************************/
//三位空间点P0到（PA-PB）组成的直线距离
template <typename TT>
TT CJW_Math<TT>::Distance_PointToLine(TT P0[3],TT PA[3],TT PB[3])
{
    TT P0A[3]={P0[0]-PA[0],P0[1]-PA[1],P0[2]-PA[2]};
    TT BA[3]={PB[0]-PA[0],PB[1]-PA[1],PB[2]-PA[2]};
    TT PABCross[3]={0}; MyVector3Cross(P0A,BA,PABCross);
    TT BANorm2=MyVectorNorm2(BA,3);
    if(BANorm2<=0) return MyVectorNorm2(P0A,3);
    return MyVectorNorm2(PABCross,3)/BANorm2;
}
//空间内两条直线距离 A1-A2和B1-B2的距离
template <typename TT>
TT CJW_Math<TT>::Distanc_LineToLine(TT A1[3],TT A2[3],TT B1[3],TT B2[3])
{
    TT SA[3]={A1[0]-A2[0],A1[1]-A2[1],A1[2]-A2[2]};
    TT SB[3]={B1[0]-B2[0],B1[1]-B2[1],B1[2]-B2[2]};
    TT AB[3]={A1[0]-B1[0],A1[1]-B1[1],A1[2]-B1[2]};
    TT dirn[3];MyVector3Cross(SA,SB,dirn);
    TT dirN0=MyVectorNorm2(dirn,3);
    if(dirN0<=0)//两直线平行
    {
        if(MyVectorNorm1(SB,3)<=0)
        {
            return Distance_PointToLine(B1,A1,A2);
        }
        else
        {
            return Distance_PointToLine(A1,B1,B2);
        }
    }
    else
    {
        return fabs(MyVector3Dot(AB,dirn)/dirN0);
    }
}
//判断平面内 点 是否在 多边形(逆时针，点表示） 内部 //返回-1表示不在，0表示在边界，1表示在内部
//特别的 pplygon 存储顶点信息 X1 Y1 Z1 X2 Y2 Z2。。。。。
template <typename TT>
int CJW_Math<TT>::Inpolygon(int type,TT P0[3],TT* polygon,int polyN)
{
    if(polyN<3) return IN_POLYGON_OUT;
    if(type==IN_POLYGON_CROSS){return Inpolygon_Cross(P0,polygon,polyN);}
    else if(type==IN_POLYGON_ANGLE){return Inpolygon_Angle(P0,polygon,polyN);}
    return IN_POLYGON_OUT;
}
//判断平面内 点 距离 多边形(逆时针，点表示） 最近距离（内部大于0，外部小于0）
//特别的 pplygon 存储顶点信息 X1 Y1 Z1 X2 Y2 Z2。。。。
template <typename TT>
TT CJW_Math<TT>::MinDis_Inpolygon(int type,TT P0[3],TT* polygon,int polyN)
{
    if(type==IN_POLYGON_CROSS){return MinDis_Inpolygon_Cross(P0,polygon,polyN);}
    else if(type==IN_POLYGON_ANGLE){return MinDis_Inpolygon_Angle(P0,polygon,polyN);}
    return 0;
}

//叉乘算法计算点是否在多边形内部及距离，只适用于凸多边形
template <typename TT>
TT CJW_Math<TT>::TriPoint_Plane_Cross(TT p0[2],TT p1[2],TT p2[2])
{
    return (p1[0]-p0[0])*(p2[1]-p0[1])-(p1[1]-p0[1])*(p2[0]-p0[0]);
}
//内部 //返回-1表示不在，0表示在边界，1表示在内部
template <typename TT>
int CJW_Math<TT>::Inpolygon_Cross(TT P0[3],TT* polygon,int pointN)
{
    int point3N=3*pointN;
    TT flag0=TriPoint_Plane_Cross(&(polygon[0]),P0,&(polygon[3]));
    //printf("the F0 is %.3f\n",flag0);
    for(int thei=3;thei<point3N;thei=thei+3)
    {
        int nexti=(thei+3)%point3N;
        TT flag=flag0*TriPoint_Plane_Cross(&(polygon[thei]),P0,&(polygon[nexti]));
        //printf("the F%d is %.3f\n",thei,flag);
        if(flag==0){return IN_POLYGON_SIDE;}
        else if(flag<0){return IN_POLYGON_OUT;}
    }
    return IN_POLYGON_IN;
}
//最近距离（内部大于0，外部小于0）
template <typename TT>
TT CJW_Math<TT>::MinDis_Inpolygon_Cross(TT P0[3],TT* polygon,int pointN)
{
    int point3N=3*pointN;
    TT flag0=TriPoint_Plane_Cross(&(polygon[0]),P0,&(polygon[3]));
    TT MinDis=Distance_PointToLine(P0,&(polygon[0]),&(polygon[3]));
    int inp=IN_POLYGON_IN;
    //printf("the F0 is %.3f\n",flag0);
    for(int thei=3;thei<point3N;thei=thei+3)
    {
        int nexti=(thei+3)%point3N;
        TT flag=flag0*TriPoint_Plane_Cross(&(polygon[thei]),P0,&(polygon[nexti]));
        //printf("the F%d is %.3f\n",thei,flag);
        if(flag==0){return IN_POLYGON_SIDE;}
        else if(flag<0){inp=IN_POLYGON_OUT;}
        MinDis=fmin(MinDis,Distance_PointToLine(P0,&(polygon[thei]),&(polygon[nexti])));
    }
    return (MinDis*inp);
}
//内角合方法计算点是否在多边形内部及距离，只适用于所有多边形
template <typename TT>
TT CJW_Math<TT>::TriPoint_Plane_Angle(TT p0[2],TT p1[2],TT p2[2])
{
    TT dp1[2]={p1[0]-p0[0],p1[1]-p0[1]};
    TT dp2[2]={p2[0]-p0[0],p2[1]-p0[1]};
    TT m=(dp1[0]*dp1[0]+dp1[1]*dp1[1])*(dp2[0]*dp2[0]+dp2[1]*dp2[1]);
    if(m<=0){return 0;}
    TT n=dp1[0]*dp2[0]+dp1[1]*dp2[1];
    TT theangle=n/sqrt(m);
    if(fabs(theangle)>=1){return 0;}
    return acos(theangle);
}
//内部 //返回-1表示不在，0表示在边界，1表示在内部
template <typename TT>
int CJW_Math<TT>::Inpolygon_Angle(TT P0[3],TT* polygon,int pointN)
{
    TT allangle=0;
    int point3N=3*pointN;
    for(int thei=0;thei<point3N;thei=thei+3)
    {
        int nexti=(thei+3)%point3N;
        TT theangle=TriPoint_Plane_Angle(P0,&(polygon[thei]),&(polygon[nexti]));
        if(theangle<=0) {return IN_POLYGON_SIDE;}
        allangle=allangle+theangle;
    }
    //printf("allangle is %.3f\n",allangle/M_PI*180);
    if(fabs(fmod(allangle,(2*M_PI)))<MATH_ZERO){return IN_POLYGON_IN;}
    return IN_POLYGON_OUT;
}
//最近距离（内部大于0，外部小于0）
template <typename TT>
TT CJW_Math<TT>::MinDis_Inpolygon_Angle(TT P0[3],TT* polygon,int pointN)
{
    TT allangle=0;
    int point3N=3*pointN;
    TT MinDis=INF;
    for(int thei=0;thei<point3N;thei=thei+3)
    {
        int nexti=(thei+3)%point3N;
        TT theangle=TriPoint_Plane_Angle(P0,&(polygon[thei]),&(polygon[nexti]));
        if(theangle<=0) {return IN_POLYGON_SIDE;}
        MinDis=fmin(MinDis,Distance_PointToLine(P0,&(polygon[thei]),&(polygon[nexti])));
        allangle=allangle+theangle;
    }
    //printf("allangle is %.3f\n",allangle/M_PI*180);
    if(fabs(fmod(allangle,(2*M_PI)))<MATH_ZERO){return (MinDis*IN_POLYGON_IN);}
    return (MinDis*IN_POLYGON_OUT);
}

/*************************************************************************/
template class CJW_Math<double>;
template class CJW_Math<float>;
/*************************************************************************/


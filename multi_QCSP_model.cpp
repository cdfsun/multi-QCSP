#include <ilcplex/ilocplex.h>
#include  <fstream>
#include  <iostream>

//#define TEST
#define ReduConst// use reduced set in constraints
//#define UserActive//activate the usercallback
#define _AFXDLL
//#define OUTBRANCH//输出分支节点处的解
//#define SubTour

#include "afxtempl.h"
#include  <afx.h>
#include <afxdb.h>
#include <time.h>
#include "stdafx.h"
#include <iomanip>


//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif



// 唯一的应用程序对象

CWinApp theApp;

using namespace std;

const   IloInt nbTask = 15;
const   IloInt nbBay=10;
const   IloInt nbCrane=2;

const   IloInt nbTime=100;

const   IloInt crane_move_time = 1;

const IloBool ORIGINAL=1;

IloInt SubCutNum;
IloNum ObjVal;
IloNum LB;
IloNum LB0;
IloNum UB;
IloNum GapF;

IloNum GapF2;

#define ModNo 1//选择要调用的模型编号

#define instance_no 1//选择测试的算例数量

//#define Bay_Task_Equal//

//#define new_constraints

typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<IloBoolArray>    BoolMatrix;
typedef IloArray<IloArray<IloBoolArray> > BoolMatrix2;
typedef IloArray<IloArray<IloNumArray> > NumMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolArray> > > BoolMatrix3;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloBoolVarArray>    BoolVarMatrix;
typedef IloArray<IloArray<IloBoolVarArray> > BoolVarMatrix2;
typedef IloArray<IloArray<IloNumVarArray> > NumVarMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolVarArray> > > BoolVarMatrix3;

ILOSTLBEGIN





//BoolMatrix B_SI;
//IloNumArray  yFsol[nbTime][nbSystem][nbBay];
//IloNumVarArray   uF[nbTime][nbSystem];
//IloNumArray   uFsol[nbTime][nbSystem];
//IloBoolVarArray   zF[nbTime][nbSystem];
//IloBoolVarArray   xF[nbTime][nbSystem];
//IloNumArray   zFsol[nbTime][nbSystem];
//IloBoolArray   B_P_KI[nbMaterial];
//IloBoolArray   B_Q_SK[nbSystem];


BOOL CBMP_R(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloBoolArray  *nbprecR, IloIntArray  nbLocation, IloIntArray  nbreadyT,
BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix zF, IloNumVar CF, IloNumVarArray QC_CF, IloNumVarArray Task_CF,
	BoolVarMatrix xF_begin, BoolVarMatrix xF_end, 
	IloNum *ObjVal);// with allocation varaible yF

BOOL CBMP(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloBoolArray  *nbprecR, IloIntArray  nbLocation, IloIntArray  nbreadyT,
	BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix zF, IloNumVar CF, IloNumVarArray QC_CF, IloNumVarArray Task_CF,
	BoolVarMatrix xF_begin, BoolVarMatrix xF_end,
	IloNum *ObjVal);//kim and park




int _tmain(int argc, char* argv[], char* envp[])
{
	double GapAve[instance_no];
	double GapAve2[instance_no];
	double DuraAve[instance_no];
	double ObjAve[instance_no];
	double IterAve[instance_no];
	for (int my = 1; my <= instance_no; my++)
	{
		cout << "data-" << my << "数据运行中..." << endl;
		char* filename;
		char dream[100] = "test/data";
		filename = dream;
		char C1[3];
		char C2[3];
		char C3[3];
		char C4[3];
		itoa(my, C1, 10);

		//此处可编辑规模，以输入
		itoa(nbBay, C2, 10);
		itoa(nbCrane, C3, 10);
		itoa(nbTask, C4, 10);

		strcpy(filename, "test/");
		//此处可编辑规模，以输入
		strcat(filename, C4);
		strcat(filename, "-");
		strcat(filename, C2);
		strcat(filename, "-");
		strcat(filename, C3);

		strcat(filename, "/data");
		strcat(filename, "-");
		strcat(filename, C1);
		strcat(filename, ".txt");
		clock_t start = 0, finish = 0;
		clock_t start_MP = 0, finish_MP = 0;
		clock_t start_SP = 0, finish_SP = 0;
		start = clock();
		double  duration;
		double  duration_MP = 0;
		double  duration_SP = 0;

		if (argc <= 1)
		{
			cerr << "Usage: " << argv[0] << " <model>" << endl;
			cerr << "  model = 0 -> convex piecewise linear model, " << endl;
			cerr << "  model = 1 -> concave piecewise linear model. [default]" << endl;
		}

		IloBool convex;
		if (argc <= 1)
			convex = IloFalse;
		else
			convex = atoi(argv[1]) == 0 ? IloTrue : IloFalse;

		//*************************//
		//     打开文件data.txt    //
		//*************************//
		ifstream fin(filename);
		IloEnv env;
		try 
		{

			//	定义原问题模型变量
			IloModel model(env);

			IloNum gap;

			//	定义临时变量
			IloInt i,k,j,t,kk; 

			//	定义连续决策变量CF
			IloNumVar   CF(env, 0, 1800);

			//IloNumArray tem_C(env, nbCrane);// 子问题计算的completion time
			IloNumVarArray QC_CF(env, nbCrane,0, IloInfinity);// each QC's completion time
			IloNumVarArray Task_CF(env, nbTask, 0, IloInfinity);// each task's completion time

			//	定义决策变量CxF,CyF,CzF
			BoolVarMatrix yF(env, nbCrane);
			for(k = 0; k < nbCrane; k++)
			{				
				yF[k]=IloBoolVarArray(env,nbTask);
			}

			BoolVarMatrix2 xF(env);
			for(k = 0; k < nbCrane; k++)
			{				
				xF.add(BoolVarMatrix(env));
				for(i = 0; i < nbTask; i++)
				{				
					xF[k].add(IloBoolVarArray(env,nbTask));
				}
			}
			BoolVarMatrix zF(env, nbTask);
			for (i = 0; i < nbTask; i++)
			{				
				zF[i]=IloBoolVarArray(env,nbTask);
			}

			BoolVarMatrix xF_begin(env, nbCrane);
			for (i = 0; i < nbCrane; i++)
			{
				xF_begin[i] = IloBoolVarArray(env, nbTask);
			}
			BoolVarMatrix xF_end(env, nbCrane);
			for (i = 0; i < nbCrane; i++)
			{
				xF_end[i] = IloBoolVarArray(env, nbTask);
			}

			//	定义保存决策变量CxF,CyF,CzF的最优值
			BoolMatrix yF_best(env, nbCrane);
			for(k = 0; k < nbCrane; k++)
			{				
				yF_best[k]=IloBoolArray(env,nbTask);
			}

			BoolMatrix2 xF_best(env);
			for(k = 0; k < nbCrane; k++)
			{
				xF_best.add(BoolMatrix(env));
				for(i = 0; i <nbTask; i++)
					xF_best[k].add(IloBoolArray(env,nbTask));
			}
			BoolMatrix zF_best(env, nbTask);
			for (i = 0; i <nbTask; i++)
			{				
				zF_best[i]=IloBoolArray(env,nbTask);
			}

			BoolMatrix xF_begin_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				xF_begin_best[k] = IloBoolArray(env, nbTask);
			}

			BoolMatrix xF_end_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				xF_end_best[k] = IloBoolArray(env, nbTask);
			}


			IloNum   CF_best;
			IloNumArray QC_CF_best(env, nbCrane);
			IloNumArray Task_CF_best(env, nbTask);




			//*******************************************************//
			//                 定义输入参数                          //
			//*******************************************************// 	 
			//约束中参数
			//	定义能力常量参数s
			IloNum   nbs;
			//定义吊机位置初始状态参数
			IloBoolArray nbb(env, nbCrane);
			//	定义任务量参数nbQ
			IloNumArray  nbQ(env, nbTask);

			IloIntArray  nbreadyT(env, nbCrane);

			//	定义任务所在贝位参数nbQ
			IloIntArray  nbLocation(env, nbTask);

			IloBoolArray  nbprecR[nbTask];//优先级关系
			for (i = 0; i < nbTask; i++)  nbprecR[i] = IloBoolArray(env, nbTask);
			for (int ai = 0; ai < nbTask; ai++)
				for (int bi = 0; bi < nbTask; bi++)
					nbprecR[ai][bi] = 0;

			//*******************************************************//
			//                 定义callback参数                      //
			//*******************************************************// 	 


	 
			//*******************************************************//
			//                 读入参数数据                          //
			//*******************************************************// 	 

			IloIntArray ls1(env, 8);

			fin >> ls1;
			//cout << ls1[2] << "  ";
			fin >> nbQ;
			fin >> nbLocation;
			fin >> nbreadyT;
			fin >> nbb;


	

			//for (i = 0; i < nbBay; i++) cout<<nbQ[i]<<"  ";
			//cout<<endl;

			////读入数据nbLocation[i]
			//for (i = 0; i < nbTask; i++)
			//{
			//	nbLocation[i]--;
			//	//cout<<nbLocation[i]<<"  ";
			//}

			////	读入数据nbb[i]
			//for (i = 0; i < nbCrane; i++)
			//{
			//	//fin >> nbb[i];
			//	cout<<nbb[i]<<"  ";
			//}
			//cout<<endl;

			IloIntArray prec(env, 2);

			//	读入数据nbprecR
			for (i = 0; i < ls1[2]; i++)
			{
				fin >> prec;

				for (int ai = 0; ai < nbTask; ai++)
					if (ai == prec[0] - 1)
						for (int bi = 0; bi < nbTask; bi++)
							if (bi == prec[1] - 1)
								nbprecR[ai][bi] = 1;


			}


			//for (int ai = 0; ai < nbTask; ai++)
			//	for (int bi = 0; bi < nbTask; bi++)
			//		if (nbprecR[ai][bi] == 1)
			//		{
			//	cout << "(" << ai + 1 << "," << bi + 1 << "), ";
			//		}cout << endl;


			nbs = 1;




			//************************************************************************************//
			//      至此已知数据输入完毕                                                          //
			//************************************************************************************//
			char* filename1;
			char dream1[100] = "result/C_data";
			filename1 = dream1;

			itoa(my, C1, 10);

			//此处可编辑规模，以输入
			itoa(nbBay, C2, 10);
			itoa(nbCrane, C3, 10);
			itoa(nbTask, C4, 10);

			strcpy(filename1, "result/");
			//此处可编辑规模，以输入
			strcat(filename1, C4);
			strcat(filename1, "-");
			strcat(filename1, C2);
			strcat(filename1, "-");
			strcat(filename1, C3);

			strcat(filename1, "/Cdata");
			strcat(filename1, "-");
			strcat(filename1, C1);
			strcat(filename1, ".txt");

			ofstream fout(filename1);

			cout<<filename1<<endl;

			//*******************************************************//
			//                 求解主问题 CBMP                       //
			//*******************************************************// 





	//主问题模型

			//CBMP_R(model, nbs, nbb, nbQ, nbprecR,  nbLocation, nbreadyT,
			//	xF, yF, zF, CF, QC_CF,Task_CF,
			//	xF_begin, xF_end, 
			//	&ObjVal);

			CBMP(model, nbs, nbb, nbQ, nbprecR, nbLocation, nbreadyT,
				xF, yF, zF, CF, QC_CF, Task_CF,
				xF_begin, xF_end,
				&ObjVal);

	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(model);

#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_FEASIBILITY);
	cplex.setParam(IloCplex::EpGap,0.0006);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim,300);
	//cplex.setParam(IloCplex::Threads, 4);

	//cplex.setParam(IloCplex::VarSel, 3);
	//cplex.setParam(IloCplex::MIPEmphasis, 3);
	//cplex.setParam(IloCplex::Probe, 3);

	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;

	//************************************************************************************//
	//      初处理过程：no re-assignment                                                  //
	//************************************************************************************//
	UB=5000;


	LB = 0;
	//cout<<H_UB<<endl;

	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		model.end();

		return FALSE;
	}











			//**********************************//
			//             输出最好解           //
			//**********************************//	

	gap = cplex.getMIPRelativeGap();
	ObjVal = cplex.getBestObjValue();
	CF_best = cplex.getValue(CF);

			//cout<<"xF_best"<<endl;
			for (k = 0; k < nbCrane; k++)
			{
				for (i = 0; i < nbTask; i++)
				{
					for (j = 0; j < nbTask; j++)
					{
						if (cplex.isExtracted(xF[k][i][j]))
						{
							if (cplex.getValue(xF[k][i][j])>0.1)
								xF_best[k][i][j] = 1;
							else
								xF_best[k][i][j] = 0;
							if (xF_best[k][i][j] == 1) cout << " QC " << k << " move from " << i << " to " << j << endl;
						}
						else
							xF_best[k][i][j] = 0;
					}
					//cout<<"    ";cout<<endl;
				}
				//cout << "    "; cout << endl;
			}
			cout << endl << endl;

			cout << "yF_best" << endl;
			for (k = 0; k < nbCrane; k++)
			{
				for (i = 0; i < nbTask; i++)
				{

					if (cplex.isExtracted(yF[k][i]))
					{
						if (cplex.getValue(yF[k][i])>0.1)
							yF_best[k][i] = 1;
						else
							yF_best[k][i] = 0;
					}
					else
						yF_best[k][i] = 0;
					cout << yF_best[k][i] << "  ";


					if (cplex.isExtracted(xF_begin[k][i]))
					{
						if (cplex.getValue(xF_begin[k][i])>0.1)
							xF_begin_best[k][i] = 1;
						else
							xF_begin_best[k][i] = 0;
						//cout<<xF_begin_best[k][i]<<"  ";
					}
					else
						xF_begin_best[k][i] = 0;

					if (cplex.isExtracted(xF_end[k][i]))
					{
						if (cplex.getValue(xF_end[k][i])>0.1)
							xF_end_best[k][i] = 1;
						else
							xF_end_best[k][i] = 0;
						//cout<<xF_end_best[k][i]<<"  ";
					}
					else
						xF_end_best[k][i] = 0;

				}
				cout << "    "; cout << endl;
			}
			cout << endl << endl;


			cout << "xF_end_best" << endl;
			for (k = 0; k < nbCrane; k++)
			{
				for (i = 0; i < nbTask; i++)
				{

					if (cplex.isExtracted(xF_end[k][i]))
					{
						if (cplex.getValue(xF_end[k][i])>0.1)
							xF_end_best[k][i] = 1;
						else
							xF_end_best[k][i] = 0;
						cout<<xF_end_best[k][i]<<"  ";
					}
					else
						xF_end_best[k][i] = 0;

				}
				cout << "    "; cout << endl;
			}
			cout << endl << endl;


			cout << "xF_begin_best" << endl;
			for (k = 0; k < nbCrane; k++)
			{
				for (i = 0; i < nbTask; i++)
				{

					if (cplex.isExtracted(xF_begin[k][i]))
					{
						if (cplex.getValue(xF_begin[k][i])>0.1)
							xF_begin_best[k][i] = 1;
						else
							xF_begin_best[k][i] = 0;
						cout<<xF_begin_best[k][i]<<"  ";
					}
					else
						xF_begin_best[k][i] = 0;

				}
				cout << "    "; cout << endl;
			}
			cout << endl << endl;

			cout<<"QC_CF_best"<<endl;
			for (k = 0; k < nbCrane; k++)
			{
				QC_CF_best[k] = cplex.getValue(QC_CF[k]);
				cout << QC_CF_best[k] << "  ";
			}
			cout << endl << endl;

			cout << "Task_CF_best" << endl;
			for (i = 0; i < nbTask; i++)
			{
				Task_CF_best[i] = cplex.getValue(Task_CF[i]);
				cout << Task_CF_best[i] << "  ";
			}
			cout << endl << endl;
				
			cout << "CF_best: " << CF_best << endl;
			cout << "ObjVal: " << ObjVal << endl;

			//**********************************//
			//             计算目标值           //
			//**********************************//	
			double a,b; 
			a=0,b=0;
			//	计算OBJ1
			//for(i = 0; i < nbMaterial; i++)for(j = 0; j < nbBay; j++) for(k = 0; k < nbBay; k++) a+=nbcF[i][j][k]*xF_best[i][j][k];			
			//for(j = 0; j < nbBay; j++) for(k = 0; k < nbBay; k++) a+=nNs*wF_best[k][k];


			//计算OBJ2
			b=CF_best;

			UB = CF_best;
			LB = ObjVal;

			GapAve[my-1]=gap;

out_end:
			//fout.precision(10);
			////fout<<"OBJ1 = "<<a<<endl;
			//fout<<"OBJ2 = "<<b<<endl;
			//fout.precision(10);
			//fout<< "\nobj1+obj2 = " << ObjVal << endl;





			fout.precision(10);
			fout<< "\LB = " << LB << endl;
			fout<< "\UB = " << UB << endl;
			//fout<< "\H_UB = " << H_UB << endl;



			finish=clock();
			duration = (double)(finish - start) / CLOCKS_PER_SEC;
			fout<<"Time="<<duration<<endl;
			cout<<"Time="<<duration<<endl;

			DuraAve[my-1]=duration;
			ObjAve[my - 1]=UB;

			//cout<<SubCutNum<<endl;
			//cout<<"lala"<<endl;

			cplex.clearModel();
			cplex.clear();
			cplex.end();		
			model.end();


		}
		catch (IloException& e) 
		{
			cerr << "ERROR: " << e.getMessage() << endl;
		}
		catch (...) 
		{
			cerr << "Error" << endl;
		}
TERMINATE : env.end();
	}
	double gapa=0;
	double gapa2 = 0;
	double durationa=0;
	double obja = 0;
	double itera = 0;
	for (int i=0;i<instance_no;i++)
	{
		gapa+=GapAve[i];
	}

	for (int i=0;i<instance_no;i++)
	{
		durationa+=DuraAve[i];
	}

	for (int i = 0; i<instance_no; i++)
	{
		obja += ObjAve[i];
	}



	durationa=durationa/instance_no;
	gapa=gapa/instance_no;
	gapa2 = gapa2 / instance_no;
	obja = obja / instance_no;
	itera = itera / instance_no;
	cout<<"average gap: "<<gapa<<endl;
	cout<<"average CPU: "<<durationa<<endl;
	cout << "average obj: " << obja << endl;


	system("pause");
	return 0;
}










BOOL CBMP_R(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloBoolArray  *nbprecR, IloIntArray  nbLocation, IloIntArray  nbreadyT,
	BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix zF, IloNumVar CF, IloNumVarArray QC_CF, IloNumVarArray Task_CF,
	BoolVarMatrix xF_begin, BoolVarMatrix xF_end, 
	IloNum *ObjVal)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;

	IloInt nbM = 1800;



	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	//**********************************//
	//            主问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj(env);
	//  建立子问题目标函数表达式 

	obj += CF;

	//for (k = 0; k < nbCrane; k++)
	//	obj +=  (1/(float)nbCrane)* QC_CF[k];

	//	将目标函数加入到主问题模型
	model.add(IloMinimize(env, obj));//

	obj.end();

	IloNumArray ctime_t[nbTask];// the travelling time of QC moving from the bay i to the bay j
	for (i = 0; i < nbTask; i++)
		ctime_t[i] = IloNumArray(env, nbTask);


	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i] == nbLocation[j]) ctime_t[i][j] = 0;
			else if (nbLocation[i] > nbLocation[j])
			{
				ctime_t[i][j] = (nbLocation[i] - nbLocation[j])*crane_move_time;
			}
			else
			{
				ctime_t[i][j] = (nbLocation[j] - nbLocation[i])*crane_move_time;
			}
		}

	}

	//**********************************//
	//            子问题 约束           //
	//**********************************//


	IloRangeArray  c1(env);
	for (k = 0; k < nbCrane; k++)
	{

		IloExpr  epa(env);
		epa += CF;
		epa -= QC_CF[k];

		c1.add(epa >= 0);
		epa.end();

	}
	model.add(c1);
	c1.end();

	//IloRangeArray  c0(env);
	//for (i = 0; i < nbTask; i++)
	//{

	//	IloExpr  epa(env);
	//	epa += CF;
	//	epa -= Task_CF[i];

	//	c0.add(epa >= 0);
	//	epa.end();

	//}
	//model.add(c0);
	//c0.end();


	//约束（2）//identify the start task of each QC
	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF_begin[k][i];
		c2.add(epa == 1);
		epa.end();
	}
	model.add(c2);
	c2.end();

	//约束（3）//identify the end task of each QC
	IloRangeArray  c3(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF_end[k][i];
		c3.add(epa == 1);
		epa.end();
	}
	model.add(c3);
	c3.end();

	//约束（4）//task allocation
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa += yF[k][i] - xF_end[k][i];
		for (j = 0; j < nbTask; j++)
		{
			//if (i != j)
				epa -= xF[k][i][j];
		}


		c4.add(epa == 0);
		epa.end();
		}
	model.add(c4);
	c4.end();

	//约束（5）//task allocation
	IloRangeArray  c5(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa += yF[k][i] - xF_begin[k][i];
		for (j = 0; j < nbTask; j++)
		{
			//if (i != j)
				epa -= xF[k][j][i];
		}
		c5.add(epa == 0);
		epa.end();
		}
	model.add(c5);
	c5.end();

	//约束（6）//each task scheduled	
	IloRangeArray  c6(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
		{
			epa += yF[k][i];
		}
		c6.add(epa == 1);
		epa.end();
	}


	model.add(c6);
	c6.end();


	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
			for (j = 0; j < nbTask; j++)
			{
				IloExpr  epa(env);
				epa += Task_CF[i] + ctime_t[i][j] + nbQ[j] - Task_CF[j];
				epa -= nbM;
				epa += nbM*xF[k][i][j];

				c7.add(epa <= 0);
				epa.end();
			}
		}
	model.add(c7);
	c7.end();

	IloRangeArray  c8(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i] != nbLocation[j])
			{
				IloExpr  epa(env);
				epa += Task_CF[i] + nbQ[j] - Task_CF[j];
				epa -= nbM*(1 - zF[i][j]);
				c8.add(epa <= 0);
				epa.end();
			}
		}
	}
	model.add(c8);
	c8.end();

	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbLocation[i] != nbLocation[j])
			{
			IloExpr  epa(env);
			epa += Task_CF[j] - nbQ[j] - Task_CF[i];
			epa -= nbM*zF[i][j];

			c9.add(epa <= 0);
			epa.end();
			}
	}

	model.add(c9);
	c9.end();

	IloRangeArray  c10(env);
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i]<nbLocation[j])
			{
				IloExpr  epa(env);
				for (int kk = 0; kk <= k; kk++)
					epa += yF[kk][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa += yF[kk][i];
				epa -= zF[i][j] + zF[j][i];
				c10.add(epa <= 1);
				epa.end();
			}
		}
		}
	model.add(c10);
	c10.end();

	IloRangeArray  c11(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbLocation[i] == nbLocation[j])
			{
			IloExpr  epa(env);
			epa += Task_CF[i] + nbQ[j] - Task_CF[j];
			epa -= nbM*(1 - zF[i][j]);
			for (k = 0; k < nbCrane; k++)
			{
				epa += crane_move_time*xF_begin[k][j];
				for (int ii = 0; ii < nbTask; ii++)
					if (nbLocation[ii] != nbLocation[i])
						epa += crane_move_time*xF[k][ii][j];
			}
			c11.add(epa <= 0);
			epa.end();
			}
	}

	model.add(c11);
	c11.end();


	IloRangeArray  c12(env);

	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbLocation[i] == nbLocation[j])
			{
			IloExpr  epa(env);
			epa += Task_CF[j] - nbQ[j] - Task_CF[i];
			epa -= nbM*zF[i][j];
			for (k = 0; k < nbCrane; k++)
			{
				epa -= crane_move_time*xF_begin[k][j];
				for (int ii = 0; ii < nbTask; ii++)
					if (nbLocation[ii] != nbLocation[i])
						epa -= crane_move_time*xF[k][ii][j];
			}
			c12.add(epa <= 0);
			epa.end();
			}
	}

	model.add(c12);
	c12.end();


	IloRangeArray  c13(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbprecR[i][j] == 1)
			{
				//IloExpr  epa(env);
				//epa += zF[i][j] + zF[j][i];

				//c13.add(epa == 1);
				//epa.end();

				c13.add(zF[i][j] == 1);
				c13.add(zF[j][i] == 0);
			}
		}
	}
	model.add(c13);
	c13.end();

	IloRangeArray  c14(env);//i<j, and both located in the same bay, then i should be processed before j
	for (i = 0; i < nbTask; i++)
	{
		for (j = i + 1; j < nbTask; j++)
		{
			if (nbLocation[i] == nbLocation[j])
			{
				c14.add(zF[i][j] == 1);
				c14.add(zF[j][i] == 0);
			}
		}
	}
	model.add(c14);
	c14.end();

	IloRangeArray  c15(env);// 没加初始移动时间
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
		{
			IloExpr  epa(env);
			epa += nbQ[j] - Task_CF[j];
			epa += nbreadyT[k];//QC ready time, already defined in the main function
			epa -= nbM*(1 - xF_begin[k][j]);

			c15.add(epa <= 0);
			epa.end();
		}

	}
	model.add(c15);
	c15.end();

	IloRangeArray  c16(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
		{
			IloExpr  epa(env);
			epa += Task_CF[j] - QC_CF[k];
			epa -= nbM*(1 - xF_end[k][j]);

			c16.add(epa <= 0);
			epa.end();
		}

	}
	model.add(c16);
	c16.end();

	//**********************************//
	//            新 约束           //
	//**********************************//
#ifdef new_constraints

	////**********************************//
	////            有效不等式          //
	////**********************************//

	//优先级约束
	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			IloExpr  epa(env);
			for (j = 0; j < nbTask; j++)
				if (nbprecR[i][j] == 1)
					epa += xF[k][j][i];
			model.add(epa <= 0);
		}
	}


	////**********************************//
	////            有效不等式   version 2        //
	////**********************************//

	BoolVarMatrix wbF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		wbF[k] = IloBoolVarArray(env, nbTask);
	}

	//only one detour task in each bay
	IloRangeArray  v6(env);
	for (k = 0; k < nbCrane; k++)
		for (int b = 0; b < nbBay; b++)
		{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			if (nbLocation[i] - 1 == b)
			{
			epa += wbF[k][i];
			}
		v6.add(epa <= 1);
		epa.end();
		}
	model.add(v6);
	v6.end();


	for (k = 0; k < nbCrane; k++)
		for (i = 1; i < nbTask; i++)
			model.add(wbF[k][i] - yF[k][i] <= 0);

	IloRangeArray  v7(env);
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		IloExpr  epa(env);
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] > nbLocation[i])
				epa += xF[k][j][i];

		//for (j = 0; j < nbTask; j++)
		//	if (nbLocation[j] < nbLocation[i])
		//		epa += xF[k][i][j];

		epa -= wbF[k][i];

		v7.add(epa == 0);
		epa.end();
		}
	model.add(v7);
	v7.end();


#endif


	return TRUE;




}

BOOL CBMP(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloBoolArray  *nbprecR, IloIntArray  nbLocation, IloIntArray  nbreadyT,
	BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix zF, IloNumVar CF, IloNumVarArray QC_CF, IloNumVarArray Task_CF,
	BoolVarMatrix xF_begin, BoolVarMatrix xF_end,
	IloNum *ObjVal)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;

	IloInt nbM = 1800;



	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	//**********************************//
	//            主问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj(env);
	//  建立子问题目标函数表达式 

	obj += CF;

	//for (k = 0; k < nbCrane; k++)
	//	obj += (float)(1 / nbCrane) * QC_CF[k];

	//	将目标函数加入到主问题模型
	model.add(IloMinimize(env, obj));//

	obj.end();

	IloNumArray ctime_t[nbTask];// the travelling time of QC moving from the bay i to the bay j
	for (i = 0; i < nbTask; i++)
		ctime_t[i] = IloNumArray(env, nbTask);


	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i] == nbLocation[j]) ctime_t[i][j] = 0;
			else if (nbLocation[i] > nbLocation[j])
			{
				ctime_t[i][j] = (nbLocation[i] - nbLocation[j])*crane_move_time;
			}
			else
			{
				ctime_t[i][j] = (nbLocation[j] - nbLocation[i])*crane_move_time;
			}
		}

	}

	//**********************************//
	//            子问题 约束           //
	//**********************************//
	//for (i = 0; i < nbTask; i++)
	//	cout << nbLocation[i] << "  ";
	//QC travel limits
	IloRangeArray  v00(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation[i] < 2 * k + 1)
				epa += yF[k][i];
			if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
				epa += yF[k][i];
		}
		v00.add(epa <= 0);
		epa.end();
	}
	model.add(v00);
	v00.end();

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{

		IloExpr  epa(env);
		epa += CF;
		epa -= QC_CF[k];

		c2.add(epa >= 0);
		epa.end();

	}
	model.add(c2);
	c2.end();

	//约束（6）//each task scheduled	
	IloRangeArray  c60(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
		{
			epa += yF[k][i];
		}
		c60.add(epa == 1);
		epa.end();
	}


	model.add(c60);
	c60.end();

	//约束（2）//identify the start task of each QC
	IloRangeArray  c3(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF_begin[k][i];
		c3.add(epa == 1);
		epa.end();
	}
	model.add(c3);
	c3.end();

	//约束（3）//identify the end task of each QC
	IloRangeArray  c4(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF_end[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	model.add(c4);
	c4.end();


	//约束（4）//task allocation
	IloRangeArray  c6(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa += yF[k][i] - xF_end[k][i];
		for (j = 0; j < nbTask; j++)
		{
			if (i != j)
				epa -= xF[k][i][j];
		}


		c6.add(epa == 0);
		epa.end();
		}
	model.add(c6);
	c6.end();

	//约束（5）//task allocation
	IloRangeArray  c5(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa += yF[k][i] - xF_begin[k][i];
		for (j = 0; j < nbTask; j++)
		{
			if (i != j)
				epa -= xF[k][j][i];
		}
		c5.add(epa == 0);
		epa.end();
		}
	model.add(c5);
	c5.end();



	IloRangeArray  c7(env);
		for (i = 0; i < nbTask; i++)
		{
		for (j = 0; j < nbTask; j++)
		{
			IloExpr  epa(env);
			epa += Task_CF[i] + ctime_t[i][j] + nbQ[j] - Task_CF[j];
			epa -= nbM;
			for (k = 0; k < nbCrane; k++)
				epa += nbM*xF[k][i][j];

			c7.add(epa <= 0);
			epa.end();
		}
		}
	model.add(c7);
	c7.end();

	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			{
				IloExpr  epa(env);
				epa += Task_CF[i] + nbQ[j] - Task_CF[j];
				epa -= nbM*(1 - zF[i][j]);
				c9.add(epa <= 0);
				epa.end();
			}
		}
	}
	model.add(c9);
	c9.end();

	IloRangeArray  c8(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
			IloExpr  epa(env);
			epa += Task_CF[i] + nbQ[j] - Task_CF[j] ;
			c8.add(epa <= 0);
			epa.end();
			}
	}
	model.add(c8);
	c8.end();

	IloRangeArray  c11(env);
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i]<nbLocation[j])
			{
				IloExpr  epa(env);
				for (int kk = 0; kk <= k; kk++)
				{
					epa += xF_begin[kk][j];
					for (int ii = 0; ii < nbTask; ii++)
						epa += xF[kk][ii][j];
				}

				for (int kk = 0; kk <= k; kk++)
				{
					epa -= xF_begin[kk][i];
					for (int ii = 0; ii < nbTask; ii++)
						epa -= xF[kk][ii][i];
				}

				epa -= nbM*(zF[i][j] + zF[j][i]);
				c11.add(epa <= 0);
				epa.end();
			}
		}
		}
	model.add(c11);
	c11.end();

	IloRangeArray  v14(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (i != j)
			{

			IloExpr  epa(env);
			epa += zF[i][j] + zF[j][i];

			v14.add(epa <= 1);
			epa.end();

			}
	}
	model.add(v14);
	v14.end();


	IloRangeArray  c10(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbprecR[i][j] == 1)
			{
				IloExpr  epa(env);
				epa += zF[i][j];

				c10.add(epa == 1);
				epa.end();
			}
		}
	}
	model.add(c10);
	c10.end();

	IloRangeArray  c14(env);//i<j, and both located in the same bay, then i should be processed before j
	for (i = 0; i < nbTask; i++)
	{
		for (j = i + 1; j < nbTask; j++)
		{
			if (nbLocation[i] == nbLocation[j])
			{
				c14.add(zF[i][j] + zF[j][i] == 1);
			}
		}
	}
	model.add(c14);
	c14.end();


	IloRangeArray  c010(env);
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i]<nbLocation[j])
			{
				IloExpr  epa(env);
				for (int kk = 0; kk <= k; kk++)
					epa += yF[kk][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa += yF[kk][i];
				epa -= zF[i][j] + zF[j][i];
				c010.add(epa <= 1);
				epa.end();
			}
		}
		}
	model.add(c010);
	c010.end();


	IloRangeArray  c15(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
		{
			IloExpr  epa(env);
			epa += nbQ[j] - Task_CF[j];
			epa += nbreadyT[k];//QC ready time, already defined in the main function
			if (nbb[k] >= nbLocation[j])
				epa += nbb[k] - nbLocation[j];
			else
				epa += nbLocation[j] - nbb[k] ;
			epa -= nbM*(1 - xF_begin[k][j]);

			c15.add(epa <= 0);
			epa.end();
		}

	}
	model.add(c15);
	c15.end();

	IloRangeArray  c16(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
		{
			IloExpr  epa(env);
			epa += Task_CF[j] - QC_CF[k];
			epa -= nbM*(1 - xF_end[k][j]);

			c16.add(epa <= 0);
			epa.end();
		}

	}
	model.add(c16);
	c16.end();




	//**********************************//
	//            测试 约束           //
	//**********************************//
	//IloExpr ddd(env);
	//for (k = 0; k < nbCrane; k++)
	//	for (i = 0; i < nbTask-1; i++)
	//		for (j = i+1; j < nbTask; j++)
	//				ddd += xF[k][j][i];

	//model.add(ddd == nbCrane);//
#ifdef new_constraints

	BoolVarMatrix wbF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		wbF[k] = IloBoolVarArray(env, nbTask);
	}

	//only one detour task in each bay
	IloRangeArray  v6(env);
	for (k = 0; k < nbCrane; k++)
		for (int b = 0; b < nbBay; b++)
		{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			if (nbLocation[i] - 1 == b)
			{
			epa += wbF[k][i];
			}
		v6.add(epa <= 1);
		epa.end();
		}
	model.add(v6);
	v6.end();

	for (k = 0; k < nbCrane; k++)
	{
		for (i = 1; i < nbTask; i++)
		{
			IloExpr  epa(env);
			for (j = 0; j < i; j++)
				epa += xF[k][i][j];
			epa -= wbF[k][i];
			model.add(epa <= 0);
		}
	}

	for (k = 0; k < nbCrane; k++)
		for (i = 1; i < nbTask; i++)
			model.add(wbF[k][i] - yF[k][i] <= 0);

	//IloRangeArray  v7(env);
	//for (k = 0; k < nbCrane; k++)
	//	for (i = 0; i < nbTask; i++)
	//	{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] > nbLocation[j])
	//			epa += xF[k][j][i];
	//	epa -= wbF[k][i];

	//	v7.add(epa == 0);
	//	epa.end();
	//	}
	//model.add(v7);
	//v7.end();
#endif

	return TRUE;




}

#include <ilcplex/ilocplex.h>
#include  <fstream>
#include  <iostream>
#include <vector>


///////////////*****调试记录********////////////

////////7-12：lift的时候，把Si0和SjT的都改变了，结果后面又用了原始的Si0和SjT，所以出错

///////////////*****调试记录********////////////
#define _AFXDLL

#include <time.h>
#include "stdafx.h"
#include <iomanip>
#include<vector>
#include <random> 

// 唯一的应用程序对象

//CWinApp theApp;

using namespace std;

#define nbTask 50
#define nbBay  15
#define nbCrane 4

#define QCmove_time 1 //t
#define safe_margin 1 //delta
#define multiple_1 1
#define solmodMP 0 //是否利用优先级限定，进行近似求解,0-最优，1-bid
#define nbM 2000//一个很大的整数，可以当初始上界
#define nbM2 2000//
#define threshold_M 0.9//发现cplex lazycallback有时候取了分数解，消除这一影响
#define DEBUGiter
#define computime_limit 18000 //in seconds

#define start_instance 1//第一个算例的编号
#define end_instance 10//测试的算例数量，算例的最大编号

//#define random_generate_data
//#define random_lowerlimit 25
//#define random_upperlimit 100

const IloBool ORIGINAL = 1;

IloInt SubCutNum;
IloNum ObjVal;
IloNum LB;
IloNum LB0;
IloNum UB;
IloNum GapF;

IloInt cbcut_count = 0;// cb cut
IloInt cbcut_single_count = 0;// cb cut
IloInt cbcut_subset_count = 0;// cb cut
IloInt previo_cut_count = 0;
IloInt subtourcut_count = 0;
IloInt SP_count = 0;//子问题个数
IloInt opt_cut_count = 0;
IloInt inf_cut_count = 0;
IloNum GapF2;
IloInt cb2016_count = 0;//combinatorial cut by Sampaio et al. (2016)

//为了子问题求解定义的全局变量
int  nbreadyT[nbCrane];
int  dT[nbCrane];
int nbs = 1;
int nbLocation_new[nbTask];
int nbQ_new[nbTask];
int nbb_new[nbCrane];
int nbprecR[nbTask][nbTask];//优先级关系
int travel_time[nbTask][nbTask];
//////********Delta_ijuv*******///////
int para_Delta[nbCrane][nbCrane][nbTask][nbTask];
int QCno_sideA = nbCrane;// the number of QCs located on the berth side A


#define subP_solve
#define BenchMark_instance
#define QCinitial_leftmost

//-------------------------------------------------------------
// 以下控制各种加速策略的使用
//--------------------------------------------------------------
#define move_inequality
#define lb_inequality
#define Interference_inequality
#define prec_inequality
#define UB_ineq//Precedence-based makespan inequalities
#define opt_property_fix

//no good cuts
#define ub_cbcut//ub based no good cb cut
#define ub_cbcut_domin // dominating cut with Si0 and SiT
#define prec_cut //precedence based no good cb cuts
#define Single_QC_prec_cut// dominating cuts for Si0 or SiT only
#define QC_collision_cut
#define preccut_tripleQCs
////multiple cuts
#define opt_property_cut
#define alter_cut_substitution
//#define alternative_cut //alternative expression of no good cb cuts
#define subset_cb_cut //subset based dominating cut
#define activ_subsetcut_0// dominating subset cut for S0i or SjT, involving variable y
#define activ_subsetcut_1// reduce the tasks in S0i or SjT, 

//////infeasibility cuts
#define subtour_cut
#define subtour_lift_ST
#define subtour_alter//S0T subtour
#define inf_route_prec_cut 1 //16年文献的cut（27），采用为1，不采用为0
#define prec_viocut// new precedence violated cut involving S0i and S0j
#define pre_vio_singleQC//precedence violation cut 27, occurred in the sequence of single QC
#define safe_margin_viocut
#define constructive_cut

//heuristic rules
#define time_alloc_heuristic
#define pre_heuristic
//#define BidSolve

////abandoned strategies
#define preced_new
//#define bay_subtour_ineq
//#define activ_subsetcut_2// reduce the tasks in both S0i and SjT
//#define uuF_sequence
//#define LazyCallBack // branch and check via callback function
//#define opt_property_fix2
//#define indented_berth
//#define insert_subtour

typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<IloBoolArray>    BoolMatrix;
typedef IloArray<IloArray<IloBoolArray>> BoolMatrix2;
typedef IloArray<IloArray<IloNumArray>> NumMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolArray>>> BoolMatrix3;
typedef IloArray<IloArray<IloArray<IloNumArray>>> NumMatrix3;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloBoolVarArray>    BoolVarMatrix;
typedef IloArray<IloArray<IloBoolVarArray> > BoolVarMatrix2;
typedef IloArray<IloArray<IloNumVarArray> > NumVarMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolVarArray> > > BoolVarMatrix3;

ILOSTLBEGIN

bool HeuristicFun(IloModel model, IloNum nbs, IloIntArray  nbb, IloIntArray  nbreadyT2,
	BoolVarMatrix yF, IloIntVar CF, IloInt* CF_best, IloNumArray2 yF_best,
	IloNum* ObjVal, IloNum* UB, IloIntArray  nbLocation, int my);

bool Check_HeuristicFun(IloModel model, IloNum nbs, IloIntArray  nbb, IloIntArray  nbreadyT, IloIntVar CF, IloNumArray2 yF_best,
	IloNum* ObjVal, IloIntArray  nbLocation, int my);


bool CBMP_R(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation,
	BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, IloNumVarArray uuF,
	BoolVarMatrix xF_begin, BoolVarMatrix xF_end, IloInt* CF_best,
	IloNum* ObjVal, int sol_mode);

bool Bid_SP_1(IloCplex cplex, IloIntVar CF, IloNumArray2 yF_best, IloNumArray2 yF_current2, IloNumArray2 yF2_current, IloNum* ObjVal);

bool AddCBCuts(IloEnv env, IloRangeArray cutPool, BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix xF_begin, BoolVarMatrix xF_end, int task_no[],
	IloNumArray2 yF_best, int record_next_task[][nbTask], int record_preceding_task[][nbTask], int first_task[], int last_task[]);

IloInt CSP(IloCplex subcplex, IloNumArray3  xF_best, IloNumArray2 yF_best, IloNumArray2 xF_begin_best, IloNumArray2 xF_end_best);

bool CBMP_bid(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation,
	IloNum* ObjVal, IloNumArray2 yF2_current, IloNumArray2 yF_current2);

bool LB_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, 
	BoolVarMatrix yF, IloIntVar CF, IloInt* CF_best, IloNumArray2 yF_best,
	IloNum* ObjVal, IloNum* UB_Fun, IloIntArray  nbLocation, int my);

int _tmain(int argc, char* argv[], char* envp[])
{
	double GapAve[end_instance - start_instance + 1];
	double GapAve2[end_instance - start_instance + 1];
	double DuraAve[end_instance - start_instance + 1];
	double DuraAve_MP[end_instance - start_instance + 1];
	double DuraAve_SP[end_instance - start_instance + 1];
	double ObjAve[end_instance - start_instance + 1];
	double LBAve[end_instance - start_instance + 1];
	double IterAve[end_instance - start_instance + 1];
	double ObjUniQCSP[end_instance - start_instance + 1];
	double ObjbiQCSP[end_instance - start_instance + 1];
	double LB0_BD[end_instance - start_instance + 1];
	int preccut_countAve[end_instance - start_instance + 1];
	int subcut_countAve[end_instance - start_instance + 1];
	int preccut_single_countAve[end_instance - start_instance + 1];
	int preccut_subset_countAve[end_instance - start_instance + 1];
	int preccut_count_2016Ave[end_instance - start_instance + 1];
	int opt_cut_countAve[end_instance - start_instance + 1];
	int inf_cut_countAve[end_instance - start_instance + 1];
	int SP_countAve[end_instance - start_instance + 1];
	int nodesNO[end_instance - start_instance + 1];

	int LB_initial[end_instance - start_instance + 1];
	int IterNo[end_instance - start_instance + 1];

#ifdef indented_berth
	cout << "indented berth, totally " << nbCrane << " QCs, please type in the number of QCs located on the A side:";
	cin >> QCno_sideA;
#endif

	for (int my = start_instance; my <= end_instance; my++)
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
		clock_t start_TA = 0, finish_TA = 0;
		start = clock();
		double  duration;
		double  duration_MP = 0;
		double  duration_SP = 0;
		double  duration_SP2 = 0;
		double  duration_TA = 0;

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

			//	定义原问题模型
			IloModel model(env);

			//	定义松弛问题模型
			IloModel RLmodel(env);

			IloNum gap;

			//	定义临时变量
			IloInt i, k, j, t, kk;

			//	定义连续决策变量CF
			IloIntVar   CF(env, 0, 1800);

			//IloNumArray tem_C(env, nbCrane);// 子问题计算的completion time
			IloNumVarArray QC_CF(env, nbCrane, 0, IloInfinity);// each QC's completion time
			IloNumVarArray Task_CF(env, nbTask, 0, IloInfinity);// each task's completion time
			IloNumVarArray uuF(env, nbTask, 0, nbTask);//auxiliary variable, for subtours eliminating
			IloNumVarArray QC_maxWait(env, nbCrane, 0, IloInfinity);//auxiliary variable, for cuts adding

			//	定义决策变量CxF,CyF,CzF
			BoolVarMatrix yF(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				yF[k] = IloBoolVarArray(env, nbTask);
			}

			BoolVarMatrix2 xF(env);
			for (k = 0; k < nbCrane; k++)
			{
				xF.add(BoolVarMatrix(env));
				for (i = 0; i < nbTask; i++)
				{
					xF[k].add(IloBoolVarArray(env, nbTask));
				}
			}
			BoolVarMatrix zF(env, nbTask);
			for (i = 0; i < nbTask; i++)
			{
				zF[i] = IloBoolVarArray(env, nbTask);
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
			//**********************************//
			//            新变量           //
			//**********************************//

			//	定义保存决策变量CxF,CyF,CzF的最优值
			//BoolMatrix yF_best(env, nbCrane);
			//for (k = 0; k < nbCrane; k++)
			//{
			//	yF_best[k] = IloBoolArray(env, nbTask);
			//}

			//BoolMatrix2 xF_best(env);
			//for (k = 0; k < nbCrane; k++)
			//{
			//	xF_best.add(BoolMatrix(env));
			//	for (i = 0; i < nbTask; i++)
			//		xF_best[k].add(IloBoolArray(env, nbTask));
			//}
			BoolMatrix zF_best(env, nbTask);
			for (i = 0; i < nbTask; i++)
			{
				zF_best[i] = IloBoolArray(env, nbTask);
			}

			//BoolMatrix xF_begin_best(env, nbCrane);
			//for (k = 0; k < nbCrane; k++)
			//{
			//	xF_begin_best[k] = IloBoolArray(env, nbTask);
			//}

			//BoolMatrix xF_end_best(env, nbCrane);
			//for (k = 0; k < nbCrane; k++)
			//{
			//	xF_end_best[k] = IloBoolArray(env, nbTask);
			//}

			// recording the current best solution for the master problem at each iteration
			IloNumArray2 yF_best(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				yF_best[i] = IloNumArray(env, nbTask);
			}
			IloNumArray2 xF_begin_best(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				xF_begin_best[i] = IloNumArray(env, nbTask);
			}
			IloNumArray2 xF_end_best(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				xF_end_best[i] = IloNumArray(env, nbTask);
			}
			IloNumArray3 xF_best(env, nbCrane);
			for (int k = 0; k < nbCrane; ++k) {
				xF_best[k] = IloNumArray2(env, nbTask);
				for (i = 0; i < nbTask; ++i) {
					xF_best[k][i] = IloNumArray(env, nbTask);
				}
			}
			IloNumArray2 yF_current2(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				yF_current2[i] = IloNumArray(env, nbTask);
			}
			IloNumArray2 yF2_current(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				yF2_current[i] = IloNumArray(env, nbTask);
			}

			// recording optimal solution for the whole problem
			IloNumArray2 yF_opt(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				yF_opt[i] = IloNumArray(env, nbTask);
			}
			IloNumArray2 xF_begin_opt(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				xF_begin_opt[i] = IloNumArray(env, nbTask);
			}
			IloNumArray2 xF_end_opt(env, nbCrane);
			for (i = 0; i < nbCrane; ++i) {
				xF_end_opt[i] = IloNumArray(env, nbTask);
			}
			IloNumArray3 xF_opt(env, nbCrane);
			for (int k = 0; k < nbCrane; ++k) {
				xF_opt[k] = IloNumArray2(env, nbTask);
				for (i = 0; i < nbTask; ++i) {
					xF_opt[k][i] = IloNumArray(env, nbTask);
				}
			}

			IloInt   CF_best;
			IloNumArray QC_CF_best(env, nbCrane);
			IloNumArray Task_CF_best(env, nbTask);
			//*******************************************************//
			//                 定义主问题输出的参数                         //
			//*******************************************************// 	 

			//约束中参数
			//	定义需代入子问题计算的work zone数目
			IloInt   nbn;
			//定义 work zone 起始 section
			IloIntArray nbn_start(env, nbBay);
			//	定义 work zone 终止 section
			IloIntArray  nbn_end(env, nbBay);
			//	定义每个work zone的crane个数
			IloIntArray  nbn_crane_no(env, nbBay);
			//	定义每个work zone的最后一个crane的编号
			IloIntArray  nbn_crane_index(env, nbBay);
			//	定义每个work zone的第一个crane的编号
			IloIntArray  nbn_crane_start(env, nbBay);



			IloInt cut_1_no = 0;
			IloInt cut_2_no = 0;





			//*******************************************************//
			//                 定义输入参数                          //
			//*******************************************************// 	 
			//约束中参数
			//	定义能力常量参数s
			//IloNum   nbs;
			//定义吊机位置初始状态参数
			IloBoolArray nbb(env, nbCrane);
			//	定义任务量参数nbQ
			IloNumArray  nbQ(env, nbTask);

			IloIntArray  nbreadyT2(env, nbCrane);
			//IloIntArray  dT(env, nbCrane);

			//	定义任务所在贝位参数nbQ
			IloIntArray  nbLocation(env, nbTask);

			//优先级关系初始化
			for (int ai = 0; ai < nbTask; ai++)
				for (int bi = 0; bi < nbTask; bi++)
					nbprecR[ai][bi] = 0;

			//*******************************************************//
			//                 定义callback参数                      //
			//*******************************************************// 	 



			//*******************************************************//
			//                 读入参数数据                          //
			//*******************************************************// 
			// 		

			IloIntArray ls1(env, 8);

			fin >> ls1;
			//cout << ls1[0] << "  ";			
			fin >> nbQ;
			fin >> nbLocation;
			fin >> nbreadyT2;
			fin >> nbb;

			//nbb[2] = 5;
			//nbb[2] = 10;
			//nbreadyT2[1] = 300;


			//for (i = 0; i < nbTask; i++) cout<<nbQ[i]<<"  ";
			//cout<<endl;

			////读入数据nbLocation[i]
			//for (i = 0; i < nbTask; i++)
			//{
			//	//fin >> nbLocation[i];
			//	cout<<nbLocation[i]<<"  ";
			//}
			////cout<<endl;

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

#ifdef	random_generate_data
			for (int ai = 0; ai < nbTask; ai++)
			{
				if (nbQ[ai] < random_lowerlimit || nbQ[ai]>random_upperlimit)
				{

					// 创建一个随机数引擎，以当前时间作为种子
					std::random_device rd;  // 用于获取随机种子
					std::mt19937 gen(rd()); // 以随机设备rd的值为种子初始化Mersenne Twister生成器

					// 定义一个分布，例如在40到99之间均匀分布的整数
					std::uniform_int_distribution<> dis(random_lowerlimit, random_upperlimit);

					// 使用分布和引擎生成随机数
					nbQ[ai] = dis(gen);
				}
			}
#endif

			for (int ai = 0; ai < nbTask; ai++)
			{
				nbLocation_new[ai] = nbLocation[ai];
				nbQ_new[ai] = nbQ[ai];
				//cout << nbLocation_new[ai] << ", ";
			}
			//cout << endl;
			for (k = 0; k < nbCrane; k++)
				nbb_new[k] = nbb[k];

			for (k = 0; k < nbCrane; k++)
				nbreadyT[k] = nbreadyT2[k];

			for (k = 0; k < nbCrane; k++)
				dT[k] = 5000;





			nbs = 1;



			//////********Delta_ijuv*******///////
			for (k = 0; k < nbCrane; k++)
				for (int k2 = 0; k2 < nbCrane; k2++)
					for (i = 0; i < nbTask; i++)
						for (int j = 0; j < nbTask; j++)
							para_Delta[k][k2][i][j] = -1;
			for (k = 0; k < nbCrane; k++)
			{
				for (int k2 = k; k2 < nbCrane; k2++)
				{
					for (i = 0; i < nbTask - 1; i++)
					{
						for (int j = i + 1; j < nbTask; j++)
						{
#ifdef indented_berth
							if ((nbLocation_new[i] >= nbLocation_new[j] - (k2 - k) * (1 + safe_margin) + 1)
								&& (k2 < QCno_sideA || k >= QCno_sideA))
#endif							
							{

								if (nbLocation_new[i] >= nbLocation_new[j] - (k2 - k) * (1 + safe_margin) + 1)

								{
									para_Delta[k][k2][i][j] = QCmove_time * (k2 - k) * (1 + safe_margin) - QCmove_time * (nbLocation_new[j] - nbLocation_new[i]);
									para_Delta[k2][k][j][i] = QCmove_time * (k2 - k) * (1 + safe_margin) - QCmove_time * (nbLocation_new[j] - nbLocation_new[i]);
								}
								if (k != k2)
								{
									para_Delta[k2][k][i][j] = QCmove_time * (k2 - k) * (1 + safe_margin) + QCmove_time * (nbLocation_new[j] - nbLocation_new[i]);
									para_Delta[k][k2][j][i] = QCmove_time * (k2 - k) * (1 + safe_margin) + QCmove_time * (nbLocation_new[j] - nbLocation_new[i]);
								}
							}
#ifdef indented_berth
							if (k < QCno_sideA && k2 >= QCno_sideA)
							{
								if (nbLocation_new[i] >= nbLocation_new[j] -  safe_margin)
								{
									para_Delta[k][k2][i][j] = para_Delta[k2][k][j][i] = para_Delta[k2][k][j][i] = para_Delta[k2][k][j][i] = QCmove_time * (1 + safe_margin) - QCmove_time * (nbLocation_new[j] - nbLocation_new[i]);
								}
							}
#endif		
							//cout << k + 1 << "-" << i + 1 << ", " << k2 + 1 << "-" << j + 1 << ": " << para_Delta[k][k2][i][j]<<";   " << k2 + 1 << "-" << i + 1 << ", " << k + 1 << "-" << j + 1 << ": " << para_Delta[k2][k][i][j] <<endl;

						}
					}
				}
			}

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

			//cout<<filename1<<endl;


			////展示task位置
			for (int b = 0; b < nbBay; b++)
			{
				cout << "bay " << b + 1 << ": ";
				fout << "bay " << b + 1 << ": ";
				for (i = 0; i < nbTask; i++)
				{
					if (nbLocation[i] == b + 1)
					{
						cout << i + 1 << "(" << nbQ_new[i] << ")  ";
						fout << i + 1 << "(" << nbQ_new[i] << ")  ";
					}

				}
				cout << endl;
				fout << endl;
			}


			//*******************************************************//
			//                 Unidirectional                       //
			//*******************************************************// 

			nbn = 0;//统计子问题个数
			UB = 3000;
			LB = 0;
			IloNum UB_temp1, UB_temp2;
			UB_temp1 = UB_temp2 = 3000;

			IloInt CF_ub;
			IloNum ObjVal2;

			int iterNo = 0;

#ifdef	pre_heuristic

			IloNum  tempbound = nbM;

			IloNum LB0;
			LB_Same_Direction(model, nbs, nbb, nbQ,
				yF, CF, &CF_best, yF_best,
				&ObjVal, &tempbound, nbLocation, my);

			LB0 = ObjVal;
			LB = LB0;


			HeuristicFun(model, nbs, nbb, nbreadyT2,
				yF, CF, &CF_ub, yF_best,
				&ObjVal2, &UB_temp1, nbLocation, my);

			//更新目前最优调度方案
			for (int k = 0; k < nbCrane; k++)
			{
				for (i = 0; i < nbTask; i++)
				{
					yF_opt[k][i] = xF_begin_opt[k][i] = xF_end_opt[k][i] = 0;
					for (j = 0; j < nbTask; j++)
						xF_opt[k][i][j] = 0;
				}

				int current_i = -1;
				for (i = 0; i < nbTask; i++)
				{
					yF_opt[k][i] = yF_best[k][i];
					if (yF_best[k][i] == 1)
					{
						if (current_i == -1)
							xF_begin_opt[k][i] = 1;
						else
							xF_opt[k][current_i][i] = 1;
						current_i = i;
					}
				}
			}

			UB = UB_temp1;



			//re-numbering the bays
			IloBoolArray nbb2(env, nbCrane);
			IloIntArray  nbLocation2(env, nbTask);
			IloIntArray  nbreadyT3(env, nbCrane);
			for (j = 0; j < nbTask; j++)
			{
				nbLocation2[j] = nbBay + 1 - nbLocation[j];
			}
			for (k = 0; k < nbCrane; k++)
				nbb2[nbCrane - 1 - k] = nbBay + 1 - nbb[k];
			for (k = 0; k < nbCrane; k++)
				nbreadyT3[nbCrane - 1 - k] = nbreadyT[k];

			HeuristicFun(model, nbs, nbb2, nbreadyT3,
				yF, CF, &CF_ub, yF_best,
				&ObjVal2, &UB_temp1, nbLocation2, my);

			//更新目前最优调度方案
			if (UB > UB_temp1)
			{
				for (int k = 0; k < nbCrane; k++)
				{
					for (i = 0; i < nbTask; i++)
					{
						yF_opt[k][i] = xF_begin_opt[k][i] = xF_end_opt[k][i] = 0;
						for (j = 0; j < nbTask; j++)
							xF_opt[k][i][j] = 0;
					}

					int current_i = -1;
					for (i = 0; i < nbTask; i++)
					{
						yF_opt[nbCrane - 1 - k][i] = yF_best[k][i];
						if (yF_best[k][i] == 1)
						{
							if (current_i == -1)
								xF_begin_opt[nbCrane - 1 - k][i] = 1;
							else
								xF_opt[nbCrane - 1 - k][current_i][i] = 1;
							current_i = i;
						}
					}
				}

				UB = UB_temp1;
			}
			fout << endl << "initial lower bound=" << LB << endl;
			fout << endl << "unidirectional makespan=" << UB << endl;
			ObjUniQCSP[my - start_instance] = UB;

			if (LB == UB)
			{
				fout << "early termination." << endl;
				goto out_end;
			}

#ifdef BidSolve
			if(nbCrane<=4)
			{
				IloModel bidMod(env);
				bool indibid;
				indibid = CBMP_bid(bidMod, nbs, nbb, nbQ, nbLocation,
					&ObjVal, yF2_current, yF_current2);
				if (indibid && ObjVal < UB)
				{
					UB = ObjVal;
					//更新目前最优调度方案，待写
					for (int k = 0; k < nbCrane; k++)
					{
						for (i = 0; i < nbTask; i++)
						{
							xF_begin_opt[k][i] = xF_end_opt[k][i] = 0;
							for (j = 0; j < nbTask; j++)
								xF_opt[k][i][j] = 0;
						}

						int current_i = -1;
						for (i = 0; i < nbTask; i++)
						{
							yF_opt[k][i] = yF_best[k][i];
							if (yF_current2[k][i] == 1)
							{
								if (current_i == -1)
									xF_begin_opt[k][i] = 1;
								else
									xF_opt[k][current_i][i] = 1;
								current_i = i;
							}
						}

						for (i = nbTask - 1; i >= 0; i--)
						{
							if (yF2_current[k][i] == 1)
							{
								if (current_i == -1)
									xF_begin_opt[k][i] = 1;
								else
									xF_opt[k][current_i][i] = 1;
								current_i = i;
							}
						}
					}
				}
				//fout << endl << "bidirectional makespan=" << UB << endl;
			}
			ObjbiQCSP[my - start_instance] = UB;
#endif
#endif



			//*******************************************************//
			//                 求解主问题 CBMP                       //
			//*******************************************************// 

			start = clock();

			//主问题模型
			CBMP_R(model, nbs, nbb, nbQ, nbLocation,
				xF, yF, zF, CF, QC_CF, uuF,
				xF_begin, xF_end, &CF_best,
				&ObjVal, solmodMP);




			//**********************************//
			//            开始求解		        //
			//**********************************//
			IloCplex cplex(env);
			cplex.extract(model);

			//IloEnv subEnv;
			IloCplex subcplex(env);

			//#ifdef UserActive
			//	cplex.use(CandyUserCallback(env));
			//	//cplex.setParam(IloCplex::MIPEmphasis, 2);
			//	//cplex.setParam(IloCplex::HeurFreq, 1);
			//	//cplex.setParam(IloCplex::ParallelMode, 1);
			//	cplex.setParam(IloCplex::Threads, 4);
			//#endif

				//cplex.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_FEASIBILITY);
			cplex.setParam(IloCplex::EpGap, 0.0006);
			//cplex.setOut(env.getNullStream());
			cplex.setParam(IloCplex::TiLim, computime_limit);
			//cplex.setParam(IloCplex::Threads, 4);

			//cplex.setParam(IloCplex::VarSel, 3);
			//cplex.setParam(IloCplex::MIPEmphasis, 3);
			//cplex.setParam(IloCplex::Probe, 3);

			//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
			//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;


#ifdef LazyCallBack
			cplex.use(BendersLazyCallback(env, xF, yF, xF_begin, xF_end, CF, cplex, subcplex));
			//cplex.use(IncumbentCallback(env, xF, yF, xF_begin, xF_end, CF, cplex, subcplex));
			//cplex.use(UserCallback(env, xF, yF, xF_begin, xF_end, CF, nbprecR, subcplex));
#endif

			fout << "算法开始：" << endl;



			//**********************************//
			//             Benders 主算法           //
			//**********************************//	
#ifdef uuF_sequence
			int checkY[nbCrane][nbTask];
			for (k = 0; k < nbCrane; k++)
				for (i = 0; i < nbTask; i++)
					checkY[k][i] = 0;
#endif

			gap = 100;
			duration_MP = 0;
			duration = 0;

#ifndef LazyCallBack
			while (LB < UB - 0.8 && duration < computime_limit)
			//while (LB < UB - 0.8 && duration < computime_limit && LB <= ObjUniQCSP[my - start_instance] - 1 - safe_margin)
#endif
			//while (iterNo <= 0 )
			{

				IloInt cbcut_count_innerloop = 0;// cb cut
				IloInt cb2016_count_innerloop = 0;
				IloInt preccut_count_2016_innerloop = 0;
				IloInt S0_cbcut_single_count_innerloop = 0;
				IloInt ST_cbcut_single_count_innerloop = 0;

				start_MP = clock();
				bool h0;

				iterNo++;
				fout << endl << iterNo << "- iteration starts. " << endl;
				//cplex.exportModel("MP.lp");
				h0 = cplex.solve();

				//fout<<"h0= "<<h0<<endl;

				finish_MP = clock();

				if (!h0)
				{
					//fout << iterNo << "- iteration MP no feasible solution. " << endl;

					cout << "\nno feasible solution has been found" << endl;
					fout << "\nno feasible solution has been found   " << iterNo << endl;
					//cplex.clearModel();
					//cplex.clear();
					//cplex.end();
					//model.end();


					finish = clock();
					duration = (double)(finish - start) / CLOCKS_PER_SEC;

					//return false;
					goto out_end;
				}


				//**********************************//
				//             记录最好解       //
				//**********************************//
				//fout << iterNo << "- MP solved" << endl;
				gap = cplex.getMIPRelativeGap();

				LB = cplex.getBestObjValue();
				CF_best = cplex.getValue(CF);

				if (CF_best - CF_best * gap < LB)
					LB = CF_best - CF_best * gap;


				if (iterNo == 1)
				{
					LB0 = LB;
					LB0_BD[my - start_instance] = LB;
				}

				double mptime;
				mptime = (double)(finish_MP - start_MP) / CLOCKS_PER_SEC;
				duration_MP += (double)(finish_MP - start_MP) / CLOCKS_PER_SEC;

				finish = clock();
				duration = (double)(finish - start) / CLOCKS_PER_SEC;

				fout << iterNo << "- MP value: " << CF_best << ", lower bound: " << LB << ", MP time: " << mptime << ", current time: " << duration << endl;

				//cout<<"yF_best"<<endl;
				for (k = 0; k < nbCrane; k++)
				{
					for (i = 0; i < nbTask; i++)
					{

						if (cplex.isExtracted(yF[k][i]))
						{
							if (cplex.getValue(yF[k][i]) > 0.1)
								yF_best[k][i] = 1;
							else
								yF_best[k][i] = 0;
						}
						else
							yF_best[k][i] = 0;
						//cout << yF_best[k][i] << "  ";

						if (cplex.isExtracted(xF_begin[k][i]))
						{
							if (cplex.getValue(xF_begin[k][i]) > 0.1)
								xF_begin_best[k][i] = 1;
							else
								xF_begin_best[k][i] = 0;
							//cout<<xF_begin_best[k][i]<<"  ";
						}
						else
							xF_begin_best[k][i] = 0;

						if (cplex.isExtracted(xF_end[k][i]))
						{
							if (cplex.getValue(xF_end[k][i]) > 0.1)
								xF_end_best[k][i] = 1;
							else
								xF_end_best[k][i] = 0;
							//cout<<xF_end_best[k][i]<<"  ";
						}
						else
							xF_end_best[k][i] = 0;

						for (int j = 0; j < nbTask; j++)
						{
							if (cplex.isExtracted(xF[k][i][j]))
							{
								if (cplex.getValue(xF[k][i][j]) > 0.1)
									xF_best[k][i][j] = 1;
								else
									xF_best[k][i][j] = 0;
								//cout<<xF_end_best[k][i]<<"  ";
							}
							else
								xF_best[k][i][j] = 0;

						}

					}
					//cout << "    "; cout << endl;
				}

#ifdef uuF_sequence
				for (i = 0; i < nbTask; i++)
				{
					for (j = 0; j < nbTask; j++)
					{

						if (i != j)
						{
							//IloExpr epa(env);
							//for (k = 0; k < nbCrane; k++)
							//	epa += xF[k][i][j];

							//model.add(uuF[j] - uuF[i] +1+ nbTask * (1 - epa) >= 0);


							for (k = 0; k < nbCrane; k++)
							{
								model.add(uuF[j] - uuF[i] - 1 + nbTask * (1 - xF[k][i][j]) >= 0);
							}
						}
					}
				}
				if (LB < UB)
				{
					for (k = 0; k < nbCrane; k++)
					{
						for (i = 0; i < nbTask - 1; i++)
							if (yF_best[k][i] == 1)
							{
								for (j = i + 1; j < nbTask; j++)
									if (yF_best[k][j] == 1)
									{
										if (checkY[k][i] == 0 || checkY[k][j] == 0)
										{
											cplex.addCut(uuF[j] - uuF[i] - 1 + nbTask * (1 - xF[k][i][j]) >= 0);											
										}

									}
							}
					}
				}
				for (k = 0; k < nbCrane; k++)
					for (i = 0; i < nbTask; i++)
						if (yF_best[k][i] == 1)
							checkY[k][i] = 1;
#endif
#ifndef LazyCallBack
				int record_next_task[nbCrane][nbTask];//建立查找字典，记录后续任务，免得后面要一遍遍查找
				int record_preceding_task[nbCrane][nbTask];//建立查找字典，记录前续任务，免得后面要一遍遍查找
				int first_task[nbCrane];//记录头部任务
				int last_task[nbCrane];//记录尾部任务


				for (int k = 0; k < nbCrane; k++)
					for (i = 0; i < nbTask; i++)
					{
						record_next_task[k][i] = -1; //初始化字典
						record_preceding_task[k][i] = -1;
					}

				int task_no[nbCrane];

				//cout << "lazy: start recording the task sequence" << endl;
				for (int k = 0; k < nbCrane; k++)
				{
					//cout << k +1<< ":  ";

					task_no[k] = 0;

					for (i = 0; i < nbTask; i++)
					{

						if (xF_begin_best[k][i] > threshold_M)
							first_task[k] = i;

						if (xF_end_best[k][i] > threshold_M)
							last_task[k] = i;

						if (yF_best[k][i] > threshold_M)
							task_no[k] = task_no[k] + 1;


						for (int j = 0; j < nbTask; j++)
							if (i != j)
							{
								if (xF_best[k][i][j] > threshold_M)
								{
									record_preceding_task[k][j] = i;
									record_next_task[k][i] = j;
									//cout << i+1 << "-" << j +1<< ",  ";
								}
							}

					}
					//cout << endl;
				}
				//cout << endl << endl;

#ifdef DEBUGiter
				for (int k = 0; k < nbCrane; k++)
				{
					int sumQ = 0;
					int sumQ2 = 0;
					fout << "QC  " << k + 1 << ":  ";
					int ci;
					ci = first_task[k];
					//sumQ += nbQ[ci];
					while (ci >= 0)
					{
						fout << ci + 1 << "-";
						ci = record_next_task[k][ci];
						//sumQ += nbQ[ci];
					}
					fout << endl;
					//fout << ", sum=" << sumQ << endl;
					fout << "QC  " << k + 1 << ":  ";
					for (j = 0; j < nbTask; j++)
						if (yF_best[k][j] == 1)
						{
							fout << j + 1 << ", ";
							sumQ2 += nbQ[j];
						}
					fout << ", sum=" << sumQ2 << endl;
				}
#endif // DEBUG

				//for (int k = 0; k < nbCrane; k++)
				//{
				//	cout << k << ":  ";
				//	int ci;
				//	ci = first_task[k];
				//	while (ci >= 0)
				//	{
				//		cout << ci << "-";
				//		ci = record_next_task[k][ci];
				//	}

				//	cout << endl;
				//}


#ifdef time_alloc_heuristic
				if (LB < UB)
				{
					start_TA = clock();
					int uni_prec_count = 0;
					int uni_prec_count2 = 0;
					int indic_uni_opt = 0;

#ifndef indented_berth
					for (k = 0; k < nbCrane; k++)
						if(nbb[k]!=1+(1+safe_margin)*k || nbreadyT[k] != 0)
							indic_uni_opt = 1;
					if (indic_uni_opt == 0)
					{
						for (i = 0; i < nbTask; i++)
						{
							for (j = 0; j < nbTask; j++)
								if (nbprecR[i][j] == 1)
								{
									int ki = -1;
									int kj = -1;
									for (k = 0; k < nbCrane; k++)
									{
										if (yF_best[k][i] == 1)
											ki = k;
										if (yF_best[k][j] == 1)
											kj = k;
									}
									if (ki < kj)
										uni_prec_count++;
									if (ki > kj)
										uni_prec_count2++;

								}
						}
						if (uni_prec_count == 0)
						{
							UB = LB;

							//更新目前最优调度方案
							for (int k = 0; k < nbCrane; k++)
							{
								for (i = 0; i < nbTask; i++)
								{
									yF_opt[k][i] = xF_begin_opt[k][i] = xF_end_opt[k][i] = 0;
									for (j = 0; j < nbTask; j++)
										xF_opt[k][i][j] = 0;
								}

								int current_i = -1;
								for (i = 0; i < nbTask; i++)
								{
									yF_opt[k][i] = yF_best[k][i];
									if (yF_best[k][i] == 1)
									{
										if (current_i == -1)
											xF_begin_opt[k][i] = 1;
										else
											xF_opt[k][current_i][i] = 1;
										current_i = i;
									}
								}
							}							

						}
					}
#endif
#ifdef indented_berth
					uni_prec_count = 1;
#endif
#ifndef BidSolve
					if (indic_uni_opt != 0 || uni_prec_count != 0)//bidirectional
					{
						bool solh;
						IloNum ObjVal_bid;
						solh = Bid_SP_1(subcplex, CF, yF_best, yF_current2, yF2_current, &ObjVal_bid);
						if (solh)
						{
							if (ObjVal_bid < UB)
							{
								UB = ObjVal_bid;
								//更新目前最优调度方案，待写
								for (int k = 0; k < nbCrane; k++)
								{
									for (i = 0; i < nbTask; i++)
									{
										xF_begin_opt[k][i] = xF_end_opt[k][i] = 0;
										for (j = 0; j < nbTask; j++)
											xF_opt[k][i][j] = 0;
									}

									int current_i = -1;
									for (i = 0; i < nbTask; i++)
									{
										yF_opt[k][i] = yF_best[k][i];
										if (yF_current2[k][i] == 1)
										{
											if (current_i == -1)
												xF_begin_opt[k][i] = 1;
											else
												xF_opt[k][current_i][i] = 1;
											current_i = i;
										}
									}

									for (i = nbTask - 1; i >= 0; i--)
									{
										if (yF2_current[k][i] == 1)
										{
											if (current_i == -1)
												xF_begin_opt[k][i] = 1;
											else
												xF_opt[k][current_i][i] = 1;
											current_i = i;
										}
									}
								}

							}
						}
					}
#endif
					finish_TA = clock();
					duration_TA += (double)(finish_TA - start_TA) / CLOCKS_PER_SEC;
				}
#endif

				IloRangeArray cutPool(env);

				bool SP_activate;

				SP_activate = AddCBCuts(env, cutPool, xF, yF, xF_begin, xF_end, task_no,
					yF_best, record_next_task, record_preceding_task, first_task, last_task);

				if (cutPool.getSize() > 0)
					for (int nn = 0; nn < cutPool.getSize(); nn++)
						cplex.addCut(cutPool[nn]);

				// //////////////////////////////////
				////****  solving subproblem  **** ///
				// //////////////////////////////////

				//当没有子环也没有cb cut的时候，可以解子问题
				if (SP_activate)
				{
					SP_count++;

					fout << iterNo << "- SP start. " << endl;
					clock_t start_SP = 0, finish_SP = 0;
					start_SP = clock();

					IloInt current_obj;

					current_obj = CSP(subcplex, xF_best, yF_best, xF_begin_best, xF_end_best);
					////双向模型的子问题，间接测试BD子问题是否正确
					//current_obj = Split_SP_Callback(subcplex, nbprecR, yF_best);

					//cout << "UB =" << UB << ", LB=" << CF_best << ",current UB=" << current_obj << endl;

					finish = clock();
					duration = (double)(finish - start) / CLOCKS_PER_SEC;

					fout << iterNo << "-current UB=" << current_obj << ", current time= " << duration << endl;

					IloExpr optcut(env);
					for (k = 0; k < nbCrane; k++)
						for (i = 0; i < nbTask; i++)
						{
							if (xF_begin_best[k][i] > threshold_M)
								optcut += 1 - xF_begin[k][i];
							if (xF_end_best[k][i] > threshold_M)
								optcut += 1 - xF_end[k][i];

							for (j = 0; j < nbTask; j++)
							{
								if (i != j && xF_best[k][i][j] > threshold_M)//precedence
									optcut += 1 - xF[k][i][j];
							}
						}
					if (current_obj < 2)
					{
						cplex.addCut(optcut >= 1);
						inf_cut_count++;
					}
					else if (current_obj - LB > 0.9)
					{
						cplex.addCut(CF + (current_obj - LB) * optcut >= current_obj);
						opt_cut_count++;
					}
					optcut.end();

					if (current_obj > 2 && current_obj <= UB)
					{
						UB = current_obj;
						for (k = 0; k < nbCrane; k++)
						{
							for (i = 0; i < nbTask; i++)
							{
								xF_begin_opt[k][i] = xF_begin_best[k][i];
								xF_end_opt[k][i] = xF_end_best[k][i];
								yF_opt[k][i] = yF_opt[k][i];
								for (j = 0; j < nbTask; j++)
									xF_opt[k][i][j] = xF_best[k][i][j];
							}
						}
					}
#ifdef Split_subP_solve
					else
					{
						//split subproblem
						if (nbCrane > 2)
						{
							for (int k = 0; k < nbCrane - 1; k++)
							{
								IloInt current_obj2;
								current_obj2 = CSP_Split(subcplex, xF_best, yF_best, xF_begin_best, xF_end_best, k);

								IloExpr optcut2(env);
								for (int k2 = k; k2 <= k + 1; k2++)
									for (i = 0; i < nbTask; i++)
									{
										if (xF_begin_best[k2][i] > threshold_M)
											optcut2 += 1 - xF_begin[k2][i];
										if (xF_end_best[k2][i] > threshold_M)
											optcut2 += 1 - xF_end[k2][i];

										for (j = 0; j < nbTask; j++)
										{
											if (i != j && xF_best[k2][i][j] > threshold_M)//precedence
												optcut2 += 1 - xF[k2][i][j];
										}
									}
								if (current_obj2 < 2)
								{
									cplex.addCut(optcut2 >= 1);
									inf_cut_count++;
								}
								else if (current_obj2 - LB > 0.9)
								{
									cplex.addCut(CF + (current_obj2 - LB) * optcut2 >= current_obj2);
									opt_cut_count++;
								}

								optcut2.end();
							}
						}
					}
#endif





					finish_SP = clock();
					duration_SP += (double)(finish_SP - start_SP) / CLOCKS_PER_SEC;

				}

				cutPool.end();

				//fout << iterNo << "- iteration ends. " << endl;

				//fout.precision(10);
				//fout << "cb cut =" << preccut_count << endl;
				//fout << "cb cut 2016 =" << cb2016_count << endl;
				//fout << "precedence cut 2016=" << previo_cut_count << endl;
				//fout << "subtour cut =" << subtourcut_count << endl;
				//fout  << "inf cut =" << inf_cut_count << endl;
				//fout << "opt cut =" << opt_cut_count << endl;
				//fout << "SP count =" << SP_count << endl;

				finish = clock();
				duration = (double)(finish - start) / CLOCKS_PER_SEC;
				fout << iterNo << "- final UB =" << UB << ", current time= " << duration << endl;
#endif		
			}//endwhile
			cplex.clearModel();
			cplex.clear();
			cplex.end();
			model.end();
			//nodesNO[my - start_instance] = cplex.getNnodes();
out_end:
				/////////算法结束/////////////

			finish = clock();
			duration = (double)(finish - start) / CLOCKS_PER_SEC;

			gap = (UB - LB) / LB;

			GapAve2[my - 1] = gap;


			//cout<<endl<<endl;


			fout << endl << "xF_opt" << endl;
			for (k = 0; k < nbCrane; k++)
			{
				for (i = 0; i < nbTask; i++)
				{
					for (j = 0; j < nbTask; j++)
					{
						if (xF_opt[k][i][j] == 1)
							fout << " QC " << k + 1 << " move from " << i + 1 << " to " << j + 1 << endl;
					}
				}
				//fout << endl;
			}
			//cout << "CF_best: " << CF_best << endl;


			//UB = ObjVal;


			//**********************************//
			//             计算目标值           //
			//**********************************//	
			double a, b;
			a = 0, b = 0;
			//	计算OBJ1
			//for(i = 0; i < nbMaterial; i++)for(j = 0; j < nbBay; j++) for(k = 0; k < nbBay; k++) a+=nbcF[i][j][k]*xF_best[i][j][k];			
			//for(j = 0; j < nbBay; j++) for(k = 0; k < nbBay; k++) a+=nNs*wF_best[k][k];


			//计算OBJ2
			b = CF_best;

			GapAve[my - 1] = gap;


			//fout.precision(10);
			////fout<<"OBJ1 = "<<a<<endl;
			//fout<<"OBJ2 = "<<b<<endl;
			//fout.precision(10);
			//fout<< "\nobj1+obj2 = " << ObjVal << endl;

			fout.precision(10);
			fout << "\LB0 = " << LB0 << endl;
			fout << "\LB = " << LB << endl;
			fout << "\UB = " << UB << endl;

			fout << endl << "cb cut =" << cbcut_count << endl;
			fout << "cb single cut =" << cbcut_single_count << endl;
			fout << "cb subset cut =" << cbcut_subset_count << endl;
			fout << "cb cut 2016 =" << cb2016_count << endl;
			fout << "precedence cut 2016=" << previo_cut_count << endl;
			fout << "subtour cut =" << subtourcut_count << endl;

			fout << endl << "inf cut =" << inf_cut_count << endl;
			fout << "opt cut =" << opt_cut_count << endl;
			fout << "SP count =" << SP_count << endl;


			IterNo[my - start_instance] = iterNo;
			preccut_countAve[my - start_instance] = cbcut_count;
			preccut_single_countAve[my - start_instance] = cbcut_single_count;
			preccut_subset_countAve[my - start_instance] = cbcut_subset_count;
			subcut_countAve[my - start_instance] = subtourcut_count;
			preccut_count_2016Ave[my - start_instance] = previo_cut_count;
			opt_cut_countAve[my - start_instance] = opt_cut_count;
			inf_cut_countAve[my - start_instance] = inf_cut_count;
			SP_countAve[my - start_instance] = SP_count;


			cbcut_count = 0;
			subtourcut_count = 0;
			previo_cut_count = 0;
			opt_cut_count = 0;
			inf_cut_count = 0;
			SP_count = 0;
			cb2016_count = 0;
			cbcut_subset_count = 0;

			LB_initial[my - start_instance] = LB0;
			GapAve[my - start_instance] = gap;
			DuraAve[my - start_instance] = duration;
			//DuraAve[my - start_instance] = duration;
			DuraAve_MP[my - start_instance] = duration_MP - duration_SP2;
			DuraAve_SP[my - start_instance] = duration_SP + duration_SP2;

			fout << "IterNo=" << iterNo << endl;
			fout << "\nGAP = " << gap << endl;

			fout << "MP_Time=" << duration_MP - duration_SP2 << endl;
			fout << "SP_Time=" << duration_SP + duration_SP2 << endl;
			fout << "Time allocation_Time=" << duration_TA << endl;
			fout << "Time=" << duration << endl;

			//finish=clock();
			//duration = (double)(finish - start) / CLOCKS_PER_SEC;
			//fout<<"Time="<<duration<<endl;
			//cout<<"Time="<<duration<<endl;

			ObjAve[my - start_instance] = UB;
			LBAve[my - start_instance] = LB;

			//fout << "Time=" << duration << endl;
			//cout << "Time=" << duration << endl;

			//duration_SP = 0;
			//cb_cut_count = 0;
			//cut_count = 0;
			//user_cut_count = 0;
			//prec_cut_count = 0;
			//ST_cb_cut_count = 0;

		}
		catch (IloException& e)
		{
			cerr << "ERROR: " << e.getMessage() << endl;
		}
		catch (...)
		{
			cerr << "Error" << endl;
		}
	TERMINATE: env.end();
	}


	double gapa = 0;
	double gapa2 = 0;
	double durationa = 0;
	double durationa_mp = 0;
	double durationa_sp = 0;
	double obja = 0;
	//double itera = 0;
	for (int i = 0; i < end_instance - start_instance + 1; i++)
	{
		gapa += GapAve[i];
		//gapa2 += GapAve2[i];
		durationa += DuraAve[i];
		durationa_mp += DuraAve_MP[i];
		durationa_sp += DuraAve_SP[i];
		obja += ObjAve[i];
		//itera += IterAve[i];
	}


	durationa = durationa / (end_instance - start_instance + 1);
	durationa_mp = durationa_mp / (end_instance - start_instance + 1);
	durationa_sp = durationa_sp / (end_instance - start_instance + 1);
	gapa = gapa / (end_instance - start_instance + 1);
	//gapa2 = gapa2 / (end_instance - start_instance + 1);
	obja = obja / (end_instance - start_instance + 1);
	//itera = itera / (end_instance - start_instance + 1);
	cout << "average gap: " << gapa << endl;
	//cout << "average gap: " << gapa2 << endl;
	cout << "average CPU: " << durationa << endl;
	cout << "average CPU of MP: " << durationa_mp << endl;
	cout << "average CPU of SP: " << durationa_sp << endl;
	cout << "average obj: " << obja << endl;
	//cout << "average iter: " << itera << endl;

	ofstream oFile, oFile1;
	string thefilename1 = "haha.csv";

	for (int i = start_instance; i <= end_instance; i++) {
		oFile1.open(thefilename1, ios::app);//| ios::trunc // 这样就很容易的输出一个需要的excel 文件  

		oFile1 << i << "," << nbTask << "," << nbBay << "," << nbCrane << ","
			<< multiple_1 * LB_initial[i - start_instance] << ","
			<< multiple_1 * LBAve[i - start_instance] << ","
			<< multiple_1 * ObjAve[i - start_instance] << ","
			<< preccut_countAve[i - start_instance] << ","
			<< preccut_single_countAve[i - start_instance] << ","
			<< preccut_subset_countAve[i - start_instance] << ","
			<< subcut_countAve[i - start_instance] << ","
			<< preccut_count_2016Ave[i - start_instance] << ","
			<< opt_cut_countAve[i - start_instance] << ","
			<< inf_cut_countAve[i - start_instance] << ","
			<< nodesNO[i - start_instance] << "," << SP_countAve[i - start_instance] << ","
			<< IterNo[i - start_instance] << ","
			<< fixed << setprecision(2) << DuraAve_MP[i - start_instance] << ","
			<< fixed << setprecision(2) << DuraAve_SP[i - start_instance] << ","
			<< fixed << setprecision(2) << DuraAve[i - start_instance] << "," 
			<< fixed << setprecision(2) << ObjUniQCSP[i - start_instance]  << "," 
			<< fixed << setprecision(2) << ObjbiQCSP[i - start_instance] << "," 
			<< fixed << setprecision(2) << LB0_BD[i - start_instance]<<"," 
			<< fixed << setprecision(2) << ObjAve[i - start_instance] << "," << endl;

		oFile1.close();
	}


	system("pause");
	return 0;
}



IloInt CSP(IloCplex subcplex, IloNumArray3  xF_best, IloNumArray2 yF_best, IloNumArray2 xF_begin_best, IloNumArray2 xF_end_best)
{
	IloEnv env = subcplex.getEnv();
	IloInt i, j, k;

	subcplex.clearModel();

	IloNumVarArray QC_CF(env, nbCrane, 0, 5000);
	IloNumVarArray Task_CF(env, nbTask, 0, 5000);
	IloIntVar CF(env);

	BoolVarMatrix zF(env, nbTask);
	for (i = 0; i < nbTask; i++)
	{
		zF[i] = IloBoolVarArray(env, nbTask);
	}


	//for (k = 0; k < nbCrane; k++)
	//	for (i = 0; i < nbTask; i++)
	//	{
	//	if (xF_begin_best[k][i] > threshold_M)
	//		xF_begin_best[k][i] = 1;
	//	else
	//		xF_begin_best[k][i] = 0;

	//	if (xF_end_best[k][i] > threshold_M)
	//		xF_end_best[k][i] = 1;
	//	else
	//		xF_end_best[k][i] = 0;

	//	for (j = 0; j < nbTask; j++)
	//	{
	//		if (xF_best[k][i][j] > threshold_M)
	//			xF_best[k][i][j] = 1;
	//		if (xF_best[k][i][j]  <= 1- threshold_M)
	//			xF_best[k][i][j] = 0;
	//	}

	//	}

	//for (int ai = 0; ai < nbTask; ai++)
	//	for (int bi = 0; bi < nbTask; bi++)
	//		if (nbprecR[ai][bi] == 1)
	//		{
	//	cout << "(" << ai + 1 << "," << bi + 1 << "), ";
	//		}cout << endl;


	//问题模型
	IloModel submodel(env);

	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	IloNumArray ctime_t[nbTask];// the travelling time of QC moving from the bay i to the bay j
	for (i = 0; i < nbTask; i++)
		ctime_t[i] = IloNumArray(env, nbTask);


	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation_new[i] == nbLocation_new[j]) ctime_t[i][j] = 0;
			else if (nbLocation_new[i] > nbLocation_new[j])
			{
				ctime_t[i][j] = (nbLocation_new[i] - nbLocation_new[j]) * QCmove_time;
			}
			else
			{
				ctime_t[i][j] = (nbLocation_new[j] - nbLocation_new[i]) * QCmove_time;
			}
		}

	}

	//**********************************//
	//            子问题 约束           //
	//**********************************//

////constraints (1r)
	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{

		IloExpr  epa(env);
		epa += CF;
		epa -= QC_CF[k];

		c2.add(epa >= 0);
		epa.end();

	}
	submodel.add(c2);
	c2.end();

	////constraints (1g)
	IloRangeArray  c7(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (i != j)
			{
				int ind1 = 0;
				for (k = 0; k < nbCrane; k++)
					if (xF_best[k][i][j] == 1)
						ind1 = 1;
				if (ind1)
				{
					IloExpr  epa(env);
					epa += Task_CF[i] + ctime_t[i][j] + nbQ_new[j] - Task_CF[j];
					c7.add(epa <= 0);
					epa.end();
				}
			}
	}
	submodel.add(c7);
	c7.end();

	////constraints (1h)
	IloRangeArray  c8(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			//if (i != j)
			if (nbLocation_new[i] != nbLocation_new[j])
			{
				IloExpr  epa(env);
				epa += Task_CF[i] + nbQ_new[j] - Task_CF[j];
				epa -= nbM * (1 - zF[i][j]);
				c8.add(epa <= 0);
				epa.end();
			}
		}
	}
	submodel.add(c8);
	c8.end();

	////constraints (9)
	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{	//if (i!=j)
			if (nbLocation_new[i] != nbLocation_new[j])
			{
				IloExpr  epa(env);
				epa += Task_CF[j] - nbQ_new[j] - Task_CF[i];
				epa -= nbM * zF[i][j];

				c9.add(epa <= 0);
				epa.end();
			}
		}
	}
	submodel.add(c9);
	c9.end();



	//////constraints (10)
	//IloRangeArray  c10(env);
	//for (i = 0; i < nbTask-1; i++)
	//{
	//	for (j = i+1; j < nbTask; j++)
	//	{
	//		if (nbLocation_new[i]<=nbLocation_new[j]+safe_margin)
	//		{
	//			int sumy = 0;
	//			for (int kk = 0; kk <= k; kk++)
	//				sumy += yF_best[kk][j];
	//			for (int kk = k; kk < nbCrane; kk++)
	//				sumy += yF_best[kk][i];

	//			if (sumy >=1.9)
	//			{
	//				IloExpr  epa(env);
	//				epa += zF[i][j] + zF[j][i];
	//				c10.add(epa >= 1);
	//				epa.end();
	//			}

	//		}
	//	}
	//}
	//submodel.add(c10);
	//c10.end();

	// new constraints (11)-(12)
	for (int k = 0; k < nbCrane - 1; k++)
	{
		for (int k2 = k + 1; k2 < nbCrane; k2++)
		{
			for (i = 0; i < nbTask - 1; i++)
			{
				for (int j = i + 1; j < nbTask; j++)
				{
					if (para_Delta[k][k2][i][j] > 0 && yF_best[k][i] > threshold_M && yF_best[k2][j] > threshold_M)
					{
						//(10)
						submodel.add(zF[i][j] + zF[j][i] == 1);

						//(11)
						submodel.add(Task_CF[i] + para_Delta[k][k2][i][j] + nbQ_new[j] - Task_CF[j] - nbM * (1 - zF[i][j]) <= 0);
						//(12)
						submodel.add(Task_CF[j] + para_Delta[k][k2][i][j] + nbQ_new[i] - Task_CF[i] - nbM * (1 - zF[j][i]) <= 0);

					}
					else if (para_Delta[k2][k][i][j] > 0 && yF_best[k2][i] == 1 && yF_best[k][j] == 1)
					{
						//(10)
						submodel.add(zF[i][j] + zF[j][i] == 1);


						//(11)
						submodel.add(Task_CF[i] + para_Delta[k2][k][i][j] + nbQ_new[j] - Task_CF[j] - nbM * (1 - zF[i][j]) <= 0);
						//(12)
						submodel.add(Task_CF[j] + para_Delta[k2][k][i][j] + nbQ_new[i] - Task_CF[i] - nbM * (1 - zF[j][i]) <= 0);

					}
				}
			}
		}
	}

	////constraints (1m)
	IloRangeArray  c13(env);
	for (i = 0; i < nbTask - 1; i++)
	{
		for (j = i + 1; j < nbTask; j++)
		{

			if (nbLocation_new[i] + safe_margin >= nbLocation_new[j] && nbLocation_new[i] <= nbLocation_new[j])
			{
				if (nbprecR[i][j] != 1)
				{
					IloExpr  epa(env);
					epa += zF[i][j] + zF[j][i];

					c13.add(epa == 1);
					epa.end();
				}
			}
			else
				break;
		}
	}
	submodel.add(c13);
	c13.end();

	////constraints (1n)
	IloRangeArray  c14(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbprecR[i][j] == 1)
			{
				c14.add(zF[i][j] == 1);
				c14.add(zF[j][i] == 0);
			}
		}
	}
	submodel.add(c14);
	c14.end();

	////constraints (1o)
	IloRangeArray  c15(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
			if (xF_begin_best[k][j] == 1)
			{
				IloExpr  epa(env);
				epa += nbQ_new[j] - Task_CF[j];

				if (nbLocation_new[j] > nbb_new[k])
					epa += (nbLocation_new[j] - nbb_new[k]) * QCmove_time;
				else
					epa -= (nbLocation_new[j] - nbb_new[k]) * QCmove_time;

				epa += nbreadyT[k];//QC ready time, already defined in the main function

				c15.add(epa <= 0);
				epa.end();
			}
	}
	submodel.add(c15);
	c15.end();

	////constraints (1p)
	IloRangeArray  c16(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (xF_end_best[k][j] == 1)
			{
				IloExpr  epa(env);
				epa += Task_CF[j] - QC_CF[k];


				c16.add(epa <= 0);
				epa.end();
			}
		}

	}
	submodel.add(c16);
	c16.end();

	////constraints (1q)
	IloRangeArray  c17(env);
	for (k = 0; k < nbCrane; k++)
	{

		IloExpr  epa(env);
		epa += dT[k];
		epa -= QC_CF[k];

		c17.add(epa >= 0);
		epa.end();

	}
	submodel.add(c17);
	c17.end();



	////  fixing (14)
	//IloRangeArray  c0(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbTask; i++)
	//	{
	//		for (j = 0; j < nbTask; j++)
	//			if(i!=j && xF_best[k][i][j]==1)
	//			{
	//				c0.add(zF[i][j] == 1);
	//				c0.add(zF[j][i] == 0);
	//			}
	//	}

	//}
	//submodel.add(c0);
	//c0.end();


		//**********************************//
		//            开始求解		        //
		//**********************************//
	subcplex.extract(submodel);




	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.01);
	//cplex.setOut(env.getNullStream());
	subcplex.setParam(IloCplex::TiLim, 3600);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	bool h0;
	h0 = subcplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		//cout<<"\nno feasible solution has been found"<<endl; return false;
		subcplex.clearModel();
		submodel.end();

		return 0;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	IloInt ObjVal = subcplex.getBestObjValue();




	//IloInt CF_best=cplex.getValue(CF);

	subcplex.clearModel();
	submodel.end();

	return ObjVal;

}


bool CBMP_R(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation,
	BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, IloNumVarArray uuF,
	BoolVarMatrix xF_begin, BoolVarMatrix xF_end, IloInt* CF_best,
	IloNum* ObjVal, int sol_mode)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;
	//BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, BoolVarMatrix uF, IloIntVarArray thetaF, IloBoolVarArray vF,

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

	//	将目标函数加入到主问题模型
	model.add(IloMinimize(env, obj));//
	//obj1.end();
	obj.end();



#ifdef	uuF_sequence
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			//if (i!=j)
			if (nbprecR[i][j]==1)
			{
				//IloExpr epa(env);
				//for (k = 0; k < nbCrane; k++)
				//	epa += xF[k][i][j];

				//model.add(uuF[j] - uuF[i] +1+ nbTask * (1 - epa) >= 0);

				model.add(uuF[j] - uuF[i] >= 1);
			}
		}
	}
#endif


	//NumMatrix travel_time(env, nbTask);
	//for (i = 0; i < nbTask; i++)
	//{
	//	travel_time[i] = IloNumArray(env, nbTask);
	//}
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{

			if (nbLocation[i] == nbLocation[j]) travel_time[i][j] = 0;
			else if (nbLocation[i] > nbLocation[j])
			{
				travel_time[i][j] = (nbLocation[i] - nbLocation[j]) * QCmove_time;
			}
			else
			{
				travel_time[i][j] = (nbLocation[j] - nbLocation[i]) * QCmove_time;
			}

		}

	}

	NumMatrix2 ctime(env);// the travelling time of crane k moving from the bay li to the bay lj and the processing time of task j
	for (k = 0; k < nbCrane; k++)
	{
		ctime.add(NumMatrix(env));
		for (i = 0; i < nbTask; i++)
		{
			ctime[k].add(IloNumArray(env, nbTask));
		}
	}

	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			for (j = 0; j < nbTask; j++)
			{
				ctime[k][i][j] = nbQ[j] + travel_time[i][j];
			}

		}
	}


	//**********************************//
	//            MP问题 约束           //
	//**********************************//

	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			model.add(xF[k][i][i] == 0);


	IloRangeArray  c1(env);
	for (k = 0; k < nbCrane; k++)
	{

		IloExpr  epa(env);
		epa += CF - nbreadyT[k];
		for (i = 0; i < nbTask; i++)
		{
			epa -= nbQ[i] * xF_begin[k][i];

			if (nbb[k] >= nbLocation[i])
				epa -= QCmove_time * (nbb[k] - nbLocation[i]) * xF_begin[k][i];
			else
				epa -= QCmove_time * (nbLocation[i] - nbb[k]) * xF_begin[k][i];

			for (j = 0; j < nbTask; j++)
				epa -= ctime[k][i][j] * xF[k][i][j];

		}
		c1.add(epa >= 0);
		epa.end();

	}
	model.add(c1);
	c1.end();

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
				if (i != j)
					epa -= xF[k][i][j];
			}
			c4.add(epa == 0);
			epa.end();
		}
	model.add(c4);
	c4.end();

	//约束（5）//task allocation
	IloRangeArray  c5(env);
	for (j = 0; j < nbTask; j++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			epa += yF[k][j] - xF_begin[k][j];
			for (i = 0; i < nbTask; i++)
			{
				if (i != j)
					epa -= xF[k][i][j];
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




	//********* available bay region of each QC ***********//

#ifndef indented_berth
	//QC travel limits
	IloRangeArray  v000(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation[i] < (1 + safe_margin) * k + 1)
				epa += yF[k][i];
			if (nbLocation[i] > nbBay - (1 + safe_margin) * (nbCrane - k - 1))
				epa += yF[k][i];
		}
		v000.add(epa == 0);
		epa.end();
	}
	model.add(v000);
	v000.end();
#endif

#ifdef indented_berth
	IloRangeArray  v00001(env);
	for (k = 0; k < nbCrane; k++)
	{
		if (k< QCno_sideA)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
			{
				if (nbLocation[i] < (1 + safe_margin) * k + 1)
					epa += yF[k][i];
				if (nbLocation[i] > nbBay - (1 + safe_margin) * (QCno_sideA - k - 1))
					epa += yF[k][i];
			}
			v00001.add(epa == 0);
			epa.end();
		}
		else
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
			{
				if (nbLocation[i] < (1 + safe_margin) * (k- QCno_sideA) + 1)
					epa += yF[k][i];
				if (nbLocation[i] > nbBay - (1 + safe_margin) * (nbCrane - k - 1))
					epa += yF[k][i];
			}
			v00001.add(epa == 0);
			epa.end();
		}
	}
	model.add(v00001);
	v00001.end();
#endif


	//********* task pairs (i1,j1),(i2,j2), Balas. el al. 1995 PCATSP*********//
	for (i = 0; i < nbTask - 2; i++)
	{
		for (j = i + 1; j < nbTask-1; j++)
		{

			if (nbprecR[i][j] == 1)
			{
				for (int i2 = i+1; i2 < nbTask - 1; i2++)
				{
					if (i2 != j)
					{
						for (int j2 = i2 + 1; j2 < nbTask; j2++)
						{
							if (nbprecR[i2][j2] == 1)
							{								
								//if (j2 != j  && nbLocation[j2] >= nbLocation[i] - 2 - safe_margin
								//	&& nbLocation[j2] <= nbLocation[i]+2+safe_margin && nbLocation[j2]!= nbLocation[i])
									for (k = 0; k < nbCrane; k++)
									{
										model.add(xF[k][i][j2] + xF[k][j][i2] <= 1);
										//enhanced inq
										//model.add(xF[k][j2][i] + xF[k][i][j2] + xF[k][j][i2] <= 1);

										if (nbLocation[j2] == nbLocation[i] + safe_margin +1)
											model.add(xF[k][j2][i] + xF[k][j][i2] <= 1);
									}
							}
							else
								break;
						}
					}
				}
			}
			else
				break;
		}

	}




	////**********************************//
	////            有效不等式          //
	////**********************************//

	BoolVarMatrix alphaF(env, nbCrane);//left
	for (k = 0; k < nbCrane; k++)
	{
		alphaF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix betaF(env, nbCrane);//right
	for (k = 0; k < nbCrane; k++)
	{
		betaF[k] = IloBoolVarArray(env, nbBay);
	}
	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
	NumVarMatrix moveF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		moveF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	//IloNumVarArray t0kF(env, nbCrane, 0, 100);
	//IloNumVarArray t0kF2(env, nbCrane, 0, 100);
	IloNumVarArray detourT(env, nbCrane, 0, 100);

	IloNumVarArray detourEndT(env, nbCrane, 0, 100);


	/////******** 1) waiting time  ************////////
#ifdef lb_inequality


	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -1000, IloInfinity);

	//计算makespan lower bound
	IloRangeArray  lb24(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF;
		//epa -= t0kF[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs) * yF[k][i];
		for (i = 0; i < nbBay; i++)
			epa -= mF[k][i];
		for (i = 0; i < nbTask; i++)
		{
			if (nbb[k] >= nbLocation[i])
				epa -= QCmove_time * (nbb[k] - nbLocation[i]) * xF_begin[k][i];
			else
				epa -= QCmove_time * (nbLocation[i] - nbb[k]) * xF_begin[k][i];

			for (j = 0; j < nbTask; j++)
				epa -= travel_time[i][j] * xF[k][i][j];

		}
		//for (i = 0; i < nbBay; i++)
		//	epa += QCmove_time * (i * alphaF[k][i] - i * betaF[k][i]);
		lb24.add(epa >= 0);
		epa.end();
	}
	model.add(lb24);
	lb24.end();

	// TF: the time spend before leaving bay i
	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			//epa += t0kF[k];
			for (int h = 0; h <= i; h++)
			{
				epa += mF[k][h];
			}
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
					epa += (nbQ[j] / nbs) * yF[k][j];		
			model.add(epa + moveF[k][i] - TF[k][i] == 0);
			epa.end();

		}
	}
	//for (k = 0; k < nbCrane; k++)
	//	if(nbb[k]==k*(1+safe_margin)+1)
	//{
	//	for (i = k * (1 + safe_margin)+1; i < nbBay; i++)
	//	{
	//		IloExpr epa(env);
	//		for (int h = 0; h <=i; h++)
	//			epa += nbM2 * betaF[k][h];
	//		model.add( QCmove_time * i - nbb[k]+1- moveF[k][i] -epa <= 0);
	//		epa.end();
	//	}
	//}
	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
				{
					if (nbb[k] - 1 <= i)
					{
						if (nbb[k] >= nbLocation[j])
							epa += (nbb[k] - nbLocation[j]) * xF_begin[k][j];
						else
							epa += (nbLocation[j] - nbb[k]) * xF_begin[k][j];
					}
					else
						epa += (i - nbLocation[j] + 1) * xF_begin[k][j];
					for (int j2 = 0; j2 < nbTask; j2++)
					{
						if (nbLocation[j2] - 1 <= i)
							epa += travel_time[j2][j] * xF[k][j2][j];
						else
							epa += QCmove_time * (i - nbLocation[j] + 1) * (xF[k][j2][j] + xF[k][j][j2]);
					}
				}
			model.add(epa - moveF[k][i] == 0);
			epa.end();
		}
	}
	//// TF: the time spend before leaving bay i,considering xijk
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		IloExpr  epa(env);
	//		for (int h = 0; h <= i; h++)
	//		{
	//			epa += nbBay * alphaF[k][h];
	//		}
	//		for (j = 0; j < nbTask; j++)
	//			if (nbLocation[j] - 1 <= i)
	//			{
	//				epa += (nbQ[j] / nbs) * yF[k][j];
	//				if (nbb[k] - 1 <= i)
	//				{
	//					if (nbb[k] >= nbLocation[j])
	//						epa += (nbb[k] - nbLocation[j]) * xF_begin[k][j];
	//					else
	//						epa += (nbLocation[j] - nbb[k]) * xF_begin[k][j];
	//				}
	//				else
	//					epa += (i - nbLocation[j] + 1) * xF_begin[k][j];
	//				for (int j2 = 0; j2 < nbTask; j2++)
	//				{
	//					if (nbLocation[j2] - 1 <= i)
	//						epa += travel_time[j2][j] * xF[k][j2][j];
	//					else
	//						epa += QCmove_time * (i - nbLocation[j] + 1) * (xF[k][j2][j] + xF[k][j][j2]);
	//				}
	//			}
	//		model.add(epa - nbBay - TF[k][i] <= 0);
	//		epa.end();
	//	}
	//}
	//约束（3）
	IloRangeArray  lb25(env);
	for (i = 0; i < nbBay - safe_margin - 1; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
#ifdef indented_berth
			if (k + 1 < QCno_sideA || k >= QCno_sideA)
#endif
			{
				IloExpr  epa(env);
				epa += TF[k][i] - TF[k + 1][i + safe_margin + 1];

				if (i >= 1)
				{
					for (int h = 0; h < i; h++)
						epa += nbM2 * betaF[k][h];
				}
				for (int h = 0; h <= i + safe_margin + 1; h++)
					epa -= nbM2 * alphaF[k + 1][h];

				lb25.add(epa + nbM2 >= 0);
				epa.end();
			}
		}
	model.add(lb25);
	lb25.end();

	IloRangeArray  lb26(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
		{
			epa += alphaF[k][i];
			epa2 += betaF[k][i];
		}
		lb26.add(epa == 1);
		epa.end();
		lb26.add(epa2 == 1);
		epa2.end();
	}
	model.add(lb26);
	lb26.end();

	IloRangeArray  lb27(env);
	for (i = 0; i < nbTask; i++)
		if (nbLocation[i] < nbBay)
		{
			for (k = 0; k < nbCrane; k++)
			{
				IloExpr  epa(env);
				epa += yF[k][i];
				for (int h = nbLocation[i]; h < nbBay; h++)
					epa += alphaF[k][h];
				lb27.add(epa <= 1);
				epa.end();
			}
		}

	model.add(lb27);
	lb27.end();

	IloRangeArray  lb28(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			epa -= (nbLocation[i] - 1) * yF[k][i];
			for (int h = 0; h < nbBay; h++)
				epa += h * betaF[k][h];
			lb28.add(epa >= 0);
			epa.end();
		}
	model.add(lb28);
	lb28.end();

	IloRangeArray  lb29(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i * alphaF[k][i] - i * betaF[k][i];
		lb29.add(epa <= 0);
		epa.end();
	}
	model.add(lb29);
	lb29.end();

	IloRangeArray  lb30(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
#ifdef indented_berth
		if (k + 1 < QCno_sideA || k >= QCno_sideA)
#endif
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * alphaF[k + 1][i] - i * alphaF[k][i];
			lb30.add(epa - safe_margin - 1 >= 0);
			epa.end();
		}
	}
	model.add(lb30);
	lb30.end();


	IloRangeArray  lb31(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
#ifdef indented_berth
		if (k + 1 < QCno_sideA || k >= QCno_sideA)
#endif
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * betaF[k + 1][i] - i * betaF[k][i];
			lb31.add(epa - safe_margin - 1 >= 0);
			epa.end();
		}
	}
	model.add(lb31);
	lb31.end();

	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa1(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa1 += QCmove_time * i * alphaF[k][i];
	//	model.add(t0kF[k] - epa1 + QCmove_time * nbb[k] - QCmove_time >= 0);
	//	epa1.end();
	//}

#endif

	////////////********** QC travel **************//////////////
#ifdef move_inequality	

	IloNumVarArray epsilonF(env, nbCrane, 0, nbBay);
	IloNumVarArray zetaF(env, nbCrane, 0, nbBay);
	IloRangeArray  c33(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		IloExpr  epaf(env);
		IloExpr  epar(env);

		for (i = 0; i < nbTask; i++)
		{
			if (nbb[k] >= nbLocation[i])
				epa += (nbb[k] - nbLocation[i]) * xF_begin[k][i];
			else
				epa += (nbLocation[i] - nbb[k]) * xF_begin[k][i];
			if (i < nbTask - 1)
				for (j = i + 1; j < nbTask; j++)
				{
					if (travel_time[i][j] > 0)
					{
						epaf += travel_time[i][j] * xF[k][i][j];
					}
					if (travel_time[j][i] > 0)
					{
						epar += travel_time[j][i] * xF[k][j][i];
					}
				}
		}
		if (nbb[k] == (1 + safe_margin) * k + 1)
		{
			IloExpr  epaf2(env);
			for (i = 0; i < nbBay; i++)
				epaf2 -= QCmove_time * i * betaF[k][i];
			epaf2 += QCmove_time * (nbb[k] - 1);
			c33.add(epa + epaf + epaf2 >= 0);
			epaf2.end();
		}
		else if (nbb[k] == nbBay - (1 + safe_margin) * (nbCrane - k - 1))
		{
			IloExpr  epar2(env);
			for (i = 0; i < nbBay; i++)
				epar2 += QCmove_time * i * alphaF[k][i];
			epar2 -= QCmove_time * (nbb[k] - 1);
			c33.add(epa + epar + epar2 >= 0);
			epar2.end();
		}
		else
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			IloExpr  epaab(env);
			for (i = 0; i < nbBay; i++)
			{
				epaab += QCmove_time * i * (betaF[k][i] - alphaF[k][i]);
				epa1 -= QCmove_time * i * betaF[k][i];
				epa2 += QCmove_time * i * alphaF[k][i];
			}
			epa1 += epsilonF[k] + QCmove_time * (nbb[k] - 1);
			epa2 += epsilonF[k] - QCmove_time * (nbb[k] - 1);
			c33.add(epa1 >= 0);
			c33.add(epa2 >= 0);
			c33.add(zetaF[k] + epsilonF[k] - epaab >= 0);
			c33.add(zetaF[k] - epsilonF[k] + epaab >= 0);
			c33.add(epa + epaf + epar - epaab - zetaF[k] >= 0);
			epa1.end();
			epa2.end();
			epaab.end();
		}
		epa.end();
		epaf.end();
		epar.end();


	}
	model.add(c33);
	c33.end();


	//考虑precedence导致的回返
	for (k = 1; k < nbCrane; k++)
	{
#ifdef indented_berth
		if (k < QCno_sideA || k - 1 >= QCno_sideA)
#endif
			for (i = 0; i < nbTask - 2; i++)
				if (nbprecR[i][i + 1] == 1)
				{
					IloExpr  epa(env);
					IloExpr  epa2(env);

					epa2 -= yF[k][i + 1] - 1;
					for (int ll = 0; ll < k; ll++)
					{
						epa2 -= yF[ll][i];
					}
					for (int ll = i + 2; ll < nbTask; ll++)
						if ((nbLocation[ll] > nbLocation[i + 1])
							|| (nbLocation[ll] == nbLocation[i + 1] && nbprecR[ll][i + 1] != 1 && nbprecR[i + 1][ll] != 1))
						{
							epa2 += xF[k][ll][i + 1];
						}

					model.add(epa2 >= 0);
					epa.end();
					epa2.end();

				}
	}


	////#ifdef  BenchMark_instance
	//	for (k = 0; k < nbCrane - 1; k++)
	//	{
	//#ifdef indented_berth
	//		if (k + 1 < QCno_sideA || k >= QCno_sideA)
	//#endif
	//		for (i = 0; i < nbTask; i++)
	//			if (nbLocation[i] <= nbBay - (1 + safe_margin) * (nbCrane - k - 1))
	//			{
	//				IloExpr  epa1(env);
	//				epa1 += yF[k][i];
	//				for (j = 0; j < nbTask; j++)
	//					if (nbLocation[j] > nbLocation[i] + safe_margin
	//						&& nbLocation[j] <= nbBay - (1 + safe_margin) * (nbCrane - k - 2))
	//						epa1 -= yF[k + 1][j];
	//				model.add(epa1 <= 0);
	//				epa1.end();
	//			}
	//	}
	////#endif

#endif

#ifdef Interference_inequality	
	IloRangeArray  lb35(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
#ifdef indented_berth
		if (k + 1 < QCno_sideA || k >= QCno_sideA)
#endif
		if (nbreadyT[k] == nbreadyT[k + 1])
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
				epa += nbLocation[i] * (xF_begin[k + 1][i] - xF_begin[k][i]);
			lb35.add(epa - safe_margin - 1 >= 0);
			epa.end();
		}
	}
	model.add(lb35);
	lb35.end();

	IloRangeArray  lb36(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
#ifdef indented_berth
		if (k + 1 < QCno_sideA || k >= QCno_sideA)
#endif
		if (dT[k] >UB && dT[k + 1]>UB)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
				epa += nbLocation[i] * (xF_end[k + 1][i] - xF_end[k][i]);
			lb36.add(epa - safe_margin - 1 >= 0);
			epa.end();
		}
	}
	model.add(lb36);
	lb36.end();
#endif
#ifdef opt_property_fix
	for (i = 0; i < nbTask - 1; i++)
	{
		for (j = i + 1; j < nbTask; j++)
		{
			if (nbprecR[i][j] == 1)
			{
				for (k = 0; k < nbCrane; k++)
					if (dT[k] >= UB)
						model.add(xF_end[k][i] == 0);
				break;
			}
		}
	}
	for (j = 1; j < nbTask; j++)
	{
		for (i = j - 1; i >= 0; i--)
		{
			if (nbprecR[i][j] == 1)
			{
				for (k = 0; k < nbCrane; k++)
					if (nbreadyT[k] == 0)
						model.add(xF_begin[k][j] == 0);
				break;
			}
		}
	}
#ifdef opt_property_fix2
	for (k = 1; k < nbCrane; k++)
#ifdef indented_berth
		if (k < QCno_sideA || k-1 >= QCno_sideA)
#endif
	{
		int task_added = -1;
		int sumQi = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation[i] == (1 + safe_margin) * k + 1 && task_added == -1)
			{
				task_added = i;
				sumQi = nbQ_new[i];

				int sumQj = 0;
				int indi_a = -1;
				for (j = 0; j < nbTask; j++)
				{
					if (nbLocation[j] == (1 + safe_margin) * (k - 1) + 1)
					{
						sumQj += nbQ_new[j];
						if (indi_a == -1 && nbreadyT[k - 1] + sumQj < nbreadyT[k] + sumQi)
						{
							model.add(xF_begin[k][i] - xF_begin[k - 1][j] <= 0);
						}
						else if (nbreadyT[k - 1] + sumQj < nbreadyT[k] + sumQi)
						{
							model.add(xF_begin[k][i] - xF[k - 1][indi_a][j] <= 0);
						}
						indi_a = j;

					}
					else if (indi_a >= 0)
						break;
				}
			}
			else if (task_added >= 0)
				break;
		}
	}

#endif
#endif


	//////////////************ precedence **************//////////////

#ifdef 	prec_inequality
	////优先级约束,同吊机，优先级后的不能比前面的先作业
	for (i = 0; i < nbTask - 1; i++)
	{
		for (j = i + 1; j < nbTask; j++)
		{

			if (nbprecR[i][j] == 1)
				for (k = 0; k < nbCrane; k++)
					model.add(xF[k][j][i] == 0);
		}
	}

	//优先级约束,优先级task之间加其他task，如果涉及bay太多加上每次MP求解时间会长很多 	
		for (i = 0; i < nbTask-1; i++)
		{			
			for (j = i+1; j < nbTask; j++)
			{				
				if (nbprecR[i][j] == 1)
				{
					for (k = 0; k < nbCrane; k++)
					{ 
						model.add(xF[k][j][i] == 0);//pre(39)
						for (int i2 = 0; i2 < nbTask; i2++)
						{
							if (i2 != i && i2 != j && (nbLocation_new[i2] == nbLocation_new[i] +1 || nbLocation_new[i2] == nbLocation_new[i] - 1))
							//if (i2 != i && i2 != j)
							{
								model.add(xF[k][j][i2] + xF[k][i2][i] <= 1);//pre(38)
							}
						}
								
					}
				} 
			}
			
		}

	//优先级约束,同吊机，中间有跨越优先级的不能连续作业 variable fixing (40)
	for (i = 0; i < nbTask - 2; i++)
	{
		int ind_KY = 0;
		for (j = i + 2; j < nbTask; j++)
		{
			for (int i1 = i + 1; i1 < j; i1++)
			{
				if (nbprecR[i][i1] == 1 && nbprecR[i1][j] == 1)
				{
					for (k = 0; k < nbCrane; k++)
					{
						model.add(xF[k][i][j] == 0);
					}
					break;
				}
			}
		}
	}


//#ifdef preced_new
//	//precedence violation elimination inequalities (41),(42)
//	for (i = 0; i < nbTask - 3; i++)
//	{
//		for (j = i + 1; j < nbTask - 2; j++)
//		{
//
//			if (nbprecR[i][j] == 1)
//			{
//				for (int i2 = j + 1; i2 < nbTask - 1; i2++)
//				{
//					if (nbLocation[i2] != nbLocation[i])
//					{
//						for (int j2 = i2 + 1; j2 < nbTask; j2++)
//						{
//							if (nbprecR[i2][j2] == 1)
//							{
//								////if (BD2016)
//								////{
//									//for (k = 0; k < nbCrane; k++)
//									//{
//									//	////tprecedence violation elimination inequalities (41)
//									//	model.add(xF[k][i][j2] + xF[k][j][i2] + xF[k][i2][j] <= 1);
//									//	//model.add(xF[k][i2][j] + xF[k][j2][i] + xF[k][i][j2] <= 1);
//
//									//	//model.add(xF[k][i][j2] + xF[k][j][i2] <= 1);
//
//									//	//model.add(xF[k][i][j2] + xF[k][i2][j] <= 1);
//
//									//}
//								////}
//
//								////precedence violation elimination inequalities (44)
//								//if (i2 != i && i2 != j && nbLocation_new[i2] <= nbLocation_new[i] + 5 && nbLocation_new[i2] >= nbLocation_new[i] - 5)
//								//{
//								//	IloExpr  epa(env);
//								//	for (k = 0; k < nbCrane; k++)
//								//	{
//								//		epa += xF[k][j][i2] + xF[k][j2][i];
//								//	}
//								//	model.add(epa <= 1);
//								//	epa.end();
//								//}
//
//								//precedence violation elimination inequalities (44)	
//								if (i2 != i && i2 != j && nbLocation_new[i2] <= nbLocation_new[i] + 5 && nbLocation_new[i2] >= nbLocation_new[i] - 5)
//								{
//									for (k = 0; k < nbCrane; k++)
//									{
//										model.add(xF[k][j][i2] + xF[k][j2][i] <= 1);
//									}
//								}
//							}
//
//
//						}
//
//					}
//				}
//			}
//
//		}
//
//	}
//
//#endif

#endif

	IloRangeArray  v00010(env);
	for (int b = 0; b < nbBay - safe_margin; b++)
	{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		IloExpr  epa4(env);
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation[i] >= b && nbLocation[i] <= b + safe_margin)
			{
				for (k = 0; k < nbCrane; k++)
					if(nbreadyT[k]==0)
						epa += xF_begin[k][i];

				for (k = 0; k < nbCrane; k++)
					if (dT[k] >= UB)
						epa2 += xF_end[k][i];
			}
		}
		v00010.add(epa <= 1);
		v00010.add(epa2 <= 1);
		epa.end();
		epa2.end();
		epa3.end();
		epa4.end();
	}
	model.add(v00010);
	v00010.end();

#ifdef safe_margin_viocut
	int bay_first_task[nbBay], bay_last_task[nbBay];
	for (int b = 0; b < nbBay; b++)
		bay_first_task[b] = bay_last_task[b] = -1;
	for (i = 0; i < nbTask; i++)
	{
		if (bay_first_task[nbLocation_new[i] - 1] == -1)
		{
			bay_first_task[nbLocation_new[i] - 1] = i;
			//if (i < nbTask - 1)
			//{
			//	for (j = i + 1; j < nbTask; j++)
			//	{
			//		if (nbprecR[i][j] == 1)
			//		{
			//			bay_first_task[nbLocation_new[i] - 1] = i;
			//		}
			//	}
			//}
		}
		else
		{
			if ((i == nbTask - 1)|| (i < nbTask - 1 && nbLocation_new[i + 1] > nbLocation_new[i]))
				bay_last_task[nbLocation_new[i] - 1] = i;
		}
	}

	for (k = 0; k < nbCrane; k++)		
	{
		if (nbreadyT[k] == 0)
		{
			for (int b = 0; b < nbBay; b++)
			{
				if (bay_first_task[b] >= 0)
				{
					IloExpr epa0(env);
					epa0 += xF_begin[k][bay_first_task[b]];
					IloExpr epa3(env);
					epa3 += xF_begin[k][bay_first_task[b]];
					if (bay_first_task[b] + 1 < nbTask)
						epa3 += xF[k][bay_first_task[b]][bay_first_task[b] + 1];
					for (i = 0; i < nbTask; i++)
					{
						if (nbLocation_new[i] - 1 >= b - safe_margin
							&& nbLocation_new[i] - 1 <= b + safe_margin 
							&& nbLocation_new[i] - 1 != b
							&& i != bay_first_task[nbLocation_new[i] - 1]
							&& bay_first_task[nbLocation_new[i] - 1] != bay_last_task[nbLocation_new[i] - 1]
							&&nbQ[bay_first_task[nbLocation_new[i] - 1]]>=nbBay)
						{
							epa0 += xF[k][bay_first_task[b]][i];
							if (bay_first_task[b] + 1 < nbTask)
								epa3 += xF[k][bay_first_task[b] + 1][i];
						}
					}
					model.add(epa0 <= 1);
					epa0.end();
					if (bay_first_task[b] + 1 < nbTask
						&& bay_first_task[nbLocation_new[i] - 1] != bay_last_task[nbLocation_new[i] - 1])
						model.add(epa3 <= 2);
					epa3.end();
				}
			}
		}
	}
	for (k = 0; k < nbCrane; k++)
	{
		if (dT[k] >= UB)
		{
			for (int b = 0; b < nbBay; b++)
			{
				if (bay_last_task[b] >= 0)
				{
					IloExpr epa0(env);
					epa0 += xF_end[k][bay_last_task[b]];
					IloExpr epa3(env);
					epa3 += xF_end[k][bay_last_task[b]];
					if(bay_last_task[b] - 1>=0)
						epa3 += xF[k][bay_last_task[b] - 1][bay_last_task[b]];
					for (i = 0; i < nbTask; i++)
					{
						if (nbLocation_new[i] - 1 >= b - safe_margin
							&& nbLocation_new[i] - 1 <= b + safe_margin 
							&& nbLocation_new[i] - 1 != b
							&& i != bay_last_task[nbLocation_new[i] - 1]
							&& bay_first_task[nbLocation_new[i] - 1] != bay_last_task[nbLocation_new[i] - 1]
							&& nbQ[bay_first_task[nbLocation_new[i] - 1]] >= nbBay)
						{
							epa0 += xF[k][i][bay_last_task[b]];
							if (bay_last_task[b] - 1 >= 0)
								epa3 += xF[k][i][bay_last_task[b] - 1];

						}
					}
					model.add(epa0 <= 1);
					epa0.end();
					if (bay_last_task[b] - 1 >= 0
						&& bay_first_task[nbLocation_new[i] - 1] != bay_last_task[nbLocation_new[i] - 1])
						model.add(epa3 <= 2);
					epa3.end();
				}
			}
		}
	}
#endif


	////**********************************//
	////            传统子环消除       //
	////**********************************//


	//两个task之间的子环去除, subtour (22)
	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask - 1; i++)
		{
			for (j = i + 1; j < nbTask; j++)
			{
				model.add(xF[k][i][j] + xF[k][j][i] <= 1);
			}

		}
	}


#ifdef bay_subtour_ineq
	//两个bay之间的子环去除
	for (k = 0; k < nbCrane; k++)
		for (int b = 0; b < nbBay - 1; b++)
		{
			IloExpr  epa(env);
			int count_task = 0;
			for (i = 0; i < nbTask; i++)
			{

				if (nbLocation[i] - 1 == b || nbLocation[i] == b + 2)
				{
					count_task++;
					if (i < nbTask - 1)
					{
						for (j = i + 1; j < nbTask; j++)
							if (j != i && (nbLocation[j] - 1 == b || nbLocation[j] == b + 2))
								epa += xF[k][i][j] + xF[k][j][i];
					}
				}
			}
			if (count_task > 1 && count_task < 6)
				model.add(epa <= count_task - 1);
			epa.end();
		}
#endif

#ifdef UB_ineq
	int maxQi0[nbTask];//minimum waiting time before the start of each task due to precedence
	int maxQjT[nbTask];//minimum time interval between completion of each task and the makespan due to precedence

	maxQi0[0] = 0;
	for (int i = 1; i < nbTask; i++)
	{
		maxQi0[i] = 0;
		for (int j = 0; j < i; j++)
		{

			if (nbprecR[j][i] == 1)
				maxQi0[i] += nbQ[j];
		}
	}
	maxQjT[nbTask - 1] = 0;
	for (int i = 0; i < nbTask - 1; i++)
	{
		maxQjT[i] = 0;
		for (int j = i + 1; j < nbTask; j++)
		{

			if (nbprecR[i][j] == 1)
				maxQjT[i] += nbQ[j];
		}
	}

	for (int k = 0; k < nbCrane; k++)//for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF;
		for (i = 0; i < nbTask; i++)
		{
			epa -= (maxQi0[i] + nbQ_new[i]) * xF_begin[k][i] + maxQjT[i] * xF_end[k][i];
			for (j = 0; j < nbTask; j++)
				epa -= ctime[k][i][j] * xF[k][i][j];
		}
		model.add(epa >= 0);
		epa.end();

		if (nbreadyT[k] != 0)
		{

			IloExpr  epa2(env);
			epa2 += CF - nbreadyT[k];
			for (i = 0; i < nbTask; i++)
			{
				epa2 -= nbQ_new[i] * xF_begin[k][i] + maxQjT[i] * xF_end[k][i];
				for (j = 0; j < nbTask; j++)
					epa2 -= ctime[k][i][j] * xF[k][i][j];
			}
			model.add(epa2 >= 0);
			epa2.end();
		}
	}
#endif


//#ifdef	BenchMark_instance
	if (sol_mode == 1)
	{
		for (i = 0; i < nbTask - 1; i++)
		{
			for (j = i + 1; j < nbTask; j++)
			{

				if (nbprecR[i][j] == 1)
				{
					for (k = 0; k < nbCrane - 1; k++)
					{

						IloExpr Bidcons_i1k1(env);
						IloExpr Bidcons_j1k1(env);
						IloExpr Bidcons_i2k1(env);
						IloExpr Bidcons_j2k1(env);
						for (int i2 = 0; i2 < nbTask; i2++)
						{
							if (nbLocation[i2] < nbLocation[j])
							{
								Bidcons_i1k1 += xF[k][i2][i];
								Bidcons_j1k1 += xF[k][i2][j];
							}
							if (nbLocation[i2] > nbLocation[i])
							{
								Bidcons_i2k1 += xF[k][i2][i];
								Bidcons_j2k1 += xF[k][i2][j];
							}
						}
						IloExpr Bidcons_i1k2(env);
						IloExpr Bidcons_j2k2(env);
						for (int k2 = k; k2 < nbCrane; k2++)
						{
							for (int i2 = 0; i2 < nbTask; i2++)
							{
								if (nbLocation[i2] >= nbLocation[j])
								{
									Bidcons_i1k2 += xF[k2][i][i2];
								}
								if (nbLocation[i2] >= nbLocation[j])
								{
									Bidcons_j2k2 += xF[k2][i2][j];
								}
							}
							if (k2 > k)
							{
								IloExpr Bidcons_j1k2(env);
								IloExpr Bidcons_i2k2(env);
								for (int i2 = 0; i2 < nbTask; i2++)
								{

									if (nbLocation[i2] < nbLocation[j])
									{
										Bidcons_j1k2 += xF[k2][i2][j];
									}
									if (nbLocation[i2] > nbLocation[i])
									{
										Bidcons_i2k2 += xF[k2][i2][i];
									}
								}
								model.add(Bidcons_i1k1 + Bidcons_j1k2 <= 1);
								model.add(Bidcons_j2k1 + Bidcons_i2k2 <= 1);
								Bidcons_j1k2.end();
								Bidcons_i2k2.end();
							}
						}
						model.add(Bidcons_j1k1 <= Bidcons_i1k2);
						model.add(Bidcons_i2k1 <= Bidcons_j2k2);
						Bidcons_j2k1.end();
						Bidcons_i1k1.end();
						Bidcons_j1k1.end();
						Bidcons_i2k1.end();
						Bidcons_i1k2.end();
						Bidcons_j2k2.end();

					}
				}
			}
		}
	}
	//if (sol_mode==2)
	//{
	//	for (i = 0; i < nbTask; i++)
	//		for (j = 0; j < nbTask; j++)
	//			if (nbprecR[i][j] == 1)
	//			{
	//				for (int l = 0; l < nbCrane - 1; l++)
	//					for (k = l + 1; k < nbCrane; k++)
	//					{
	//						IloExpr  epa(env);
	//						epa += yF[l][i] + yF[k][j];
	//						model.add(epa <= 1);
	//					}
	//			}
	//}
//#endif

	//for (i = 0; i < nbTask - 1; i++)
	//{
	//	if (nbQ_new[i] == 1)
	//	{
	//		if (i > 0 && nbprecR[i - 1][i] == 1)
	//			for (k = 0; k < nbCrane; k++)
	//				model.add(yF[k][i]-yF[k][i-1]>=0);
	//		else
	//			for (k = 0; k < nbCrane; k++)
	//				model.add(yF[k][i] - yF[k][i + 1] >= 0);
	//	}
	//}


	return true;

}

bool Bid_SP_1(IloCplex cplex, IloIntVar CF, IloNumArray2 yF_best, IloNumArray2 yF_current2, IloNumArray2 yF2_current, IloNum* ObjVal)
{

	//xF: forward trip, yF: retracing trip


	IloEnv env = cplex.getEnv();
	IloInt i, j, k;

	//问题模型
	IloModel submodel(env);

	//问题模型
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//

	IloNumVarArray QC_CF(env, nbCrane, 0, 3000);

	IloNumVarArray ckF(env, nbCrane, 0, 3000);//
	IloNumVarArray chaF(env, nbCrane, 0, 3000);//

	IloIntVarArray thetaF(env, nbCrane, 0, 3000);
	IloBoolVarArray vF(env, nbCrane);

	NumVarMatrix wF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		wF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CwF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CwF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix zF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		zF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix endCzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		endCzF[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix xF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		xF[k] = IloBoolVarArray(env, nbTask);
	}

	BoolVarMatrix yF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		yF[k] = IloBoolVarArray(env, nbTask);
	}


	//**********************************//
	//           objective 原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//for (k = 0; k < nbCrane; k++)
	//	obj2 += chaF[k];

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();

	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
			if (yF_best[k][i] > threshold_M)
				yF_best[k][i] = 1;
			else
				yF_best[k][i] = 0;
		}


	//**********************************//
	//            variable fixing           //
	//**********************************//
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
				int i_k, j_k;
				for (k = 0; k < nbCrane; k++)
				{
					if (yF_best[k][i] > threshold_M)
						i_k = k;
					if (yF_best[k][j] > threshold_M)
						j_k = k;
				}

				if (i_k > j_k)
				{
					submodel.add(xF[i_k][i] == 1);
					submodel.add(yF[i_k][i] == 0);
				}
				else if (i_k < j_k)
				{
					submodel.add(xF[j_k][j] == 0);
					submodel.add(yF[j_k][j] == 1);
				}
			}


	//**********************************//
	//            Constraints           //
	//**********************************//
	//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

	//变量固定
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			submodel.add(xF[k][i] + yF[k][i] == yF_best[k][i]);

	for (k = 0; k < nbCrane; k++)
		submodel.add(CF - ckF[k] >= 0);

	//for (k = 0; k < nbCrane; k++)
	//	submodel.add(chaF[k] - ckF[k] + UB>= 0);

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		//epa += CF - t0kF[k] - gammaF[k];

		epa += ckF[k] - t0kF[k] - gammaF[k];

		for (i = 0; i < nbTask; i++)
			epa -= (nbQ_new[i] / nbs) * (xF[k][i] + yF[k][i]);
		for (j = 0; j < nbBay; j++)
		{
			epa -= QCmove_time * j * endCzF[k][j];
			epa -= wF[k][j] + CwF[k][j];
		}
		//epa -= thetaF[k];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)
			epa += QCmove_time * j * zF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	// thetaF 取值
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
			//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += QCmove_time * i * endCzF[k][i];
			submodel.add(epa - QCmove_time * (nbLocation_new[j] - 1) * xF[k][j] >= 0);
			submodel.add(epa - QCmove_time * nbLocation_new[j] * yF[k][j] >= 0);
			epa.end();
		}
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*endCzF[k][i];
	//	model.add(epa - thetaF[k] == 0);
	//	epa.end();
	//}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i * zF[k][i];
		submodel.add(epa - nbb_new[k] + 1 <= 0);
		epa.end();
	}


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
			submodel.add(yF[k][i] - vF[k] <= 0);
			//model.add(nbLocation[i] * yF[k][i] - thetaF[k] <= 0);
		}
	for (k = 0; k < nbCrane - 1; k++)
	{
		//model.add(thetaF[k + 1] - thetaF[k] >= 2);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i * endCzF[k + 1][i] - i * endCzF[k][i];
		submodel.add(epa >= 1 + safe_margin);
		epa.end();
	}


	//约束（4）// 所有任务被分配到QC上
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += xF[k][i] + yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	////约束（5）//endCzF 取在最后一个xF处
	//IloRangeArray  c5(env);
	//for (i = 0; i < nbBay; i++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] - 1 == i)
	//			epa += xF[k][j];
	//	epa -= endCzF[k][i];
	//	c5.add(epa >= 0);
	//	epa.end();
	//	}
	//model.add(c5);
	//c5.end();

	//约束（6）// zF 小于最小的xF的bay
	IloRangeArray  c6(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			IloExpr  epa2(env);
			for (j = 0; j < nbTask; j++)
				if (nbLocation_new[j] - 1 < i)
				{
					epa += yF[k][j];
					epa2 += xF[k][j];
				}
			for (j = i; j < nbBay; j++)
			{
				epa += 100 * CzF[k][j];
				epa2 += 100 * zF[k][j];
			}


			c6.add(epa <= 100);
			c6.add(epa2 <= 100);
			epa.end();
			epa2.end();
		}
	submodel.add(c6);
	c6.end();



	//约束（7） zF unique
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		for (i = 0; i < nbBay; i++)
		{
			epa += zF[k][i];
			epa2 += CzF[k][i];
			epa3 += endCzF[k][i];
		}
		epa2 -= vF[k];
		c7.add(epa == 1);
		c7.add(epa2 == 0);
		c7.add(epa3 == 1);
		epa.end();
		epa2.end();
		epa3.end();
	}
	submodel.add(c7);
	c7.end();


	//	建立约束(8)  z 和 z之间隔 /delta +1
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i * zF[k + 1][i] - i * zF[k][i];
		c8.add(epa >= 1 + safe_margin);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += i * CzF[k + 1][i] - i * CzF[k][i];
		epa2 -= nbBay * (vF[k] + vF[k + 1] - 2);
		c8.add(epa2 >= 1 + safe_margin);
		epa2.end();

	}
	submodel.add(c8);
	c8.end();

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time * i * zF[k][i];
		submodel.add(epa1 + epa2 - QCmove_time * (nbb_new[k] - 1) >= 0);
		submodel.add(epa1 - epa2 + QCmove_time * (nbb_new[k] - 1) >= 0);
		epa1.end();
		epa2.end();
	}


	//QC travel limits
	IloRangeArray  v00(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation_new[i] < 2 * k + 1)
				epa += xF[k][i] + yF[k][i];
			if (nbLocation_new[i] > nbBay - 2 * (nbCrane - k - 1))
				epa += xF[k][i] + yF[k][i];
		}
		v00.add(epa <= 0);
		epa.end();
	}
	submodel.add(v00);
	v00.end();



	////////////////sub-problem

	// vF 取值(2)
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += yF[k][i];
		epa -= vF[k];
		submodel.add(epa >= 0);
		epa.end();
	}
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			submodel.add(yF[k][i] - vF[k] <= 0);

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF[k][i];
		submodel.add(epa >= 1);
		epa.end();
	}
	// vF 取值(3) vF与violation关系
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
				for (k = 0; k < nbCrane; k++)
				{
					IloExpr epa(env);
					IloExpr epa2(env);

					epa += vF[k] + 1;
					epa2 += vF[k] + 1;

					if (k < nbCrane - 1)
					{
						epa -= xF[k][i] + yF[k][i];
						for (int kk = k + 1; kk < nbCrane; kk++)
							epa -= xF[kk][j] + yF[kk][j];
						submodel.add(epa >= 0);
					}

					if (k > 0)
					{
						epa2 -= xF[k][j] + yF[k][j];
						for (int kk = 0; kk < k; kk++)
							epa2 -= xF[kk][i] + yF[kk][i];
						submodel.add(epa2 >= 0);
					}
					epa.end();
					epa2.end();
				}
			}


	//precedence
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
				for (k = 0; k < nbCrane; k++)
				{

					{
						IloExpr epa(env);
						epa += yF[k][i];
						for (int kk = k; kk < nbCrane; kk++)
							epa -= yF[kk][j];
						submodel.add(epa <= 0);
						epa.end();
					}

					{
						IloExpr epa2(env);
						epa2 += xF[k][j];
						for (int kk = k; kk < nbCrane; kk++)
							epa2 -= xF[kk][i];
						submodel.add(epa2 <= 0);
						epa2.end();
					}

					if (k < nbCrane - 1)
					{
						IloExpr epa3(env);
						epa3 += xF[k][i];
						for (int kk = k + 1; kk < nbCrane; kk++)
							epa3 += xF[kk][j];
						//epa3 += xF[kk][j] + yF[kk][j];
						submodel.add(epa3 <= 1);
						epa3.end();
					}


					if (k < nbCrane - 1)
					{
						IloExpr epa4(env);
						epa4 += yF[k][j];

						//for (int kk = k + 1; kk < nbCrane; kk++)
						//	epa4 += yF[kk][i];
						//epa4 -= 1;


						for (int kk = k + 1; kk < nbCrane; kk++)
							epa4 -= xF[kk][i];
						for (int kk = 0; kk <= k; kk++)
							epa4 -= xF[kk][i] + yF[kk][i];
						submodel.add(epa4 <= 0);
						epa4.end();
					}

					for (int kk = 0; kk <= k; kk++)
					{
						submodel.add(yF[k][j] + xF[kk][i] - vF[kk] <= 1);
					}

				}
			}

	////back time gammaF
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += gammaF[k];
		for (i = 0; i < nbBay; i++)
			epa += i * CzF[k][i] - i * endCzF[k][i];
		epa += nbBay - nbBay * vF[k];
		submodel.add(epa >= 0);
		epa.end();
	}

	//////////////////////////////////////
	//////waiting time ///////
	/////////////////////////////////////

	///////forward trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += t0kF[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation_new[j] - 1 <= i)
					epa += (nbQ_new[j] / nbs) * xF[k][j];
			for (j = 0; j <= i; j++)
			{

				epa += wF[k][j];
				epa -= QCmove_time * j * zF[k][j];
				//epa -= 1000 * CzF[k][j];
			}
			submodel.add(epa + QCmove_time * i - TF[k][i] == 0);
			epa.end();

		}

	/////约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 2; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
			IloExpr  epa(env);

			epa += TF[k][i] - TF[k + 1][i + 2];
			epa += 1600;
			for (j = 0; j <= i + 2; j++)
				epa -= 1600 * zF[k + 1][j];
			for (j = 0; j < i; j++)
				epa += 1600 * endCzF[k][j];

			c3.add(epa >= 0);
			epa.end();
		}
	submodel.add(c3);
	c3.end();

	//////retracing trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += t0kF[k];
			for (j = 0; j < nbTask; j++)
				epa += (nbQ_new[j] / nbs) * xF[k][j];
			for (j = 0; j < nbBay; j++)
			{
				epa += QCmove_time * 2 * j * endCzF[k][j] - QCmove_time * j * zF[k][j];
				epa += wF[k][j];
			}


			for (j = 0; j < nbTask; j++)
				if (nbLocation_new[j] - 1 >= i)
					epa += (nbQ_new[j] / nbs) * yF[k][j];
			for (j = i; j < nbBay; j++)
			{
				epa += CwF[k][j];
			}

			submodel.add(epa - QCmove_time * i - CTF[k][i] == 0);
			epa.end();

		}

	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 2; i++)
		{
			IloExpr  epa(env);
			epa += CTF[k + 1][i + 2] - CTF[k][i];

			epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];
			for (j = 0; j <= i + 2; j++)
				epa -= 2000 * CzF[k + 1][j];

			submodel.add(epa >= 0);
			epa.end();
		}

	// forward and retrace
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbBay; j++)
			epa += j * CzF[k + 1][j] - j * endCzF[k][j];
		epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
		submodel.add(epa >= 2);
		epa.end();

	}

	//**********************************//
	//            开始求解		        //
	//**********************************//
	//IloCplex cplex(env);
	cplex.extract(submodel);

#ifdef UserActive
	//cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.01);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 720);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;

	bool h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found for bidQCSP" << endl;
		cplex.clearModel();
		cplex.clear();
		//cplex.end();
		submodel.end();

		return false;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*ObjVal = cplex.getBestObjValue();
	*ObjVal = cplex.getValue(CF);

	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.getValue(xF[k][i]) > 0.9)
				yF_current2[k][i] = 1;
			else
				yF_current2[k][i] = 0;

			if (cplex.getValue(yF[k][i]) > 0.9)
				yF2_current[k][i] = 1;
			else
				yF2_current[k][i] = 0;
		}


	cplex.clearModel();
	cplex.clear();
	//cplex.end();
	submodel.end();

	return true;

}

bool HeuristicFun(IloModel model, IloNum nbs, IloIntArray  nbb, IloIntArray  nbreadyT2,
	BoolVarMatrix yF, IloIntVar CF, IloInt* CF_best, IloNumArray2 yF_best,
	IloNum* ObjVal, IloNum* UB, IloIntArray  nbLocation, int my)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;


	//问题模型
	IloModel submodel(env);
	BoolVarMatrix CuF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CuF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix CvF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CvF[k] = IloBoolVarArray(env, nbBay);
	}
	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	IloIntVarArray startF(env, nbCrane, 0, 100);
	IloIntVarArray endF(env, nbCrane, 0, 100);

	IloNumVarArray t0kF(env, nbCrane, 0, 100);

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);

	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//

	//IloRangeArray  c0(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] >(nbBay - 2 * (nbCrane - k - 1)))
	//			epa += yF[k][j];
	//	c0.add(epa == 0);
	//	epa.end();
	//}
	//submodel.add(c0);
	//c0.end();


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += nbreadyT2[k] + t0kF[k] - QCmove_time * startF[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
					epa += (nbQ_new[j] / nbs) * yF[k][j];
			for (j = 0; j <= i; j++)
			{

				epa += mF[k][j];
			}
			submodel.add(epa + QCmove_time * i - TF[k][i] == 0);
			epa.end();

		}

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - nbreadyT2[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ_new[i] / nbs) * yF[k][i];
		for (i = 0; i < nbBay; i++)
			epa -= mF[k][i];


		epa += QCmove_time * startF[k];
		epa -= QCmove_time * endF[k];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 1-safe_margin; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
			IloExpr  epa(env);
			epa += TF[k][i] - TF[k + 1][i + 1+safe_margin];

			epa += 2000 * (CvF[k][i] + CuF[k + 1][i + 1 + safe_margin]);
			//epa += 2000 * (CuF[k][i] + CuF[k + 1][i + 2]);

			c3.add(epa >= 0);
			epa.end();
		}
	submodel.add(c3);
	c3.end();

	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	//约束（40）
	IloRangeArray  c40(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += yF[k][i];
		c40.add(epa >= 1);
		epa.end();
	}
	submodel.add(c40);
	c40.end();

	//约束（5）
	IloRangeArray  c5(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			epa += (nbLocation[i] - 1) * yF[k][i];
			epa -= startF[k];
			epa += nbBay * (1 - yF[k][i]);
			c5.add(epa >= 0);
			epa.end();
		}
	submodel.add(c5);
	c5.end();

	//约束（6）
	IloRangeArray  c6(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			epa -= (nbLocation[i] - 1) * yF[k][i];
			epa += endF[k];
			c6.add(epa >= 0);
			epa.end();
		}
	submodel.add(c6);
	c6.end();



	//	建立约束(9)
	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
		for (j = 0; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
				for (int l = 0; l < nbCrane - 1; l++)
					for (k = l + 1; k < nbCrane; k++)
					{
						IloExpr  epa(env);
						epa += yF[l][i] + yF[k][j];
						c9.add(epa <= 1);
					}

			}
	submodel.add(c9);
	c9.end();

	//约束（5）
	IloRangeArray  c10a(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			IloExpr  epa3(env);
			epa1 += i - startF[k];
			epa2 += nbBay * (1 - CuF[k][i]);
			epa3 += nbBay * CuF[k][i];
			c10a.add(epa1 - epa2 + 0.1 <= 0);
			c10a.add(epa3 + epa1 >= 0);
			epa1.end();
			epa2.end();
			epa3.end();
		}
	submodel.add(c10a);
	c10a.end();


	//约束（5）
	IloRangeArray  c10a2(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			IloExpr  epa3(env);
			epa1 += endF[k] - i;
			epa2 += nbBay * (1 - CvF[k][i]);
			epa3 += nbBay * CvF[k][i];
			c10a2.add(epa1 - epa2 + 0.1 <= 0);
			c10a2.add(epa3 + epa1 >= 0);
			epa1.end();
			epa2.end();
			epa3.end();
		}
	submodel.add(c10a2);
	c10a2.end();

	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += startF[k + 1] - startF[k];
		epa2 += endF[k + 1] - endF[k];
		submodel.add(epa1 >= 1 + safe_margin);
		submodel.add(epa2 >= 1 + safe_margin);
		epa1.end();
		epa2.end();
	}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		epa2 += QCmove_time * startF[k] - QCmove_time * (nbb[k] - 1);
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		//submodel.add(startF[k] <= nbb[k] - 1);
		//submodel.add(mF[k][0] == 0);
		epa1.end();
		epa2.end();
	}

	//submodel.add(yF[3][27]==1);


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);
	//cplex.exportModel("LP_format_of_CB.LP");
	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/CB_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/CB_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);
#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 1800);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	bool h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return false;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return false;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	*UB = cplex.getValue(CF);
	*CF_best = cplex.getValue(CF);
	//*UB=(int)(*UB)+1; 

	//cout<<"yF_best"<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		//int sumQK = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(yF[k][i]))
			{
				if (cplex.getValue(yF[k][i]) > 0.1)
				{
					yF_best[k][i] = 1;
					//cout << i + 1 << "  ";
				}

				else
					yF_best[k][i] = 0;
				//cout << yF_best[k][i] << "  ";

				//sumQK += nbQ[i] * yF_best[k][i];
			}
			else
				yF_best[k][i] = 0;
		}//cout << sumQK<<"    "<<endl;
	}
	//cout << endl << endl;

	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "uF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CuF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "vF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CvF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << cplex.getValue(startF[k]) << "    "<<cplex.getValue(endF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//	cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//for (i = 0; i < nbBay; i++)
	//{
	//	cout << i+1 << ":  ";
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j]==i+1)
	//		{
	//			cout <<j+1 << "  ";
	//		}
	//	cout << endl;
	//}cout << "    "; cout << endl;


	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return true;



}

bool Check_HeuristicFun(IloModel model, IloNum nbs, IloIntArray  nbb, IloIntArray  nbreadyT, IloIntVar CF, IloNumArray2 yF_best,
	IloNum* ObjVal, IloIntArray  nbLocation, int my)

{
	IloEnv env = model.getEnv();
	IloInt i, j, k;


	//问题模型
	IloModel submodel(env);
	BoolVarMatrix CuF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CuF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix CvF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CvF[k] = IloBoolVarArray(env, nbBay);
	}
	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	IloIntVarArray startF(env, nbCrane, 0, 100);
	IloIntVarArray endF(env, nbCrane, 0, 100);

	IloNumVarArray t0kF(env, nbCrane, 0, 100);

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);

	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//

	//IloRangeArray  c0(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] >(nbBay - 2 * (nbCrane - k - 1)))
	//			epa += yF[k][j];
	//	c0.add(epa == 0);
	//	epa.end();
	//}
	//submodel.add(c0);
	//c0.end();


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += nbreadyT[k] + t0kF[k] - QCmove_time * startF[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i && yF_best[k][i] == 1)
					epa += nbQ_new[j] / nbs;
			for (j = 0; j <= i; j++)
			{

				epa += mF[k][j];
			}
			submodel.add(epa + QCmove_time * i - TF[k][i] == 0);
			epa.end();

		}

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - nbreadyT[k];
		for (i = 0; i < nbTask; i++)
			if (yF_best[k][i] == 1)
				epa -= nbQ_new[i] / nbs;
		for (i = 0; i < nbBay; i++)
			epa -= mF[k][i];


		epa += QCmove_time * startF[k];
		epa -= QCmove_time * endF[k];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 2; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
			IloExpr  epa(env);
			epa += TF[k][i] - TF[k + 1][i + 2];

			epa += 2000 * (CvF[k][i] + CuF[k + 1][i + 2]);
			//epa += 2000 * (CuF[k][i] + CuF[k + 1][i + 2]);

			c3.add(epa >= 0);
			epa.end();
		}
	submodel.add(c3);
	c3.end();



	//约束（5）
	IloRangeArray  c5(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
			if (yF_best[k][i] == 1)
			{
				IloExpr  epa(env);
				epa += (nbLocation[i] - 1);
				epa -= startF[k];
				c5.add(epa >= 0);
				epa.end();
			}
	submodel.add(c5);
	c5.end();

	//约束（6）
	IloRangeArray  c6(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
			if (yF_best[k][i] == 1)
			{
				IloExpr  epa(env);
				epa -= (nbLocation[i] - 1);
				epa += endF[k];
				c6.add(epa >= 0);
				epa.end();
			}
	submodel.add(c6);
	c6.end();




	//约束（5）
	IloRangeArray  c10a(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			IloExpr  epa3(env);
			epa1 += i - startF[k];
			epa2 += nbBay * (1 - CuF[k][i]);
			epa3 += nbBay * CuF[k][i];
			c10a.add(epa1 - epa2 + 0.1 <= 0);
			c10a.add(epa3 + epa1 >= 0);
			epa1.end();
			epa2.end();
			epa3.end();
		}
	submodel.add(c10a);
	c10a.end();


	//约束（5）
	IloRangeArray  c10a2(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			IloExpr  epa3(env);
			epa1 += endF[k] - i;
			epa2 += nbBay * (1 - CvF[k][i]);
			epa3 += nbBay * CvF[k][i];
			c10a2.add(epa1 - epa2 + 0.1 <= 0);
			c10a2.add(epa3 + epa1 >= 0);
			epa1.end();
			epa2.end();
			epa3.end();
		}
	submodel.add(c10a2);
	c10a2.end();

	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += startF[k + 1] - startF[k];
		epa2 += endF[k + 1] - endF[k];
		submodel.add(epa1 >= 2);
		submodel.add(epa2 >= 2);
		epa1.end();
		epa2.end();
	}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		epa2 += QCmove_time * startF[k] - QCmove_time * (nbb[k] - 1);
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		//submodel.add(startF[k] <= nbb[k] - 1);
		//submodel.add(mF[k][0] == 0);
		epa1.end();
		epa2.end();
	}

	//submodel.add(yF[3][27]==1);


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);
	//cplex.exportModel("LP_format_of_CB.LP");
	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/CB_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/CB_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);
#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 1800);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	bool h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return false;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return false;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	//*UB = cplex.getValue(CF);
	//*CF_best = cplex.getValue(CF);
	//*UB=(int)(*UB)+1; 


	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "uF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CuF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "vF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CvF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << cplex.getValue(startF[k]) << "    "<<cplex.getValue(endF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//	cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//for (i = 0; i < nbBay; i++)
	//{
	//	cout << i+1 << ":  ";
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j]==i+1)
	//		{
	//			cout <<j+1 << "  ";
	//		}
	//	cout << endl;
	//}cout << "    "; cout << endl;


	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return true;



}


bool LB_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ,
	BoolVarMatrix yF, IloIntVar CF, IloInt* CF_best, IloNumArray2 yF_best,
	IloNum* ObjVal, IloNum* UB_Fun, IloIntArray  nbLocation, int my)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloModel submodel(env);
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	//IloNumVarArray betaF(env, nbCrane, 0, 100);

	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix CzF2(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF2[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);


	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//




	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
		//for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		epa += CF - nbreadyT[k];
#ifdef Benchmark_M2011
		epa -= t0kF[k];
#endif
		for (int ii = 0; ii < nbTask; ii++)
			epa -= (nbQ[ii] / nbs) * yF[k][ii];
		//for (j = 0; j < nbLocation[i]-1; j++)	
		for (j = 0; j < nbBay; j++)
			epa -= mF[k][j];
		//epa -= (nbLocation[i] - 1)*yF[k][i];

		//epa -= betaF[k];

		for (j = 0; j < nbBay; j++)//for (j = 2 * k; j < nbBay - 2 * (nbCrane - k - 1); j++)
			epa -= QCmove_time * j * CzF2[k][j];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)//for (j = 2 * k; j < nbBay - 2 * (nbCrane - k - 1); j++)
			epa += QCmove_time * j * CzF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] < nbBay)
			{
				IloExpr  epa(env);

				for (i = nbLocation[j]; i < nbBay; i++)
					epa += CzF[k][i];
				//epa += betaF[k];
				epa += yF[k][j];

				submodel.add(epa <= 1);
				epa.end();
			}


	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
			//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * CzF2[k][i];
			submodel.add(epa - (nbLocation[j] - 1) * yF[k][j] >= 0);
			epa.end();
		}

	//约束（7）
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += CzF[k][i];
		c7.add(epa == 1);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += CzF2[k][i];
		c7.add(epa2 == 1);
		epa2.end();
	}
	submodel.add(c7);
	c7.end();

	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa1(env);
	//	IloExpr  epa2(env);
	//	epa1 += t0kF[k];
	//	for (i = 0; i < nbBay; i++)
	//		epa2 += QCmove_time*i*CzF[k][i];
	//	epa2 -= QCmove_time*nbb[k] - QCmove_time;
	//	submodel.add(epa1 - epa2 >= 0);
	//	submodel.add(epa1 + epa2 >= 0);
	//	epa1.end();
	//	epa2.end();
	//}

	//13o
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i * CzF[k][i] - i * CzF2[k][i];
		submodel.add(epa <= 0);
		epa.end();
	}


	//	建立约束(8)
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		if ((1 + safe_margin) * (k + 1) + 1 <= nbBay)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * CzF[k + 1][i] - i * CzF[k][i];
			c8.add(epa >= 1 + safe_margin);
			epa.end();
		}

		if ((1 + safe_margin) * (k + 2) <= nbBay)
		{
			IloExpr  epa2(env);
			for (i = 0; i < nbBay; i++)
				epa2 += i * CzF2[k + 1][i] - i * CzF2[k][i];
			c8.add(epa2 >= 1 + safe_margin);
			epa2.end();
		}
		//IloExpr  epa3(env);
		//for (i = 0; i < nbBay; i++)
		//	epa3 += i*CzF2[k][i] - i*CzF[k][i];
		//c8.add(epa3 >= 0);
		//epa3.end();

	}
	submodel.add(c8);
	c8.end();

	//for (k = 0; k < nbCrane; k++)
	//{
	//	if ((1 + safe_margin)*k + 1 > nbBay)
	//	{
	//		for (i = 0; i < nbTask; i++)
	//			submodel.add(yF[k][i] == 0);
	//	}
	//}

	IloRangeArray  c010(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] > (nbBay - 2 * (nbCrane - k - 1)))
				epa += yF[k][j];
		c010.add(epa == 0);
		epa.end();
	}
	model.add(c010);
	c010.end();


	////	建立约束(9)
	//IloRangeArray  c9(env);
	//for (i = 0; i < nbTask; i++)
	//	for (j = 0; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	for (int l = 0; l < nbCrane - 1; l++)
	//		for (k = l + 1; k < nbCrane; k++)
	//		{
	//		IloExpr  epa(env);
	//		epa += yF[l][i] + yF[k][j];
	//		c9.add(epa <= 1);
	//		}

	//		}
	//submodel.add(c9);
	//c9.end();




	//定义TF
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)//for (i = (1 + safe_margin) * k; i < nbBay - (1 + safe_margin) * (nbCrane - k - 1); i++)
		{
			IloExpr  epa(env);
#ifdef Benchmark_M2011
			epa += t0kF[k] + nbreadyT[k];
#endif
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)//if (nbLocation[j] - 1 <= i &&nbLocation[j] - 1 >= 2 * k)
					epa += (nbQ[j] / nbs) * yF[k][j];

			for (j = 0; j <= i; j++)//for (j = 2 * k; j <= i; j++)
			{
				epa += mF[k][j];
				epa -= QCmove_time * j * CzF[k][j];
			}

			submodel.add(epa + QCmove_time * i - TF[k][i] == 0);
			epa.end();

		}

	//约束（3）
	IloRangeArray  c3(env);
	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 1 - safe_margin; i++)//for (i = 2 * k; i < nbBay - 2 * (nbCrane - k - 1); i++)
		{
			IloExpr  epa(env);

			epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];

			for (j = 0; j <= i + 1 + safe_margin; j++)
			{
				epa -= nbM * CzF[k + 1][j];
			}
			for (j = 0; j < i; j++)
			{
				epa += nbM * CzF2[k][j];
			}

			epa += nbM;

			c3.add(epa >= 0);
			epa.end();
		}
	submodel.add(c3);
	c3.end();

#ifdef Benchmark_M2011 ////constraints 22,23

	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa1(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa1 += QCmove_time * i * CzF[k][i];
	//	model.add(t0kF[k] + epa1 - QCmove_time * nbb[k] + QCmove_time >= 0);
	//	model.add(t0kF[k] - epa1 + QCmove_time * nbb[k] - QCmove_time >= 0);
	//	epa1.end();
	//}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time * i * CzF[k][i];
		epa2 -= QCmove_time * nbb[k] - QCmove_time;
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		epa1.end();
		epa2.end();
	}
#endif 


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);

	//cplex.exportModel("LP_F.LP");

	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/F_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/F_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);

#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	//cplex.setParam(IloCplex::TiLim, time_limit_LB);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	bool h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return false;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return false;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	*UB_Fun = cplex.getValue(CF);
	*CF_best = cplex.getValue(CF);
	//*UB_Fun=(int)(*UB_Fun)+1; 


	//cout<<"yF_best"<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		//int sumQK = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(yF[k][i]))
			{
				if (cplex.getValue(yF[k][i]) > 0.1)
				{
					yF_best[k][i] = 1;
					//cout << i + 1 << "  ";
				}

				else
					yF_best[k][i] = 0;
				//cout << yF_best[k][i] << "  ";

				//sumQK += nbQ[i] * yF_best[k][i];
			}
			else
				yF_best[k][i] = 0;
		}//cout << sumQK<<"    "<<endl;
	}
	//cout << endl << endl;

	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//		if (cplex.getValue(CzF[k][i])>0.1)
	//			cout << i << "    " << k; 
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "betaF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//		if (cplex.getValue(CzF2[k][i])>0.1)
	//			cout << i << "    " << k; 
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return true;


}

bool AddCBCuts(IloEnv env, IloRangeArray cutPool, BoolVarMatrix2 xF, BoolVarMatrix yF, BoolVarMatrix xF_begin, BoolVarMatrix xF_end, int task_no[],
	IloNumArray2 yF_best, int record_next_task[][nbTask], int record_preceding_task[][nbTask], int first_task[], int last_task[])
{
	IloInt cbcut_count_innerloop = 0;// cb cut
	//IloInt cb2016_count_innerloop = 0;
	//IloInt preccut_count_2016_innerloop = 0;
	IloInt S0_cbcut_single_count_innerloop = 0;
	IloInt ST_cbcut_single_count_innerloop = 0;

	//IloEnv env = cplex.getEnv();

	// /////////////////////////////////////////
////****  recording the operation of each QC  **** ///
// //////////////////////////////////////////

	vector<vector<int>>  S_0T;//存储各个QC的从0到T的任务链条
	vector<vector<int>>  Sk_0T;//辅助任务链条，用于insert subtour形成新的链条
	vector<vector<int>>  Sk_subtour;//辅助任务链条，用于insert subtour形成新的链条

	int sumQall = 0;//所有任务的时间加和
	int sumQk[nbCrane];//每个QC所有任务时间加和
	int task_QC[nbTask];//k_i
	for (int i = 0; i < nbTask; i++)
	{
		task_QC[i] = -1;
		for (int k = 0; k < nbCrane; k++)
		{
			if (yF_best[k][i] >= threshold_M)
				task_QC[i] = k;
		}
		//cout << i + 1 << " task,QC= " << task_QC[i]+1 << endl;
	}

	for (int k = 0; k < nbCrane; k++)
	{
		vector<int> S_0;//存储从0到T的任务链条
		int current_i;
		current_i = first_task[k];

		S_0.push_back(current_i);

		//int reclast = -1;

		//while (current_i != last_task[k] && current_i>=0)
		while (current_i >= 0)
		{
			//if(record_next_task[k][current_i]<0)
			//	reclast = current_i;
			current_i = record_next_task[k][current_i];
			if (current_i >= 0)
			{
				S_0.push_back(current_i);
			}

		}

		S_0T.push_back(S_0);
		Sk_0T.push_back(S_0);

		//cout << k+1<<" QC: ";
		//for (int i2 = 0; i2 < Sk_0T[k].size(); i2++)
		//	cout << Sk_0T[k][i2]+1 << ", ";
		//cout << endl;

		S_0.clear();
		vector<int>().swap(S_0);

		vector<int> S_00;//给subtour存入一个-1元素进行初始化
		S_00.push_back(-1);
		Sk_subtour.push_back(S_00);
		S_00.clear();
		vector<int>().swap(S_00);

		sumQk[k] = 0;
		for (int ii2 = 0; ii2 < nbTask; ii2++)
		{
			if (yF_best[k][ii2] >= threshold_M)
				sumQk[k] += nbQ_new[ii2];
		}

		sumQall += sumQk[k];

	}


	// /////////////////////////////////////////
	////****  subtour cut  **** ///
	// //////////////////////////////////////////
		////////2020-5-14备注：识别出现子环的QC，从而能针对并没有子环的两个QC生成prec—cut
	int subtour_count_sum = 0;
	int subtour_count[nbCrane];

	int task_in_subtour[nbTask];//记录task是否在subtour里，1-在，0-否
	for (int j = 0; j < nbTask; j++)
		task_in_subtour[j] = 1;
	for (int k = 0; k < nbCrane; k++)
		for (int i2 = 0; i2 < Sk_0T[k].size(); i2++)
			task_in_subtour[Sk_0T[k][i2]] = 0;

	int task_in_subtour2[nbTask];//记录task是否在subtour里，1-在，0-否
	for (int j = 0; j < nbTask; j++)
		task_in_subtour2[j] = task_in_subtour[j];

	int sub_tour_set_no[nbCrane];

	for (int k = 0; k < nbCrane; k++)
	{
		subtour_count[k] = 0;
		sub_tour_set_no[k] = 0;
	}

	//cout << "subtour" << endl;
#ifdef subtour_cut
	///加子环
		////subtour cut
	for (int j = 0; j < nbTask; j++)
	{
		if (task_in_subtour[j] == 1)
		{
			for (int k = 0; k < nbCrane; k++)
			{
				if (yF_best[k][j] > threshold_M)
				{
					IloExpr subtour(env);
					IloExpr subtour2(env);
					vector<int> S_j;//存储回环

					S_j.push_back(j);
					int current_task2;
					current_task2 = record_next_task[k][j];
					//由于两个任务组成的子环已经被提前去除
					//cout << "subtour cut involve QC " << k + 1 << ", cut: ";
					while (current_task2 != j && current_task2 >= 0)
					{

						for (int n2 = 0; n2 < S_j.size(); n2++)
						{
							subtour += xF[k][current_task2][S_j[n2]] + xF[k][S_j[n2]][current_task2];
							//cout << current_task2+1 << "+" << S_j[n2] +1<< ", ";
						}

						S_j.push_back(current_task2);

						current_task2 = record_next_task[k][current_task2];
					}
					//fout << endl;

					int county;
					county = S_j.size();

					//cout << county<<endl;

					if (county > 1)
					{
						cutPool.add(subtour <= county - 1);//cut 22
						//for (int n2 = 0; n2 < S_j.size(); n2++)
						//	subtour2 += (1 - 1 / county) * yF[k][S_j[n2]];
						//cutPool.add(subtour - subtour2 <= 0);
						subtour_count[k] = subtour_count[k] + 1;
						subtourcut_count++;
						//cbcut_count_innerloop++;
						//if (sub_tour_set_no[k] < county)
						//	sub_tour_set_no[k] = county;
						//cout << "subtour cut involve QC " << k + 1 << ", " << county << " tasks from " << j + 1 << " to " << current_task2 + 1 << ": ";
						//for (int n2 = 0; n2 < S_j.size(); n2++)
						//	cout << S_j[n2] + 1 << ", ";
					}

					//task that has been added into subtour cut, should be removed from task_in_subtour
					for (int n2 = 0; n2 < S_j.size(); n2++)
					{
						task_in_subtour[S_j[n2]] = 0;
						Sk_subtour[k].push_back(S_j[n2]);//add into the vector of subtour
					}

					S_j.clear();
					vector<int>().swap(S_j);
					subtour.end();
					subtour2.end();
				}
			}
		}
	}

		

#ifdef subtour_alter ////Sk_0T本身，26
	for (int k = 0; k < nbCrane; k++)
	{
		if (Sk_0T[k].size() < task_no[k])//不全在S_0里，有子环
		{
			//cout << k +1<< " QC: last task=" <<last_task[k]+1<<endl;
			//cout <<"S0="<< Sk_0T[k].size() << "，including: " ;
			//for (int i2 = 0; i2 < Sk_0T[k].size(); i2++)
			//	cout << Sk_0T[k][i2] +1<< ", ";
			//cout << endl;
			//cout << k+1 << " QC: yk==1 totally " << task_no[k] << " tasks, including: ";
			//for (j = 0; j < nbTask; j++)
			//	if (yF_best[k][j] > threshold_M)
			//		cout << j+1 << ", ";
			//cout << endl;
			int sum_avg = 0;
			int sum_avg2 = 0;
			int sum_ub = 0;
			for (int i2 = 0; i2 < Sk_0T[k].size(); i2++)
				sum_avg += nbQ_new[Sk_0T[k][i2]];
			for (int j = 0; j < nbTask; j++)
				sum_avg2 += nbQ_new[j];

			for (int kk3 = 0; kk3 < nbCrane; kk3++)
			{
				if (kk3 != k)
				{
					if (UB <= dT[kk3])
						sum_ub += UB - nbreadyT[kk3];
					else
						sum_ub += dT[kk3] - nbreadyT[kk3];
				}
			}

			if (sum_avg2 - sum_ub >= sum_avg)//Sk_0T之外的必须分配给k，可以加 cut
			{
				if (S_0T[k].size() >= 1)
				{
					IloExpr epast1(env);
					IloExpr epast2(env);
					int s0Tsize;
					s0Tsize = S_0T[k].size();
					for (int i2 = 0; i2 < S_0T[k].size(); i2++)
					{
						epast1 += xF_begin[k][S_0T[k][i2]] + xF_end[k][S_0T[k][i2]];
						if (i2 < S_0T[k].size() - 1)
							for (int i3 = i2 + 1; i3 < S_0T[k].size(); i3++)
								epast1 += xF[k][S_0T[k][i2]][S_0T[k][i3]] + xF[k][S_0T[k][i3]][S_0T[k][i2]];
						epast2 += yF[k][S_0T[k][i2]];
					}
					//cutPool.add(epast1 <= s0Tsize);//cut
					cutPool.add(epast1 - epast2 <= 0);//cut
					subtourcut_count++;
					//cbcut_count_innerloop++;
					epast1.end();
					epast2.end();
					//for (int i2 = 0; i2 < S_0T[k].size(); i2++)
					//	cout << S_0T[k][i2]+1 << ", ";
					//cout << endl;
					//cout << s0Tsize << endl;
					//cout << endl;
				}
			}
		}
	}
#endif

#endif
	//fout << "1a" << endl;

	for (int k = 0; k < nbCrane; k++)
	{
		sub_tour_set_no[k] = Sk_subtour[k].size()-1;
		subtour_count_sum += subtour_count[k];
	}


	//// /////////////////////////////////////////
	////****  precedence violation cut  **** /////
	// //////////////////////////////////////////
	// 
	cout << "precedence violation" << endl;
	//#ifndef insert_subtour////前面加subtour cut时，已经把子环加进Sk_0T了//回正
	//	for (int k = 0; k < nbCrane; k++)
	//	{
	//		Sk_0T[k].clear();
	//		vector<int>().swap(Sk_0T[k]);
	//		for (int ii2 = 0; ii2 < S_0T[k].size(); ii2++)
	//			Sk_0T[k].push_back(S_0T[k][ii2]);
	//	}
	//#endif
 
	//precedence violation cut for sequence of single QC
#ifdef pre_vio_singleQC//precedence violation cut for single QC
	for (int k = 0; k < nbCrane; k++)
	{
		//Sk_0T序列中的precedence violation
		if (Sk_0T[k].size() > 2)
		{
			for (int jj = 0; jj < Sk_0T[k].size() - 2; jj++)
			{
				for (int ii = jj + 2; ii < Sk_0T[k].size(); ii++)
				{
					if (nbprecR[Sk_0T[k][ii]][Sk_0T[k][jj]] == 1)
					{
						IloExpr previo40(env);
						for (int i3 = jj + 1; i3 <= ii; i3++)
						{
							previo40 += xF[k][Sk_0T[k][jj]][Sk_0T[k][i3]];
							if (i3 < ii)
							{
								for (int i4 = i3 + 1; i4 <= ii; i4++)
									previo40 += xF[k][Sk_0T[k][i3]][Sk_0T[k][i4]] + xF[k][Sk_0T[k][i4]][Sk_0T[k][i3]];
							}
						}
						cutPool.add(previo40 <= ii - jj - 1);//Precedence violation elimination cut 55
						previo40.end();
						//cbcut_count_innerloop++;

						//precVioCut_count_27++;
						previo_cut_count++;
					}
				}
			}
		}
		//subtour里的precedence violation
		if (Sk_subtour[k].size() > 3)
		{
			for (int jj = 0; jj < Sk_subtour[k].size(); jj++)
			{
				IloExpr pre_vio_cut(env);
				vector<int> S_j;//存储

				int current_task2;
				current_task2 = record_next_task[k][Sk_subtour[k][jj]];
				while (current_task2 != Sk_subtour[k][jj] && current_task2 >= 0)
				{
					pre_vio_cut += xF[k][Sk_subtour[k][jj]][current_task2];
					if(S_j.size()>=1)
						for (int n2 = 0; n2 < S_j.size(); n2++)
							pre_vio_cut += xF[k][current_task2][S_j[n2]] + xF[k][S_j[n2]][current_task2];
					S_j.push_back(current_task2);
					current_task2 = record_next_task[k][current_task2];

					if (nbprecR[current_task2][Sk_subtour[k][jj]] == 1)
					{
						int county;
						county = S_j.size();

						//cout << county<<endl;

						if (county >= 1)
						{
							cutPool.add(pre_vio_cut <= county - 1);//cut 22;
						}
					}
				}
				//fout << endl;

				S_j.clear();
				vector<int>().swap(S_j);
				pre_vio_cut.end();
			}
		}
	}

	////sumSj0+yik<=|Sj0|+1,lifted cut: 
	for (int j = 1; j < nbTask; j++)
	{
		for (int i = j - 1; i < nbTask - 1; i++)
		{
			if (nbprecR[i][j] == 1)
			{
				for (int k = 0; k < nbCrane; k++)
				{
					//if (yF_best[k][j] > threshold_M)
					if (yF_best[k][i] > threshold_M && yF_best[k][j] > threshold_M)
					{
						int indic = 0;//记录：在sequence里，i是否在j之前
						//int indicT = 0;
						int position_j = -1;
						for (int i2 = 0; i2 < Sk_0T[k].size(); i2++)
						{
							if (Sk_0T[k][i2] == j)
							{
								position_j = i2;
								for (int j2 = 0; j2 <= i2; j2++)
								{
									if (Sk_0T[k][j2] == i)
										indic = 1;
								}
							}
						}
						if (position_j >= 0 && indic == 0)
						{
							IloExpr stour23(env);
							for (int i3 = 0; i3 <= position_j; i3++)
							{
								stour23 += xF_begin[k][Sk_0T[k][i3]];
								if (i3 < position_j)
								{
									stour23 -= yF[k][Sk_0T[k][i3]];
									for (int i4 = i3 + 1; i4 <= position_j; i4++)
										stour23 += xF[k][Sk_0T[k][i3]][Sk_0T[k][i4]] + xF[k][Sk_0T[k][i4]][Sk_0T[k][i3]];
								}
							}
							stour23 += yF[k][i];
							cutPool.add(stour23 <= 1);//cut23
							stour23.end();

							previo_cut_count++;
						}

					}
				}
			}
		}
	}

#endif

	//记录是否有涉及两个 precedence pair 的 infeasible route
	int prec_viola_count[nbCrane][nbCrane];
	for (int k = 0; k < nbCrane; k++)
		for (int k2 = 0; k2 < nbCrane; k2++)
			prec_viola_count[k][k2] = 0;
#ifdef inf_route_prec_cut//根据16年文献(27)
	//if (LB < UB - 0.8)
	{
		////////// precedence violation
		for (int i = 0; i < nbTask - 3; i++)
		{
			for (int j = i + 1; j < nbTask - 2; j++)
				if (nbprecR[i][j] == 1)//precedence
				{
					int k_Si;
					int k_Sj;
					for (int k = 0; k < nbCrane; k++)
					{
						if (yF_best[k][i] >= threshold_M)
							k_Si = k;
						if (yF_best[k][j] >= threshold_M)
							k_Sj = k;
					}
					//if (k_Si != k_Sj && subtour_count[k_Si] == 0 && subtour_count[k_Sj] == 0)//不涉及子环
					if (k_Si != k_Sj)
					{

						for (int i2 = j + 1; i2 < nbTask - 1; i2++)
							if (nbLocation_new[i] != nbLocation_new[i2])
							{
								for (int j2 = i2 + 1; j2 < nbTask; j2++)
								{
									if (nbprecR[i2][j2] == 1)//precedence
									{
										int k_Si2;
										int k_Sj2;


										for (int k = 0; k < nbCrane; k++)
										{
											if (yF_best[k][i2] >= threshold_M)
												k_Si2 = k;
											if (yF_best[k][j2] >= threshold_M)
												k_Sj2 = k;
										}
										//两个 prec pair
										if (k_Si == k_Sj2 && k_Sj == k_Si2)//满足条件
										{

											//if (inf_route_prec_cut)
											{
												IloExpr BDcut_Sj2i(env);
												IloExpr BDcut_Sji2(env);
												IloExpr BDcut_Sij2(env);

												vector<int> S_j2i;

												int check_j2_i = 0;//Sj2i是否存在
												int count_Sj2i;
												int task_i0, task_i2;
												int current_i;

												int ind_check_j2i = 0;//看是否有回环
												int ind_check_ji2 = 0;
												int ind_check_ij2 = 0;


												//S_j2i
												current_i = i;
												S_j2i.push_back(current_i);

												while (current_i != j2 && current_i >= 0 && ind_check_j2i == 0)
												{
													current_i = record_preceding_task[k_Si][current_i];

													if (current_i == j2)
														check_j2_i = 1;
													else if (current_i >= 0)
													{
														for (int n2 = 0; n2 < S_j2i.size(); n2++)
														{
															BDcut_Sj2i += xF[k_Si][current_i][S_j2i[n2]] + xF[k_Si][S_j2i[n2]][current_i];
														}
														S_j2i.push_back(current_i);
													}

													if (current_i == i)
													{//有回环
														ind_check_j2i = 1;
													}

												}
												for (int n2 = 0; n2 < S_j2i.size(); n2++)
												{
													BDcut_Sj2i += xF[k_Si][j2][S_j2i[n2]];
												}

												count_Sj2i = S_j2i.size();


												//S_ji2
												vector<int> S_ji2;
												int current_i2;
												int check_j_i2 = 0;
												int count_Sji2;
												current_i2 = i2;
												S_ji2.push_back(current_i2);
												while (current_i2 != j && current_i2 >= 0 && ind_check_ji2 == 0)
												{
													current_i2 = record_preceding_task[k_Sj][current_i2];

													if (current_i2 == j)
														check_j_i2 = 1;
													else if (current_i2 >= 0)
													{
														for (int n2 = 0; n2 < S_ji2.size(); n2++)
														{
															BDcut_Sji2 += xF[k_Sj][current_i2][S_ji2[n2]] + xF[k_Sj][S_ji2[n2]][current_i2];

														}
														S_ji2.push_back(current_i2);
													}
													if (current_i2 == i2)
													{//有回环
														ind_check_ji2 = 1;
													}

												}
												for (int n2 = 0; n2 < S_ji2.size(); n2++)
												{
													BDcut_Sji2 += xF[k_Sj][j][S_ji2[n2]];

												}
												count_Sji2 = S_ji2.size();

												if (check_j_i2 == 1 && check_j2_i == 1 && ind_check_j2i == 0 && ind_check_ji2 == 0)
												{
													//precVioCut_count_2016++;
													prec_viola_count[k_Si][k_Sj] = 1;
													cutPool.add(BDcut_Sji2 + BDcut_Sj2i <= count_Sj2i + count_Sji2 - 1);

													//cbcut_count_innerloop++;

													previo_cut_count++;

												}
												BDcut_Sj2i.end();
												BDcut_Sji2.end();
												S_j2i.clear();
												vector<int>().swap(S_j2i);
												S_ji2.clear();
												vector<int>().swap(S_ji2);

											}
										}
										////#ifdef precut_3//三个 prec pair
										//if (j2 < nbTask - 2)
										//{
										//	for (int i3 = j2 + 1; i3 < nbTask - 1; i3++)
										//	{
										//		if (nbLocation_new[i2] != nbLocation_new[i3])
										//		{
										//			for (int j3 = i3 + 1; j3 < nbTask; j3++)
										//			{
										//				if (nbprecR[i3][j3] == 1)//precedence
										//				{
										//					int k_Si3;
										//					int k_Sj3;
										//					for (int k = 0; k < nbCrane; k++)
										//					{
										//						if (yF_best[k][i3] >= threshold_M)
										//							k_Si3 = k;
										//						if (yF_best[k][j3] >= threshold_M)
										//							k_Sj3 = k;
										//					}
										//					if (k_Si3 == k_Sj2 && k_Si == k_Sj3 && k_Si2 == k_Sj && k_Si != k_Si2 && k_Si != k_Si3 && k_Si2 != k_Si3)//满足条件
										//					{
										//						//if (inf_route_prec_cut)
										//						{
										//							//IloExpr BDcut_Sj2i(env);
										//							IloExpr BDcut_Sji2(env);
										//							IloExpr BDcut_Sj2i3(env);
										//							IloExpr BDcut_Sj3i(env);
										//							vector<int> S_j3i;
										//							int check_j3_i = 0;
										//							int count_Sj3i;
										//							int task_i0, task_i3;
										//							int current_i;
										//							int ind_check_j3i = 0;//看是否有回环
										//							int ind_check_ji2 = 0;
										//							int ind_check_j2i3 = 0;
										//							//S_j3i
										//							current_i = i;
										//							S_j3i.push_back(current_i);
										//							while (current_i != j3 && current_i >= 0 && ind_check_j3i == 0)
										//							{
										//								current_i = record_preceding_task[k_Si][current_i];
										//								if (current_i == j3)
										//									check_j3_i = 1;
										//								else if (current_i >= 0)
										//								{
										//									for (int n2 = 0; n2 < S_j3i.size(); n2++)
										//									{
										//										BDcut_Sj3i += xF[k_Si][current_i][S_j3i[n2]] + xF[k_Si][S_j3i[n2]][current_i];
										//									}
										//									S_j3i.push_back(current_i);
										//								}
										//								if (current_i == i)
										//								{//有回环
										//									ind_check_j3i = 1;
										//								}
										//							}
										//							for (int n2 = 0; n2 < S_j3i.size(); n2++)
										//							{
										//								BDcut_Sj3i += xF[k_Si][j3][S_j3i[n2]];
										//							}
										//							count_Sj3i = S_j3i.size();
										//							//S_ji2
										//							vector<int> S_ji2;
										//							int current_i2;
										//							int check_j_i2 = 0;
										//							int count_Sji2;
										//							current_i2 = i2;
										//							S_ji2.push_back(current_i2);
										//							while (current_i2 != j && current_i2 >= 0 && ind_check_ji2 == 0)
										//							{
										//								current_i2 = record_preceding_task[k_Sj][current_i2];
										//								if (current_i2 == j)
										//									check_j_i2 = 1;
										//								else if (current_i2 >= 0)
										//								{
										//									for (int n2 = 0; n2 < S_ji2.size(); n2++)
										//									{
										//										BDcut_Sji2 += xF[k_Sj][current_i2][S_ji2[n2]] + xF[k_Sj][S_ji2[n2]][current_i2];
										//									}
										//									S_ji2.push_back(current_i2);
										//								}
										//								if (current_i2 == i2)
										//								{//有回环
										//									ind_check_ji2 = 1;
										//								}
										//							}
										//							for (int n2 = 0; n2 < S_ji2.size(); n2++)
										//							{
										//								BDcut_Sji2 += xF[k_Sj][j][S_ji2[n2]];
										//							}
										//							count_Sji2 = S_ji2.size();
										//							//S_j2i3
										//							vector<int> S_j2i3;
										//							int current_i3;
										//							int check_j2_i3 = 0;
										//							int count_Sj2i3;
										//							current_i3 = i3;
										//							S_j2i3.push_back(current_i3);
										//							while (current_i3 != j2 && current_i3 >= 0 && ind_check_j2i3 == 0)
										//							{
										//								current_i3 = record_preceding_task[k_Sj2][current_i3];
										//								if (current_i3 == j2)
										//									check_j2_i3 = 1;
										//								else if (current_i3 >= 0)
										//								{
										//									for (int n2 = 0; n2 < S_j2i3.size(); n2++)
										//									{
										//										BDcut_Sj2i3 += xF[k_Sj2][current_i3][S_j2i3[n2]] + xF[k_Sj2][S_j2i3[n2]][current_i3];
										//									}
										//									S_j2i3.push_back(current_i3);
										//								}
										//								if (current_i3 == i3)
										//								{//有回环
										//									ind_check_j2i3 = 1;
										//								}
										//							}
										//							for (int n2 = 0; n2 < S_j2i3.size(); n2++)
										//							{
										//								BDcut_Sj2i3 += xF[k_Sj2][j2][S_j2i3[n2]];
										//							}
										//							count_Sj2i3 = S_j2i3.size();
										//							if (check_j_i2 == 1 && check_j3_i == 1 && check_j2_i3 == 1 && ind_check_j2i3 == 0 && ind_check_ji2 == 0 && ind_check_j3i == 0)
										//							{
										//								//precVioCut_count_2016++;
										//								//prec_viola_count[k_Si][k_Sj] = 1;
										//								cutPool.add(BDcut_Sji2 + BDcut_Sj3i + BDcut_Sj2i3 <= count_Sj3i + count_Sj2i3 + count_Sji2 - 1);
										//								//preccut_count_2016_innerloop++;
										//								//cbcut_count_innerloop++;
										//								previo_cut_count++;
										//							}
										//							BDcut_Sj2i3.end();
										//							BDcut_Sji2.end();
										//							BDcut_Sj3i.end();
										//							S_j2i3.clear();
										//							vector<int>().swap(S_j2i3);
										//							S_ji2.clear();
										//							vector<int>().swap(S_ji2);
										//							S_j3i.clear();
										//							vector<int>().swap(S_j3i);
										//						}
										//					}
										//				}
										//			}
										//		}
										//	}
										//}//end 三个prec  pair
	//#endif									
									}
								}
							}
					}
				}
		}
	}
#endif

	// /////////////////////////////////////////
////****  recording Si0 and SiT of each task i  **** ///
// /////////////////////////////////////////


	float QC_ub[nbCrane];//记录QC的makespan的当前上界
	for (int k = 0; k < nbCrane; k++)
	{
		QC_ub[k] = UB;
		if (dT[k] < UB)
			QC_ub[k] = dT[k];
		//#ifndef  insert_subtour
		//#ifndef  constructive_cut
				////前面加subtour cut时，已经把子环加进Sk_0T了
				////回正
		//Sk_0T[k].clear();
		//vector<int>().swap(Sk_0T[k]);
		//for (int ii2 = 0; ii2 < S_0T[k].size(); ii2++)
		//	Sk_0T[k].push_back(S_0T[k][ii2]);
		//#endif
				//cout << "QC " << k << ": ";
				//for (int ii2 = 0; ii2 < Sk_0T[k].size(); ii2++)
				//	cout << Sk_0T[k][ii2] +1<< "-";
				//cout << endl;
	}



	IloExprArray BDcut_S0(env, nbTask);//S_i^0						
	IloExprArray BDcut_ST(env, nbTask);//S_i^T

	IloExprArray BDcut_S0_pure(env, nbTask);//S_i^0						
	IloExprArray BDcut_ST_pure(env, nbTask);//S_i^T
	IloExprArray BDcut_S0_x(env, nbTask);//S_i^0						
	IloExprArray BDcut_ST_x(env, nbTask);//S_i^T

	IloExprArray BDcut_S0T_pure(env, nbCrane);//S_0T						
	IloExprArray BDcut_S0T_0(env, nbCrane);//S_0T_T
	IloExprArray BDcut_S0T_T(env, nbCrane);//S_0T_T
	int Ctime_S0T[nbCrane];//sum Q
	int count_S0T[nbCrane];
	int Ctime_subtour[nbCrane];//sum Q


	int Ctime_S0[nbTask];//sum Q
	int Ctime_ST[nbTask];//sum Q
	int count_S0[nbTask];
	int count_ST[nbTask];

	int task_position[nbTask];//记录i在k_i的operation sequence里的位置
	int precQ0[nbTask];//sum Q of all tasks have higher precedence with current task
	int precQT[nbTask];//sum Q of all tasks have lower precedence with current task

	int count_domin_cut_0[nbTask];// 记录是否加过dominated cut,避免重复加cut
	int count_domin_cut_T[nbTask];

	for (int k = 0; k < nbCrane; k++)
	{
		BDcut_S0T_0[k] = IloExpr(env);
		BDcut_S0T_T[k] = IloExpr(env);
		BDcut_S0T_pure[k] = IloExpr(env);
		Ctime_S0T[k] = count_S0T[k] = Ctime_subtour[k] = 0;
		if (S_0T[k].size() >= 2)
			for (int jj = 0; jj < S_0T[k].size() - 1; jj++)
				for (int jj2 = jj + 1; jj2 < S_0T[k].size(); jj2++)
					BDcut_S0T_pure[k] += xF[k][S_0T[k][jj]][S_0T[k][jj2]] + xF[k][S_0T[k][jj2]][S_0T[k][jj]];
		for (int jj2 = 0; jj2 < S_0T[k].size(); jj2++)
		{
			BDcut_S0T_T[k] += xF_end[k][S_0T[k][jj2]];
			BDcut_S0T_0[k] += xF_begin[k][S_0T[k][jj2]];
			Ctime_S0T[k] += nbQ_new[S_0T[k][jj2]];
			count_S0T[k]++;
		}
		if (Sk_subtour[k].size() > 1)
			for (int jj2 = 1; jj2 < Sk_subtour[k].size(); jj2++)
				Ctime_subtour[k] += nbQ_new[Sk_subtour[k][jj2]];
	}

	// ////////////////////////////////////////////
	////**** build S0 and ST sets, with nogood cut and opt cut for single QC  **** /////
	// ///////////////////////////////////////////
	for (int i = 0; i < nbTask; i++)
	{
		BDcut_S0[i] = IloExpr(env);
		BDcut_ST[i] = IloExpr(env);
		BDcut_S0_pure[i] = IloExpr(env);
		BDcut_ST_pure[i] = IloExpr(env);
		BDcut_S0_x[i] = IloExpr(env);
		BDcut_ST_x[i] = IloExpr(env);

		Ctime_S0[i] = Ctime_ST[i] = count_S0[i] = count_ST[i] = 0;


		count_domin_cut_0[i] = count_domin_cut_T[i] = -1;

		precQ0[i] = precQT[i] = 0;
		for (int ii = 0; ii < nbTask; ii++)
		{
			if (nbprecR[ii][i] == 1 && ii != i)
				precQ0[i] += nbQ_new[ii];
			if (nbprecR[i][ii] == 1 && ii != i)
				precQT[i] += nbQ_new[ii];
		}

		int k_Si;
		k_Si = task_QC[i];

		task_position[i] = -1;

		for (int jj = 0; jj < S_0T[k_Si].size(); jj++)
			if (S_0T[k_Si][jj] == i)
				task_position[i] = jj;
		//count_Si = record_i;//不含task 0

		//S_i0 and S_iT for each task i
		if (task_position[i] >= 0)
		{
			//S_i0
			if (task_position[i] >= 1)
			{
				for (int n2 = 0; n2 < task_position[i]; n2++)
				{
					count_S0[i]++;
					Ctime_S0[i] += nbQ_new[S_0T[k_Si][n2]];
					BDcut_S0[i] += xF_begin[k_Si][S_0T[k_Si][n2]];
					BDcut_S0_x[i] += xF[k_Si][S_0T[k_Si][n2]][i];
					if (n2 < task_position[i] - 1)
					{
						for (int n3 = n2 + 1; n3 < task_position[i]; n3++)
						{
							BDcut_S0[i] += xF[k_Si][S_0T[k_Si][n2]][S_0T[k_Si][n3]] + xF[k_Si][S_0T[k_Si][n3]][S_0T[k_Si][n2]];
							BDcut_S0_pure[i] += xF[k_Si][S_0T[k_Si][n2]][S_0T[k_Si][n3]] + xF[k_Si][S_0T[k_Si][n3]][S_0T[k_Si][n2]];
						}
					}
				}
				count_S0[i]++;//包括task 0
			}
			else
			{
				BDcut_S0_x[i] += xF_begin[k_Si][i];
				count_S0[i]++;//包括task 0
			}

			//S_iT
			if (task_position[i] < Sk_0T[k_Si].size() - 1)
			{
				for (int n2 = task_position[i] + 1; n2 < Sk_0T[k_Si].size(); n2++)
				{
					count_ST[i]++;
					Ctime_ST[i] += nbQ_new[S_0T[k_Si][n2]];
					BDcut_ST[i] += xF_end[k_Si][S_0T[k_Si][n2]];
					BDcut_ST_x[i] += xF[k_Si][i][S_0T[k_Si][n2]];
					if (n2 < S_0T[k_Si].size() - 1)
					{
						for (int n3 = n2 + 1; n3 < S_0T[k_Si].size(); n3++)
						{
							BDcut_ST[i] += xF[k_Si][S_0T[k_Si][n2]][S_0T[k_Si][n3]] + xF[k_Si][S_0T[k_Si][n3]][S_0T[k_Si][n2]];
							BDcut_ST_pure[i] += xF[k_Si][S_0T[k_Si][n2]][S_0T[k_Si][n3]] + xF[k_Si][S_0T[k_Si][n3]][S_0T[k_Si][n2]];
						}
					}
				}
				count_ST[i]++;// 包括task T
			}
			else
			{
				count_ST[i]++;// 包括task T
				BDcut_ST_x[i] += xF_end[k_Si][i];
			}

			//cout << i+1 << ", i=" << Sk_0T[k_Si][task_position[i]]+1 << ", record_i=" << task_position[i]
			//	<< ", count_Si0=" << count_S0[i] << ", count_SiT=" << count_ST[i] <<endl;
			//for (int ii2 = 0; ii2 < Sk_0T[k_Si].size(); ii2++)
			//	cout << Sk_0T[k_Si][ii2]+1 << "-";
			//cout << endl;
		}

		//加cut：single QC dominating cut
		if (task_position[i] >= 0)//不在子环里，可以加cut
		{
#ifdef Single_QC_prec_cut//dominated cut on single QC
			///// single QC precedence cut -针对SjT单独违反prec cut-去头

			if (precQ0[i] + nbQ_new[i] + Ctime_ST[i] > QC_ub[k_Si])
			{
				int inic0 = 0;//指示是否加过sebset cut
				count_domin_cut_0[i] = 1;
				cbcut_single_count++;
				cbcut_count_innerloop++;
				ST_cbcut_single_count_innerloop++;

#ifdef  activ_subsetcut_0
				////第一种	按xij原来顺序
				if (count_ST[i] > 2)
				{
					int difference_UB;
					difference_UB = precQ0[i] + nbQ_new[i] + Ctime_ST[i] - QC_ub[k_Si];
					int count_subset = 0;
					IloExpr epa2(env);
					for (int n2 = task_position[i] + 1; n2 < Sk_0T[k_Si].size(); n2++)
					{
						if (nbQ_new[Sk_0T[k_Si][n2]] < difference_UB)
						{
							difference_UB -= nbQ_new[Sk_0T[k_Si][n2]];
							count_subset++;
							epa2 += yF[k_Si][Sk_0T[k_Si][n2]];
						}
					}

					if (count_subset > 0)
					{
						cutPool.add(BDcut_ST[i] - epa2 + yF[k_Si][i] <= count_ST[i] - 1- count_subset);
						inic0 = 1;
						//Alternative cuts based on immediate processing
#ifdef alternative_cut
						IloExpr epa(env);
						int indic2 = 0;
						for (int j0 = 0; j0 < nbTask; j0++)
						{
							int indic = 1;//是否包含在SiT里
							for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
								if (S_0T[k_Si][i2] == j0)
									indic = 0;
							if (indic == 1 && precQ0[j0] + nbQ_new[j0] + Ctime_ST[i] > QC_ub[k_Si])
							{
								for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
								{
									epa += xF[k_Si][j0][S_0T[k_Si][i2]];
								}
								indic2++;
							}
						}
						if (indic2)
						{
							cutPool.add(BDcut_ST[i] - epa2 + epa <= count_ST[i] - 1 - count_subset);
						}
						epa.end();
						//cbcut_multicut_count++;

#endif
//#ifdef alternative_cut2
//						IloExpr epa(env);
//						int indic2 = 0;
//						for (int j0 = 0; j0 < nbTask; j0++)
//						{
//							int indic = 1;//是否包含在SiT里
//							for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
//								if (S_0T[k_Si][i2] == j0)
//									indic = 0;
//							if (indic == 1 && precQ0[j0] + nbQ_new[j0] + Ctime_ST[i] > QC_ub[k_Si])
//							{
//								for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
//								{
//									epa += xF[k_Si][j0][S_0T[k_Si][i2]];									
//								}
//								indic2++;
//							}
//						}
//						int lowind =  nbTask;
//						int upind = 0;
//						for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
//						{
//							if (lowind > S_0T[k_Si][i2])
//								lowind = S_0T[k_Si][i2];
//							if ( upind < S_0T[k_Si][i2])
//								upind = S_0T[k_Si][i2];
//						}
//						IloExpr epa3(env);
//						for (int i3 = lowind; i3 <= upind; i3++)
//							if(yF_best[k_Si][i3]==0)
//							{ 
//								epa3 += xF_end[k_Si][i3] - yF[k_Si][i3];
//								for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
//									epa3 += xF[k_Si][i3][S_0T[k_Si][i2]] + xF[k_Si][S_0T[k_Si][i2]][i3];
//							}
//
//						if (indic2)
//						{
//							cutPool.add(BDcut_ST[i] - epa2 + epa + epa3<= count_ST[i] - 1 - count_subset);
//						}
//						epa.end();
//						epa3.end();
//						//cbcut_multicut_count++;
//						
//#endif
					}
					epa2.end();

				}
#endif
				if(inic0 == 0)//未加sebset cut
					cutPool.add(BDcut_ST[i] + yF[k_Si][i] <= count_ST[i] - 1);
////纯subset，没有y
//#ifdef  activ_subsetcut_1
//
//				////第一种	按xij原来顺序，依次删
//				if (count_ST[i] > 2)
//				{
//					vector<int> S_j22;
//					for (int n2 = task_position[i] + 1; n2 < Sk_0T[k_Si].size(); n2++)
//						S_j22.push_back(Sk_0T[k_Si][n2]);
//					int difference_UB;
//					difference_UB = precQ0[i] + nbQ_new[i] + Ctime_ST[i] - QC_ub[k_Si];
//					int count_Sj22_subset, count_Sj22_subset0;
//					count_Sj22_subset = S_j22.size();
//					count_Sj22_subset0 = S_j22.size();
//					IloExpr BDcut_newj22(env);
//					if (count_Sj22_subset >= 2)
//					{
//						int ind_lin = 0;
//						while (count_Sj22_subset >= 2 && ind_lin == 0)
//						{
//							if (nbQ_new[S_j22[count_Sj22_subset0 - count_Sj22_subset]] < difference_UB)
//							{
//								difference_UB -= nbQ_new[S_j22[count_Sj22_subset0 - count_Sj22_subset]];
//								count_Sj22_subset--;
//							}
//							else
//								ind_lin = 1;
//						}
//
//						for (int n1 = count_Sj22_subset0 - count_Sj22_subset; n1 < count_Sj22_subset0 - 1; n1++)
//						{
//							BDcut_newj22 += xF_end[k_Si][S_j22[n1]];
//							for (int n2 = n1 + 1; n2 < count_Sj22_subset0; n2++)
//								BDcut_newj22 += xF[k_Si][S_j22[n1]][S_j22[n2]] + xF[k_Si][S_j22[n2]][S_j22[n1]];
//						}
//						BDcut_newj22 += xF_end[k_Si][S_j22[count_Sj22_subset0 - 1]];
//					}
//					else if (count_Sj22_subset == 1)
//						BDcut_newj22 += xF_end[k_Si][S_j22[count_Sj22_subset0 - 1]];
//
//					if (count_Sj22_subset < S_j22.size())
//					{
//						cutPool.add(BDcut_newj22 + yF[k_Si][i] <= count_Sj22_subset);
//						//cbcut_multicut_count++;
//					}
//
//					BDcut_newj22.end();
//					S_j22.clear();
//					vector<int>().swap(S_j22);
//				}
//#endif
			}
#ifdef  constructive_cut
			else if (Sk_subtour[k_Si].size() >= 3 && count_ST[i] > 1
				&& precQ0[i] + nbQ_new[i] + Ctime_ST[i] + Ctime_subtour[k_Si] > QC_ub[k_Si])
			{
				{
					IloExpr epa(env);
					for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
					{
						//epa += xF_end[k_Si][S_0T[k_Si][i2]];
						for (int i3 = 1; i3 < Sk_subtour[k_Si].size(); i3++)
						{
							epa += xF[k_Si][S_0T[k_Si][i2]][Sk_subtour[k_Si][i3]] + xF[k_Si][Sk_subtour[k_Si][i3]][S_0T[k_Si][i2]];
						}
					}
					int coutsubtour = 0;
					coutsubtour = count_ST[i];
					int difference_UB;
					difference_UB = precQ0[i] + nbQ_new[i] + Ctime_ST[i] - QC_ub[k_Si];
					int indi_UB = 0;
					for (int i2 = 1; i2 < Sk_subtour[k_Si].size(); i2++)
					{
						if (indi_UB)
							epa -= yF[k_Si][Sk_subtour[k_Si][i2]];
						else
							coutsubtour++;
						epa += xF_end[k_Si][Sk_subtour[k_Si][i2]];
						if (i2 < Sk_subtour[k_Si].size() - 1)
							for (int i3 = i2 + 1; i3 < Sk_subtour[k_Si].size(); i3++)
							{
								epa += xF[k_Si][Sk_subtour[k_Si][i2]][Sk_subtour[k_Si][i3]] + xF[k_Si][Sk_subtour[k_Si][i3]][Sk_subtour[k_Si][i2]];
							}
						difference_UB += nbQ_new[Sk_subtour[k_Si][i2]];
						if(difference_UB>0)
							indi_UB = 1;
					}
					cutPool.add(BDcut_ST[i] + yF[k_Si][i] + epa <= coutsubtour - 1);
					epa.end();
				}
			}
#endif
			///// single QC precedence cut - 针对Si0单独违反prec cut- 去尾
			if (nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + precQT[i] > UB)
			{
				int inic0 = 0;//指示是否加过sebset cut
				count_domin_cut_T[i] = 1;
				cbcut_single_count++;
				cbcut_count_innerloop++;
				S0_cbcut_single_count_innerloop++;
#ifdef  activ_subsetcut_0
				//第一种，yik
				if (count_S0[i] > 2)
				{

					//按照si0中任务xij顺序删
					int difference_UB;
					difference_UB = nbreadyT[k_Si]+ Ctime_S0[i]+ nbQ_new[i] + precQT[i] - UB;

					int count_subset = 0;
					IloExpr BDcut_new2_2(env);
					for (int n1 = task_position[i]-1; n1 >=0; n1--)
					{
						if (nbQ_new[Sk_0T[k_Si][n1]] < difference_UB)
						{
							difference_UB -= nbQ_new[Sk_0T[k_Si][n1]];
							count_subset++;
							BDcut_new2_2 += yF[k_Si][Sk_0T[k_Si][n1]];
						}
					}

					if (count_subset > 0)
					{
						cutPool.add(BDcut_S0[i] + yF[k_Si][i] - BDcut_new2_2 <= count_S0[i] - 1- count_subset);
						inic0 = 1;
						//////Alternative cuts based on immediate processing
//#ifdef alternative_cut
//						IloExpr epa(env);
//						int indic2 = 0;
//						for (int j0 = 0; j0 < nbTask; j0++)
//						{
//
//							int indic = 1;//是否包含在Si0里
//							for (int i2 = 0; i2 < task_position[i]; i2++)
//								if (Sk_0T[k_Si][i2] == j0)
//									indic = 0;
//
//							if (indic && nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[j0] + precQT[j0] > UB)
//							{
//								for (int i2 = 0; i2 < task_position[i]; i2++)
//								{
//									epa += xF[k_Si][Sk_0T[k_Si][i2]][j0];
//								}
//								indic2++;
//							}
//						}
//						if (indic2)
//						{
//							cutPool.add(BDcut_S0[i] + epa - BDcut_new2_2 <= count_S0[i] - 1 - count_subset);
//						}
//						epa.end();
//#endif
//#ifdef alternative_cut2
//						IloExpr epa(env);
//						int indic2 = 0;
//						for (int j0 = 0; j0 < nbTask; j0++)
//						{
//
//							int indic = 1;//是否包含在Si0里
//							for (int i2 = 0; i2 < task_position[i]; i2++)
//								if (Sk_0T[k_Si][i2] == j0)
//									indic = 0;
//
//							if (indic && nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[j0] + precQT[j0] > UB)
//							{
//								for (int i2 = 0; i2 < task_position[i]; i2++)
//								{
//									epa += xF[k_Si][Sk_0T[k_Si][i2]][j0];								
//								}
//								indic2++;
//							}
//						}
//						int lowind = nbTask;
//						int upind = 0;
//						for (int i2 = 0; i2 < task_position[i]; i2++)
//						{
//							if (lowind > S_0T[k_Si][i2])
//								lowind = S_0T[k_Si][i2];
//							if (upind < S_0T[k_Si][i2])
//								upind = S_0T[k_Si][i2];
//						}
//						IloExpr epa3(env);
//						for (int i3 = lowind; i3 <= upind; i3++)
//							if (yF_best[k_Si][i3] ==0)
//							{
//								epa3 += xF_begin[k_Si][i3] - yF[k_Si][i3];
//								for (int i2 = 0; i2 < task_position[i]; i2++)
//									epa3 += xF[k_Si][i3][S_0T[k_Si][i2]] + xF[k_Si][S_0T[k_Si][i2]][i3];
//							}
//						if (indic2)
//						{
//							cutPool.add(BDcut_S0[i] + epa - BDcut_new2_2 + epa3<= count_S0[i] - 1 - count_subset);
//						}
//						epa.end();
//						epa3.end();
//#endif

					}

					BDcut_new2_2.end();

				}
#endif
				if (inic0 == 0)//未加sebset cut
					cutPool.add(BDcut_S0[i] + yF[k_Si][i] <= count_S0[i] - 1);
//纯subset，没有y
//#ifdef  activ_subsetcut_1
//				if (count_S0[i] > 2)
//				{
//
//					//按照si0中任务xij顺序删
//					int count_Si2_subset;
//					vector<int> S_i2;
//					for (int n1 = 0; n1 < task_position[i]; n1++)
//						S_i2.push_back(S_0T[k_Si][n1]);
//					count_Si2_subset = S_i2.size();
//					int difference_UB;
//					difference_UB = nbreadyT[k_Si] + precQT[i] + nbQ_new[i] + Ctime_S0[i] - UB;
//
//					IloExpr BDcut_new2(env);
//					if (count_Si2_subset >= 2)
//					{
//						int ind_lin = 0;
//						while (count_Si2_subset >= 2 && ind_lin == 0)
//						{
//							if (nbQ_new[S_i2[count_Si2_subset - 1]] < difference_UB)
//							{
//								difference_UB -= nbQ_new[S_i2[count_Si2_subset - 1]];
//								count_Si2_subset--;
//							}
//							else
//								ind_lin = 1;
//						}
//
//						for (int n1 = 0; n1 < count_Si2_subset; n1++)
//						{
//							BDcut_new2 += xF_begin[k_Si][S_i2[n1]];
//							if (n1 < count_Si2_subset - 1)
//								for (int n2 = n1 + 1; n2 < count_Si2_subset; n2++)
//									BDcut_new2 += xF[k_Si][S_i2[n1]][S_i2[n2]] + xF[k_Si][S_i2[n2]][S_i2[n1]];
//						}
//					}
//					else if (count_Si2_subset == 1)
//						BDcut_new2 += xF_begin[k_Si][S_i2[0]];
//					if (count_Si2_subset < S_i2.size())
//					{
//						cutPool.add(BDcut_new2 + yF[k_Si][i] <= count_Si2_subset);
//						//cbcut_multicut_count++;
//					}
//					BDcut_new2.end();
//					S_i2.clear();
//					vector<int>().swap(S_i2);
//				}
//#endif

			}
#ifdef  constructive_cut
			else if (Sk_subtour[k_Si].size() >= 3 && count_ST[i] > 1
				&& nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + precQT[i] + Ctime_subtour[k_Si] > UB)
			{
				IloExpr epa(env);
				if (task_position[i] >= 1)
					for (int i2 = 0; i2 < task_position[i]; i2++)
					{
						for (int i3 = 1; i3 < Sk_subtour[k_Si].size(); i3++)
						{
							epa += xF[k_Si][S_0T[k_Si][i2]][Sk_subtour[k_Si][i3]] + xF[k_Si][Sk_subtour[k_Si][i3]][S_0T[k_Si][i2]];
						}
					}
				int coutsubtour = 0;
				coutsubtour = count_S0[i];
				int difference_UB;
				difference_UB = nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + precQT[i] - UB;
				int indic_ub=0;
				for (int i2 = 1; i2 < Sk_subtour[k_Si].size(); i2++)
				{
					if(indic_ub)
						epa -= yF[k_Si][Sk_subtour[k_Si][i2]];
					else
						coutsubtour++;
					epa += xF_begin[k_Si][Sk_subtour[k_Si][i2]];
					if (i2 < Sk_subtour[k_Si].size() - 1)
						for (int i3 = i2 + 1; i3 < Sk_subtour[k_Si].size(); i3++)
						{
							epa += xF[k_Si][Sk_subtour[k_Si][i2]][Sk_subtour[k_Si][i3]] + xF[k_Si][Sk_subtour[k_Si][i3]][Sk_subtour[k_Si][i2]];
						}
					difference_UB += nbQ_new[Sk_subtour[k_Si][i2]];
					if (difference_UB > 0)
						indic_ub = 1;
				}

				cutPool.add(BDcut_S0[i] + yF[k_Si][i] + epa <= coutsubtour - 1);
				epa.end();
			}
#endif
			//cout << "prec_cut added!" << endl;
#endif
//Alternative cuts by substitution
#ifdef alter_cut_substitution// reverse no good cut
			int DiffValue = 0;
			DiffValue += sumQall;//pik need modification
			for (int k2 = 0; k2 < nbCrane; k2++)
				if(k2!=k_Si)
					DiffValue -= QC_ub[k2] - nbreadyT[k2];
			//if (count_domin_cut_0[i] < 1 && precQ0[i] + DiffValue - Ctime_S0[i] > 0 )
			if (precQ0[i] + DiffValue - Ctime_S0[i] > QC_ub[k_Si] )
			{
				//cutPool.add(BDcut_S0[i] + BDcut_S0_x[i] + xF_begin[k_Si][i] <= count_S0[i] - 1);

				IloExpr epa(env);
				for (int i2 = 0; i2 < task_position[i]; i2++)
					epa -= yF[k_Si][S_0T[k_Si][i2]];
				for (int j0 = 0; j0 < nbTask; j0++)
				{
					int indic = 1;//是否包含在Si0里
					for (int i2 = 0; i2 < task_position[i]; i2++)
						if (Sk_0T[k_Si][i2] == j0)
							indic = 0;
					if (indic && precQ0[j0] + DiffValue - Ctime_S0[i] > QC_ub[k_Si])
					{
						epa += xF_begin[k_Si][j0];
						for (int i2 = 0; i2 < task_position[i]; i2++)
						{
							epa += xF[k_Si][Sk_0T[k_Si][i2]][j0];
						}
					}
				}
				cutPool.add(BDcut_S0[i] + epa <= 0);
				epa.end();
			}
			//if (count_domin_cut_T[i] < 1 && nbreadyT[k_Si] + DiffValue - Ctime_ST[i] + precQT[i] > 0 )
			if (nbreadyT[k_Si] + DiffValue - Ctime_ST[i] + precQT[i] > UB )
			{
				//cutPool.add(BDcut_ST[i] + BDcut_ST_x[i] + xF_end[k_Si][i] <= count_ST[i] - 1);
				IloExpr epa(env);
				for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
					epa -= yF[k_Si][S_0T[k_Si][i2]];
				for (int j0 = 0; j0 < nbTask; j0++)
				{
					int indic = 1;//是否包含在Si0里
					for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
						if (Sk_0T[k_Si][i2] == j0)
							indic = 0;
					if (indic && nbreadyT[k_Si] + DiffValue - Ctime_ST[i] + precQT[j0] > UB)
					{
						epa += xF_end[k_Si][j0];
						for (int i2 = task_position[i] + 1; i2 < S_0T[k_Si].size(); i2++)
						{
							epa += xF[k_Si][j0][Sk_0T[k_Si][i2]];
						}
					}
				}
				cutPool.add(BDcut_ST[i] + epa <= 0);
				epa.end();
			}
#endif

		}
//#ifdef  constructive_cut
//		else
//		{
//			if (count_domin_cut_0[S_0T[k_Si][count_S0T[k_Si]-1]] < 1 && precQ0[i] + nbQ_new[i] + Ctime_S0T[k_Si] > QC_ub[k_Si])
//			{
//				cutPool.add(BDcut_S0T_pure[k_Si] + BDcut_S0T_T[k_Si] + yF[k_Si][i] <= count_S0T[k_Si]);
//			}
//			if (count_domin_cut_0[S_0T[k_Si][0]] < 1 && precQT[i] + nbQ_new[i] + Ctime_S0T[k_Si] > UB)
//			{
//				cutPool.add(BDcut_S0T_pure[k_Si] + BDcut_S0T_0[k_Si] + yF[k_Si][i] <= count_S0T[k_Si]);
//			}
//		}
//#endif
	}

	// ////////////////////////////////////////////
	////**** no good combinatorial cut  **** /////
	// ///////////////////////////////////////////

	for (int i = 0; i < nbTask; i++)
	{
		if (task_position[i] >= 0)
		{
			int k_Si;
			k_Si = task_QC[i];

			if (i < nbTask - 1)
				for (int j = i + 1; j < nbTask; j++)
				{
					if (task_position[j] >= 0 && k_Si != task_QC[j])
					{
						int k_Sj;
						k_Sj = task_QC[j];
						//precedence no good cut
						if (nbprecR[i][j] == 1)//precedence									
						{
							if (count_domin_cut_0[i] < 1 && count_domin_cut_T[j] < 1 
								&& nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] > QC_ub[k_Sj])//precedence									
							{
								int inic0 = 0;//指示是否加过sebset cut
								//cout << "UB: " << Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] << "," << UB << endl;

								cbcut_count++;
								cbcut_count_innerloop++;

//#ifdef  activ_subsetcut_0
//								//第一种，yik
//								if (count_S0[i] > 2)
//								{
//
//									//按照si0中任务xij顺序删
//									int difference_UB;
//									difference_UB = nbreadyT[k_Si] + nbQ_new[i] + Ctime_S0[i] + nbQ_new[j] + Ctime_ST[j] - QC_ub[k_Sj];
//
//									int count_subset = 0;
//									IloExpr BDcut_new2_2(env);
//									for (int n1 = task_position[i] - 1; n1 >= 0; n1--)
//									{
//										if (nbQ_new[Sk_0T[k_Si][n1]] < difference_UB)
//										{
//											difference_UB -= nbQ_new[Sk_0T[k_Si][n1]];
//											count_subset++;
//											BDcut_new2_2 += yF[k_Si][Sk_0T[k_Si][n1]];
//										}
//									}
//
//									if (count_subset > 0)
//									{
//										cutPool.add(BDcut_S0[i]- BDcut_new2_2 + BDcut_ST[j] + yF[k_Si][i] + yF[k_Sj][j] <= count_S0[i] + count_ST[j] - 1- count_subset);
//										inic0 = 1;
//									}
//
//									BDcut_new2_2.end();
//								}
//
//#endif
//#ifdef  activ_subsetcut_0
//								////第一种	按xij原来顺序
//								if (count_ST[j] > 2)
//								{
//									int difference_UB;
//									difference_UB = nbreadyT[k_Si] + nbQ_new[i] + Ctime_S0[i] + nbQ_new[j] + Ctime_ST[j] - QC_ub[k_Sj];
//									int count_subset = 0;
//									IloExpr epa2(env);
//									for (int n2 = task_position[j] + 1; n2 < Sk_0T[k_Sj].size(); n2++)
//									{
//										if (nbQ_new[Sk_0T[k_Sj][n2]] < difference_UB)
//										{
//											difference_UB -= nbQ_new[Sk_0T[k_Sj][n2]];
//											count_subset++;
//											epa2 += yF[k_Sj][Sk_0T[k_Sj][n2]];
//										}
//									}
//									if (count_subset > 0)
//									{
//										cutPool.add(BDcut_S0[i] + yF[k_Si][i] + BDcut_ST[j] + yF[k_Sj][j] - epa2 <= count_S0[i] + count_ST[j] - 1 - count_subset);
//										//cbcut_multicut_count++;
//										inic0 = 1;
//									}
//									epa2.end();
//								}
//#endif
								if(inic0 == 0)
									cutPool.add(BDcut_S0[i] + BDcut_ST[j] + yF[k_Si][i] + yF[k_Sj][j] <= count_S0[i] + count_ST[j] - 1);
//#ifdef  activ_subsetcut_1
//								if (count_S0[i] > 2 && count_S0[i] >= count_ST[j])
//								{
//
//									//第一种，按照si0中任务xij顺序删
//									int count_Si2_subset;
//									vector<int> S_i2;
//									for (int n1 = 0; n1 < task_position[i]; n1++)
//										S_i2.push_back(S_0T[k_Si][n1]);
//									count_Si2_subset = S_i2.size();
//									int difference_UB;
//									difference_UB = nbreadyT[k_Si] + nbQ_new[i] + Ctime_S0[i] + nbQ_new[j] + Ctime_ST[j] - QC_ub[k_Sj];
//
//									IloExpr BDcut_new2(env);
//									IloExpr BDcut_new2_2(env);
//									if (count_Si2_subset >= 2)
//									{
//										int ind_lin = 0;
//										while (count_Si2_subset >= 2 && ind_lin == 0)
//										{
//											if (nbQ_new[S_i2[count_Si2_subset - 1]] < difference_UB)
//											{
//												difference_UB -= nbQ_new[S_i2[count_Si2_subset - 1]];
//												count_Si2_subset--;
//												BDcut_new2_2 += yF[k_Si][S_i2[count_Si2_subset - 1]];
//											}
//											else
//												ind_lin = 1;
//										}
//
//										for (int n1 = 0; n1 < count_Si2_subset - 1; n1++)
//										{
//											BDcut_new2 += xF_begin[k_Si][S_i2[n1]];
//											for (int n2 = n1 + 1; n2 < count_Si2_subset; n2++)
//												BDcut_new2 += xF[k_Si][S_i2[n1]][S_i2[n2]] + xF[k_Si][S_i2[n2]][S_i2[n1]];
//										}
//										BDcut_new2 += xF_begin[k_Si][S_i2[count_Si2_subset - 1]];
//									}
//									else if (count_Si2_subset == 1)
//										BDcut_new2 += xF_begin[k_Si][S_i2[0]];
//									if (count_Si2_subset < S_i2.size())
//									{
//										//cutPool.add(BDcut_new2 + yF[k_Si][i] + BDcut_ST[j] + yF[k_Sj][j] <= count_Si2_subset + count_ST[j]);
//										cutPool.add(BDcut_S0[i] - BDcut_new2_2 + yF[k_Si][i] + BDcut_ST[j] + yF[k_Sj][j] <= count_Si2_subset + count_ST[j]);
//									}
//									BDcut_new2.end();
//									BDcut_new2_2.end();
//									S_i2.clear();
//									vector<int>().swap(S_i2);
//								}
//#endif
//#ifdef  activ_subsetcut_1
//
//								////第一种	按xij原来顺序
//								if (count_ST[j] > 2 && count_S0[i] <= count_ST[j])
//								{
//									vector<int> S_j22;
//									for (int n2 = task_position[j] + 1; n2 < Sk_0T[k_Sj].size(); n2++)
//										S_j22.push_back(Sk_0T[k_Sj][n2]);
//									int difference_UB;
//									difference_UB = nbreadyT[k_Si] + nbQ_new[i] + Ctime_S0[i] + nbQ_new[j] + Ctime_ST[j] - QC_ub[k_Sj];
//									int count_Sj22_subset, count_Sj22_subset0;
//									count_Sj22_subset = S_j22.size();
//									count_Sj22_subset0 = S_j22.size();
//									IloExpr BDcut_newj22(env);
//									IloExpr BDcut_newj22_2(env);
//									if (count_Sj22_subset >= 2)
//									{
//										int ind_lin = 0;
//										while (count_Sj22_subset >= 2 && ind_lin == 0)
//										{
//											if (nbQ_new[S_j22[count_Sj22_subset0 - count_Sj22_subset]] < difference_UB)
//											{
//												difference_UB -= nbQ_new[S_j22[count_Sj22_subset0 - count_Sj22_subset]];
//												count_Sj22_subset--;
//												BDcut_newj22_2 += yF[k_Sj][S_j22[count_Sj22_subset0 - count_Sj22_subset]];
//											}
//											else
//												ind_lin = 1;
//										}
//
//										for (int n1 = count_Sj22_subset0 - count_Sj22_subset; n1 < count_Sj22_subset0 - 1; n1++)
//										{
//											BDcut_newj22 += xF_end[k_Sj][S_j22[n1]];
//											for (int n2 = n1 + 1; n2 < count_Sj22_subset0; n2++)
//												BDcut_newj22 += xF[k_Sj][S_j22[n1]][S_j22[n2]] + xF[k_Sj][S_j22[n2]][S_j22[n1]];
//										}
//										BDcut_newj22 += xF_end[k_Sj][S_j22[count_Sj22_subset0 - 1]];
//									}
//									else if (count_Sj22_subset == 1)
//										BDcut_newj22 += xF_end[k_Sj][S_j22[count_Sj22_subset0 - 1]];
//
//									if (count_Sj22_subset < S_j22.size())
//									{
//										//cutPool.add(BDcut_S0[i] + yF[k_Si][i] + BDcut_newj22 + yF[k_Sj][j] <= count_S0[i] + count_Sj22_subset);
//										cutPool.add(BDcut_S0[i] + yF[k_Si][i] + BDcut_ST[j] + yF[k_Sj][j] - BDcut_newj22_2 <= count_S0[i] + count_Sj22_subset);
//										//cbcut_multicut_count++;
//									}
//									BDcut_newj22.end();
//									BDcut_newj22_2.end();
//									S_j22.clear();
//									vector<int>().swap(S_j22);
//								}
//#endif
							}
#ifdef ub_cbcut_domin
							int indi_Sj0 = 1;
							int indi_SiT = 1;
							for (int n2 = task_position[i] + 1; n2 < S_0T[k_Si].size(); n2++)
							{
								if (para_Delta[k_Si][k_Sj][S_0T[k_Si][n2]][j] < 0)
									indi_SiT = 0;
							}
							for (int n2 = 0; n2 < task_position[j]; n2++)
							{
								if (para_Delta[k_Sj][k_Si][S_0T[k_Sj][n2]][i] < 0)
									indi_Sj0 = 0;
							}

							//第一种，Si0 and SiT
							if (indi_SiT && count_domin_cut_0[i] < 1 && count_domin_cut_T[i] < 1 &&
								nbreadyT[k_Si] + Ctime_S0[i] + Ctime_ST[i] + nbQ_new[i] + nbQ_new[j] > UB)
								{
									if (count_ST[i] >= 2)
									{
										int totalsum = 0;
										int sumk[nbCrane];//记录可加入cut的QC，1为可加入
										for (int k = 0; k < nbCrane; k++)
										{
											sumk[k] = 0;
											if (k != k_Si)
											{
												int indiSiT = 1;
												for (int n2 = task_position[i] + 1; n2 < Sk_0T[k_Si].size(); n2++)
												{
													if (para_Delta[k_Si][k][Sk_0T[k_Si][n2]][j] < 0)
														indiSiT = 0;
												}
												if (indiSiT)
													sumk[k] = 1;

												totalsum += sumk[k];
											}
										}

										if (totalsum)
										{
											IloExpr cuty(env);
											for (int k = 0; k < nbCrane; k++)
												if (k != k_Si && sumk[k] == 1)
													cuty += yF[k][j];
											int inic0 = 0;//指示是否加过sebset cut
											
//#ifdef  activ_subsetcut_0
//											//第一种，yik
//											if (count_S0[i] > 2)
//											{
//
//												//按照si0中任务xij顺序删
//												int difference_UB;
//												difference_UB = nbreadyT[k_Si] + nbQ_new[i] + Ctime_S0[i] + nbQ_new[j] + Ctime_ST[i] - UB;
//
//												int count_subset = 0;
//												IloExpr BDcut_new2_2(env);
//												for (int n1 = task_position[i] - 1; n1 >= 0; n1--)
//												{
//													if (nbQ_new[Sk_0T[k_Si][n1]] < difference_UB)
//													{
//														difference_UB -= nbQ_new[Sk_0T[k_Si][n1]];
//														count_subset++;
//														BDcut_new2_2 += yF[k_Si][Sk_0T[k_Si][n1]];
//													}
//												}
//
//												if (count_subset > 0)
//												{
//													cutPool.add(BDcut_S0[i] - BDcut_new2_2+ BDcut_ST[i] + yF[k_Si][i] + cuty <= count_S0[i] + count_ST[i] - 1- count_subset);
//													inic0 = 1;
//												}
//
//												BDcut_new2_2.end();
//											}
//
//#endif
#ifdef  activ_subsetcut_0
											////第一种	按xij原来顺序
											if (count_ST[i] > 2)
											{
												int difference_UB;
												difference_UB = nbreadyT[k_Si] + nbQ_new[i] + Ctime_S0[i] + nbQ_new[j] + Ctime_ST[i] - UB;
												int count_subset=0;
												IloExpr epa2(env);
												if (count_ST[i] >= 2)
												{
													for (int n2 = task_position[i] + 1; n2 < Sk_0T[k_Si].size(); n2++)
													{
														if (nbQ_new[Sk_0T[k_Si][n2]] < difference_UB)
														{
															difference_UB -= nbQ_new[Sk_0T[k_Si][n2]];
															count_subset++;
															epa2 += yF[k_Si][Sk_0T[k_Si][n2]];
														}
													}
												}

												if (count_subset>0)
												{
													cutPool.add(BDcut_S0[i] + yF[k_Si][i] + BDcut_ST[i] -epa2 + cuty <= count_S0[i] + count_ST[i] -1 - count_subset);
													//cbcut_multicut_count++;
													inic0 = 1;

#ifdef alter_cut_substitution
													int DiffValue = 0;
													DiffValue += sumQall;//pik need modification
													for (int k2 = 0; k2 < nbCrane; k2++)
															DiffValue -= QC_ub[k2]- nbreadyT[k2];
													if (nbreadyT[k_Si] + nbQ_new[j] + DiffValue > 0)
														cutPool.add(BDcut_ST[i] + BDcut_ST_x[i] + xF_end[k_Si][i] - epa2 + cuty <= count_ST[i]- count_subset);

#endif
												}
												epa2.end();

											}
#endif
											if(inic0 == 0)
												cutPool.add(BDcut_S0[i] + BDcut_ST[i] + yF[k_Si][i] + cuty <= count_S0[i] + count_ST[i] - 1);
	
											//cbcut_count++;
											//cbcut_count_innerloop++;
											cuty.end();
										}
									}
									else if (task_position[i] == Sk_0T[k_Si].size() - 1)
									{
										IloExpr cuty(env);
										for (int k = 0; k < nbCrane; k++)
											if (k != k_Si)
												cuty += yF[k][j];
										cutPool.add(BDcut_S0[i] + xF_end[k_Si][i] + cuty <= count_S0[i]);
										//cbcut_count++;
										//cbcut_count_innerloop++;
										cuty.end();
									}

								}
							//第二种，Sj0 and SjT
							if (indi_Sj0 && count_domin_cut_0[j] < 1 && count_domin_cut_T[j] < 1
								&& nbreadyT[k_Sj] + Ctime_S0[j] + Ctime_ST[j] + nbQ_new[i] + nbQ_new[j] > QC_ub[k_Sj])
								{
									if (count_S0[j] >= 2)
									{
										int totalsum = 0;
										int sumk[nbCrane];//记录可加入cut的QC，1为可加入
										for (int k = 0; k < nbCrane; k++)
										{
											sumk[k] = 0;
											if (k != k_Sj)
											{
												int indiSjT = 1;
												for (int n2 = 0; n2 < task_position[j]; n2++)
												{
													if (para_Delta[k_Sj][k][Sk_0T[k_Sj][n2]][i] < 0)
														indiSjT = 0;
												}
												if (indiSjT)
													sumk[k] = 1;

												totalsum += sumk[k];
											}
										}

										if (totalsum)
										{
											int inic0 = 0;//指示是否加过sebset cut
											IloExpr cuty(env);
											for (int k = 0; k < nbCrane; k++)
												if (k != k_Sj && sumk[k] == 1)
													cuty += yF[k][i];
											
#ifdef  activ_subsetcut_0
											//第一种，yik
											if (count_S0[i] > 2)
											{

												//按照si0中任务xij顺序删
												int difference_UB;
												difference_UB = nbreadyT[k_Sj] + Ctime_S0[j] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] - QC_ub[k_Sj];

												int count_subset = 0;
												IloExpr BDcut_new2_2(env);
												for (int n1 = task_position[j] - 1; n1 >= 0; n1--)
												{
													if (nbQ_new[Sk_0T[k_Sj][n1]] < difference_UB)
													{
														difference_UB -= nbQ_new[Sk_0T[k_Sj][n1]];
														count_subset++;
														BDcut_new2_2 += yF[k_Sj][Sk_0T[k_Sj][n1]];
													}
												}

												if (count_subset > 0)
												{
													cutPool.add(BDcut_S0[j]  - BDcut_new2_2+ BDcut_ST[j] + yF[k_Sj][j] + cuty <= count_S0[j] + count_ST[j] - 1- count_subset);
													inic0 = 1;
#ifdef alter_cut_substitution
													int DiffValue = 0;
													DiffValue += sumQall;//pik need modification
													for (int k2 = 0; k2 < nbCrane; k2++)
														DiffValue -= QC_ub[k2]- nbreadyT[k2];
													if ( nbreadyT[k_Sj] + nbQ_new[i] + DiffValue > 0)							
														cutPool.add(BDcut_S0[j] + BDcut_S0_x[j] + xF_begin[k_Sj][j] - BDcut_new2_2  + cuty <= count_S0[j]- count_subset);
#endif
												}

												BDcut_new2_2.end();
											}

#endif
//#ifdef  activ_subsetcut_0
//////第一种	按xij原来顺序
//											if (count_ST[j] > 2)
//											{
//												int difference_UB;
//												difference_UB = nbreadyT[k_Sj] + Ctime_S0[j] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] - QC_ub[k_Sj];
//												int count_subset = 0;
//												IloExpr epa2(env);
//												for (int n2 = task_position[j] + 1; n2 < Sk_0T[k_Sj].size(); n2++)
//												{
//													if (nbQ_new[Sk_0T[k_Sj][n2]] < difference_UB)
//													{
//														difference_UB -= nbQ_new[Sk_0T[k_Sj][n2]];
//														count_subset++;
//														epa2 += yF[k_Sj][Sk_0T[k_Sj][n2]];
//													}
//												}
//												if (count_subset > 0)
//												{
//													cutPool.add(BDcut_S0[j]  + BDcut_ST[j] + yF[k_Sj][j] + cuty - epa2 <= count_S0[j] + count_ST[j] - 1 - count_subset);
//													//cbcut_multicut_count++;
//													inic0 = 1;
//												}
//												epa2.end();
//											}
//#endif
											if (inic0 == 0)
												cutPool.add(BDcut_S0[j] + BDcut_ST[j] + yF[k_Sj][j] + cuty <= count_S0[j] + count_ST[j] - 1);
											cuty.end();
										}
									}
								}
#endif

								//#ifdef alter_cut_substitution // for si0 and sjT
								//								int DiffValue2 = 0;
								//								DiffValue2 += sumQall;//pik need modification
								//								for (int k2 = 0; k2 < nbCrane; k2++)
								//										DiffValue2 -= QC_ub[k2];
								//								if (count_domin_cut_0[i] < 1 && count_domin_cut_0[j] < 1 &&
								//									nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + DiffValue2 - Ctime_S0[j]>0)
								//								{
								//									cutPool.add(BDcut_S0[j] + BDcut_S0_x[j] + xF_begin[k_Sj][j] + yF[k_Si][i] + BDcut_S0[i] <= count_S0[i] + count_S0[j] - 1);
								//								}
								//#endif

//#ifdef alter_cut_substitution // for si0 and sjT
//								//int DiffValue2 = 0;
//								//DiffValue2 += sumQall;//pik need modification
//								//for (int k2 = 0; k2 < nbCrane; k2++)
//								//		DiffValue2 -= QC_ub[k2];
//								//if (count_domin_cut_0[i] < 1 && count_domin_cut_0[j] < 1 &&
//								//	nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + DiffValue2 - Ctime_S0[j]>0)
//								if (count_domin_cut_0[i] < 1 && count_domin_cut_T[j] < 1 && count_S0[j]< count_ST[j]
//									&& nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] > QC_ub[k_Sj])//precedence
//								{
//									int count_subset = 0;
//									IloExpr BDcut_new2_2(env);
//									for (int n1 = task_position[j] - 1; n1 >= 0; n1--)
//									{
//											count_subset++;
//											BDcut_new2_2 += yF[k_Sj][Sk_0T[k_Sj][n1]];									
//									}
//
//									cutPool.add(BDcut_S0[j] - BDcut_new2_2  + BDcut_S0_x[j] + xF_begin[k_Sj][j] + yF[k_Si][i] + BDcut_S0[i] <= count_S0[i] + count_S0[j] - count_subset - 1);
//									BDcut_new2_2.end();
//								}
//#endif
						}
						//else if (para_Delta[k_Si][k_Sj][i][j] > 0 && count_domin_cut_0[i] < 1 && count_domin_cut_T[i] < 1
						//	&& count_domin_cut_0[j] < 1 && count_domin_cut_T[j] < 1)
						//	{
						//		//							if (nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] > QC_ub[k_Sj]
						//		//								&& nbreadyT[k_Sj] + Ctime_S0[j] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[i] > QC_ub[k_Si])
						//		//							{
						//		//#ifdef ub_cbcut
						//		//								//basic no good cut
						//		//								cutPool.add(BDcut_S0[i] + BDcut_ST[i] + BDcut_S0[j] + BDcut_ST[j] + yF[k_Si][i] + yF[k_Sj][j]
						//		//									<= count_S0[i] + count_ST[i] + count_S0[j] + count_ST[j] - 3);
						//		//#endif
						//		//								//cbcut_count++;
						//		//								//cbcut_count_innerloop++;
						//		//							}
						//		int DiffValuei, DiffValuej;
						//		DiffValuei = DiffValuej = 0;
						//		DiffValuej += sumQall;//pik need modification
						//		for (int k2 = 0; k2 < nbCrane; k2++)
						//		{
						//			if (k2 != k_Si)
						//				DiffValuei -= QC_ub[k2];
						//			if (k2 != k_Sj)
						//				DiffValuej -= QC_ub[k2];
						//		}
						//		int indi_Si = 1;
						//		for (int n1 = 0; n1 < task_position[i]; n1++)
						//		{
						//			for (int n2 = 0; n2 < task_position[j]; n2++)
						//				if (para_Delta[k_Si][k_Sj][S_0T[k_Si][n1]][S_0T[k_Sj][n2]] < 0 || para_Delta[k_Sj][k_Si][S_0T[k_Sj][n2]][S_0T[k_Si][n1]] < 0)
						//					indi_Si = 0;
						//		}
						//		int count_subset2 = 0;
						//		IloExpr BDcut_new2_2(env);
						//		for (int n1 = task_position[j] - 1; n1 >= 0; n1--)
						//		{
						//			count_subset2++;
						//			BDcut_new2_2 += yF[k_Sj][Sk_0T[k_Sj][n1]];
						//		}
						//		int count_subset = 0;
						//		IloExpr BDcut_new2_1(env);
						//		for (int n1 = task_position[i] - 1; n1 >= 0; n1--)
						//		{
						//			count_subset++;
						//			BDcut_new2_1 += yF[k_Si][Sk_0T[k_Si][n1]];
						//		}
						//		if (indi_Si && nbreadyT[k_Si] + Ctime_S0[i] + DiffValuej > UB && nbreadyT[k_Si]==0 && nbreadyT[k_Sj] == 0
						//			&& nbreadyT[k_Sj] + Ctime_S0[j] + DiffValuei > UB)
						//		{
						//			cutPool.add(BDcut_S0[i] + yF[k_Si][i] + yF[k_Si][j] + BDcut_S0[j] - BDcut_new2_1  - BDcut_new2_2 <= count_S0[i] + count_S0[j] - 1- count_subset- count_subset2);
						//		}
						//		BDcut_new2_1.end();
						//		BDcut_new2_2.end();
						//}
					}

#ifndef prec_cut
					if (nbprecR[i][j] == 1 && task_position[j] >= 0 && k_Si != task_QC[j]&& count_Si0 >= 0 && count_SjT >= 0)
					{
						cb2016_count++;
						//cb2016_count_innerloop++;
						cbcut_count_innerloop++;
						//if (activated_2016 == 1 && record_preceding_task[k_Si][i] >= 0 && record_next_task[k_Sj][j] >= 0)
						cutPool.add(BDcut_Si0 + BDcut_Si0_x + xF_begin[k_Si][i]
							+ xF_begin[k_Sj][j] + BDcut_SjT_x + BDcut_SjT <= count_Si0 + count_SjT - 1);
					}

#endif

#ifdef safe_margin_viocut
				//#ifdef BenchMark_instance
				//					if (para_Delta[k_Si][k_Sj][i][j] > 0)
				//#endif
				//#ifndef BenchMark_instance
					if (nbprecR[i][j] == 1 && task_position[j] >= 0 && k_Si != task_QC[j])
						//#endif
					{
						int k_Sj;
						k_Sj = task_QC[j];
						////case 1
						//IloExpr epals1(env);
						//for (int n1 = task_position[i] + 1; n1 < Sk_0T[k_Si].size(); n1++)
						//	epals1 += yF[k_Si][Sk_0T[k_Si][n1]];
						//int indicSi = 1;
						//if (task_position[i] < Sk_0T[k_Si].size() - 1)
						//	for (int n3 = task_position[i] + 1; n3 < Sk_0T[k_Si].size(); n3++)
						//	{
						//		if ((Sk_0T[k_Si][n3] < j && para_Delta[k_Si][k_Sj][Sk_0T[k_Si][n3]][j] < 0)
						//			|| (Sk_0T[k_Si][n3] > j && para_Delta[k_Sj][k_Si][j][Sk_0T[k_Si][n3]] < 0))
						//			indicSi = 0;
						//	}
						//if (indicSi && task_position[i] < Sk_0T[k_Si].size() - 1 && count_domin_cut_T[i] < 1 && dT[k_Si]>UB)
						//{
						//	//cutPool.add(BDcut_ST[i] + BDcut_ST_x[i] + yF[k_Sj][j] <= count_ST[i]);
						//	cutPool.add(BDcut_ST[i] + BDcut_ST_x[i] + yF[k_Sj][j] - epals1 <= 1);
						//	cbcut_count++;
						//}
						//epals1.end();
						//case 2
						IloExpr epals(env);
						for (int n1 = 0; n1 < task_position[j]; n1++)
							epals += yF[k_Sj][Sk_0T[k_Sj][n1]];
						int indicSj = 1;
						if (task_position[j] > 0)
							for (int n3 = 0; n3 < task_position[j]; n3++)
							{
								if ((Sk_0T[k_Sj][n3] > i && para_Delta[k_Si][k_Sj][i][Sk_0T[k_Sj][n3]] < 0)
									|| (Sk_0T[k_Sj][n3] < i && para_Delta[k_Sj][k_Si][Sk_0T[k_Sj][n3]][i] < 0))
									indicSj = 0;
							}
						if (indicSj && task_position[j] > 0 && count_domin_cut_0[j] < 1 && nbreadyT[k_Sj]==0 && nbQ_new[i]>=5)
						{
							//cutPool.add(BDcut_S0[i] + BDcut_S0_x[i] + yF[k_Si][i] <= count_S0[i]);
							cutPool.add(BDcut_S0[j] + BDcut_S0_x[j] + yF[k_Si][i] -epals<= 1);
							cbcut_count++;
							//cbcut_count_innerloop++;
						}
						epals.end();
						//case 3
						if (j < nbTask - 1)
						{
							for (int iij = j + 1; iij < nbTask; iij++)
							{
								if (nbprecR[j][iij] == 1 && yF_best[k_Si][iij] == 1)
								{
									int iij_found = -1;
									for (int n3 = task_position[i] + 1; n3 < Sk_0T[k_Si].size(); n3++)
										if (Sk_0T[k_Si][n3] == iij)
											iij_found = n3;
									if (iij_found >= task_position[i] + 1)
									{
										int indicSi2 = 1;
										for (int n3 = task_position[i] + 1; n3 <= iij_found; n3++)
											if (para_Delta[k_Si][k_Sj][Sk_0T[k_Si][n3]][j] <= 0)
												indicSi2 = 0;
										if (indicSi2)
										{
											IloExpr BDcut_safemargin(env);
											for (int n3 = task_position[i]; n3 < iij_found; n3++)
											{
												for (int n4 = n3 + 1; n4 <= iij_found; n4++)
													BDcut_safemargin += xF[k_Si][Sk_0T[k_Si][n3]][Sk_0T[k_Si][n4]] + xF[k_Si][Sk_0T[k_Si][n4]][Sk_0T[k_Si][n3]];
											}
											cutPool.add(BDcut_safemargin + yF[k_Sj][j] <= iij_found - task_position[i]);
											cbcut_count++;
											//cbcut_count_innerloop++;
											BDcut_safemargin.end();
										}
									}

									//if (dT[k_Si] > UB && k_Si > k_Sj)
									//{
									//	int indopt = 1;
									//	for (int taskindi = j; taskindi < nbTask; taskindi++)
									//		if (nbLocation_new[taskindi] > nbLocation_new[i] + safe_margin
									//			&& yF_best[k_Si][taskindi] > threshold_M)
									//			indopt = 0;
									//	//int indopt2 = 0;
									//	//for (int taskindi = j; taskindi < nbTask; taskindi++)
									//	//	if (nbLocation_new[taskindi] <= nbLocation_new[i] + safe_margin
									//	//		&& nbLocation_new[taskindi] >= nbLocation_new[i]
									//	//		&& yF_best[k_Si][taskindi] > threshold_M)
									//	//		indopt2 = 1;
									//	//if (indopt && indopt2)
									//	if (indopt)//全部右侧task都在safetymargin里，可以加cut
									//	{
									//		IloExpr BDcut_safemargin2(env);
									//		BDcut_safemargin2 += yF[k_Si][i] + yF[k_Si][iij] - 2;

									//		for (int ksf = 0; ksf < k_Si; ksf++)
									//			BDcut_safemargin2 += yF[ksf][j];

									//		for (int taskindi = j; taskindi < nbTask; taskindi++)
									//			if (nbLocation_new[taskindi] > nbLocation_new[i] + safe_margin)
									//				BDcut_safemargin2 -= yF[k_Si][taskindi];
									//		cutPool.add(BDcut_safemargin2 <= 0);
									//		BDcut_safemargin2.end();
									//	}
									//}
								}
							}
						}
					}



#endif

#ifdef preccut_tripleQCs
					if (nbprecR[i][j] == 1 && count_domin_cut_0[i] < 1)//precedence multiple QCs
					{
						if (i < nbTask - 3 && j < nbTask - 2)//precedence
						{
							int k_Sj;
							for (int k = 0; k < nbCrane; k++)
							{
								if (yF_best[k][j] >= threshold_M)
									k_Sj = k;
							}
							//if (k_Si != k_Sj && subtour_count[k_Si] == 0 && subtour_count[k_Sj] == 0)//不涉及子环
							if (k_Si != k_Sj)
							{
								for (int i2 = j + 1; i2 < nbTask - 1; i2++)
									if (nbLocation_new[i] != nbLocation_new[i2])
									{
										for (int j2 = i2 + 1; j2 < nbTask; j2++)
										{
											if (nbprecR[i2][j2] == 1 && count_domin_cut_T[j2] < 1)//precedence
											{
												int k_Si2;
												int k_Sj2;


												for (int k = 0; k < nbCrane; k++)
												{
													if (yF_best[k][i2] >= threshold_M)
														k_Si2 = k;
													if (yF_best[k][j2] >= threshold_M)
														k_Sj2 = k;
												}
												//两个 prec pair
												if (k_Sj == k_Si2 && k_Sj2 != k_Si2)//满足条件
												{

													//if (inf_route_prec_cut)
													{
														IloExpr BDcut_Sji2(env);
														int ind_check_ji2 = 0;

														//S_i2
														vector<int> S_ji2;
														int Ctime_ji2 = 0;
														int current_i2;
														int check_j_i2 = 0;
														int count_Sji2 = 0;
														current_i2 = i2;
														while (current_i2 != j && current_i2 >= 0 && ind_check_ji2 == 0)
														{
															current_i2 = record_preceding_task[k_Sj][current_i2];

															if (current_i2 == j)
																check_j_i2 = 1;

															else if (current_i2 >= 0)
															{
																if (S_ji2.size() > 0)
																	for (int n2 = 0; n2 < S_ji2.size(); n2++)
																	{
																		BDcut_Sji2 += xF[k_Sj][current_i2][S_ji2[n2]] + xF[k_Sj][S_ji2[n2]][current_i2];

																	}
																S_ji2.push_back(current_i2);
																BDcut_Sji2 += xF[k_Sj][current_i2][i2] + xF[k_Sj][j][current_i2];
																Ctime_ji2 += nbQ_new[current_i2];

															}
															if (current_i2 == i2)
															{//有回环
																ind_check_ji2 = 1;
															}

														}
														count_Sji2 = S_ji2.size();

														if (check_j_i2 == 1 && ind_check_ji2 == 0
															&& nbreadyT[k_Si] + Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ji2 + nbQ_new[i2] + nbQ_new[j2] + Ctime_ST[j2] > QC_ub[k_Sj2])
														{
															cutPool.add(BDcut_S0[i] + BDcut_ST[j2] + yF[k_Si][i] + BDcut_Sji2 + yF[k_Sj][j2] <= count_S0[i] + count_ST[j2] + count_Sji2);
															//cbcut_count_innerloop++;
														}

														BDcut_Sji2.end();
														S_ji2.clear();
														vector<int>().swap(S_ji2);

													}

												}
											}
										}
									}
							}
						}
					}
#endif
				}
		}
//		else
//		{
//			int k_Si;
//			k_Si = task_QC[i];
//
//			if (i < nbTask - 1)
//				for (int j = i + 1; j < nbTask; j++)
//				{
//
//					if (count_domin_cut_0[i] < 1 && count_domin_cut_T[j] < 1 && task_position[j] >= 0 && k_Si != task_QC[j])
//					{
//						int k_Sj;
//						k_Sj = task_QC[j];
//
//						if (nbprecR[i][j] == 1)//precedence									
//						{
//							if (nbreadyT[k_Si] + Ctime_S0T[k_Si] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] > QC_ub[k_Sj])//precedence									
//							{
//								//cout << "UB: " << Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] << "," << UB << endl;
//#ifdef prec_cut
//								cutPool.add(BDcut_S0T_pure[k_Si] + BDcut_S0T_0[k_Si] + BDcut_ST[j] + yF[k_Si][i] + yF[k_Sj][j] <= count_S0T[k_Si] + count_ST[j]);
//#endif
//								cbcut_count++;
//								//cbcut_count_innerloop++;
//
//							}
//						}
//
//					}
////					else if (count_domin_cut_0[i] < 1 && count_domin_cut_T[j] < 1 && k_Si != task_QC[j])
////					{
////						int k_Sj;
////						k_Sj = task_QC[j];
////						if (nbprecR[i][j] == 1)//precedence									
////						{
////							if (nbreadyT[k_Si] + Ctime_S0T[k_Si] + nbQ_new[i] + nbQ_new[j] + Ctime_S0T[k_Sj] > QC_ub[k_Sj])//precedence									
////							{
////								//cout << "UB: " << Ctime_S0[i] + nbQ_new[i] + nbQ_new[j] + Ctime_ST[j] << "," << UB << endl;
////#ifdef prec_cut
////								cutPool.add(BDcut_S0T_pure[k_Si] + BDcut_S0T_0[k_Si] + BDcut_S0T_pure[k_Sj] + BDcut_S0T_T[k_Sj] + yF[k_Si][i] + yF[k_Sj][j] <= count_S0T[k_Si] + count_S0T[k_Sj] + 1);
////#endif
////								cbcut_count++;
////								cbcut_count_innerloop++;
////
////							}
////						}
////
////					}
//				}
//		}
	}


	// ////////////////////////////////////////////

#ifdef opt_property_cut
	for (int k = 1; k < nbCrane; k++)
	{
		int task_added = -1;
		int sumQi = 0;
		for (int i = 0; i < S_0T[k].size(); i++)
		{
			if (nbLocation_new[i] == (1 + safe_margin) * k + 1)
				task_added = i;
			else
				break;
		}

		if (task_added >= 0)//可加cut
		{
			IloExpr epa(env);
			for (int i = 0; i <= task_added; i++)
			{
				epa += xF_begin[k][i];
				sumQi += nbQ_new[i];
				if (i < task_added)
					for (int j = i + 1; j <= task_added; j++)
						epa += xF[k][j][i] + xF[k][i][j];
			}
			int sumQj = 0;
			int indi_a = -1;
			for (int j = 0; j < nbTask; j++)
			{
				if (nbLocation_new[j] == (1 + safe_margin) * (k - 1) + 1)
				{
					sumQj += nbQ_new[j];
					if (indi_a == -1 && sumQj < sumQi && first_task[k - 1] != j)
					{
						cutPool.add(epa - task_added - xF_begin[k - 1][j] <= 0);
					}
					else if (sumQj < sumQi && record_next_task[k - 1][indi_a] != j)
					{
						cutPool.add(epa - task_added - xF[k - 1][indi_a][j] <= 0);
					}
					indi_a = j;
				}
				else if (indi_a >= 0)
					break;
			}

			epa.end();
		}
	}
#endif

	cout << "cut adding end." << endl;


	BDcut_S0.end();
	BDcut_ST.end();
	BDcut_S0_pure.end();
	BDcut_ST_pure.end();
	BDcut_S0_x.end();
	BDcut_ST_x.end();
	BDcut_S0T_0.end();
	BDcut_S0T_T.end();
	BDcut_S0T_pure.end();

	//释放内层
	for (auto& i : S_0T) {
		vector<int>().swap(i);
	}
	// 释放外层
	vector<vector<int>>().swap(S_0T);

	//释放内层
	for (auto& i : Sk_0T) {
		vector<int>().swap(i);
	}
	// 释放外层
	vector<vector<int>>().swap(Sk_0T);

	//释放内层
	for (auto& i : Sk_subtour) {
		vector<int>().swap(i);
	}
	// 释放外层
	vector<vector<int>>().swap(Sk_subtour);

	if (cbcut_count_innerloop == 0
		&& S0_cbcut_single_count_innerloop == 0
		&& ST_cbcut_single_count_innerloop == 0
		//&& preccut_count_2016_innerloop == 0
		//&& cb2016_count_innerloop == 0
		&& subtour_count_sum == 0)//当没有子环也没有cb cut的时候，可以解子问题
		return true;
	else return false;
}

bool CBMP_bid(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation,
	IloNum* ObjVal, IloNumArray2 yF2_current, IloNumArray2 yF_current2)
{

	//xF: forward trip, yF: retracing trip


	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//retracing time
	IloNumVarArray thetaF2(env, nbCrane, 0, 100);


	BoolVarMatrix xF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		xF[k] = IloBoolVarArray(env, nbTask);
	}

	BoolVarMatrix yF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		yF[k] = IloBoolVarArray(env, nbTask);
	}

	BoolVarMatrix uF(env, nbTask);
	for (i = 0; i < nbTask; i++)
	{
		uF[i] = IloBoolVarArray(env, nbTask);
	}

	IloBoolVarArray vF(env, nbTask);
	IloIntVarArray thetaF(env, nbCrane, 0, safe_margin * nbBay);

	BoolVarMatrix zF(env, nbCrane);
	for (i = 0; i < nbCrane; i++)
	{
		zF[i] = IloBoolVarArray(env, nbBay);
	}
	IloIntVar   CF(env, 0, 2000);
	IloNumVarArray QC_CF(env, nbCrane, 0, IloInfinity);


	NumVarMatrix wF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		wF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CwF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CwF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix endCzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		endCzF[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);


	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//			obj2 += uF[i][j];

	//	将目标函数加入到原问题模型
	model.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//
	//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

	////////////////测试////////////////
	//for (k = 0; k < nbCrane; k++)
	//	for (i = 5*(k+1); i < nbBay; i++)
	//		if (i < nbBay)
	//		{
	//			model.add(zF[k][i] == 0);
	//			for (j = 0; j < nbTask; j++)
	//				if (nbLocation[j]==i+1)
	//					model.add(xF[k][j] + yF[k][j] == 0);
	//		}


	////////////////问题约束////////////////
	//13c
	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - gammaF[k] - nbreadyT[k] - thetaF2[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs) * (xF[k][i] + yF[k][i]);
		for (j = 0; j < nbBay; j++)
		{
			//epa -= QCmove_time*j*endCzF[k][j];
			epa -= wF[k][j] + CwF[k][j];
		}

		//for (j = 0; j < nbBay; j++)
		//	epa += QCmove_time*j*zF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	model.add(c2);
	c2.end();

	//约束 13d // 所有任务被分配到QC上
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += xF[k][i] + yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	model.add(c4);
	c4.end();

	//13e
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			model.add(yF[k][i] - vF[k] <= 0);


	////13f 13g
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] < nbBay)
			{
				IloExpr  epa(env);

				for (i = nbLocation[j]; i < nbBay; i++)
					epa += zF[k][i];
				epa += xF[k][j];

				model.add(epa <= 1);
				epa.end();

				IloExpr  epa2(env);

				for (i = nbLocation[j]; i < nbBay; i++)
					epa2 += CzF[k][i];
				epa2 += yF[k][j];

				model.add(epa2 <= 1);
				epa2.end();

				//IloExpr  epa3(env);
				//for (i = 0; i < nbLocation[j] - 1; i++)
				//	epa3 += endCzF[k][i];
				//epa3 += xF[k][j];

				//model.add(epa3 <= 1);
				//epa3.end();

			}



	//13h, 13i,ending bay
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
			//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * endCzF[k][i];
			model.add(epa - (nbLocation[j] - 1) * (xF[k][j] + yF[k][j]) >= 0);
			model.add(epa - nbLocation[j] * yF[k][j] >= 0);
			epa.end();
		}


	//约束 13j,13k : zF unique
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		for (i = 0; i < nbBay; i++)
		{
			epa += zF[k][i];
			epa2 += CzF[k][i];
			epa3 += endCzF[k][i];
		}
		epa2 -= vF[k];
		c7.add(epa == 1);
		c7.add(epa2 == 0);
		c7.add(epa3 == 1);
		epa.end();
		epa2.end();
		epa3.end();
	}
	model.add(c7);
	c7.end();

	//13l
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time * i * zF[k][i];
		model.add(epa1 + epa2 - QCmove_time * nbb[k] + QCmove_time >= 0);
		model.add(epa1 - epa2 + QCmove_time * nbb[k] - QCmove_time >= 0);
		epa1.end();
		epa2.end();
	}


	//13m
	IloRangeArray  c13m(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += thetaF2[k];
		for (j = 0; j < nbBay; j++)
		{
			epa -= QCmove_time * j * endCzF[k][j] - QCmove_time * j * zF[k][j];
		}
		c13m.add(epa >= 0);
		epa.end();
	}
	model.add(c13m);
	c13m.end();

	////13n : back time gammaF
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += gammaF[k];
		for (i = 0; i < nbBay; i++)
			epa += QCmove_time * i * CzF[k][i] - QCmove_time * i * endCzF[k][i];
		epa += QCmove_time * nbBay - QCmove_time * nbBay * vF[k];
		model.add(epa >= 0);
		epa.end();
	}

	//13o
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i * zF[k][i] - i * endCzF[k][i];
		model.add(epa <= 0);
		epa.end();
	}

	//	建立约束 13p,13q:  z 和 z之间隔 /delta +1
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
		if ((1 + safe_margin) * (k + 1) + 1 <= nbBay)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * zF[k + 1][i] - i * zF[k][i];
			c8.add(epa >= 1 + safe_margin);
			epa.end();

			IloExpr  epa2(env);
			for (i = 0; i < nbBay; i++)
				epa2 += i * CzF[k + 1][i] - i * CzF[k][i];
			epa2 += nbBay * (1 - vF[k + 1]);
			c8.add(epa2 >= 1 + safe_margin);
			epa2.end();

		}
	model.add(c8);
	c8.end();

	//13r
	for (k = 0; k < nbCrane - 1; k++)
		if ((1 + safe_margin) * (k + 2) <= nbBay)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i * endCzF[k + 1][i] - i * endCzF[k][i];
			model.add(epa >= 1 + safe_margin);
			epa.end();
		}

	for (k = 0; k < nbCrane; k++)
	{
		if ((1 + safe_margin) * k + 1 > nbBay)
		{
			for (i = 0; i < nbTask; i++)
				model.add(xF[k][i] + yF[k][i] == 0);
		}

		//else if ((1 + safe_margin)*(k+1) + 1 > nbBay && k+1< nbCrane)
		//{
		//	for (i = 0; i < nbTask; i++)
		//		if (nbLocation[i]<=(1 + safe_margin)*k)
		//			model.add(xF[k+1][i] + yF[k+1][i] == 0);
		//}
	}

	//////1+delta 版本有错误
	////	建立约束 13p,13q:  z 和 z之间隔 /delta +1
	//IloRangeArray  c8(env);
	//for (k = 0; k < nbCrane - 1; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*zF[k + 1][i] - i*zF[k][i];
	//	c8.add(epa >= 1 );
	//	epa.end();

	//	IloExpr  epa2(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
	//	epa2 += nbBay * (1 - vF[k + 1]);
	//	c8.add(epa2 >= 1);
	//	epa2.end();

	//}
	//model.add(c8);
	//c8.end();

	////13r
	//for (k = 0; k < nbCrane - 1; k++)
	//	{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*endCzF[k + 1][i] - i*endCzF[k][i];
	//	model.add(epa >= 1 );
	//	epa.end();
	//	}




	//////约束（5）//endCzF 取在最后一个xF处之后
	//////约束（6）// zF 小于最小的xF的bay
	//for (j = 0; j < nbTask; j++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);

	//	epa += yF[k][j];
	//	for (i = 0; i < nbLocation[j]; i++)
	//		epa -= CzF[k][i];

	//	model.add(epa <= 0);
	//	epa.end();

	//	IloExpr  epa2(env);
	//	epa2 += xF[k][j];
	//	for (i = 0; i < nbLocation[j]; i++)
	//		epa2 -= zF[k][i];
	//	model.add(epa2 <= 0);
	//	epa2.end();

	//	IloExpr  epa3(env);
	//	epa3 += xF[k][j];
	//	for (i = nbLocation[j] - 1; i < nbBay; i++)
	//		epa3 -= endCzF[k][i];
	//	model.add(epa3 <= 0);
	//	epa3.end();

	//	}


	////QC travel limits
	//IloRangeArray  v00(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//	{
	//		if (nbLocation[i] < 2 * k + 1)
	//			epa += xF[k][i] + yF[k][i];
	//		if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
	//			epa += xF[k][i] + yF[k][i];
	//	}
	//	v00.add(epa <= 0);
	//	epa.end();

	//	if (k < nbCrane)
	//		for (j = nbBay - 2 * (nbCrane - k - 1); j < nbBay; j++)
	//			model.add(endCzF[k][j] == 0);
	//	if (k > 0)
	//		for (j = 0; j < 2 * k; j++)
	//			model.add(zF[k][j]+CzF[k][j] == 0);

	//}
	//model.add(v00);
	//v00.end();



////////////////sub-problem

	//// vF 取值(2)
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//			epa += yF[k][i];
	//	epa -= vF[k];
	//	model.add(epa >= 0);
	//	epa.end();
	//}


	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//		epa += xF[k][i];
	//	model.add(epa >= 1);
	//	epa.end();
	//}

	//// vF 取值(3) vF与violation关系
	//for (i = 0; i < nbTask; i++)
	//	for (j = 0; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{ 
	//		for (k = 0; k < nbCrane; k++)
	//		{
	//			IloExpr epa(env);
	//			IloExpr epa2(env);
	//			
	//			epa += vF[k] +1;
	//			epa2 += vF[k] + 1;

	//			if (k < nbCrane - 1)
	//			{
	//				epa -= xF[k][i] + yF[k][i];
	//				for (int kk = k+1; kk < nbCrane; kk++)
	//					epa -= xF[kk][j] + yF[kk][j];
	//				model.add(epa >= 0);
	//			}
	//			
	//			if (k > 0)
	//			{
	//				epa2 -= xF[k][j] + yF[k][j];
	//				for (int kk = 0; kk < k; kk++)
	//					epa2 -= xF[kk][i]+yF[kk][i];
	//				model.add(epa2 >= 0);
	//			}
	//			epa.end();
	//			epa2.end();
	//		}
	//		}


	//precedence
	for (i = 0; i < nbTask; i++)
		for (j = 0; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
				for (k = 0; k < nbCrane; k++)
				{

					{
						IloExpr epa(env);
						epa += yF[k][i];
						for (int kk = k; kk < nbCrane; kk++)
							epa -= yF[kk][j];
						model.add(epa <= 0);
						epa.end();
					}

					{
						IloExpr epa2(env);
						epa2 += xF[k][j];
						for (int kk = k; kk < nbCrane; kk++)
							epa2 -= xF[kk][i];
						model.add(epa2 <= 0);
						epa2.end();
					}

					if (k < nbCrane - 1)
					{
						IloExpr epa3(env);
						epa3 += xF[k][i];
						for (int kk = k + 1; kk < nbCrane; kk++)
							epa3 += xF[kk][j];
						//epa3 += xF[kk][j] + yF[kk][j];
						model.add(epa3 <= 1);
						epa3.end();
					}

					if (k < nbCrane - 1)
					{
						/*for (int kk = k + 1; kk < nbCrane; kk++)
						{
						IloExpr epa3(env);
						IloExpr epa4(env);
						epa3 += xF[k][i] + xF[kk][j];
						epa4 += yF[k][j] + yF[kk][i];
						model.add(epa3 <= 1);
						model.add(epa4 <= 1);
						epa3.end();
						epa4.end();
						}*/
						IloExpr epa3(env);
						IloExpr epa4(env);
						epa3 += xF[k][i];
						epa4 += yF[k][j];
						for (int kk = k + 1; kk < nbCrane; kk++)
						{
							epa3 += xF[kk][j];
							epa4 += yF[kk][i];
						}
						model.add(epa3 <= 1);
						model.add(epa4 <= 1);
						epa3.end();
						epa4.end();

					}

					//if (k < nbCrane - 1)
					//{
					//	IloExpr epa4(env);
					//	epa4 += yF[k][j];

					//	//for (int kk = k + 1; kk < nbCrane; kk++)
					//	//	epa4 += yF[kk][i];
					//	//epa4 -= 1;
					//	

					//	for (int kk = k + 1; kk < nbCrane; kk++)
					//		epa4 -= xF[kk][i];
					//	for (int kk = 0; kk <= k; kk++)
					//		epa4 -= xF[kk][i] + yF[kk][i];
					//	model.add(epa4 <= 0);
					//	epa4.end();
					//}

					//for (int kk = 0; kk <= k; kk++)
					//{
					//	model.add(yF[k][j] + xF[kk][i]-vF[kk]<=1);
					//}

				}
			}



	//////////////////////////////////////
		//////waiting time ///////
	/////////////////////////////////////


	///////forward trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += t0kF[k] + nbreadyT[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
					epa += (nbQ[j] / nbs) * xF[k][j];
			for (j = 0; j <= i; j++)
			{

				epa += wF[k][j];
				epa -= QCmove_time * j * zF[k][j];
				//epa -= 1000 * CzF[k][j];
			}
			model.add(epa + QCmove_time * i - TF[k][i] == 0);
			epa.end();

		}

	/////约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 1 - safe_margin; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
			IloExpr  epa(env);

			epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];
			epa += 2000;
			for (j = 0; j <= i + 1 + safe_margin; j++)
				epa -= 2000 * zF[k + 1][j];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];

			c3.add(epa >= 0);
			epa.end();
		}
	model.add(c3);
	c3.end();

	//////retracing trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += t0kF[k] + nbreadyT[k];
			for (j = 0; j < nbTask; j++)
				epa += (nbQ[j] / nbs) * xF[k][j];
			for (j = 0; j < nbBay; j++)
			{
				epa += 2 * QCmove_time * j * endCzF[k][j] - QCmove_time * j * zF[k][j];
				epa += wF[k][j];
			}


			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 >= i)
					epa += (nbQ[j] / nbs) * yF[k][j];
			for (j = i; j < nbBay; j++)
			{
				epa += CwF[k][j];
			}

			model.add(epa - QCmove_time * i - CTF[k][i] == 0);
			epa.end();

		}

	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 1 - safe_margin; i++)
		{
			IloExpr  epa(env);
			epa += CTF[k + 1][i + 1 + safe_margin] - CTF[k][i];

			epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];
			for (j = 0; j <= i + 1 + safe_margin; j++)
				epa -= 2000 * CzF[k + 1][j];

			model.add(epa >= 0);
			epa.end();
		}

	// forward and retrace
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbBay; j++)
			epa += j * CzF[k + 1][j] - j * endCzF[k][j];
		epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
		model.add(epa >= 1 + safe_margin);
		epa.end();

	}



	IloCplex cplex(env);
	cplex.extract(model);

	bool h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found for bidQCSP" << endl;
		cplex.clearModel();
		cplex.clear();
		//cplex.end();
		model.end();

		return false;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	*ObjVal = cplex.getBestObjValue();
	//*ObjVal = cplex.getValue(CF);

	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.getValue(xF[k][i]) > 0.9)
				yF_current2[k][i] = 1;
			else
				yF_current2[k][i] = 0;

			if (cplex.getValue(yF[k][i]) > 0.9)
				yF2_current[k][i] = 1;
			else
				yF2_current[k][i] = 0;
		}


	cplex.clearModel();
	cplex.clear();
	//cplex.end();
	model.end();

	return true;

}
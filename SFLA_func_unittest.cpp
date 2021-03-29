#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include <google/protobuf/util/json_util.h>

#include "common/file_system/file_system.h"
#include "common/file_system/file_system_helper.h"
#include "common/utility/string_helper.h"
#include "common/trash/file_system_utility.h"
#include "common/utility/time_utility.h"
#include "common/config/json_config_helper.h"

#include <random>  
#include <math.h>       /* sin */

#include<iostream>  
#include<time.h>  
#include<stdlib.h>  
#include<cmath>  
#include<fstream>  
#include<iomanip>  
#include<vector>
#include<string>

using namespace ficus;
using namespace std;

DEFINE_string(logDir, "./SFLA_func_unittest_log", "");

DEFINE_string(config_path, "/mnt/WXRG0274_ssd/yyang4/ficus2/my/", "");
DEFINE_string(output, "", "");

DEFINE_int32(MCN, 20, "混合迭代次数");
DEFINE_int32(miter, 2, "族群内更新次数");

DEFINE_int32(M, 20, "族群数");

DEFINE_int32(I, 10, "每个族群中的个体数");
DEFINE_int32(P, 20, "个体总数"); // = M*I

DEFINE_int32(SD, 20, "解空间对应的维度"); // 比如其中一种组合方案为Xi, Xi为一个SD维的向量

DEFINE_double(step, 2, "蛙跳的步长, [-step, step]");

/**
 * 0:ABC(标准人工峰群算法); 1:EABC算法; 2:GABC算法;　
 * 3:ABC/best/1算法; 4:DE(差分进化算法); 5:本文算法
 */
DEFINE_int32(algorithm_type, 6, "算法类型:ABC, EABC, LSABC, DEABC, HSABC(本文算法), IEABC, SFLA");

DEFINE_double(lb, 0, "搜索空间下界"); // 解 Xi 是 SD 维的向量,向量每个元素的取值范围 [lb, up]

DEFINE_double(up, 1, "搜索空间上界");
DEFINE_int32(loop, 1, "循环次数");

DEFINE_int32(func_param_count, 1, "");
DEFINE_int32(func, 0, "目标函数");

const int PI = 3.14159265359;

// 算法枚举到名称的映射, 需要保存测试信息, for 画图
const static std::map<int, std::string> gAlg2Name = {
    {0, "ABC"}, {1, "EABC"}, {2, "LSABC"}, {3, "DEABC"}, {4, "HSABC"}, {5, "IEABC"}, {6, "SFLA"}
};

// 每个算法循环执行 loop 次,保存每次的结果,最终的评价指标为平均值,最优质,最差值,标准差
int gLoop = 0;
std::map<std::string, std::vector<double>> gAlg2LoopResut;

// 算法枚举到名称的映射, 需要保存测试信息, for 画图
const static std::map<int, std::string> gFunc2Name= {
    {1, "f1"}, {2, "f2"}, {3, "f3"}, {4, "f4"}, {5, "f5"}, {6, "f6"}, {7, "f7"}, {8, "f8"}
};

class SflaAlgorithm
{
public:
    // 种群的定义
    struct Individal
    {
        Individal()
        {
            solution.resize(FLAGS_SD, 0);
            rsolution.resize(FLAGS_SD, 0);
            viriable.resize(FLAGS_func_param_count, 0.0);
        }

        // 记录真实的最小值, 目标函数的值
        double trueFit;
        // 适应度值
        double fitness;
        // 解决方案, 记录的是每个子任务对应的原子服务下标
        std::vector<double> solution;
        // 反向解
        std::vector<double> rsolution;

        // 函数的自变量个数为 FLAGS_func_param_count, 用指定f(x)求解时需要
        std::vector<double> viriable;
    };

    int initializeData()
    {
        if (FLAGS_P != FLAGS_M * FLAGS_I){
            LOG(ERROR) << "initializeData FLAGS_P:" << FLAGS_P << " M:" << FLAGS_M
                << " I:" << FLAGS_I;
            return FICUS_INVALID_ARGUMENT_ERROR;
        }
        LOG(ERROR) << "initializeData FLAGS_P:" << FLAGS_P << " M:" << FLAGS_M << " MCN:" << FLAGS_MCN
            << " miter:" << FLAGS_miter << " I:" << FLAGS_I
            << " SD:" << FLAGS_SD << " step:" << FLAGS_step << " lb:" << FLAGS_lb << " up:"
            << FLAGS_up << " algorithm_type:" << FLAGS_algorithm_type;
        pw.resize(FLAGS_M);
        pb.resize(FLAGS_M);
        group.resize(FLAGS_P);
        memeplex.resize(FLAGS_M, std::vector<Individal>(FLAGS_I, Individal()));
        temp.resize(FLAGS_M);
        result.resize(FLAGS_MCN, 0.0);
        timeCost.resize(FLAGS_MCN, 0.0);
        return FICUS_SUCC;
    }

    int start()
    {
        srand((unsigned)time(NULL));

        auto funcName = gFunc2Name.at(FLAGS_func);
        JsonConfigHelper result;
        CHECK_RTN_LOGE_OF_FUNC(result.Load(FLAGS_config_path + "/" + FLAGS_output));
        const auto& algName = gAlg2Name.at(FLAGS_algorithm_type);
#if 1
        ofstream resultof;
        resultof.open(funcName + "_result.txt", ios::app);
#endif

        // 初始化种群,
        CHECK_RTN_LOGE_OF_FUNC(initilize());

        uint64_t totalCostTime = TimeUtility::GetTimeStamp();

        int iter = 0;
        while(iter < FLAGS_MCN) {
            partition();
            VLOG(1) << "start partition finish iter:" << iter << " MCN:" << FLAGS_MCN;

            memetic();
            VLOG(1) << "start memetic finish iter:" << iter << " MCN:" << FLAGS_MCN;

            renew();
            VLOG(1) << "start renew finish iter:" << iter << " MCN:" << FLAGS_MCN;

            getGlobalBest();

            if (gLoop == FLAGS_loop){
                result[algName]["time"].append(int(TimeUtility::GetTimeStamp() - totalCostTime));
                result[algName]["fitness"].append(globalBest.trueFit);
            }

            ++iter;
            if (iter == FLAGS_MCN) {
                // 保存每次的结果
                auto& loopResult = gAlg2LoopResut[algName];
                loopResult.push_back(globalBest.trueFit);
                double actualFx = calculationTruefit(globalBest);

                if (gLoop == FLAGS_loop){
                    result[algName]["min_fitness"] = globalBest.trueFit;
                    auto& solution = globalBest.solution;
                    for (int i = 0; i < solution.size(); ++i) {
                        result[algName]["solution"].append(solution[i]);
                    }
                    result[algName]["actual_min_fitness"] = actualFx;
                    // 升序列
                    std::sort(loopResult.begin(), loopResult.end());
                    double sum = std::accumulate(std::begin(loopResult), std::end(loopResult), 0.0);
                    double mean = sum / loopResult.size(); // 均值

                    double accum  = 0.0;
                    std::for_each(std::begin(loopResult), std::end(loopResult), [&](const double d) {
                        accum += (d-mean)*(d-mean);
                    });
                    double stdev = sqrt(accum/(loopResult.size())); // 方差

                    LOG(INFO) << "start result output:" << FLAGS_output 
                        << " alg:" << algName << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                        << " infer_fx:" << globalBest.trueFit 
                        << " actual_min_fx:" << actualFx << " mean:" << mean << " stdev:" << stdev << " better:" 
                        << loopResult.front() << " worst:" << loopResult.back()
                        << " SD:" << FLAGS_SD << " miter:" << FLAGS_miter
                        << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                        << " step:" << FLAGS_step << " M:" << FLAGS_M << " I:" << FLAGS_I << " P:" << FLAGS_P;
                } else {
                    LOG(INFO) << "start result output:" << FLAGS_output 
                        << " alg:" << algName << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                        << " infer_fx:" << globalBest.trueFit 
                        << " actual_min_fx:" << actualFx 
                        << " SD:" << FLAGS_SD << " miter:" << FLAGS_miter
                        << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                        << " step:" << FLAGS_step << " M:" << FLAGS_M << " I:" << FLAGS_I << " P:" << FLAGS_P;
                }

                resultof << setprecision(100) << globalBest.trueFit << " actual_fx:" << actualFx
                    << " output:" << FLAGS_output
                    << " algName:" << algName << " iter:" << iter << " mcn:"
                    << FLAGS_MCN << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                    << " SD:" << FLAGS_SD << " miter:" << FLAGS_miter
                    << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                    << " step:" << FLAGS_step << " M:" << FLAGS_M << " I:" << FLAGS_I << " P:" << FLAGS_P << std::endl; 
            }
        }
        resultof.close();
        CHECK_RTN_LOGE_OF_FUNC(result.Save(FLAGS_config_path + "/" + FLAGS_output));
        return FICUS_SUCC;
    }

private:
    int GetAllLineOfFile(const std::string &file, std::vector<std::string> &lines)
    {
        std::string oneLine;
        std::ifstream list(file.c_str());
        while (getline(list, oneLine))
        {
            lines.push_back(oneLine);
        }
        return FICUS_SUCC;
    }

    // 随机产生区间内的随机数
    double random(double start, double end)
    {
        double randValue = start+(end-start)*rand()/(RAND_MAX + 1.0);
        return randValue;
    }

    int randomSolution(std::vector<double> &solution, std::vector<double>& rsolution)
    {
        CHECK_ASSERT_RTN((solution.size() == FLAGS_SD));
        std::vector<double> chk(FLAGS_SD, 0.0);
        std::set<double> forbid{0.0, 1.0};
        chk[0] = random(0, 1);
        while(forbid.count(chk[0])){
            chk[0] = random(0, 1);
        }
        for (int i = 0; i < solution.size(); ++i) {
            solution[i] = random(FLAGS_lb, FLAGS_up);
            if (FLAGS_algorithm_type == 3){
                rsolution[i] = random(0, 1) * (FLAGS_up - FLAGS_lb) - solution[i];
            } else if (FLAGS_algorithm_type == 4 || FLAGS_algorithm_type == 5 || FLAGS_algorithm_type == 6) {
                if (i < FLAGS_SD - 1) {
                    chk[i+1] = sin(PI*chk[i]);
                }
                solution[i] = FLAGS_lb + chk[i]*(FLAGS_up - FLAGS_lb);
                rsolution[i] = FLAGS_lb + FLAGS_up - solution[i];
            }
            if (rsolution[i] < FLAGS_lb) {
                rsolution[i] = FLAGS_lb;
            }   
            if (rsolution[i] > FLAGS_up){
                rsolution[i] = FLAGS_up;
            }
            if (solution[i] < FLAGS_lb) {
                solution[i] = FLAGS_lb;
            }
            if (solution[i] > FLAGS_up){
                solution[i] = FLAGS_up;
            }
            VLOG(1) << "randomSolution number:" << solution[i] << " rnumber:" << rsolution[i];
        }
        return FICUS_SUCC;
    }

    // 初始化参数
    int initilize()
    {
        std::set<std::vector<double>> solotionSet;
        for (auto i = 0; i < FLAGS_P; ++i){
            // 随机生成一个解, 如果该解已经存在, 则重新生成
            auto& solution = group[i].solution;
            auto& rsolution = group[i].rsolution;
            CHECK_RTN_LOGE_OF_FUNC(randomSolution(solution, rsolution));
            while (solotionSet.count(solution)) {
                CHECK_RTN_LOGE_OF_FUNC(randomSolution(solution, rsolution));
            }
            solotionSet.insert(solution);
        }

        CHECK_RTN_LOGE_OF_FUNC(topNSolution());

        for (int idx = 0; idx < FLAGS_P; ++idx){
            group[idx].trueFit = calculationTruefit(group[idx]);  
            group[idx].fitness = calculationFitness(group[idx].trueFit);  
        }

        CHECK_RTN_LOGE_OF_FUNC(getGlobalBest());
        return FICUS_SUCC;
    }

    int topNSolution()
    {
        struct Data {
            std::vector<double> solution;
            double trueFit;
            double fitNess;
        };
        int idx = -1;
        std::vector<Data> dataList(FLAGS_P * 2);
        for (int i = 0; i < FLAGS_P; ++i) {
            ++idx;
            dataList[idx].solution = group[i].solution;
            dataList[idx].trueFit = calculationTruefit(group[i]);
            dataList[idx].fitNess = calculationFitness(dataList[idx].trueFit);

            ++idx;
            dataList[idx].solution = group[i].rsolution;
            group[i].solution = group[i].rsolution;
            dataList[idx].trueFit = calculationTruefit(group[i]);
            dataList[idx].fitNess = calculationFitness(dataList[idx].trueFit);
        }

        // 贪婪选择topN
        std::sort(dataList.begin(), dataList.end(),
            [this](const Data& l, const Data& r) -> bool 
            {
                return l.trueFit < r.trueFit;
            });

        for(int i = 0; i < FLAGS_P; ++i) {
            group[i].solution = dataList[i].solution;
        }
        return FICUS_SUCC;
    }

    // 分组
    void partition()
    {
        // 升序
        std::sort(group.begin(), group.end(),
            [this](const Individal& l, const Individal& r) -> bool 
            {
                return l.fitness < r.fitness; // 最糟糕的放在前面
            });

        /* 分组 */
        int k = 0;
        for (int i = 0; i < FLAGS_I; i++) {
            for (int j = 0; j < FLAGS_M; j++) {
                memeplex[j][i] = group[k];
                ++k;
            }
        }
        globalBest = group[FLAGS_P - 1];
        for (int i = 0; i < FLAGS_M; i++) {
            pw[i] = memeplex[i][0];
            pb[i] = memeplex[i][FLAGS_I - 1];
        }
    }

    void memetic()
    {
        int i, j, k, l, n;
        double a,b;
        std::vector<double> c(FLAGS_SD, 0.0);
        std::set<double> forbid{0.0, 0.25, 0.50, 0.75, 1.0};
        double miu = 4.0;
        for (n = 0; n<FLAGS_miter; n++) {
            for (i = 0; i<FLAGS_M; i++) {
                temp[i] = memeplex[i][0];
                for (j = 0; j<FLAGS_SD; j++) {  
                    temp[i].solution[j] = random(0, 1)*(pb[i].solution[j] - pw[i].solution[j]);
                    temp[i].solution[j] = temp[i].solution[j] > FLAGS_step ? FLAGS_step : temp[i].solution[j];
                    temp[i].solution[j] = temp[i].solution[j] < -FLAGS_step ? -FLAGS_step : temp[i].solution[j];
                    temp[i].solution[j] += pw[i].solution[j];

                    if (temp[i].solution[j] < FLAGS_lb) {
                        temp[i].solution[j] = FLAGS_lb;
                    }
                    if (temp[i].solution[j] > FLAGS_up){
                        temp[i].solution[j] = FLAGS_up;
                    }
                }
                temp[i].trueFit = calculationTruefit(temp[i]);
                temp[i].fitness = calculationFitness(temp[i].trueFit);
                if (temp[i].trueFit<pw[i].trueFit) { // 子群最好的, 目标函数最小化, 适应度值则是最大化
                    memeplex[i][0] = temp[i];
                } else {
                    for (k = 0; k<FLAGS_SD; k++) {
                        temp[i].solution[k] = random(0, 1)*(globalBest.solution[k] - pw[i].solution[k]);
                        temp[i].solution[k] = temp[i].solution[k] > FLAGS_step ? FLAGS_step : temp[i].solution[k];
                        temp[i].solution[k] = temp[i].solution[k] < -FLAGS_step ? -FLAGS_step : temp[i].solution[k];
                        temp[i].solution[k] += pw[i].solution[k];

                        if (temp[i].solution[k] < FLAGS_lb) {
                            temp[i].solution[k] = FLAGS_lb;
                        }
                        if (temp[i].solution[k] > FLAGS_up){
                            temp[i].solution[k] = FLAGS_up;
                        }
                    }
                    temp[i].trueFit = calculationTruefit(temp[i]);
                    temp[i].fitness = calculationFitness(temp[i].trueFit);
                    if (temp[i].trueFit<pw[i].trueFit) { // 整体最好的, 目标函数最小化, 适应度值则是最大化
                        memeplex[i][0] = temp[i];
                    } else { // 随机解
                        filterSolution.insert(pw[i].solution); // 被放弃的方案
                        auto& solution = memeplex[i][0].solution;

                        c[0] = random(0,1);
                        while(forbid.count(c[0])) {
                            c[0] = random(0,1);
                        }
                        for (int j = 0; j< FLAGS_SD; ++j) {
                            if (j != FLAGS_SD - 1){
                                c[j+1] = miu * c[j] * (1.0 - c[j]);
                            }
                            solution[j]=FLAGS_lb+c[j]*(FLAGS_up-FLAGS_lb);
                        }
                        // 重新生成解
                        while (filterSolution.count(solution)) {
                            c[0] = random(0,1);
                            while(forbid.count(c[0])) {
                                c[0] = random(0,1);
                            }
                            for (int j = 0; j< FLAGS_SD; ++j) {
                                if (j != FLAGS_SD - 1){
                                    c[j+1] = miu * c[j] * (1.0 - c[j]);
                                }
                                solution[j]=FLAGS_lb+c[j]*(FLAGS_up-FLAGS_lb);
                            }
                        }

                        memeplex[i][0].trueFit = calculationTruefit(memeplex[i][0]);
                        memeplex[i][0].fitness = calculationFitness(memeplex[i][0].trueFit);
                    }
                }
                // 升序
                std::sort(memeplex[i].begin(), memeplex[i].end(),
                    [this](const Individal& l, const Individal& r) -> bool 
                    {
                        return l.fitness < r.fitness; // 好的放后面
                    });
                pw[i] = memeplex[i][0];
                pb[i] = memeplex[i][FLAGS_I - 1];
                VLOG(1) << "memetic M:" << FLAGS_M << " i:" << i << " trueFit:" << memeplex[i][0].trueFit
                    << " fitness:" << memeplex[i][0].fitness;
            }
        }
    }

    void renew()
    {
        int i, j, k;
        i = 0;
        for (j = 0; j<FLAGS_M; j++)
        {
            for (k = 0; k<FLAGS_I; k++)
            {
                group[i] = memeplex[j][k];
                i++;
            }
        }
    }

    int getGlobalBest()
    {
        for (int i = 0; i < FLAGS_P; ++i) {
            if (group[i].trueFit < globalBest.trueFit) {
                globalBest = group[i];
            }
        }
        return FICUS_SUCC;
    }

    // 计算真实的函数值 f(x)
    double calculationTruefit(Individal& bee)
    {
        auto& solution = bee.solution;
        CHECK_ASSERT_RTN(solution.size() == FLAGS_SD);

        double trueFit = 0.0;
        double xi,xj,tmp1,tmp2,tmp;
        int n = solution.size();
        switch(FLAGS_func){
            case 1:
                for (int i = 0; i < solution.size(); ++i){
                    xi = solution[i];
                    trueFit += xi*xi;
                }
                break;
            case 2:
                for (int i = 0; i < solution.size(); ++i){
                    xi = solution[i];
                    trueFit += xi*xi*(i+1);
                }
                break;
            case 3:
                for (int i = 0; i < solution.size(); ++i){
                    xi = solution[i];
                    trueFit += (xi + 0.5)*(xi + 0.5);
                }
                break;
            case 4:
                for (int i = 0; i < solution.size() - 1; ++i){
                    xi = solution[i];
                    xj = solution[i+1];
                    trueFit += 100*(xj - xi*xi)*(xj - xi*xi) + (xi - 1.0)*(xi - 1.0);
                }
                break;
            case 5:
                for (int i = 0; i < solution.size(); ++i){
                    xi = solution[i];
                    tmp = xi*xi - 10*cos(2*PI*xi) + 10.0;
                    trueFit += tmp;
                }   
                break;
            case 6:
                tmp1 = 0.0;
                tmp2 = 1.0;
                for (int i = 0; i < solution.size(); ++i){
                    xi = solution[i];
                    tmp1 += xi*xi;
                    tmp2 *= cos(xi/(sqrt(i+1)));
                }
                trueFit = tmp1/4000 - tmp2 + 1.0;
                break;
            case 7:
                tmp1 = 0.0;
                tmp2 = 0.0;
                for (int i = 0; i < solution.size(); ++i){
                    xi = solution[i];
                    tmp1 += xi*xi;
                    tmp2 += cos(2*PI*xi);
                }
                tmp1 = sqrt(tmp1/n) * -0.2;
                tmp2 = tmp2 / n;
                trueFit = -20*exp(tmp1) - exp(tmp2) + 20.0 + exp(1);
                break;
            case 8:
                for (int i = 0; i < solution.size() - 1; ++i){
                    xi = solution[i];
                    tmp = xi*sin(xi) + 0.1*xi;
                    trueFit += fabs(tmp);
                }
                break;
            default:
                break;
        }
        return trueFit;
    }

    // 计算适应值 fitness(i)
    double calculationFitness(double truefit)
    {  
        double fitnessResult=0;  
        if (truefit>=0) {
            fitnessResult=1/(truefit+1);
        } else {
            fitnessResult=1+abs(truefit);
        }
        return fitnessResult;
    }  

    int lbUpCheck(Individal& individial, int param2change)
    {
        /*******判断是否越界*******/
        if (individial.solution[param2change]<FLAGS_lb) {
            individial.solution[param2change]=FLAGS_lb;
        }
        if (individial.solution[param2change]>FLAGS_up) {
            individial.solution[param2change]=FLAGS_up;
        }
        return FICUS_SUCC;
    }

    // https://bbs.csdn.net/topics/200067894 返回一个符合正态分布(均值为d, 标准差为d)的随机数
    double Normal_rand(double u, double d, int nnum_i)
    {//生成一个以u为均值,d为均方差的正态随机数,大数取nnum_i
    
        if(d<=0)return(u);
    
        int nnum=nnum_i;
        double sum_nnum=0.0;
        int i;
        //srand((unsigned)time(NULL));
        for(i=0;i<nnum;i++)
        {   
            sum_nnum+=(double)rand()/RAND_MAX;
        }
        ///////////////////////////
        //以下采取了中心极限定理,
        //由多个均匀分布产生服从正态的随机数,
        //均匀分布均值是0.5,方差是1/12,
        //n个均匀分布的随机变量和减(n*0.5)
        //再除以根号下n*(1/12)应当近似服从N(0,1),
        //在此n取nnum
        //////////////////////////
        return(u+d*(sum_nnum-nnum/2)/(double)sqrt(nnum/12));
    }
    // https://bbs.csdn.net/topics/200067894 返回一个符合正态分布(均值为d, 标准差为d)的随机数
    double Normal_rand(double u,double d)
    {//生成一个以u为均值,d为均方差的正态随机数x,采用反函数法   
        double u1,u2,z,x;

        if(d<=0)return(u);
    
        u1=(double)rand()/(double)RAND_MAX;
        u2=(double)rand()/(double)RAND_MAX;
    
        if(u1>0.0000000000000000)
        z=sqrt(-2*log(u1))*sin(2*PI*u2);
        else z=0;
        x=u+d*z;
    
        return(x);
    }

    std::vector<Individal> pw;  // 族群中个体最差位置
    std::vector<Individal> pb;  // 族群中个体最好位置
    Individal globalBest;       // 全局最优个体

    std::vector<Individal> group; // 全部个体
    std::vector<std::vector<Individal>> memeplex; // 排序后的群组

    std::vector<Individal> temp;

    std::set<std::vector<double>> filterSolution;
    std::vector<double> result;    // 截止到当前迭代轮数 对应的最优解
    std::vector<uint64_t> timeCost;  // 截止到当前迭代轮数 对应的总耗时
};


int calcAverageStandard()
{
    JsonConfigHelper result;
    CHECK_RTN_LOGE_OF_FUNC(result.Load(FLAGS_config_path + "/" + FLAGS_output));

    for (auto& it: gAlg2LoopResut) {
        auto& algName = it.first;
        auto& loopResult = it.second;
        CHECK_ASSERT_RTN(FLAGS_loop == loopResult.size());
        // 升序列
        std::sort(loopResult.begin(), loopResult.end());

        double sum = std::accumulate(std::begin(loopResult), std::end(loopResult), 0.0);
        double mean = sum / loopResult.size(); // 均值

        double accum  = 0.0;
        std::for_each(std::begin(loopResult), std::end(loopResult), [&](const double d) {
            accum += (d-mean)*(d-mean);
        });
        double stdev = sqrt(accum/(loopResult.size())); //方差

        for (const auto& fx: loopResult){
            result[algName]["loop"].append(fx);
        }
        result[algName]["worst"] = loopResult.back();
        result[algName]["better"] = loopResult.front();
        result[algName]["mean"] = mean;
        result[algName]["stdev"] = stdev;
    }
    CHECK_RTN_LOGE_OF_FUNC(result.Save(FLAGS_config_path + "/" + FLAGS_output));
    return FICUS_SUCC;
}

/**
 * @brief 行为识别举横幅端到端视频 recal 对齐
 */
TEST(SflaAlgorithm, sfla)
{
    int rtn = FICUS_SUCC;
    gLoop = 1;
    while (gLoop < FLAGS_loop + 1){
        SflaAlgorithm sflaAlg;
        rtn = sflaAlg.initializeData();
        ASSERT_EQ(rtn, FICUS_SUCC);
        rtn = sflaAlg.start();
        ASSERT_EQ(rtn, FICUS_SUCC);
        ++gLoop;
    }
    rtn = calcAverageStandard();
    ASSERT_EQ(rtn, FICUS_SUCC);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    google::ParseCommandLineFlags(&argc, &argv, false);

    INIT_LOG_2(argv[0], FLAGS_logDir);

    int rtn = RUN_ALL_TESTS();
    return rtn;
}


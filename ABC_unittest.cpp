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

DEFINE_string(logDir, "./ABC_unittest_log", "");

DEFINE_string(config_path, "/mnt/WXRG0274_ssd/yyang4/ficus2/my/", "");
DEFINE_string(test_data, "qws2.txt", "");
DEFINE_string(output, "", "");

DEFINE_int32(FoodNumber, 20, "食物的数量，为采蜜蜂的数量, 子任务数");

DEFINE_int32(D, 0, "每个子任务对应的原子服务数, D维空间");

DEFINE_int32(K, 4, "每个原子服务对应的Qos维度, 时间/可靠性/价格/可用性");

DEFINE_int32(SD, 20, "解空间对应的维度"); // 比如其中一种组合方案为Xi, Xi为一个SD维的向量

DEFINE_int32(limit, 20, "限度,超过这个限度没有更新采蜜蜂变成侦查蜂");
DEFINE_int32(loop, 1, "循环次数");

DEFINE_int32(MCN, 0, "最大迭代次数");

DEFINE_int32(A, 1, "EABC, LSABC需要用到的参数");
DEFINE_double(beta_average, 0.3, "EABC, LSABC需要用到的参数");
DEFINE_double(beta_standard, 0.3, "EABC, LSABC需要用到的参数");


/**
 * 0:ABC(标准人工峰群算法); 1:EABC算法; 2:GABC算法;　
 * 3:ABC/best/1算法; 4:DE(差分进化算法); 5:本文算法
 */
DEFINE_int32(algorithm_type, 0, "算法类型:ABC, EABC, LSABC, DEABC, HSABC(本文算法)");

DEFINE_int32(lb, 0, "搜索空间下界"); // = 0

DEFINE_int32(up, 1, "搜索空间上界"); // = FLAGS_D
DEFINE_int32(step, 1, "EABC, LSABC需要用到的参数");

DEFINE_int32(func_param_count, 1, "函数的参数个数");

const int PI = 3.14159265359;

// 算法枚举到名称的映射, 需要保存测试信息, for 画图
const static std::map<int, std::string> gAlg2Name = {
    {0, "ABC"}, {1, "EABC"}, {2, "LSABC"}, {3, "DEABC"}, {4, "HSABC"}, {5, "IEABC"}
};

// 每个算法循环执行 loop 次,保存每次的结果,最终的评价指标为平均值,最优质,最差值,标准差
int gLoop = 0;
std::map<std::string, std::vector<double>> gAlg2LoopResut;


class AbcAlgorithm
{
public:
    // 种群的定义
    struct BeeGroup
    {
        BeeGroup()
        {
            services.resize(FLAGS_D, AtomService());
            minMaxJ.resize(FLAGS_K, std::vector<double>(2, 0));
            solution.resize(FLAGS_SD, 0);
            rsolution.resize(FLAGS_SD, 0);
            viriable.resize(FLAGS_func_param_count, 0.0);
        }
        // 原子服务结构体定义
        struct AtomService
        {
            AtomService() { qos.resize(FLAGS_K, 0.0); }
            // QWS数据集没有该属性,分别是:时间,价格(随机生成,[1, 50]),可用性,可靠性,
            std::vector<double> qos;
        };

        // 每个原子服务对应的Qos
        std::vector<AtomService> services;
        // 当前子任务中,所有原子服务Qos各维度的最小值,最大值
        std::vector<std::vector<double>> minMaxJ;
        // 记录真实的最小值, 目标函数的值
        double trueFit;
        // 适应度值
        double fitness;
        // 跟随概率
        double rfitness;
        // 表示实验的次数,用于与limit作比较
        int trail;
        // 解决方案, 记录的是每个子任务对应的原子服务下标
        std::vector<int> solution;
        // 反向解
        std::vector<int> rsolution;
        // 超过limit没有被更新的方案将存入这里面作为禁忌表,侦查蜂重新生成方案需要不包含在这里面
        //std::set<std::vector<int>> filterSolution;
        // 禁忌表, set 存储被放弃的索引,比如set(j1,j2),表示ij1, ij2被放弃
        std::set<int> giveUpTable;
        // 最近一次邻域探索的 j,当超过limit次没被更新时,把 param2change 存入禁忌表 giveUpTable
        int param2change = -1;

        // 函数的自变量个数为 FLAGS_func_param_count, 用指定f(x)求解时需要
        std::vector<double> viriable;
    };  

    int initializeData()
    {
        NectarSource.resize(FLAGS_FoodNumber);
        EmployedBee.resize(FLAGS_FoodNumber);
        OnLooker.resize(FLAGS_FoodNumber);
        weight.resize(FLAGS_K, 0.0);
        result.resize(FLAGS_MCN, 0.0);
        timeCost.resize(FLAGS_MCN, 0.0);
        return FICUS_SUCC;
    }

    int start()
    {
        srand((unsigned)time(NULL));

        JsonConfigHelper result;
        CHECK_RTN_LOGE_OF_FUNC(result.Load(FLAGS_config_path + "/" + FLAGS_output));
        const auto& algName = gAlg2Name.at(FLAGS_algorithm_type);

#if 1
        std::string suffix = FilePathUtility::GetOnlyFileNameWithOutExt(FLAGS_test_data);
        ofstream resultof;
        resultof.open(suffix + "_result.txt");
#endif

        // 初始化种群,
        CHECK_RTN_LOGE_OF_FUNC(initilize());
        // 保存最好的蜜源
        getBestSource();

        uint64_t totalCostTime = TimeUtility::GetTimeStamp();

        while(iter < FLAGS_MCN) {
            // 雇佣峰
            sendEmployedBees();
            VLOG(1) << "start sendEmployedBees finish iter:" << iter << " MCN:" << FLAGS_MCN;

            // 计算更随概率
            CalculateProbabilities();
            VLOG(1) << "start CalculateProbabilities finish iter:" << iter << " MCN:" << FLAGS_MCN;

            // 观察峰
            sendOnlookerBees();
            VLOG(1) << "start sendOnlookerBees finish iter:" << iter << " MCN:" << FLAGS_MCN;

            // 记录最优蜜源
            getBestSource();
            VLOG(1) << "start getBestSource finish iter:" << iter << " MCN:" << FLAGS_MCN;

            // 判断是否存在 limit 次未得到更新,如果存在,侦查蜂出现
            sendScoutBees();
            VLOG(1) << "start sendScoutBees finish iter:" << iter << " MCN:" << FLAGS_MCN;

            // 记录最优蜜源
            getBestSource();
            VLOG(1) << "start getBestSource finish iter:" << iter << " MCN:" << FLAGS_MCN;
            if (gLoop == FLAGS_loop){
                result[algName]["time"].append(int(TimeUtility::GetTimeStamp() - totalCostTime));
                result[algName]["fitness"].append(BestSource.trueFit);
            }


            ++iter;
            if (iter == FLAGS_MCN) {
                // 保存每次的结果
                auto& loopResult = gAlg2LoopResut[algName];
                loopResult.push_back(BestSource.trueFit);

                if (gLoop == FLAGS_loop)
                    result[algName]["min_fiteness"] = BestSource.trueFit;
                auto& solution = BestSource.solution;
                std::vector<double> qosSum(FLAGS_K, 0.0);
                for (int i = 0; i < solution.size(); ++i) {
                    auto dIdx = solution[i];
                    if (gLoop == FLAGS_loop)
                        result[algName]["solution"].append(dIdx);
                    // 时间,价格(随机生成,[1, 50]),可用性,可靠性,
                    auto& qos = NectarSource[i].services[dIdx].qos;
                    if (gLoop == FLAGS_loop){
                        result[algName]["qos"]["time"].append(qos[0]);
                        result[algName]["qos"]["price"].append(qos[1]);
                        result[algName]["qos"]["avaible"].append(qos[2]);
                        result[algName]["qos"]["reliable"].append(qos[3]);
                    }
                    qosSum[0] += qos[0];
                    qosSum[1] += qos[1];
                    qosSum[2] *= qos[2];
                    qosSum[3] *= qos[3];
                }
                double actualFx = 0.0;
                for (int i = 0; i < weight.size(); ++i) {
                    if (gLoop == FLAGS_loop)
                        result[algName]["weight"].append(weight[i]);
                    actualFx += weight[i] * qosSum[i];
                }
                if (gLoop == FLAGS_loop)
                    result[algName]["actual_min_fitness"] = actualFx;
                if (gLoop == FLAGS_loop) {
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
                        << " infer_fx:" << BestSource.trueFit 
                        << " actual_min_fx:" << actualFx << " mean:" << mean << " stdev:" << stdev << " better:" 
                        << loopResult.front() << " worst:" << loopResult.back()
                        << " FoodNumber:" << FLAGS_FoodNumber
                        << " D:" << FLAGS_D << " K:" << FLAGS_K << " SD:" << FLAGS_SD << " limit:" << FLAGS_limit
                        << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                        << " step:" << FLAGS_step;
                } else {
                    LOG(INFO) << "start result output:" << FLAGS_output 
                        << " alg:" << algName << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                        << " infer_fx:" << BestSource.trueFit 
                        << " actual_min_fx:" << actualFx << " FoodNumber:" << FLAGS_FoodNumber
                        << " D:" << FLAGS_D << " K:" << FLAGS_K << " SD:" << FLAGS_SD << " limit:" << FLAGS_limit
                        << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                        << " step:" << FLAGS_step;
                }

                resultof << setprecision(100) << BestSource.trueFit << " actual_fx:" << actualFx
                    << " output:" << FLAGS_output
                    << " algName:" << algName << " iter:" << iter << " mcn:"
                    << FLAGS_MCN << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                    << " FoodNumber:" << FLAGS_FoodNumber << " SD:" << FLAGS_SD << " limit:" << FLAGS_limit
                    << " lb:" << FLAGS_lb << " up:" << FLAGS_up << " step:" << FLAGS_step << std::endl;   
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
    double random(double start, double end, std::set<int> filter = {})
    {
        double randValue = start+(end-start)*rand()/(RAND_MAX + 1.0);
        if (filter.empty()){
            return randValue;
        }
        // 如果已经存在, 则循环, 禁忌表功能才有这个逻辑
        while(filter.count(int(randValue))){
            randValue = start+(end-start)*rand()/(RAND_MAX + 1.0);
        }
        return randValue;
    }

    int randomSolution(std::vector<int> &solution, std::vector<int>& rsolution)
    {
        CHECK_ASSERT_RTN((solution.size() == FLAGS_SD));
        std::vector<double> chk(FLAGS_SD, 0.0);
        std::set<double> forbid{0.0, 1.0};
        chk[0] = random(0, 1);
        while(forbid.count(chk[0])){
            chk[0] = random(0, 1);
        }
        for (int i = 0; i < solution.size(); ++i) {
            solution[i] = int(random(FLAGS_lb, FLAGS_up));
            if (FLAGS_algorithm_type == 3){
                rsolution[i] = int(random(0, 1) * (FLAGS_up - FLAGS_lb) - solution[i]);
            } else if (FLAGS_algorithm_type == 4 || FLAGS_algorithm_type == 5) {
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
        std::vector<std::string> lines;
        CHECK_RTN_LOGE_OF_FUNC(GetAllLineOfFile(FLAGS_config_path + FLAGS_test_data, lines));

        // 最后一行是权重
        std::vector<std::string> tokens = StringSplit(lines.back(), ",");
        CHECK_ASSERT_RTN((tokens.size() == FLAGS_K));
        for (int j = 0; j < tokens.size(); ++j) {
            weight[j] = std::stod(tokens[j]);
        }

        std::set<std::vector<int>> solotionSet;
        int foodIdx = -1;
        int jdx = -1;
        for (int k = 0; k < lines.size() - 1; ++k) {
            std::vector<std::string> tokens = StringSplit(lines[k], ",");
            if (k / FLAGS_D >= FLAGS_FoodNumber) {
                break;
            }
            foodIdx = k / FLAGS_D;
            jdx = k % FLAGS_D;
            VLOG(1) << "initilize lines size:" << lines.size() << " D:" << FLAGS_D
                << " FLAGS_FoodNumber:" << FLAGS_FoodNumber << " K:" << FLAGS_K
                << " foodIdx:" << foodIdx << " jdx:" << jdx << " k:" << k
                << " NectarSource size:" << NectarSource.size() << " tokens size:"
                << tokens.size() << " line:" << lines[k];

            // 初始化 QoS 非功能属性值
            auto& qos = NectarSource[foodIdx].services[jdx].qos;
            qos[0] = std::stod(tokens[0]); // 时间
            qos[1] = random(1, 50); // 价格
            qos[2] = std::stod(tokens[1]); // 可用性
            qos[3] = std::stod(tokens[4]); // 可靠性

            // 随机生成一个解, 如果该解已经存在, 则重新生成
            auto& solution = NectarSource[foodIdx].solution;
            auto& rsolution = NectarSource[foodIdx].rsolution;
            CHECK_RTN_LOGE_OF_FUNC(randomSolution(solution, rsolution));
            while (solotionSet.count(solution)) {
                CHECK_RTN_LOGE_OF_FUNC(randomSolution(solution, rsolution));
            }
            solotionSet.insert(solution);
        }
        CHECK_ASSERT_RTN((foodIdx == FLAGS_FoodNumber - 1));
        CHECK_ASSERT_RTN((jdx == FLAGS_D - 1));

        CHECK_RTN_LOGE_OF_FUNC(qosPreProcess());

        CHECK_RTN_LOGE_OF_FUNC(topNSolution());

        for (int foodIdx = 0; foodIdx < FLAGS_FoodNumber; ++foodIdx){
            /****蜜源的初始化*****/
            NectarSource[foodIdx].trueFit = calculationTruefit(NectarSource[foodIdx]);  
            NectarSource[foodIdx].fitness = calculationFitness(NectarSource[foodIdx].trueFit);  
            NectarSource[foodIdx].rfitness = 0;
            NectarSource[foodIdx].trail = 0;
            /****采蜜蜂的初始化*****/
            EmployedBee[foodIdx].trueFit = NectarSource[foodIdx].trueFit;
            EmployedBee[foodIdx].fitness = NectarSource[foodIdx].fitness;
            EmployedBee[foodIdx].rfitness = NectarSource[foodIdx].rfitness;
            EmployedBee[foodIdx].trail = NectarSource[foodIdx].trail;
            /****观察蜂的初始化****/
            OnLooker[foodIdx].trueFit = NectarSource[foodIdx].trueFit;
            OnLooker[foodIdx].fitness = NectarSource[foodIdx].fitness;
            OnLooker[foodIdx].rfitness = NectarSource[foodIdx].rfitness;
            OnLooker[foodIdx].trail = NectarSource[foodIdx].trail;
        }
        /*****最优蜜源的初始化*****/
        BestSource.trueFit = NectarSource[0].trueFit;  
        BestSource.fitness = NectarSource[0].fitness;  
        BestSource.rfitness = NectarSource[0].rfitness;  
        BestSource.trail = NectarSource[0].trail;
        BestIndex = 0;
        return FICUS_SUCC;
    }

    int topNSolution()
    {
        if (FLAGS_algorithm_type != 3 && FLAGS_algorithm_type != 4 && FLAGS_algorithm_type != 5) {
            return FICUS_SUCC;
        }
        struct Data {
            std::vector<int> solution;
            double trueFit;
            double fitNess;
        };
        int idx = -1;
        std::vector<Data> dataList(FLAGS_FoodNumber * 2);
        for (int i = 0; i < FLAGS_FoodNumber; ++i) {
            ++idx;
            dataList[idx].solution = NectarSource[i].solution;
            dataList[idx].trueFit = calculationTruefit(NectarSource[i]);
            dataList[idx].fitNess = calculationFitness(dataList[idx].trueFit);

            ++idx;
            dataList[idx].solution = NectarSource[i].rsolution;
            NectarSource[i].solution = NectarSource[i].rsolution;
            dataList[idx].trueFit = calculationTruefit(NectarSource[i]);
            dataList[idx].fitNess = calculationFitness(dataList[idx].trueFit);
        }

        // 贪婪选择topN
        std::sort(dataList.begin(), dataList.end(),
            [this](const Data& l, const Data& r) -> bool 
            {
                return l.trueFit > r.trueFit;
            });

        for(int i = 0; i < FLAGS_FoodNumber; ++i) {
            NectarSource[i].solution = dataList[i].solution;
        }
        return FICUS_SUCC;
    }

    // 数据预处理, 归一化
    int qosPreProcess()
    {
        // 计算每个子任务内所有原子服务QoS的最小和最大值
        for (int i = 0; i < FLAGS_FoodNumber; ++i) {
            auto& minMaxJ = NectarSource[i].minMaxJ;
            for (int k = 0; k < FLAGS_K; ++k) {
                minMaxJ[k][0] = NectarSource[i].services[0].qos[k];
                minMaxJ[k][1] = NectarSource[i].services[0].qos[k];
            }
            for (int j = 0; j < FLAGS_D; ++j){
                auto& qos = NectarSource[i].services[j].qos;
                for (int k = 0; k < FLAGS_K; ++k) {
                    if (qos[k] > minMaxJ[k][1]) {
                        minMaxJ[k][1] = qos[k];
                    }
                    if (qos[k] < minMaxJ[k][0]) {
                        minMaxJ[k][0] = qos[k];
                    }
                    VLOG(1) << "qosPreProcess i:" << i << " j:" << j << " k:" << k
                        << " qos:" << qos[k] << " FLAGS_D:" << FLAGS_D << " FLAGS_K:"
                        << FLAGS_K << " qos size:" << qos.size();
                }
            }
        }

        // 归一化, 时间/价格:越小越好; 可靠性/可用性:越大越好
        for (int i = 0; i < FLAGS_FoodNumber; ++i) {
            auto& minMaxJ = NectarSource[i].minMaxJ;
            for (int j = 0; j < FLAGS_D; ++j) {
                auto& qos = NectarSource[i].services[j].qos;
                for (int k = 0; k < FLAGS_K; ++k) {
                    VLOG(1) << "qosPreProcess i:" << i << " j:" << j << " k:" << k
                    << " qos:" << qos[k];
                    auto& min = minMaxJ[k][0];
                    auto& max = minMaxJ[k][1];
                    auto diff = max - min;
                    if (k < 2) { // 时间/价格
                        qos[k] = (max - qos[k]) / diff;
                    } else {    // 可用性/可靠性
                        qos[k] = (qos[k] - min) / diff;
                    }
                    VLOG(1) << "qosPreProcess i:" << i << " j:" << j << " k:" << k
                        << " qos:" << qos[k] << " min:" << min << " max:" << max;
                }
            }
        }
        return FICUS_SUCC;
    }

    // 计算真实的函数值 f(x)
    double calculationTruefit(BeeGroup& bee)
    {
        auto& solution = bee.solution;
        CHECK_ASSERT_RTN(solution.size() == FLAGS_SD);

        // 时间/价格累加; 可用性/可靠性累乘;
        std::vector<double> qosSum = NectarSource[0].services[solution.front()].qos;
        for (int i = 1; i < solution.size(); ++i) {
            // 当前方案 solution 下子任务 i 的原子服务下标为 dIdx
            auto& dIdx = solution[i];
            auto& qos = NectarSource[i].services[dIdx].qos;
            for (int k = 0; k < qos.size(); ++k){
                VLOG(1) << "i:" << i << " k:" << k << " qosSum size:"
                    << qosSum.size() << " qos size:" << qos.size() << " dIdx:" << dIdx;
                if (k < 2) {
                    qosSum[k] += qos[k];
                } else {
                    qosSum[k] *= qos[k];
                }
            }
        }

        // QoS 各维度权重计算
        double trueFit = 0.0;
        for (int k = 0; k < FLAGS_K; ++k) {
            trueFit += weight[k] * qosSum[k];
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

    // 修改采蜜蜂的函数
    int sendEmployedBees()
    {
        int i,j,k,r1,r2;
        int best = BestIndex;  
        int param2change;   // 需要改变的维数
        double Rij; // [-1,1]之间的随机数
        for (i = 0;i < FLAGS_FoodNumber; ++i) {
            param2change=(int)random(0,FLAGS_SD);  // 随机选取需要改变的维数

            for (j=0;j<FLAGS_SD;j++) { // 解空间Xi是一个SD维向量
                EmployedBee[i].solution[j]=NectarSource[i].solution[j];
            }

            /******* 采蜜蜂去更新信息 算法类型:ABC, EABC, LSABC, DEABC, HSABC(本文算法) *******/
            std::set<int> filter;
            double alpha, beta;
            double w, c1, c2;
            double fx1, fx2, fx3;
            switch(FLAGS_algorithm_type){
                case 0:
                    while (1){
                        k=(int)random(0,FLAGS_FoodNumber);
                        if (k!=i) {
                            break;
                        }
                    }
                    Rij=random(-1,1);
                    EmployedBee[i].solution[param2change]=NectarSource[i].solution[param2change]
                        +Rij*(NectarSource[i].solution[param2change]-NectarSource[k].solution[param2change]);
                    lbUpCheck(EmployedBee[i], param2change);
                    EmployedBee[i].trueFit=calculationTruefit(EmployedBee[i]);
                    break;
                case 1:
                case 5:
                    filter.insert(i);
                    filter.insert(best);
                    while (1){
                        r1=(int)random(0,FLAGS_FoodNumber);
                        r2=(int)random(0,FLAGS_FoodNumber);
                        if (!filter.count(r1) && !filter.count(r2)) {
                            break;
                        }
                    }
                    alpha = random(0,FLAGS_A);
                    beta = Normal_rand(FLAGS_beta_average, FLAGS_beta_standard) * random(0, 1);
                    EmployedBee[i].solution[param2change]=NectarSource[r1].solution[param2change]
                        +alpha*(NectarSource[best].solution[param2change]-NectarSource[r1].solution[param2change])
                        +beta*(NectarSource[r1].solution[param2change]-NectarSource[r2].solution[param2change]);
                    lbUpCheck(EmployedBee[i], param2change);
                    EmployedBee[i].trueFit=calculationTruefit(EmployedBee[i]);
                    break;
                case 2:
                    filter.insert(i);
                    //filter.insert(best);
                    while (1){
                        k=(int)random(0,FLAGS_FoodNumber);
                        if (k != i){
                            r1=(int)random(0,FLAGS_FoodNumber);
                            r2=(int)random(0,FLAGS_FoodNumber);
                            if (r1 != r2 && r2 != i) {
                                break;
                            }
                        }
                    }
                    // 搜索公式1
                    for (j=0;j<FLAGS_SD;j++) { // 解空间Xi是一个SD维向量
                        tempSource.solution[j]=NectarSource[i].solution[j];
                    }
                    Rij=random(-1,1);
                    c1 = c2 = 2.0;
                    w = 1.0 - (iter * 1.0 / FLAGS_MCN) * (1 - 0.5);
                    tempSource.solution[param2change]=w*NectarSource[i].solution[param2change]
                        +c1*Rij*(NectarSource[best].solution[param2change]-NectarSource[i].solution[param2change])
                        +c2*Rij*(NectarSource[k].solution[param2change]-NectarSource[i].solution[param2change]);
                    lbUpCheck(tempSource, param2change);
                    fx1 = calculationTruefit(tempSource);

                    // 搜索公式2
                    alpha = random(0,FLAGS_A);
                    beta = Normal_rand(FLAGS_beta_average, FLAGS_beta_standard) * random(0, 1);
                    EmployedBee[i].solution[param2change]=NectarSource[i].solution[param2change]
                        +alpha*(NectarSource[best].solution[param2change]-NectarSource[r1].solution[param2change])
                        +beta*(NectarSource[r1].solution[param2change]-NectarSource[r2].solution[param2change]);
                    lbUpCheck(EmployedBee[i], param2change);
                    EmployedBee[i].trueFit = calculationTruefit(EmployedBee[i]);
                    if (fx1 > EmployedBee[i].trueFit) {
                        EmployedBee[i].trueFit = fx1;
                        EmployedBee[i].solution[param2change] = tempSource.solution[param2change];
                    }
                    break;
                case 3:
                    while (1){
                        k=(int)random(0,FLAGS_FoodNumber);
                        if (k!=i) {
                            break;
                        }
                    }
                    Rij=random(-1,1);
                    EmployedBee[i].solution[param2change]=(iter*1.0/FLAGS_MCN)*NectarSource[best].solution[param2change]
                        +(1.0-(iter*1.0)/FLAGS_MCN)*Rij*(NectarSource[i].solution[param2change]-NectarSource[k].solution[param2change]);
                    lbUpCheck(EmployedBee[i], param2change);
                    EmployedBee[i].trueFit=calculationTruefit(EmployedBee[i]);
                    break;
                case 4:
                    while (1){
                        k=(int)random(0,FLAGS_FoodNumber);
                        if (k != i){
                            break;
                        }
                    }
                    // 搜索公式1
                    for (j=0;j<FLAGS_SD;j++) { // 解空间Xi是一个SD维向量
                        tempSource.solution[j]=NectarSource[i].solution[j];
                    }
                    Rij=random(-1,1);
                    tempSource.solution[param2change]=(iter*1.0/FLAGS_MCN)*NectarSource[best].solution[param2change]
                        +(1.0-(iter*1.0)/FLAGS_MCN)*Rij*(NectarSource[i].solution[param2change]-NectarSource[k].solution[param2change]);
                    lbUpCheck(tempSource, param2change);
                    fx1 = calculationTruefit(tempSource);

                    // 搜索公式2
                    EmployedBee[i].solution[param2change]=NectarSource[i].solution[param2change]
                        +(random(0, 1)-0.5)*2.0*(NectarSource[i].solution[param2change]-NectarSource[k].solution[param2change])
                        +random(0, 1)*2.0*(NectarSource[best].solution[param2change]-NectarSource[i].solution[param2change]);
                    lbUpCheck(EmployedBee[i], param2change);
                    EmployedBee[i].trueFit = calculationTruefit(EmployedBee[i]);
                    if (fx1 > EmployedBee[i].trueFit) {
                        EmployedBee[i].trueFit = fx1;
                        EmployedBee[i].solution[param2change] = tempSource.solution[param2change];
                    }
                    break;
                default:
                    break;
            }
            EmployedBee[i].fitness = calculationFitness(EmployedBee[i].trueFit);

            /******贪婪选择策略*******/
            if (EmployedBee[i].trueFit>NectarSource[i].trueFit) {
                for (j=0;j<FLAGS_SD;j++) {
                    NectarSource[i].solution[j]=EmployedBee[i].solution[j];
                }
                NectarSource[i].trail=0;
                NectarSource[i].trueFit=EmployedBee[i].trueFit;
                NectarSource[i].fitness=EmployedBee[i].fitness;
                NectarSource[i].param2change = param2change;
            } else {
                NectarSource[i].trail++;
            }
        }
        return FICUS_SUCC;
    }

    // 计算跟随概率
    int CalculateProbabilities()
    {
        double sumFit = 0;
        for (int i=0;i < FLAGS_FoodNumber; ++i) {
            sumFit += NectarSource[i].trueFit; // 目标函数最大化, 
        }

        for (int i = 0;i < FLAGS_FoodNumber; ++i) {
            NectarSource[i].rfitness=NectarSource[i].trueFit/sumFit;
        }
        return FICUS_SUCC;
    }

    int lbUpCheck(BeeGroup& bee, int param2change)
    {
        /*******判断是否越界*******/
        if (bee.solution[param2change]<FLAGS_lb) {
            bee.solution[param2change]=FLAGS_lb;
        }
        if (bee.solution[param2change]>FLAGS_up) {
            bee.solution[param2change]=FLAGS_up;
        }
        return FICUS_SUCC;
    }

    // 采蜜蜂与观察蜂交流信息，观察蜂 根据跟随概率 更改信息
    int sendOnlookerBees()
    {  
        int i,j,t,k,r1,r2;
        int best = BestIndex;
        double R_choosed;   //被选中的概率
        int param2change;   //需要被改变的维数
        double Rij; //[-1,1]之间的随机数
        i=0; 
        t=0;
        while(t<FLAGS_FoodNumber) {
            R_choosed=random(0,1);
            VLOG(1) << "sendOnlookerBees 00 i:" << i << " t:" << t << " R_choosed:" << R_choosed
                << " rfitness:" << NectarSource[i].rfitness;
            // 根据被选择的概率选择
            if(R_choosed<NectarSource[i].rfitness) {
                ++t;
                param2change=(int)random(0,FLAGS_SD);
                VLOG(1) << "sendOnlookerBees 11 i:" << i << " t:" << t << " R_choosed:" << R_choosed
                    << " rfitness:" << NectarSource[i].rfitness;
                for(j=0;j<FLAGS_SD;j++) {
                    OnLooker[i].solution[j]=NectarSource[i].solution[j];
                }

                /******* 采蜜蜂去更新信息 算法类型:ABC, EABC, LSABC, DEABC, HSABC(本文算法) *******/
                std::set<int> filter;
                double alpha, beta;
                double w, c1, c2;
                double fx1, fx2, fx3;
                double cr = 0.5;
                switch(FLAGS_algorithm_type){
                    case 0:
                    case 2:
                    case 3:
                        while (1) {
                            k=(int)random(0,FLAGS_FoodNumber);
                            if (k!=i) {
                                break;
                            }
                        }
                        Rij=random(-1,1);
                        OnLooker[i].solution[param2change]=NectarSource[i].solution[param2change]
                            +Rij*(NectarSource[i].solution[param2change]-NectarSource[k].solution[param2change]);
                        lbUpCheck(OnLooker[i], param2change);
                        OnLooker[i].trueFit=calculationTruefit(OnLooker[i]);
                        break;
                    case 1:
                    case 5:
                        filter.insert(i);
                        filter.insert(best);
                        while (1){
                            r1=(int)random(0,FLAGS_FoodNumber);
                            r2=(int)random(0,FLAGS_FoodNumber);
                            if (!filter.count(r1) && !filter.count(r2)) {
                                break;
                            }
                        }
                        alpha = random(0,FLAGS_A);
                        beta = Normal_rand(FLAGS_beta_average, FLAGS_beta_standard) * random(0, 1);
                        OnLooker[i].solution[param2change]=NectarSource[r1].solution[param2change]
                            +alpha*(NectarSource[best].solution[param2change]-NectarSource[r1].solution[param2change])
                            +beta*(NectarSource[r1].solution[param2change]-NectarSource[best].solution[param2change]);
                        lbUpCheck(OnLooker[i], param2change);
                        OnLooker[i].trueFit=calculationTruefit(OnLooker[i]);
                        break;
                    case 4:
                        while (1){
                            k=(int)random(0,FLAGS_FoodNumber);
                            r1=(int)random(0,FLAGS_FoodNumber);
                            r2=(int)random(0,FLAGS_FoodNumber);
                            if (k != i){
                                if (r1 != r2 && r1 != i){
                                    break;
                                }
                            }
                        }
                        // 搜索公式1
                        if (random(0, 1) < cr){
                            OnLooker[i].solution[param2change]=NectarSource[i].solution[param2change]
                                +(random(0,1)-0.5)*2.0*(NectarSource[i].solution[param2change]-NectarSource[k].solution[param2change]);
                            lbUpCheck(OnLooker[i], param2change);
                            OnLooker[i].trueFit = calculationTruefit(OnLooker[i]);
                        } else {
                            // 搜索公式2
                            OnLooker[i].solution[param2change]=NectarSource[best].solution[param2change]
                                +(random(0, 1)-0.5)*2.0*(NectarSource[r1].solution[param2change]-NectarSource[r2].solution[param2change]);
                            lbUpCheck(OnLooker[i], param2change);
                            OnLooker[i].trueFit = calculationTruefit(OnLooker[i]);
                        }
                        break;
                    default:
                        break;
                }

                VLOG(1) << "sendOnlookerBees 22 i:" << i << " t:" << t << " R_choosed:" << R_choosed
                    << " rfitness:" << NectarSource[i].rfitness;

                OnLooker[i].fitness=calculationFitness(OnLooker[i].trueFit);

                /****贪婪选择策略******/
                if (OnLooker[i].trueFit>NectarSource[i].trueFit) {
                    for (j=0;j<FLAGS_SD;j++) {
                        NectarSource[i].solution[j]=OnLooker[i].solution[j];
                    }
                    NectarSource[i].trail=0;  
                    NectarSource[i].trueFit=OnLooker[i].trueFit;
                    NectarSource[i].fitness=OnLooker[i].fitness;
                    NectarSource[i].param2change = param2change;
                } else {
                    ++NectarSource[i].trail;
                }

            }
            VLOG(1) << "sendOnlookerBees 33 i:" << i << " t:" << t << " R_choosed:" << R_choosed
                << " rfitness:" << NectarSource[i].rfitness;
            ++i;
            if (i==FLAGS_FoodNumber) {
                i = 0;
            }
        }
        return FICUS_SUCC;
    }

    // 判断是否有侦查蜂的出现，有则重新生成蜜源
    int sendScoutBees()    
    {
        double R;   // [0,1]之间的随机数
        for (int i = 0; i < FLAGS_FoodNumber; ++i) {
            if(NectarSource[i].trail < FLAGS_limit) {
                continue;
            }

            auto& solution = NectarSource[i].solution;
            std::vector<double> c(FLAGS_SD, 0.0);
            std::set<double> forbid{0.0, 0.25, 0.50, 0.75, 1.0};
            double miu = 4.0;
            switch(FLAGS_algorithm_type){
                case 0:
                case 1:
                    for (int j = 0; j< FLAGS_SD; ++j) {
                        R=random(0,1);
                        NectarSource[i].solution[j]=FLAGS_lb+R*(FLAGS_up-FLAGS_lb);
                    }
                    NectarSource[i].trail=0;
                    NectarSource[i].trueFit=calculationTruefit(NectarSource[i]);
                    break;
                case 2:
                    filterSolution.insert(solution); // 将局部最优解存入禁忌表
                    for (int j = 0; j< FLAGS_SD; ++j) {
                        R=random(0,1);
                        solution[j]=FLAGS_lb+R*(FLAGS_up-FLAGS_lb);
                    }
                    // 重新生成解
                    while (filterSolution.count(solution)) {
                        for (int j = 0; j< FLAGS_SD; ++j) {
                            R=random(0,1);
                            solution[j]=FLAGS_lb+R*(FLAGS_up-FLAGS_lb);
                        }
                    }
                    NectarSource[i].trail=0;
                    NectarSource[i].trueFit=calculationTruefit(NectarSource[i]);
                    break;
                case 3:
                    for (int j = 0; j< FLAGS_SD; ++j) {
                        if (j != FLAGS_SD - 1){
                            c[j+1] = miu * c[j] * (1.0 - c[j]);
                        }
                        solution[j]=FLAGS_lb+c[j]*(FLAGS_up-FLAGS_lb);
                    }
                    NectarSource[i].trail=0;
                    NectarSource[i].trueFit=calculationTruefit(NectarSource[i]);
                    break;
                case 4:
                    filterSolution.insert(solution); // 将局部最优解存入禁忌表
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
                    #if 0
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
                    #endif
                    NectarSource[i].trail=0;
                    NectarSource[i].trueFit=calculationTruefit(NectarSource[i]);
                    break;
                case 5:
                    filterSolution.insert(solution); // 将局部最优解存入禁忌表
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
                    #if 1
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
                    #endif
                    NectarSource[i].trail=0;
                    NectarSource[i].trueFit=calculationTruefit(NectarSource[i]);
                    break;
                default:
                    break;
            }

            NectarSource[i].fitness=calculationFitness(NectarSource[i].trueFit);
        }

        if (FLAGS_algorithm_type == 4) {
            for (int i = 0; i < 10; ++i){
                sfla();
            }
        }
        return FICUS_SUCC;
    }

    int sfla()
    {
        getBestSource();
        struct Data{
            double trueFit;
            int idx;
        };
        std::vector<Data> vec(FLAGS_FoodNumber, Data());
        for (int i = 0; i < NectarSource.size(); ++i) {
            auto& data = vec[i];
            data.trueFit = NectarSource[i].trueFit;
            data.idx = i;
        }
        std::sort(vec.begin(), vec.end(), [](const Data& l, const Data& r){
                return l.trueFit > r.trueFit; // 目标函数最大化, 差的放后面, 好的放前面
            });

        std::vector<double> c(FLAGS_SD, 0.0);
        std::set<double> forbid{0.0, 0.25, 0.50, 0.75, 1.0};
        double miu = 4.0;

        int cnt = FLAGS_FoodNumber / 10;
        for (int k = 0; k < cnt; ++k) {
            int start = k;
            int end = FLAGS_FoodNumber - 1 - start;
            auto idx = vec[start].idx;  // better
            auto jdx = vec[end].idx;    // worst
            auto& better = NectarSource[idx];
            auto& worst = NectarSource[jdx];

            for (int j = 0; j<FLAGS_SD; j++) {  
                tempSource.solution[j] = random(0, 1)*(better.solution[j] - worst.solution[j]);
                tempSource.solution[j] = tempSource.solution[j] > FLAGS_step ? FLAGS_step : tempSource.solution[j];
                tempSource.solution[j] = tempSource.solution[j] < -FLAGS_step ? -FLAGS_step : tempSource.solution[j];
                tempSource.solution[j] += worst.solution[j];

                if (tempSource.solution[j] < FLAGS_lb) {
                    tempSource.solution[j] = FLAGS_lb;
                }
                if (tempSource.solution[j] > FLAGS_up){
                    tempSource.solution[j] = FLAGS_up;
                }
            }
            tempSource.trueFit = calculationTruefit(tempSource);
            tempSource.fitness = calculationFitness(tempSource.trueFit);
            if (tempSource.trueFit>worst.trueFit) { // 子群最好的
                worst.solution = tempSource.solution;
                worst.trueFit = tempSource.trueFit;
                worst.fitness = tempSource.fitness;
            } else {
                for (int j = 0; j<FLAGS_SD; j++) {  
                    tempSource.solution[j] = random(0, 1)*(BestSource.solution[j] - worst.solution[j]);
                    tempSource.solution[j] = tempSource.solution[j] > FLAGS_step ? FLAGS_step : tempSource.solution[j];
                    tempSource.solution[j] = tempSource.solution[j] < -FLAGS_step ? -FLAGS_step : tempSource.solution[j];
                    tempSource.solution[j] += worst.solution[j];

                    if (tempSource.solution[k] < FLAGS_lb) {
                        tempSource.solution[k] = FLAGS_lb;
                    }
                    if (tempSource.solution[k] > FLAGS_up){
                        tempSource.solution[k] = FLAGS_up;
                    }
                }
                tempSource.trueFit = calculationTruefit(tempSource);
                tempSource.fitness = calculationFitness(tempSource.trueFit);
                if (tempSource.trueFit>worst.trueFit) { // 子群最好的
                    worst.solution = tempSource.solution;
                    worst.trueFit = tempSource.trueFit;
                    worst.fitness = tempSource.fitness;
                } else {// 随机解
                    filterSolution.insert(worst.solution); // 被放弃的方案
                    auto& solution = worst.solution;

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
                    worst.trueFit = calculationTruefit(worst);
                    worst.fitness = calculationFitness(worst.trueFit);
                }
            }
        }
        return FICUS_SUCC;
    }

    // 保存最优的蜜源
    int getBestSource()
    {
        for (int i = 0; i < FLAGS_FoodNumber; ++i) {
            if (NectarSource[i].trueFit > BestSource.trueFit) {
                for (int j = 0; j < FLAGS_SD; ++j) {
                    BestSource.solution[j] = NectarSource[i].solution[j];
                }
                BestSource.trueFit = NectarSource[i].trueFit;
                BestIndex = i;
            }
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

    int iter = 0;
    std::vector<BeeGroup> NectarSource;     //蜜源，注意：一切的修改都是针对蜜源而言的  
    std::vector<BeeGroup> EmployedBee;      //采蜜蜂
    std::vector<BeeGroup> OnLooker;         //观察蜂
    std::vector<float> weight;              //Qos各维度权重 
    BeeGroup BestSource;                    //记录最好蜜源
    int BestIndex;
    BeeGroup tempSource;
    std::set<std::vector<int>> filterSolution;
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
TEST(AbcAlgorithmTest, ABC)
{
    int rtn = FICUS_SUCC;
    gLoop = 1;
    while (gLoop < FLAGS_loop + 1){
        AbcAlgorithm abc;
        rtn = abc.initializeData();
        ASSERT_EQ(rtn, FICUS_SUCC);
        rtn = abc.start();
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

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

DEFINE_string(logDir, "./SFLA_unittest_log", "");

DEFINE_string(config_path, "/mnt/WXRG0274_ssd/yyang4/ficus2/my/", "");
DEFINE_string(test_data, "qws2.txt", "");
DEFINE_string(output, "", "");

DEFINE_int32(MCN, 20, "混合迭代次数");
DEFINE_int32(miter, 2, "族群内更新次数");

DEFINE_int32(M, 20, "族群数");

DEFINE_int32(I, 20, "每个族群中的个体数");
DEFINE_int32(P, 20, "个体总数"); // = M*I

DEFINE_int32(D, 0, "每个子任务对应的原子服务数, D维空间");
DEFINE_int32(K, 4, "每个原子服务对应的Qos维度, 时间/可靠性/价格/可用性");
DEFINE_int32(SD, 20, "解空间对应的维度"); // 比如其中一种组合方案为Xi, Xi为一个SD维的向量

DEFINE_int32(step, 2, "蛙跳的步长, [-step, step]");

/**
 * 0:ABC(标准人工峰群算法); 1:EABC算法; 2:GABC算法;　
 * 3:ABC/best/1算法; 4:DE(差分进化算法); 5:本文算法
 */
DEFINE_int32(algorithm_type, 6, "算法类型:ABC, EABC, LSABC, DEABC, HSABC(本文算法), IEABC, SFLA");

DEFINE_int32(lb, 0, "搜索空间下界"); // 解 Xi 是 SD 维的向量,向量每个元素的取值范围 [lb, up]

DEFINE_int32(up, 1, "搜索空间上界");
DEFINE_int32(loop, 1, "循环次数");

DEFINE_int32(func_param_count, 1, "");

const int PI = 3.14159265359;

// 算法枚举到名称的映射, 需要保存测试信息, for 画图
const static std::map<int, std::string> gAlg2Name = {
    {0, "ABC"}, {1, "EABC"}, {2, "LSABC"}, {3, "DEABC"}, {4, "HSABC"}, {5, "IEABC"}, {6, "SFLA"}
};

// 每个算法循环执行 loop 次,保存每次的结果,最终的评价指标为平均值,最优质,最差值,标准差
int gLoop = 0;
std::map<std::string, std::vector<double>> gAlg2LoopResut;

class SflaAlgorithm
{
public:
    // 种群的定义
    struct Individal
    {
        Individal()
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
        // 解决方案, 记录的是每个子任务对应的原子服务下标
        std::vector<int> solution;
        // 反向解
        std::vector<int> rsolution;

        int taskIdx;
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
            << " miter:" << FLAGS_miter << " I:" << FLAGS_I << " D:" << FLAGS_D << " K:" << FLAGS_K
            << " SD:" << FLAGS_SD << " step:" << FLAGS_step << " lb:" << FLAGS_lb << " up:"
            << FLAGS_up << " algorithm_type:" << FLAGS_algorithm_type;
        pw.resize(FLAGS_M);
        pb.resize(FLAGS_M);
        group.resize(FLAGS_P);
        memeplex.resize(FLAGS_M, std::vector<Individal>(FLAGS_I, Individal()));
        weight.resize(FLAGS_K, 0.0);
        temp.resize(FLAGS_M);
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

                if (gLoop == FLAGS_loop)
                    result[algName]["min_fitness"] = globalBest.trueFit;
                auto& solution = globalBest.solution;

                std::vector<double> qosSum{0, 0, 1.0, 1.0};
                for (int i = 0; i < solution.size(); ++i) {
                    auto dIdx = solution[i];
                    if (gLoop == FLAGS_loop)
                        result[algName]["solution"].append(dIdx);
                    // 时间,价格(随机生成,[1, 50]),可用性,可靠性,
                    auto& qos = group[taskId2GroupIdx.at(i)].services[dIdx].qos;
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
                for (int i = 0; i < weight.size(); ++i){
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
                        << " infer_fx:" << globalBest.trueFit 
                        << " actual_min_fx:" << actualFx << " mean:" << mean << " stdev:" << stdev << " better:" 
                        << loopResult.front() << " worst:" << loopResult.back()
                        << " D:" << FLAGS_D << " K:" << FLAGS_K << " SD:" << FLAGS_SD << " miter:" << FLAGS_miter
                        << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                        << " step:" << FLAGS_step << " M:" << FLAGS_M << " I:" << FLAGS_I << " P:" << FLAGS_P;
                } else {
                    LOG(INFO) << "start result output:" << FLAGS_output 
                        << " alg:" << algName << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                        << " infer_fx:" << globalBest.trueFit 
                        << " actual_min_fx:" << actualFx 
                        << " D:" << FLAGS_D << " K:" << FLAGS_K << " SD:" << FLAGS_SD << " miter:" << FLAGS_miter
                        << " MCN:" << FLAGS_MCN << " lb:" << FLAGS_lb << " up:" << FLAGS_up
                        << " step:" << FLAGS_step << " M:" << FLAGS_M << " I:" << FLAGS_I << " P:" << FLAGS_P;
                }

                resultof << setprecision(100) << globalBest.trueFit << " actual_fx:" << actualFx
                    << " output:" << FLAGS_output
                    << " algName:" << algName << " iter:" << iter << " mcn:"
                    << FLAGS_MCN << " loop:" << FLAGS_loop << " loopIdx:" << gLoop
                    << " D:" << FLAGS_D << " K:" << FLAGS_K << " SD:" << FLAGS_SD << " miter:" << FLAGS_miter
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
        std::vector<std::string> lines;
        CHECK_RTN_LOGE_OF_FUNC(GetAllLineOfFile(FLAGS_config_path + FLAGS_test_data, lines));

        // 最后一行是权重
        std::vector<std::string> tokens = StringSplit(lines.back(), ",");
        CHECK_ASSERT_RTN((tokens.size() == FLAGS_K));
        for (int j = 0; j < tokens.size(); ++j) {
            weight[j] = std::stod(tokens[j]);
        }

        std::set<std::vector<int>> solotionSet;
        int idx = -1;
        int jdx = -1;
        for (int k = 0; k < lines.size() - 1; ++k) {
            std::vector<std::string> tokens = StringSplit(lines[k], ",");
            if (k / FLAGS_D >= FLAGS_P) {
                break;
            }
            idx = k / FLAGS_D;
            jdx = k % FLAGS_D;
            VLOG(1) << "initilize lines size:" << lines.size() << " D:" << FLAGS_D
                << " FLAGS_P:" << FLAGS_P << " K:" << FLAGS_K
                << " idx:" << idx << " jdx:" << jdx << " k:" << k
                << " group size:" << group.size() << " tokens size:"
                << tokens.size() << " line:" << lines[k];

            // 初始化 QoS 非功能属性值
            auto& qos = group[idx].services[jdx].qos;
            qos[0] = std::stod(tokens[0]); // 时间
            qos[1] = random(1, 50); // 价格
            qos[2] = std::stod(tokens[1]); // 可用性
            qos[3] = std::stod(tokens[4]); // 可靠性

            // 随机生成一个解, 如果该解已经存在, 则重新生成
            auto& solution = group[idx].solution;
            auto& rsolution = group[idx].rsolution;
            CHECK_RTN_LOGE_OF_FUNC(randomSolution(solution, rsolution));
            while (solotionSet.count(solution)) {
                CHECK_RTN_LOGE_OF_FUNC(randomSolution(solution, rsolution));
            }
            solotionSet.insert(solution);
        }
        CHECK_ASSERT_RTN((idx == FLAGS_P - 1));
        CHECK_ASSERT_RTN((jdx == FLAGS_D - 1));

        CHECK_RTN_LOGE_OF_FUNC(qosPreProcess());

        //CHECK_RTN_LOGE_OF_FUNC(topNSolution());

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
            std::vector<int> solution;
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
                return l.trueFit > r.trueFit; // 目标函数最大化
            });

        for(int i = 0; i < FLAGS_P; ++i) {
            group[i].solution = dataList[i].solution;
            group[i].trueFit = dataList[i].trueFit;
            group[i].fitness = dataList[i].fitNess;
            group[i].taskIdx = i;
            taskId2GroupIdx[i] = i;
        }
        return FICUS_SUCC;
    }

    // 数据预处理, 归一化
    int qosPreProcess()
    {
        // 计算每个子任务内所有原子服务QoS的最小和最大值
        for (int i = 0; i < FLAGS_P; ++i) {
            group[i].taskIdx = i;
            taskId2GroupIdx[i] = i;
            auto& minMaxJ = group[i].minMaxJ;
            for (int k = 0; k < FLAGS_K; ++k) {
                minMaxJ[k][0] = group[i].services[0].qos[k];
                minMaxJ[k][1] = group[i].services[0].qos[k];
            }
            for (int j = 0; j < FLAGS_D; ++j){
                auto& qos = group[i].services[j].qos;
                for (int k = 0; k < FLAGS_K; ++k) {
                    if (qos[k] > minMaxJ[k][1]) {
                        minMaxJ[k][1] = qos[k];
                    }
                    if (qos[k] < minMaxJ[k][0]) {
                        minMaxJ[k][0] = qos[k];
                    }
                    VLOG(2) << "qosPreProcess i:" << i << " j:" << j << " k:" << k
                        << " qos:" << qos[k] << " FLAGS_D:" << FLAGS_D << " FLAGS_K:"
                        << FLAGS_K << " qos size:" << qos.size();
                }
            }
        }

        // 归一化, 时间/价格:越小越好; 可靠性/可用性:越大越好
        for (int i = 0; i < FLAGS_P; ++i) {
            auto& minMaxJ = group[i].minMaxJ;
            for (int j = 0; j < FLAGS_D; ++j) {
                auto& qos = group[i].services[j].qos;
                for (int k = 0; k < FLAGS_K; ++k) {
                    auto& min = minMaxJ[k][0];
                    auto& max = minMaxJ[k][1];
                    auto diff = max - min;
                    if (k < 2) { // 时间/价格
                        qos[k] = (max - qos[k]) / diff;
                    } else {    // 可用性/可靠性
                        qos[k] = (qos[k] - min) / diff;
                    }
                    VLOG(2) << "qosPreProcess i:" << i << " j:" << j << " k:" << k
                        << " qos:" << qos[k] << " min:" << min << " max:" << max;
                }
            }
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
                return l.fitness > r.fitness;  // fitness越大, 越差, 差的放前面
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
                if (temp[i].fitness<pw[i].fitness) { // 子群最好的
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
                    if (temp[i].fitness<pw[i].fitness) { // 整体最好的
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
                        return l.fitness > r.fitness;
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
                auto taskIdx = memeplex[j][k].taskIdx;
                taskId2GroupIdx[taskIdx] = i;
                i++;
            }
        }
    }

    int getGlobalBest()
    {
        for (int i = 0; i < FLAGS_P; ++i) {
            if (group[i].trueFit > globalBest.trueFit) {
                globalBest = group[i];
                VLOG(1) <<"getGlobalBest trueFit1:" << group[i].trueFit << " trueFit2:"
                    << globalBest.trueFit << " fitness1:" << group[i].fitness << " fitness2:" << globalBest.fitness;
            }
        }
        return FICUS_SUCC;
    }

    // 计算真实的函数值 f(x)
    double calculationTruefit(Individal& individal)
    {
        auto& solution = individal.solution;
        CHECK_ASSERT_RTN(solution.size() == FLAGS_SD);

        // 时间/价格累加; 可用性/可靠性累乘;
        std::vector<double> qosSum = group[taskId2GroupIdx.at(0)].services[solution.front()].qos;
        for (int i = 1; i < solution.size(); ++i) {
            // 当前方案 solution 下子任务 i 的原子服务下标为 dIdx
            auto& dIdx = solution[i];
            auto& qos = group[taskId2GroupIdx.at(i)].services[dIdx].qos;
            for (int k = 0; k < qos.size(); ++k){
                VLOG(2) << "i:" << i << " k:" << k << " qosSum size:"
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
    std::vector<float> weight;  // Qos各维度权重

    std::vector<Individal> temp;

    std::set<std::vector<int>> filterSolution;
    std::vector<double> result;    // 截止到当前迭代轮数 对应的最优解
    std::vector<uint64_t> timeCost;  // 截止到当前迭代轮数 对应的总耗时

    std::map<int, int> taskId2GroupIdx;
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


#include <bits/stdc++.h>
#include "Fixed.h"
#include "Procesing.h"
#include "ProcesingDynamic.h"

const std::string input = "../input.txt";

template<int N, int M, class P_Type, class V_Type, class Flow_Candidate, class ...Args>
std::shared_ptr<ProcesingParent> getFlowType(const std::string &flow_type, char** field_) {
    //std::cout << "Checking Flow_Candidate: " << getName<Flow_Candidate>() << " against " << flow_type << std::endl;
    if (getName<Flow_Candidate>() == flow_type) {
        return std::make_shared<Procesing<P_Type, V_Type, Flow_Candidate, N, M>>();
    }
    else if constexpr (sizeof...(Args) == 0) {
        throw std::runtime_error("Invalid type Flow: " + flow_type);
    }
    else {
        return getFlowType<N, M, P_Type, V_Type, Args...>(flow_type, field_);
    }
}

template<int N, int M, class P_Type, class V_Candidate, class ...Args>
std::shared_ptr<ProcesingParent> getVType(const std::string &type_v, const std::string &flow_type, char** field_) {
    if constexpr (sizeof...(Args) == 0) {
        throw std::runtime_error("Bad V Type");
    }
    else if (getName<V_Candidate>() == type_v) {
        return getFlowType<N, M, P_Type, V_Candidate, Args...>(flow_type, field_);
    }
    else{
        return getVType<N, M, P_Type, Args...>(type_v, flow_type, field_);
    }
}

template<int N, int M, class P_Candidate, class ...Args>
std::shared_ptr<ProcesingParent> getPType(const std::string &p_type, const std::string &type_v, const std::string &flow_type, char** field_) {
    if constexpr (sizeof...(Args) == 0) {
        throw std::runtime_error("Bad P type");
    }
    else if (getName<P_Candidate>() == p_type) {
        return getVType<N, M, P_Candidate, Args...>(type_v, flow_type, field_);
    }
    else{
        return getPType<N, M, Args...>(p_type, type_v, flow_type, field_);
    }
}
#define S(a, b) {{a, b}, {getPType<a, b, TYPES, TYPES, TYPES>(comand_arg["p-type"], comand_arg["v-type"], comand_arg["flow-type"], actualfield_)}}


#ifndef SIZES
#define SIZES S(1920,1080),S(36,84),S(66, 66)
#endif

#ifndef TYPES
#define TYPES Fixed<54, 16, false>,Fixed<32, 16, false>,float,double
#endif

// пример использования ./fluid  --file=../input.txt --p-type="FIXED(32, 16)" --v-type=DOUBLE --flow-type="FIXED(32, 16)"

int main(int argc, char* argv[]) {
    std::map<std::string, std::string> comand_arg;
    std::regex pattern("--(p-type|v-type|flow-type|file)=(.+)");

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        std::smatch matchSearch;
        if (std::regex_match(arg, matchSearch, pattern)) {
            comand_arg[matchSearch[1]] = matchSearch[2];
        }
    }

    if (comand_arg.find("file") == comand_arg.end()) {
        std::cerr << "Usage: " << argv[0] << " --file=<filename> --p-type=<p_type> --v-type=<type_v> --flow-type=<flow_type>\n";
        return 1;
    }

    std::cout << "Parsed parameters:\n";
    for (const auto& [key, value] : comand_arg) {
        std::cout << key << ": " << value << std::endl;
    }

    int realN, realM;
    std::ifstream file(comand_arg["file"]);
    if (!file) {
        std::cerr << "Can't open file";
        return 1;
    }
    file >> realN >> realM;
    
    char** actualfield_ = new char*[realN];
    std::string line;
    std::getline(file, line);  

    for (int i = 0; i < realN; ++i) {
        actualfield_[i] = new char[realM + 1];
        std::getline(file, line);
        std::strncpy(actualfield_[i], line.c_str(), realM);
        actualfield_[i][realM] = '\0';
    }

    std::map<std::pair<int, int>, std::shared_ptr<ProcesingParent>> allFluids = {SIZES};
    if (allFluids.find({realN, realM}) == allFluids.end()) {
        auto new_fluid = Procesing<double, double, double, 0, 0>(realN, realM);
        new_fluid.start(actualfield_);
    } else {
        auto fluid = allFluids[{realN, realM}];
        fluid->start(actualfield_);
    }

    for (int i = 0; i < realN; ++i) {
        delete[] actualfield_[i];
    }
    delete[] actualfield_;

    return 0;
}
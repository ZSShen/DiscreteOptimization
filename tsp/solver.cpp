#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


class TwoOpt {
public:
    TwoOpt(const std::string& fileName);
    ~TwoOpt() = default;

    TwoOpt(const TwoOpt&) = delete;
    TwoOpt(TwoOpt&&) = delete;
    TwoOpt& operator=(const TwoOpt&) = delete;
    TwoOpt& operator=(TwoOpt&&) = delete;

public:
    void run();
    std::pair<double, std::vector<int>> getAnswer() {
        return make_pair(optTourLen_, optTour_);
    }

private:
    void swapEdge(int k, int h);

private:
    int numNode_;

    double optTourLen_;
    double curTourLen_;

    std::vector<std::pair<double, double>> map_;
    std::vector<std::vector<double>> distance_;

    std::vector<int> optTour_;
    std::vector<int> curTour_;
};


TwoOpt::TwoOpt(const std::string& fileName) {

    std::ifstream input(fileName);
    input >> numNode_;

    map_.resize(numNode_);
    distance_.resize(numNode_);
    for (int i = 0 ; i < numNode_ ; ++i) {
        distance_[i].resize(numNode_);
        input >> map_[i].first >> map_[i].second;
    }

    for (int i = 0 ; i < numNode_ ; ++i) {
        for (int j = i + 1 ; j < numNode_ ; ++j) {
            double xSqr = pow(map_[j].first - map_[i].first, 2);
            double ySqr = pow(map_[j].second - map_[i].second, 2);
            distance_[i][j] = distance_[j][i] = sqrt(xSqr + ySqr);
        }
    }

    // Construct the initial tour
    optTourLen_ = 0;
    curTour_.resize(numNode_ + 1);
    optTour_.resize(numNode_ + 1);
    optTour_[0] = 0;

    for (int i = 1 ; i < numNode_ ; ++i) {
        optTour_[i] = i;
        optTourLen_ += distance_[i][i - 1];
    }

    optTour_[numNode_] = 0;
    optTourLen_ += distance_[0][numNode_ - 1];
}

void TwoOpt::run() {

    int trial = 0;
    while (trial < numNode_) {

        for (int i = 1 ; i < numNode_ - 3 ; ++i) {
            for (int j = i + 1 ; j < numNode_ ; ++j) {
                swapEdge(i, j);
                if (curTourLen_ < optTourLen_) {
                    optTourLen_ = curTourLen_;
                    optTour_ = curTour_;
                    trial = 0;
                }
            }
        }

        ++trial;
    }
}

void TwoOpt::swapEdge(int k, int h) {

    for (int i = 0 ; i < k ; ++i) {
        curTour_[i] = optTour_[i];
    }

    for (int i = k, j = h ; i <= h ; ++i, --j) {
        curTour_[i] = optTour_[j];
    }

    for (int i = h + 1 ; i <= numNode_ ; ++i) {
        curTour_[i] = optTour_[i];
    }

    curTourLen_ = 0;
    for (int i = 1 ; i <= numNode_ ; ++i) {
        int src = curTour_[i - 1];
        int dst = curTour_[i];
        curTourLen_ += distance_[src][dst];
    }
}


int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Please specify the input data path" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    TwoOpt solver(fileName);
    solver.run();

    auto pair = solver.getAnswer();
    printf("%lf 0\n", pair.first);
    for (int i = 0 ; i < pair.second.size() - 1; ++i) {
        std::cout << pair.second[i] << ' ';
    }
    std::cout << std::endl;

    return 0;
}